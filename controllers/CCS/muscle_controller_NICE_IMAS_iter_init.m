%% CCS magnetic controller: MUSCLE3 actor init and Simulink run.
% Based on work of R. Nouailletas, T. Ravensbergen, M. Schneider and G. Gros.
%
% 1. creates the MUSCLE3 instance and receives the reference equilibrium (the inverse run)
%    and the reference coil currents;
% 2. builds the references and the parameters read by the Simulink model (base workspace);
% 3. inits the PCSSP modules and runs pcssp_controller_and_ps, which exchanges with NICE at
%    every step through muscle_nice2control / muscle_control2nice.

fprintf('Controller Initialization\n')

DEBUG_LOG = true;   % writes simulink_log.csv in the run directory (columns: sig, below)

pcs_path = getenv('PCS_PATH');
if isempty(pcs_path)
    error('PCS_PATH environment variable not set');
end
run('pcssp_pcs_add_paths')

logger = py.logging.getLogger();
setdefault(py.os.environ,'MUSCLE_INSTANCE','macro');
py.logging.basicConfig( ...
    pyargs( ...
        "format", "%(asctime)s - %(name)s - %(levelname)s - %(message)s", ...
        "level", py.logging.INFO ...
    ) ...
)

%% MUSCLE3 instance: the ports must match the workflow ymmsl
logger.info("Starting Simulink actor")
ports = py.dict();
ports{py.getattr(@py.ymmsl.Operator,"F_INIT")} = py.list({"equilibrium_in_f", "pf_active_in_f"});
ports{py.getattr(@py.ymmsl.Operator,"S")} = py.list({"equilibrium_in_s", "pf_active_in_s"});
ports{py.getattr(@py.ymmsl.Operator,"O_I")} = py.list({"pf_active_out_i"});

coder.extrinsic('py.libmuscle.Instance');
% the Simulink blocks reach the instance with evalin('base', ...)
if ~exist('instance','var') == 1
    instance = py.libmuscle.Instance(ports);
end

if instance.reuse_instance()
    %% reference equilibrium and coil currents (F_INIT)
    msg_eq = instance.receive("equilibrium_in_f");
    equilibrium = imas_deserialize(uint8(msg_eq.data), 'equilibrium');
    msg_pfa = instance.receive("pf_active_in_f");
    pf_active = imas_deserialize(uint8(msg_pfa.data), 'pf_active');

    t_start = equilibrium.time(1);
    t_max = equilibrium.time(end);

    %% settings (magnetic_controller.* in the ymmsl)
    enable_KCURR = double(instance.get_setting('KCURR', pyargs('default', 1)));
    IpControlMode_KCURR = double(instance.get_setting('control_Ip', pyargs('default', 1)));
    enable_RZIp = double(instance.get_setting('RZIp', pyargs('default', 0)));
    enable_VS = double(instance.get_setting('VS', pyargs('default', 1)));
    VSControlMode = double(instance.get_setting('VS_control_mode', pyargs('default', 3)));    % 0 off, 1 VS1, 2 VS3, 3 VS1+VS3
    ShapeIpControlMode = double(instance.get_setting('shape_control', pyargs('default', 0)));  % 2 isoflux limiter, 3 isoflux X-point, 4 gaps
    usePSModel = double(instance.get_setting('PSModel', pyargs('default', 1)));
    in_out_dt = double(instance.get_setting('in_out_msg_dt', pyargs('default', 0.001)));      % [s] exchange period with NICE = sample time of the CODAC blocks
    zdot_disturbance = double(instance.get_setting('zdot_disturbance', pyargs('default', 0)));  % noise on zdot (controller_input_preprocess)
    zdot_filter = double(instance.get_setting('zdot_filter', pyargs('default', 0)));            % low pass on zdot
    nice_signs = double(instance.get_setting('nice_signs', pyargs('default', 1)));              % 1 with NICE evo_rd: voltages sent with the opposite sign

    %% controlled gaps and isoflux control points
    GAP_LIST = [3 8 15 21 25 48];   % the gaps used by shapeIp_CCS (nonzero columns of its Proj)
    N_GAPS = 62;
    N_FLUX = 62;                    % flux inputs of shapeIp_CCS; only the first numel(RfluxCTRpts) are used
    % PLACEHOLDER: control points of the CREATE test case TC_CCS_01 (3 MA limited plasma),
    % not designed for this shot. Replace them before using mode 2 or 3.
    RfluxCTRpts = [4.362071 4.1486759 5.1894636 7.0161489 8.137644 7.0624581];   % [m]
    ZfluxCTRpts = [-1.5 1.5492 3.2382093 2.6846939 0.62839253 -2.0424846];
    R_lim = 4.1053;                 % limiter point [m]: where the plasma touches the wall in mode 2
    Z_lim = 0.6039;
    FLUX_CTR_PLACEHOLDER = true;

    %% references: one value per slice of the reference equilibrium
    n_slices = length(equilibrium.time_slice);
    ip_ref = zeros(1,n_slices);
    rpl_ref = zeros(1,n_slices);
    zpl_ref = zeros(1,n_slices);
    rxp_ref = nan(1,n_slices);      % NaN where the slice has no X-point
    zxp_ref = nan(1,n_slices);
    psixp_ref = nan(1,n_slices);
    gaps_ref = zeros(N_GAPS,n_slices);
    psi_pts_ref = zeros(numel(RfluxCTRpts), n_slices);
    psi_lim_ref = zeros(1, n_slices);

    reference_current = [];
    for i = 1:length(pf_active.coil)
        reference_current = [reference_current, pf_active.coil{i}.current.data];
    end
    for i = 1:n_slices
        ts = equilibrium.time_slice{i};
        cc = ts.global_quantities.current_centre;
        ip_ref(i) = ts.global_quantities.ip;
        rpl_ref(i) = cc.r;
        zpl_ref(i) = cc.z;
        % gaps and X-point come from the diagnostics actor: without it in the inverse workflow
        % the gaps are missing and the X-point is a mesh node (jumps of 1-3 cm)
        gap = ts.boundary.gap;
        for k = GAP_LIST
            gaps_ref(k, i) = NaN;
            if numel(gap) >= k && ~isempty(gap{k}.value) && abs(gap{k}.value) < 1e30
                gaps_ref(k, i) = gap{k}.value;
            end
        end
        for k = 1:numel(ts.contour_tree.node)
            node = ts.contour_tree.node{k};
            if node.critical_type == 1
                rxp_ref(i) = node.r;
                zxp_ref(i) = node.z;
                psixp_ref(i) = node.psi;
                break
            end
        end
        % psi at the isoflux points, interpolated on NICE's mesh
        g = ts.ggd{1};
        F = scatteredInterpolant(g.r{1}.values(:), g.z{1}.values(:), g.psi{1}.values(:));
        psi_pts_ref(:, i) = F(RfluxCTRpts(:), ZfluxCTRpts(:));
        psi_lim_ref(i) = F(R_lim, Z_lim);
    end

    % slices without separatrix or X-point take the value of the nearest slice
    missing = isnan(gaps_ref(GAP_LIST, :));
    if any(missing(:))
        logger.warning(sprintf('gap references: %d values not found, taken from the nearest slice', nnz(missing)));
        gaps_ref(GAP_LIST, :) = fillmissing(gaps_ref(GAP_LIST, :), 'nearest', 2);
    end
    xp_missing = isnan(rxp_ref);
    if any(xp_missing)
        logger.warning(sprintf('X-point references: none in %d slices, taken from the nearest slice', nnz(xp_missing)));
        rxp_ref = fillmissing(rxp_ref, 'nearest');
        zxp_ref = fillmissing(zxp_ref, 'nearest');
        psixp_ref = fillmissing(psixp_ref, 'nearest');
    end

    % isoflux reference: psi(control point) - psi(boundary), where the boundary flux is the
    % limiter point's in mode 2 and the X-point's in mode 3; 0 when the point is on the boundary
    if FLUX_CTR_PLACEHOLDER && any(ShapeIpControlMode == [2 3])
        error('isoflux control points are placeholders: set the right ones before mode 2 or 3');
    end
    if ShapeIpControlMode == 3 && all(xp_missing)
        error('mode 3: no X-point in the reference equilibrium');
    end
    n_ctr = numel(RfluxCTRpts);
    dpsi_ref = zeros(N_FLUX, n_slices);
    if ShapeIpControlMode == 2
        dpsi_ref(1:n_ctr, :) = psi_pts_ref - psi_lim_ref;
    elseif ShapeIpControlMode == 3
        dpsi_ref(1:n_ctr, :) = psi_pts_ref - psixp_ref;
    end

    % read by the From Workspace blocks of the model
    Ipl_ref = timeseries(ip_ref, equilibrium.time);
    Rpl_ref = timeseries(rpl_ref, equilibrium.time);
    Zpl_ref = timeseries(zpl_ref, equilibrium.time);
    Rxp_ref = timeseries(rxp_ref, equilibrium.time);
    Zxp_ref = timeseries(zxp_ref, equilibrium.time);
    Gaps_ref = timeseries(gaps_ref, equilibrium.time);
    dPsi_ref = timeseries(dpsi_ref, equilibrium.time);
    CSPF_curr_ref = timeseries(reference_current, equilibrium.time);

    %% fixed controller parameters (Constant blocks of the model)
    schedulingVar_KCURR = 0;
    IpControlMode_RZIp = 0;         % schedulingVar_RZIp = [t, li3, elongation] comes from the measurements
    VSIpScaling = 1;                % VS action scaled with: 1 nominal Ip, 2 estimated Ip
    SchedulingVar_ShapeContr = VSControlMode;   % gain region of shapeIp_CCS (its two regions are identical today)
    SCIpScaling = 1;                % shape action scaled with: 0 nominal Ip, 1 Ip_ref, 2 Ip_est
    IpControlMode = 0;              % Ip control inside the CCS: 0 off, 1 with KCURR, 2 with the shape controller
    AllocatorMode = 1;

    logger.info(sprintf('exchange period in_out_dt = %g s', in_out_dt));

    %% PCSSP modules
    obj = pcssp_KCURR_PFPO1_obj();
    obj.init;
    obj.setup;

    obj2 = RZIp_CCS_wrapper_obj();
    obj2.init;
    obj2.setup;

    obj3 = pcssp_VS_obj();
    obj3.init;
    % VS gains from FC016 instead of the module's defaults
    S = load(fullfile(getenv('CONTROLLER_PATH'), 'params/FC016_Config_VS.mat'), 'Config_VS');
    obj3.update_nominal_param_value('VS_fp', ccs_vs_fp(obj3, S.Config_VS));
    fprintf('VS gains from FC016_Config_VS.mat (KZ_VS1 = %g)\n', S.Config_VS.Gains.KZ_VS1);
    obj3.setup;

    obj4 = shapeIp_CCS_obj();
    obj4.init;
    obj4.setup;

    obj5 = VS1_allocator_obj();
    obj5.init;
    obj5.setup;

    obj6 = pcssp_PowerSupplies_obj();
    obj6.init;
    obj6.setup;

    % data dictionary of the top model: links the modules' dictionaries and the base workspace
    top_dd_file = fullfile(fileparts(which('pcssp_controller_and_ps.slx')), 'pcssp_controller_and_ps.sldd');
    if isfile(top_dd_file)
        top_dd = Simulink.data.dictionary.open(top_dd_file);
    else
        top_dd = Simulink.data.dictionary.create(top_dd_file);
    end
    if ~top_dd.EnableAccessToBaseWorkspace
        top_dd.EnableAccessToBaseWorkspace = true;
    end
    for dd_source = {'pcssp_KCURR_PFPO1.sldd', 'pcssp_PowerSupplies.sldd'}
        if ~any(strcmp(top_dd.DataSources, dd_source{1}))
            addDataSource(top_dd, dd_source{1});
        end
    end
    if top_dd.HasUnsavedChanges
        saveChanges(top_dd);
    end
    close(top_dd);
    load_system('pcssp_controller_and_ps');
    set_param('pcssp_controller_and_ps', 'DataDictionary', 'pcssp_controller_and_ps.sldd');

    %% simulation
    in = Simulink.SimulationInput('pcssp_controller_and_ps');
    in = in.setModelParameter('StartTime', num2str(t_start), 'StopTime', num2str(t_max));
    in = in.setModelParameter('SaveOutput', 'on');

    if DEBUG_LOG
        % header of simulink_log.csv: the inputs of the logger Mux, in its order, with their widths
        csv_file = fullfile(getenv('RUN_DIR'), 'simulink_log.csv');
        if exist(csv_file, 'file'), delete(csv_file); end
        sig = {
            'KCURR_CSPF_ctrl_e', 11;
            'KCURR_Ip_ctrl_e', 1;
            'KCURR_residual_ctrl_e', 11;
            'KCURR_actual_ref', 11;
            'RZIp_CSPF_curr_cmd', 11;
            'RZIp_VS1_cmd', 1;
            'VS_VS1_cmd', 1;
            'VS_VS3_cmd', 1;
            'Shape_CSPF_curr_cmd', 11;
            'Shape_VS_decoupling', 11;
            'Output_total', 14;
            'CSPF_curr_ref', 11;
            'Zdot', 1;
            'Ip_ref', 1;
            'Z_ref', 1;
            'R_ref', 1;
            'ip', 1;
            'zpl', 1;
            'gaps', numel(GAP_LIST);
            'gaps_ref', numel(GAP_LIST);
            'rxp', 1;
            'zxp', 1;
            'rxp_ref', 1;
            'zxp_ref', 1;
            };
        h = {'t'};
        for k = 1:size(sig,1)
            if sig{k,2} == 1
                h{end+1} = sig{k,1};
            else
                for j = 1:sig{k,2}
                    h{end+1} = sprintf('%s_%d', sig{k,1}, j);
                end
            end
        end
        f = fopen(csv_file, 'w');  fprintf(f, '%s\n', strjoin(h, ','));  fclose(f);
    end

    fprintf('Start closed loop\n')
    tic;
    logs = sim(in);
    close_system('pcssp_controller_and_ps', 0);
    simulation_time = toc;
    plasma_duration = t_max - t_start;
    fprintf(['\nControl over! ' num2str(simulation_time) 's of simulation for ' num2str(plasma_duration) 's of plasma.\n']);
end
instance.reuse_instance();
