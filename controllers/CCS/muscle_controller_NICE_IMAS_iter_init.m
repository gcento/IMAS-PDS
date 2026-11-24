%% CCS magnetic controller: MUSCLE3 actor that runs the PCSSP controller model against NICE.
% Based on work of R. Nouailletas, T. Ravensbergen, M. Schneider, G. Gros and G. Centonze
%
% 1. receives the reference equilibrium and coil currents (F_INIT, from the inverse run);
% 2. writes to the base workspace what the model reads: settings, references, fixed parameters;
% 3. inits the PCSSP modules and runs the model, which exchanges with NICE at every step through its CODAC
%    blocks (muscle_nice2control.m, muscle_control2nice.m). The model runs from t = 0 on NICE's steps:
%    plasma time = Simulink time + t_start.
%
% Options: psc_enabled = 1 configures the PSC block (pulse supervision) from psc_config_dir. 
%          mcls_enabled = 1 (needs psc_enabled) runs pcssp_controller_and_ps_mcls: the same model plus the MCLS

fprintf('Controller Initialization\n')
DEBUG_LOG = true;       % writes simulink_log.csv in the run directory (columns: sig, below)

if isempty(getenv('PCS_PATH'))
    error('PCS_PATH environment variable not set');
end
run('pcssp_pcs_add_paths')

logger = py.logging.getLogger();
setdefault(py.os.environ, 'MUSCLE_INSTANCE', 'macro');
py.logging.basicConfig(pyargs("format", "%(asctime)s - %(name)s - %(levelname)s - %(message)s", "level", py.logging.INFO));

%% MUSCLE3 instance: ports as in the workflow ymmsl; the Simulink blocks reach it with evalin('base', 'instance')
logger.info("Starting Simulink actor")
ports = py.dict();
ports{py.getattr(@py.ymmsl.Operator, "F_INIT")} = py.list({"equilibrium_in_f", "pf_active_in_f"});
ports{py.getattr(@py.ymmsl.Operator, "S")} = py.list({"equilibrium_in_s", "pf_active_in_s"});
ports{py.getattr(@py.ymmsl.Operator, "O_I")} = py.list({"pf_active_out_i"});
if ~exist('instance', 'var')
    instance = py.libmuscle.Instance(ports);
end

if instance.reuse_instance()
    %% reference equilibrium and coil currents
    msg_eq = instance.receive("equilibrium_in_f");
    equilibrium = imas_deserialize(uint8(msg_eq.data), 'equilibrium');
    msg_pfa = instance.receive("pf_active_in_f");
    pf_active = imas_deserialize(uint8(msg_pfa.data), 'pf_active');
    t_start = equilibrium.time(1);  % Simulink time 0 (log_row.m, muscle_nice2control.m)
    t_max = equilibrium.time(end);  % stop time
    t_offset = t_start;             % the same, for mcls_log_row.m and the MCLS clock (MCLS_clock_offset)

    %% settings (magnetic_controller.* in the ymmsl)
    setting = @(name, default) double(instance.get_setting(name, pyargs('default', default)));
    enable_KCURR = setting('KCURR', 1);
    IpControlMode_KCURR = setting('control_Ip', 1);     % 1: KCURR gains with Ip control, 0: without
    enable_RZIp = setting('RZIp', 0);
    enable_VS = setting('VS', 1);
    VSControlMode = setting('VS_control_mode', 3);      % 0 off, 1 VS1, 2 VS3, 3 VS1+VS3
    ShapeIpControlMode = setting('shape_control', 0);   % 2 isoflux limiter, 3 isoflux X-point, 4 gaps
    usePSModel = setting('PSModel', 1);                 % 1: commands through the power supply model
    in_out_dt = setting('in_out_msg_dt', 0.001);        % [s] exchange period with NICE = sample time of the CODAC blocks
    zdot_disturbance = setting('zdot_disturbance', 0);  % noise on zdot (controller_input_preprocess.m)
    zdot_filter = setting('zdot_filter', 0);            % low pass on zdot (controller_input_preprocess.m)
    nice_signs = setting('nice_signs', 1);              % 1 with NICE evo_rd: voltages sent with the opposite sign
    psc_enabled = setting('psc_enabled', 0);            % 1: PSC block in the loop, Ipl_ref from the PSC
    if psc_enabled
        psc_config_dir = char(instance.get_setting('psc_config_dir'));   % folder with the four PSC files
    end
    mcls_enabled = setting('mcls_enabled', 0);
    logger.info(sprintf('exchange period in_out_dt = %g s', in_out_dt));

    controller_model = 'pcssp_controller_and_ps';
    if mcls_enabled
        if ~psc_enabled
            error('muscle_controller_NICE_IMAS_iter_init:mcls', 'mcls_enabled = 1 needs psc_enabled = 1');
        end
        controller_model = 'pcssp_controller_and_ps_mcls';
        mcls_dir = fullfile(getenv('CONTROLLER_PATH'), 'mcls');
        addpath(mcls_dir, fullfile(mcls_dir, 'third_party'));
        mcls_config = char(instance.get_setting('mcls_config'));    % mode-table YAML (required)
        mcls_shape_enable = setting('mcls_shape_enable', 0);       % 1: ShapeIp modes from the MCLS
        % mode table = the four Parameter data of the MCLS block
        [MCMode, StartingMode, NAlarms, NEnablingConditions, mcls_info] = mcls_load_table(mcls_config, ...
            fullfile(mcls_dir, 'mcls_library_defaults.yaml'), fullfile(getenv('RUN_DIR'), 'mcls_diary.txt'));
        logger.info(sprintf('MCLS table %s: %d mode(s) [%s], StartingMode %d, shape guard enable %d; model %s', ...
            mcls_config, numel(MCMode), strjoin(strtrim(mcls_info.MCModeNames), ', '), StartingMode, mcls_shape_enable, controller_model));
    end

    %% shape control: controlled gaps, isoflux control points, limiter point
    GAP_LIST = [3 8 15 21 25 48];   % gaps used by shapeIp_CCS (nonzero columns of its Proj)
    N_GAPS = 62;                    % gap inputs of shapeIp_CCS
    N_FLUX = 62;                    % flux inputs of shapeIp_CCS, the first numel(RfluxCTRpts) used
    % default points: RR3/ZZ3 of the pcs XPT02 case (copy_setup_pcs_params_pcs_expert.m:466-494) without its 9th
    % point, the X-point: here the mode-3 reference flux is the measured psi_xp (muscle_nice2control.m)
    RfluxCTRpts = [4.3863 4.2220 5.1385 7.2929 8.1866 7.0307 4.6170 5.2594];   % [m]
    ZfluxCTRpts = [-1.5 1.5492 3.5094 2.9244 0.6292 -2.0128 -3.6807 -3.9861];
    if ShapeIpControlMode == 2      % limited plasma: the first six points and the touch point (RR2/ZZ2)
        RfluxCTRpts = [RfluxCTRpts(1:6) 4.1000];
        ZfluxCTRpts = [ZfluxCTRpts(1:6) 0.7500];
    end
    R_lim = 4.1053;                 % [m] limiter point: PsiLimiter outside mode 2 (muscle_nice2control.m)
    Z_lim = 0.6039;
    % the settings isoflux_r/isoflux_z replace the control points; limiter = [R Z] replaces the limiter point,
    % where mode 2 then takes the boundary flux instead of psi_boundary
    isoflux_r = setting('isoflux_r', 0);
    isoflux_z = setting('isoflux_z', 0);
    if numel(isoflux_r) > 1 || any(isoflux_r ~= 0)
        if numel(isoflux_r) ~= numel(isoflux_z) || numel(isoflux_r) > N_FLUX
            error('isoflux_r and isoflux_z: %d and %d points, they need the same number (at most %d)', ...
                  numel(isoflux_r), numel(isoflux_z), N_FLUX);
        end
        RfluxCTRpts = isoflux_r(:)';
        ZfluxCTRpts = isoflux_z(:)';
        logger.info(sprintf('isoflux control points from the settings isoflux_r/isoflux_z: %d points', numel(RfluxCTRpts)));
    end
    % shapeIp_CCS's Proj uses the flux columns 1-6 in mode 2, 1-8 in mode 3: a missing point would see the
    % boundary flux itself as error (tens of Wb)
    n_min = 2*ShapeIpControlMode + 2;
    if any(ShapeIpControlMode == [2 3]) && numel(RfluxCTRpts) < n_min
        error('shape mode %d needs at least %d isoflux points, %d given', ShapeIpControlMode, n_min, numel(RfluxCTRpts));
    end
    limiter = setting('limiter', 0);
    lim_point = numel(limiter) >= 2;    % true: the mode-2 boundary flux is psi at (R_lim, Z_lim)
    if lim_point
        R_lim = limiter(1);
        Z_lim = limiter(2);
        logger.info(sprintf('limiter point from the settings: R = %g m, Z = %g m (mode 2 boundary flux)', R_lim, Z_lim));
    end

    %% references: one value per slice of the reference equilibrium, in the controller frame
    n_slices = length(equilibrium.time_slice);
    ip_ref = zeros(1, n_slices);
    rpl_ref = zeros(1, n_slices);
    zpl_ref = zeros(1, n_slices);
    rxp_ref = nan(1, n_slices);     % NaN where the slice has no X-point
    zxp_ref = nan(1, n_slices);
    psixp_ref = nan(1, n_slices);
    gaps_ref = zeros(N_GAPS, n_slices);
    psi_pts_ref = zeros(numel(RfluxCTRpts), n_slices);
    psi_lim_ref = zeros(1, n_slices);

    ip_sign = sign(equilibrium.time_slice{1}.global_quantities.ip);
    reference_current = [];
    for i = 1:length(pf_active.coil)
        reference_current = [reference_current, pf_active.coil{i}.current.data]; %#ok<AGROW>
    end
    reference_current = ip_sign * reference_current;

    for i = 1:n_slices
        ts = equilibrium.time_slice{i};
        cc = ts.global_quantities.current_centre;
        ip_ref(i) = abs(ts.global_quantities.ip);
        rpl_ref(i) = cc.r;
        zpl_ref(i) = cc.z;
        % gaps and X-point from the diagnostics actor, which writes every gap in every slice (NaN where not found)
        for k = GAP_LIST
            gaps_ref(k, i) = ts.boundary.gap{k}.value;
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
        % psi at the control points and at the limiter point, interpolated on NICE's mesh
        g = ts.ggd{1};
        F = scatteredInterpolant(g.r{1}.values(:), g.z{1}.values(:), g.psi{1}.values(:));
        psi_pts_ref(:, i) = F(RfluxCTRpts(:), ZfluxCTRpts(:));
        if lim_point
            psi_lim_ref(i) = F(R_lim, Z_lim);
        else
            psi_lim_ref(i) = ts.global_quantities.psi_boundary;   % limited plasma: touch-point flux = boundary flux
        end
    end

    % slices without X-point take the nearest slice's; with no X-point in any slice (a limited window) they are
    % 0, since the From Workspace blocks reject NaN and only mode 3 needs them
    xp_missing = isnan(rxp_ref);
    if all(xp_missing)
        if ShapeIpControlMode == 3
            error('mode 3: no X-point in the reference equilibrium');
        end
        logger.warning(sprintf('X-point references: no X-point in the reference, set to 0 (not used in mode %d)', ...
                               ShapeIpControlMode));
        rxp_ref(:) = 0;
        zxp_ref(:) = 0;
        psixp_ref(:) = 0;
    elseif any(xp_missing)
        logger.warning(sprintf('X-point references: none in %d slices, taken from the nearest slice', nnz(xp_missing)));
        rxp_ref = fillmissing(rxp_ref, 'nearest');
        zxp_ref = fillmissing(zxp_ref, 'nearest');
        psixp_ref = fillmissing(psixp_ref, 'nearest');
    end

    % isoflux reference: psi(control point) - boundary flux, 0 when the point is on the boundary. Boundary flux:
    % psi_lim_ref in mode 2, the X-point's in mode 3
    n_ctr = numel(RfluxCTRpts);
    psi_pts_ref = ip_sign * psi_pts_ref;
    psi_lim_ref = ip_sign * psi_lim_ref;
    psixp_ref = ip_sign * psixp_ref;
    dpsi_ref = zeros(N_FLUX, n_slices);
    if ShapeIpControlMode == 2
        dpsi_ref(1:n_ctr, :) = psi_pts_ref - psi_lim_ref;
    elseif ShapeIpControlMode == 3
        dpsi_ref(1:n_ctr, :) = psi_pts_ref - psixp_ref;
    end

    % From Workspace blocks of the model, on the Simulink time (they reject NaN: gaps not found are 0)
    t_ref = equilibrium.time - t_start;
    gaps_ref(isnan(gaps_ref)) = 0;
    Ipl_ref = timeseries(ip_ref, t_ref);
    Rpl_ref = timeseries(rpl_ref, t_ref);
    Zpl_ref = timeseries(zpl_ref, t_ref);
    Rxp_ref = timeseries(rxp_ref, t_ref);
    Zxp_ref = timeseries(zxp_ref, t_ref);
    Gaps_ref = timeseries(gaps_ref, t_ref);
    dPsi_ref = timeseries(dpsi_ref, t_ref);
    CSPF_curr_ref = timeseries(reference_current, t_ref);

    %% fixed parameters (Constant blocks on the module inputs)
    schedulingVar_KCURR = 0;                    % KCURR gain region
    IpControlMode_RZIp = 0;                     % RZIp without Ip control (its scheduling variable is measured)
    RZIp_R_sign = setting('RZIp_R_sign', 1);    % signs of the R and Z errors fed to RZIp_CCS, 1 in the controller frame
    RZIp_Z_sign = setting('RZIp_Z_sign', RZIp_R_sign);
    VSIpScaling = 1;                            % VS action scaled with 1 Ip_ref, 2 Ip_est
    SchedulingVar_ShapeContr = VSControlMode;   % gain region of shapeIp_CCS (its two regions are identical today)
    SCIpScaling = setting('SC_Ip_scaling', 1);  % shape action scaled with 0 nominal Ip, 1 Ip_ref, 2 Ip_est
    xp_flux_at_ref = setting('xp_flux_at_ref', 0);   % 1: mode-3 psi_xp at the reference X-point (muscle_nice2control.m)
    if xp_flux_at_ref && ShapeIpControlMode == 3
        logger.info(['xp_flux_at_ref = 1: psi_xp read at the reference X-point; without a saddle in NICE, ' ...
                     'the reference X-point also stands in for the measured position']);
    end
    IpControlMode = 1;                          % Ip control (shapeIp_CCS input): 0 off, 1 with KCURR, 2 with the shape controller
    AllocatorMode = 1;                          % VS1_allocator: 1 VS1 command onto PF2-PF5, 0 CS/PF commands unchanged

    % for the analysis of a run: isoflux settings and references; isoflux_log.csv gets one row per step
    % (muscle_nice2control.m)
    save(fullfile(getenv('RUN_DIR'), 'isoflux_refs.mat'), 'RfluxCTRpts', 'ZfluxCTRpts', 'dpsi_ref', 'rxp_ref', 'zxp_ref', ...
         't_ref', 't_start', 'ip_sign', 'SCIpScaling', 'xp_flux_at_ref', 'RZIp_R_sign', 'RZIp_Z_sign', 'R_lim', 'Z_lim', 'lim_point');
    f = fopen(fullfile(getenv('RUN_DIR'), 'isoflux_log.csv'), 'w');
    fprintf(f, 't,%s,psi_xp,r_xp,z_xp,zpl,ip,psi_lim\n', strjoin(compose('flux_%d', 1:8), ','));
    fclose(f);

    %% PCSSP modules: init creates the module's data dictionary with its nominal parameters, setup links the
    % dictionaries of its referenced models
    obj = pcssp_KCURR_PFPO1_obj();      obj.init;   obj.setup;      % coil currents (and Ip)
    obj2 = RZIp_CCS_wrapper_obj();      obj2.init;  obj2.setup;     % plasma position and current
    obj3 = pcssp_VS_obj();              obj3.init;                  % vertical stabilisation, with the FC016 gains:
    S = load(fullfile(getenv('CONTROLLER_PATH'), 'params/FC016_Config_VS.mat'), 'Config_VS');
    fp = ccs_vs_fp(obj3, S.Config_VS);
    if isprop(obj3, 'ddHdl') && ~isempty(obj3.ddHdl) && ismethod(obj3.ddHdl, 'updateField')
        obj3.ddHdl.updateField('VS_fp', fp);    % pcssp-nightly: update_nominal_param_value calls the missing ddHdl.updateValue
    else
        obj3.update_nominal_param_value('VS_fp', fp);
    end
    fprintf('VS gains from FC016_Config_VS.mat (KZ_VS1 = %g)\n', S.Config_VS.Gains.KZ_VS1);
    obj3.setup;
    obj4 = shapeIp_CCS_obj();           obj4.init;  obj4.setup;     % plasma shape: isoflux or gaps
    obj5 = VS1_allocator_obj();         obj5.init;  obj5.setup;     % VS1 command onto PF2-PF5
    obj6 = pcssp_PowerSupplies_obj();   obj6.init;  obj6.setup;     % power supply model
    SCDconf_setConf('pcssp_Simulation', 'configurations_container_pcssp.sldd');   % 'configurationSettings' for the RZIp wrapper

    %% model: its data dictionary links the modules' dictionaries and the base workspace
    top_dd_file = fullfile(fileparts(which([controller_model '.slx'])), [controller_model '.sldd']);
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
    set_param(0, 'ErrorIfLoadNewModel', 'off');
    load_system(controller_model);
    set_param(controller_model, 'DataDictionary', [controller_model '.sldd']);

    if psc_enabled      % the PSC block (Subsystem Reference to PSC) gets its four files
        psc_blk = find_system(controller_model, 'LookUnderMasks', 'all', 'ReferencedSubsystem', 'PSC');
        if isempty(psc_blk)
            error('psc_enabled = 1 but the model has no PSC block (Subsystem Reference to PSC)');
        end
        psc_setup_block(psc_blk{1}, psc_config_dir);
        logger.info(sprintf('PSC configured from %s', psc_config_dir));
    end

    %% log headers: the rows are appended during the run
    if mcls_enabled     % mcls_log.csv (mcls_log_row.m): the MCLS_mux signals of the model, then alarms and enabling conditions
        mcls_csv = fullfile(getenv('RUN_DIR'), 'mcls_log.csv');
        if exist(mcls_csv, 'file'), delete(mcls_csv); end
        hdr = [{'t_abs','t_sim','ControlModeIndex','ShapeIpControlMode','IpControlMode','VSControlMode','SCIpScaling', ...
                'VSIpScaling','SNUControl','CurrentLimitAvoidance','VS3_kick_volt_on','RefSmoothTrans','HL_contr_on', ...
                'LH_contr_on','FWHLP_contr_on','EmShutDown_contr_on','ForceControlMode','Ipl_est_pcs','mcls_time'}, ...
               arrayfun(@(k) sprintf('alarm_%d', k), 1:NAlarms, 'UniformOutput', false), ...
               arrayfun(@(k) sprintf('cond_%d', k), 1:NEnablingConditions, 'UniformOutput', false)];
        f = fopen(mcls_csv, 'w');  fprintf(f, '%s\n', strjoin(hdr, ','));  fclose(f);
    end
    if DEBUG_LOG        % simulink_log.csv (log_row.m): the inputs of the logger Mux, in its order, with their widths
        csv_file = fullfile(getenv('RUN_DIR'), 'simulink_log.csv');
        if exist(csv_file, 'file'), delete(csv_file); end
        sig = {'KCURR_CSPF_ctrl_e', 11; 'KCURR_Ip_ctrl_e', 1; 'KCURR_residual_ctrl_e', 11; 'KCURR_actual_ref', 11
               'RZIp_CSPF_curr_cmd', 11; 'RZIp_VS1_cmd', 1; 'VS_VS1_cmd', 1; 'VS_VS3_cmd', 1
               'Shape_CSPF_curr_cmd', 11; 'Shape_VS_decoupling', 11; 'Output_total', 14; 'CSPF_curr_ref', 11
               'Zdot', 1; 'Ip_ref', 1; 'Z_ref', 1; 'R_ref', 1; 'ip', 1; 'zpl', 1
               'gaps', numel(GAP_LIST); 'gaps_ref', numel(GAP_LIST); 'rxp', 1; 'zxp', 1; 'rxp_ref', 1; 'zxp_ref', 1
               'flux_ctrl', 8; 'flux_lim', 1; 'flux_xp', 1
               'psc_segment', 1; 'psc_ip_too_high', 1; 'psc_end_of_flux', 1; 'psc_z_excursion', 1};
        h = cell(1, size(sig, 1));      % columns of each signal: its name, or name_1 ... name_n
        for k = 1:size(sig, 1)
            h{k} = sig(k, 1);
            if sig{k,2} > 1
                h{k} = compose([sig{k,1} '_%d'], 1:sig{k,2});
            end
        end
        h = [{'t'}, h{:}];
        f = fopen(csv_file, 'w');  fprintf(f, '%s\n', strjoin(h, ','));  fclose(f);
    end

    %% simulation
    in = Simulink.SimulationInput(controller_model);
    in = in.setModelParameter('StartTime', '0', 'StopTime', num2str(t_max - t_start, 12));
    in = in.setModelParameter('SaveOutput', 'on');
    fprintf('Start closed loop\n')
    tic;
    logs = sim(in);
    close_system(controller_model, 0);
    simulation_time = toc;
    plasma_duration = t_max - t_start;
    fprintf(['\nControl over! ' num2str(simulation_time) 's of simulation for ' num2str(plasma_duration) 's of plasma.\n']);
end
instance.reuse_instance();      % MUSCLE3: this second call returns false and ends the instance
