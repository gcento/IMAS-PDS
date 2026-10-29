%% muscle controller initialization script. 
% 
% based on work of R. Nouailletas, T. Ravensbergen et M. Schneider, G. Gros

% This script initializes MUSCLE3 (logging, ports),
% initializes the controller params, and then launches a simulation.
% The controller receives a line_avg_density via MUSCLE3, and fires a
% pellet when the integral of the requested particle flow is larger than
% the size of a single pellet

fprintf('Controller Initialization\n')

DEBUG_LOG = true;
% run('../../initNice4matlab.m')
% run('../../../../pcs/pcssp_pcs_add_paths.m')
% get PCS stuff
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

%% Create libmuscle instance, specify all Ports
% Ports need to match the MUSCLE3 yaml description on the other side
logger.info("Starting Simulink actor")
ports = py.dict();
ports{py.getattr(@py.ymmsl.Operator,"F_INIT")} = py.list({"equilibrium_in_f", "pf_active_in_f"});
ports{py.getattr(@py.ymmsl.Operator,"S")} = py.list({"equilibrium_in_s", "pf_active_in_s"});
ports{py.getattr(@py.ymmsl.Operator,"O_I")} = py.list({"pf_active_out_i"});

%% declare muscle3 as non-codegen compatible
coder.extrinsic('py.libmuscle.Instance');

% write instance object to matlab workspace. Currently, there is no other
% way than to retrieve it from there within simulink with nasty 'evalin'
% statements.
% only populate when not yet there
if ~exist('instance','var') == 1
    instance = py.libmuscle.Instance(ports);
end

% This workflow's F_INIT ports (equilibrium_in_f/pf_active_in_f) are fed by
% waveform_editor exactly once, so there is never a second cycle. A `while`
% here calls reuse_instance() again after the one real cycle, trying to
% start a new F_INIT round against peers (nice_evo_rd) that have already
% exited -- it hangs retrying the dead connection for ~5 minutes, then
% crashes with "OSError: Bad file descriptor" instead of exiting cleanly
% (confirmed live). `if` runs the one cycle that exists and stops there.
if instance.reuse_instance()
    %% Prepare IDS Python object
    % equilibrium_python=py.imas.equilibrium();
    equilibrium_python=ids_init('equilibrium');
    % equilibrium_python=ids_gen_allocate(equilibrium_python, 'equilibrium', '');
    % equilibrium_python=py.imas.IDSFactory().equilibrium();


    %% Init NICE parameters
    msg_eq = instance.receive("equilibrium_in_f");
    equilibrium_serial=uint8(msg_eq.data);
    equilibrium = imas_deserialize(equilibrium_serial,'equilibrium');

    msg_pfa = instance.receive("pf_active_in_f");
    pfa_serial=uint8(msg_pfa.data);
    pf_active = imas_deserialize(pfa_serial, 'pf_active');

    t_start = equilibrium.time(1);
    t_max = equilibrium.time(end);

    %% call simulink model inits
    reference_current=[];
    for i=1:length(pf_active.coil)
        reference_current=[reference_current, pf_active.coil{i}.current.data];
    end

    %Init Simulink object

    % Ipl_ref/Rpl_ref/Zpl_ref track the F_INIT equilibrium (the un-controlled
    % target trajectory) over the whole run, same as CSPF_curr_ref above --
    % fed to the model via FromWorkspace blocks (not Constant), so they vary
    % over the simulation instead of pinning to a single scalar.
    n_slices = length(equilibrium.time_slice);
    ip_ref = zeros(1,n_slices);
    rgeo_ref = zeros(1,n_slices);
    zgeo_ref = zeros(1,n_slices);
    
    N_GAPS = 62;
    gaps_ref = zeros(N_GAPS,n_slices); %TODO: update gaps
    
    % Plasma position read by the controller, for the R/Z references here and the
    % measurement in muscle_IDS_NICE_output: setting plasma_position (in a case override:
    % magnetic_controller.plasma_position), 'geometric_axis' (boundary geometric axis,
    % default) or 'magnetic_axis' (global_quantities.magnetic_axis: the linear plant models
    % move the magnetic axis, while their geometric axis is the reference passed through).
    plasma_position = char(instance.get_setting('plasma_position', pyargs('default', 'geometric_axis')));
    if ~any(strcmp(plasma_position, {'geometric_axis','magnetic_axis'}))
        error('muscle_controller_NICE_IMAS_iter_init:plasma_position', ...
            'plasma_position must be ''geometric_axis'' or ''magnetic_axis'' (got ''%s'').', plasma_position);
    end
    logger.info(sprintf('plasma position = %s', plasma_position));

    for i=1:n_slices
        ip_ref(i)   = equilibrium.time_slice{i}.global_quantities.ip;
        if strcmp(plasma_position, 'magnetic_axis')
            rgeo_ref(i) = equilibrium.time_slice{i}.global_quantities.magnetic_axis.r;
            zgeo_ref(i) = equilibrium.time_slice{i}.global_quantities.magnetic_axis.z;
        else
            rgeo_ref(i) = equilibrium.time_slice{i}.boundary.geometric_axis.r;
            zgeo_ref(i) = equilibrium.time_slice{i}.boundary.geometric_axis.z;
            rgeo_ref(i) = geo_ref_with_fallback(rgeo_ref(i), equilibrium.time_slice{i}, 'r', logger);
            zgeo_ref(i) = geo_ref_with_fallback(zgeo_ref(i), equilibrium.time_slice{i}, 'z', logger);
        end
    end

    Ipl_ref =  timeseries(ip_ref,equilibrium.time);
    Rpl_ref = timeseries(rgeo_ref,equilibrium.time);
    Zpl_ref = timeseries(zgeo_ref,equilibrium.time);
    Gaps_ref = timeseries(gaps_ref, equilibrium.time);
    CSPF_curr_ref = timeseries(reference_current,equilibrium.time);

    % setting hold_references (0 = off, default): every reference (CS/PF currents, Ip, R, Z,
    % gaps) stays at the first slice, the operating point of the linear plant models
    if double(instance.get_setting('hold_references', pyargs('default', 0)))
        t_hold = [t_start, max(t_max, t_start) + 1];
        Ipl_ref = timeseries(ip_ref([1 1]), t_hold);
        Rpl_ref = timeseries(rgeo_ref([1 1]), t_hold);
        Zpl_ref = timeseries(zgeo_ref([1 1]), t_hold);
        Gaps_ref = timeseries(gaps_ref(:, [1 1]), t_hold);
        CSPF_curr_ref = timeseries(reference_current([1 1], :), t_hold);

    end

    % Test reference, e.g. for the power-supply model: setting ref_trapezoid_amplitude [A]
    % (0 = off, default). The CS/PF references stay at the first slice (t0, the
    % linearisation point) and every CS/PF coil (the 12 before VS3U/VS3L) gets the same
    % trapezoid from t0: ramp 0.5 s, flat 0.5 s, ramp down 0.5 s.

        %t_trap = t_start + [0 1.5 3 4.5];
        %t_trap(end+1) = max(t_max, t_trap(end)) + 1;          % flat after the trapezoid
        %trap_current = repmat(reference_current(1,:), numel(t_trap), 1);
        %trap_current(2:3, 1:12) = trap_current(2:3, 1:12) + 500;
        %CSPF_curr_ref = timeseries(trap_current, t_trap);

    % step test: all CS/PF references (first slice) jump by +200 A at t0 + 0.1 s
    %t_step = t_start + 0.1;
    %step_current = repmat(reference_current(1,:), 4, 1);
    %step_current(3:4, 1:12) = step_current(3:4, 1:12) + 1000;
    %CSPF_curr_ref = timeseries(step_current, [t_start, t_step, t_step + 1e-6, max(t_max, t_step) + 1]);





    

    enable_KCURR= double(instance.get_setting('KCURR', pyargs('default', 0))); %Enable CSPF current controller KCURR (=1)
    IpControlMode_KCURR= double(instance.get_setting('control_Ip', pyargs('default', 0))); %Activation of Ip Control if == 1
    schedulingVar_KCURR=0; %Ip Vector of variables for gain scheduling current controller

    enable_RZIp=double(instance.get_setting('RZIp', pyargs('default', 0)));
    IpControlMode_RZIp=0; % Activation of Ip Control if == 1
    %schedulingVar_RZIp=0; %Ip Vector of parameters used for gain scheduling in RZIp control
    %schedulingVar_RZIp =[50, 1, 1.2]';
    
    enable_VS = double(instance.get_setting('VS', pyargs('default', 0))); % Enable VS control (=1 enabled)
    VSControlMode = double(instance.get_setting('VS_control_mode', pyargs('default', 2))); % VS controller Mode {'Off (0)'; 'VS1 (1)'; 'VS3 (2)'; 'VS1+VS3 (3)'}
    %schedulingVar_VS = 0.9; %Vector of variables for gain scheduling VS controller "li3"
    VSIpScaling = 1;% VS control action scaled with (nominal or estimated) Ip {'On with nominal Ip (1)'; 'On with estimated Ip (2)'}

    SchedulingVar_GapContr = 0; % Vector of variables for gain scheduling gap control
    ShapeIpControlMode = double(instance.get_setting('Shape_control', pyargs('default', 0))); %only gap control Shape-Ip control mode {'Off (0)'; 'PositionControl Rpl and Zpl (1)'; 'Isoflux ShapeControl (2)'; 'Isoflux Shape Control with Xpoint (3)'; 'Gap Control (4)'; 'Only CSPF Currents (5)'}
    SCIpScaling = 0; %Shape control action scaled with (nominal or estimated) Ip {'Off (0)'; 'On with nominal Ip (1)'; 'On with estimated Ip (2)'}
    IpControlMode = 0; % Activation of Ip Control (In case of CCS, 0 non active, 1 if integrated with KCURR, 2 if integrated with Shape control, In case of VCS 0 non active, 1 active).

    AllocatorMode = double(instance.get_setting('Allocator', pyargs('default', 0))); ;
    usePSModel = double(instance.get_setting('PSModel', pyargs('default', 0))); 
    % Exchange period with the plant (s), the sample time of the CODAC_input, MATLAB Function and Output blocks of the Simulink model:
    in_out_dt = double(instance.get_setting('in_out_msg_dt', pyargs('default', 0.001)));
    zdot_disturbance = double(instance.get_setting('zdot_disturbance', pyargs('default', 0)));
    zdot_filter = double(instance.get_setting('zdot_filter', pyargs('default', 0)));


    logger.info(sprintf('exchange period in_out_dt = %g s', in_out_dt));
    %%  Simulink controller

    obj = pcssp_KCURR_PFPO1_obj();
    obj.init;
    obj.setup;

    obj2 = RZIp_CCS_wrapper_obj();
    %obj2 = pcssp_RZIp_CCS_obj();
    obj2.init;
    obj2.setup;

    obj3 = pcssp_VS_obj();
    obj3.init;
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

    % Data dictionary of the top model: the KCURR one it used so far, the power-supply one
    % (bus objects pcssp_PS_* of the PowerSupplies ports, created by obj6.init) and the base
    % workspace (the variables of this script). Created next to the model the first time,
    % then only checked; the model is linked to it here, before the simulation.
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

    %% create Reference for Simulink
    in = Simulink.SimulationInput('pcssp_controller_and_ps');
    in = in.setModelParameter('StartTime',num2str(t_start),'StopTime',num2str(t_max));
    in = in.setModelParameter('SaveOutput', 'on');

    fprintf('Start closed loop\n')
    if DEBUG_LOG
        csv_file = fullfile(getenv('RUN_DIR'), 'simulink_log.csv');
        if exist(csv_file, 'file'), delete(csv_file); end
            sig = { 
                'KCURR_CSPF_ctrl_e',11;
                'KCURR_Ip_ctrl_e',1;  
                'KCURR_residual_ctrl_e',11; 
                'KCURR_actual_ref',11;
                    'RZIp_CSPF_curr_cmd',11;
                    'RZIp_VS1_cmd',1;
                'VS_VS1_cmd',1;   
                'VS_VS3_cmd',1;
                    'Shape_CSPF_curr_cmd',11;
                    'Shape_VS_decoupling',11;
                'Output_total',14;
                    'CSPF_curr_ref',11;
                    'Zdot',1;
                    'Ip_ref',1;
                    'Z_ref',1;
                    'R_ref',1;
                'ip',1;
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

    tic;

    % disturbance test: voltage kick added to one supply after the controller (14 = VS3,
    % 12 = VS1); not filtered by the PS model. voltage = 0 turns it off
    %kick = struct('supply', 14, 'voltage', 200, 't_on', t_start + 0.05, 'width', 0.01);


    logs = sim(in);
    close_system('pcssp_controller_and_ps',0);
    simulation_time=toc;
    plasma_duration=t_max-t_start;

    fprintf(['\nControl over! ' num2str(simulation_time) 's of simulation for ' num2str(plasma_duration) 's of plasma.\n']);
end

instance.reuse_instance();

function value = geo_ref_with_fallback(value, ts, field, logger)
    % Recompute boundary.geometric_axis.(field) from the boundary outline
    % midpoint when the primary value is empty or carries IMAS's
    % empty-float sentinel (abs > 1e30). If the outline is itself empty,
    % fall back further to global_quantities.magnetic_axis.(field). Warns
    % only the first time any fallback is used across the whole loop.
    persistent warned
    if isempty(warned)
        warned = false;
    end
    if isempty(value) || abs(value) > 1e30
        outline_field = ts.boundary.outline.(field);
        if ~isempty(outline_field)
            value = (max(outline_field) + min(outline_field)) / 2;
            fallback = 'outline midpoint';
        else
            value = ts.global_quantities.magnetic_axis.(field);
            fallback = 'magnetic_axis';
        end
        if ~warned
            logger.warning(sprintf( ...
                ['muscle_controller_NICE_IMAS_iter_init: boundary.geometric_axis.%s ' ...
                 'was empty/sentinel, falling back to %s'], field, fallback));
            warned = true;
        end
    end
end
