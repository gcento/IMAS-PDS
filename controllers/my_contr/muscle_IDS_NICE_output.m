function [coils_current,rgeo,zgeo,zdot,ip,gaps,li3,schedulingVar_RZIp,t_out,stop_signal] = muscle_IDS_NICE_output()

global z_cur

% Receive the line average density package from MUSCLE3 via the instance
% py-class living in the base WS of matlab. Matlab automatically converts
% the data type into doubles for use in Simulink

% retrieve instance object from base ws
instance = evalin('base','instance');
logger = evalin('base','logger');

% Holds the last accepted measurement across calls, for the output_flag=-1
% case below (persists for the lifetime of the MATLAB session/model).
persistent last_coils_current last_rgeo last_zgeo last_ip last_t_out last_gaps last_zdot last_zdot_t last_li3 last_schedulingVar_RZIp
% Set once nice_evo_rd sends its final message (next_timestamp=None). The
% Stop Simulation block this feeds only halts at the end of the current
% step, so this same function still gets called again on a later step --
% without this guard, that call tries to receive on the now-closed port and
% crashes with "Port ... was closed while trying to receive on it" instead
% of the simulation stopping cleanly (confirmed live).
persistent finished
N_GAPS=62;
li3=0;
t_out=0;
stop_signal=0;
rgeo=0;
zgeo=0;
zdot=0;
ip=0;
gaps = zeros(N_GAPS,1);
coils_current=zeros(14,1);
schedulingVar_RZIp=zeros(3,1); %time,li,elongation

if isequal(finished, true)
    coils_current = last_coils_current;
    rgeo = last_rgeo;
    zgeo = last_zgeo;
    zdot = last_zdot;
    ip = last_ip;
    gaps = last_gaps;
    schedulingVar_RZIp = last_schedulingVar_RZIp;
    li3 = last_li3;
    t_out = last_t_out;
    stop_signal = 1;
    return;
end

% nice_evo_rd closes its output port directly once it reaches its own
% t_end, rather than sending one more message flagged next_timestamp=None
% -- so the receive() for that (never-coming) message throws
% "Port ... was closed while trying to receive on it" instead of the
% t_next-is-None branch below ever running. Treat that the same way: stop
% cleanly on the last accepted measurement (confirmed live).
try
    msg_eq = instance.receive("equilibrium_in_s");
    msg_pfa = instance.receive("pf_active_in_s");
catch ME
    if contains(ME.message, 'was closed while trying to receive')
        logger.warning('muscle_IDS_NICE_output: peer port closed without a next_timestamp=None sentinel, stopping cleanly');
        coils_current = last_coils_current;
        rgeo = last_rgeo;
        zgeo = last_zgeo;
        zdot = last_zdot;
        ip = last_ip;
        gaps = last_gaps;
        t_out = last_t_out;
        li3 = last_li3;
        schedulingVar_RZIp = last_schedulingVar_RZIp;
        stop_signal = 1;
        finished = true;
        return;
    end
    rethrow(ME);
end

t_out = msg_eq.timestamp;
t_next = msg_eq.next_timestamp;
last_t_out = t_out;

if isequal(t_next, py.None)
    stop_signal  = 1;
    finished = true;
end

equilibrium_serial=uint8(msg_eq.data);
equilibrium = imas_deserialize(equilibrium_serial,'equilibrium');

pf_active_serial=uint8(msg_pfa.data);
pf_active=imas_deserialize(pf_active_serial,'pf_active');

eq_flag = double(equilibrium.code.output_flag);
pfa_flag = double(pf_active.code.output_flag);
is_stale = (~isempty(eq_flag) && eq_flag(1) == -1) || (~isempty(pfa_flag) && pfa_flag(1) == -1);

if is_stale && ~isempty(last_ip)
    logger.warning(sprintf('muscle_IDS_NICE_output: received output_flag=-1 at t=%g, holding last valid measurement', t_out));
    coils_current = last_coils_current;
    rgeo = last_rgeo;
    zgeo = last_zgeo;
    zdot = last_zdot;
    ip = last_ip;
    gaps = last_gaps;
    li3 = last_li3;
    schedulingVar_RZIp = last_schedulingVar_RZIp;
else
    ip= equilibrium.time_slice{1}.global_quantities.ip; 
    rgeo=equilibrium.time_slice{1}.boundary.geometric_axis.r;
    zgeo=equilibrium.time_slice{1}.boundary.geometric_axis.z;
    if isempty(last_zgeo) || isempty(last_zdot_t)
        zdot = 0;                      
    else
        dt = t_out - last_zdot_t;
        if dt > 0
            a = (2*0.0075 - dt) / (2*0.0075 + dt);
            b =  2 / (2*0.0075 + dt);
            zdot = a*last_zdot + b*(zgeo - last_zgeo);
            %zdot=(zgeo - last_zgeo)/dt
        else
            zdot = last_zdot;     
        end
    end
    for ii=1:length(pf_active.coil)
       coil=struct(pf_active.coil{ii});
       coils_current(ii)= coil.current.data;
    end
    coils_current(3) = 0.5 * (coils_current(3) + coils_current(4));
    li3 = equilibrium.time_slice{1}.global_quantities.li_3;
    kappa = equilibrium.time_slice{1}.boundary.elongation;
    schedulingVar_RZIp = [t_out; li3; kappa];

    last_coils_current = coils_current;
    last_rgeo = rgeo;
    last_zgeo = zgeo;
    last_zdot = zdot;
    last_ip = ip;
    last_gaps = gaps;
    last_zdot_t = t_out;
    last_li3 = li3;
    last_schedulingVar_RZIp = schedulingVar_RZIp;
end

z_cur = zgeo;

