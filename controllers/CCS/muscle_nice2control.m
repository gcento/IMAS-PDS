function [coils_current,rpl,zpl,zdot,ip,r_xp,z_xp,psi_xp,PsiLimiter,flux_contr,gaps,li3,schedulingVar_RZIp,t_out,stop_signal] = muscle_nice2control()
% Receives the NICE equilibrium and pf_active from MUSCLE3 and returns the controller
% measurements. Called by the CODAC_input block at every step.

instance = evalin('base','instance');
logger = evalin('base','logger');

persistent last         % the last valid measurements, in the order of the outputs
persistent last_t_out finished
persistent F_flux       % psi interpolant on NICE's mesh, built once

%% initialization
coils_current = zeros(14,1);
rpl = 0;
zpl = 0;
zdot = 0;
ip = 0;
r_xp = NaN;             % NaN until the first X-point is found
z_xp = NaN;
psi_xp = NaN;
PsiLimiter = 0;
flux_contr = zeros(62,1);
gaps = zeros(62,1);
li3 = 0;
schedulingVar_RZIp = zeros(3,1);    % [t, li3, elongation]
t_out = 0;
stop_signal = 0;

% start from the last measurements: what is not found below keeps its last value
if ~isempty(last)
    [coils_current,rpl,zpl,zdot,ip,r_xp,z_xp,psi_xp,PsiLimiter,flux_contr,gaps,li3,schedulingVar_RZIp] = last{:};
    t_out = last_t_out;
end

%% receive
% after the last message: keep returning the last measurements with stop_signal = 1
if isequal(finished, true)
    stop_signal = 1;
    return
end

% nice_evo_rd closes its port at its t_end instead of sending a next_timestamp=None
% message: stop cleanly on the last measurements
try
    msg_eq = instance.receive("equilibrium_in_s");
    msg_pfa = instance.receive("pf_active_in_s");
catch ME
    if ~contains(ME.message, 'was closed while trying to receive')
        rethrow(ME);
    end
    logger.warning('muscle_nice2control: peer port closed, stopping cleanly');
    finished = true;
    stop_signal = 1;
    return
end

t_out = msg_eq.timestamp;
last_t_out = t_out;
if isequal(msg_eq.next_timestamp, py.None)    % last message
    finished = true;
    stop_signal = 1;
end

equilibrium = imas_deserialize(uint8(msg_eq.data), 'equilibrium');
pf_active = imas_deserialize(uint8(msg_pfa.data), 'pf_active');

% a rejected NICE step (output_flag = -1) repeats the previous state: keep the last measurements
eq_flag = double(equilibrium.code.output_flag);
pfa_flag = double(pf_active.code.output_flag);
rejected = (~isempty(eq_flag) && eq_flag(1) == -1) || (~isempty(pfa_flag) && pfa_flag(1) == -1);
if rejected && ~isempty(last)
    logger.warning(sprintf('muscle_nice2control: output_flag=-1 at t=%g, holding the last measurements', t_out));
    return
end

%% measurements
ts = equilibrium.time_slice{1};
gq = ts.global_quantities;
cc = gq.current_centre;
bd = ts.boundary;

% unset IDS floats arrive as -9e40, unset arrays as []
vals = [gq.ip, cc.r, cc.z, cc.velocity_z, gq.li_3, bd.elongation];
if numel(vals) ~= 6 || any(abs(vals) > 1e30)
    error('muscle_nice2control: ip, current centre (with velocity_z), li_3 or elongation missing at t=%g', t_out);
end
ip = gq.ip;
rpl = cc.r;
zpl = cc.z;
zdot = cc.velocity_z;
li3 = gq.li_3;
schedulingVar_RZIp = [t_out; li3; bd.elongation];

for k = 1:numel(pf_active.coil)
    coil = struct(pf_active.coil{k});
    coils_current(k) = coil.current.data;
end

% gaps, written by the diagnostics actor in every slice, NaN where not found
GAP_LIST = evalin('base', 'GAP_LIST');
for k = GAP_LIST
    gaps(k) = bd.gap{k}.value;
end

% fluxes in the Ip > 0 frame shapeIp_CCS is designed for: psi * sign(Ip), sign of the reference set by the init
ip_sign = evalin('base', 'ip_sign');

% X-point: the saddle node of contour_tree; none in a limited plasma (last value kept). psi_xp is the mode-3
% reference flux: at the saddle, or with xp_flux_at_ref = 1 at the reference X-point (below, on the mesh),
% which then also stands in for the measured position when NICE finds no saddle
mode3 = evalin('base', 'ShapeIpControlMode') == 3;
xp_at_ref = mode3 && evalin('base', 'xp_flux_at_ref');
xp_found = false;
for k = 1:numel(ts.contour_tree.node)
    node = ts.contour_tree.node{k};
    if node.critical_type == 1
        r_xp = node.r;
        z_xp = node.z;
        psi_xp = ip_sign * node.psi;
        xp_found = true;
        break
    end
end
if isnan(r_xp) && mode3 && ~xp_at_ref
    error('muscle_nice2control: mode 3 needs a diverted plasma, no X-point at t=%g', t_out);
end

% isoflux: psi at the control points, interpolated on NICE's mesh, and the boundary flux (mode 2: psi_boundary, or psi at
% the limiter point with the setting limiter; else the fixed limiter point).
% The mesh is the same at every step, so the interpolant is built once.
if numel(ts.ggd) >= 1 && numel(ts.ggd{1}.psi) >= 1 && ~isempty(ts.ggd{1}.psi{1}.values)
    g = ts.ggd{1};
    psi_nodes = g.psi{1}.values(:);
    if isempty(F_flux) || numel(F_flux.Values) ~= numel(psi_nodes)
        F_flux = scatteredInterpolant(g.r{1}.values(:), g.z{1}.values(:), psi_nodes);
    else
        F_flux.Values = psi_nodes;
    end
    R_pts = evalin('base', 'RfluxCTRpts');
    Z_pts = evalin('base', 'ZfluxCTRpts');
    flux_contr(1:numel(R_pts)) = ip_sign * F_flux(R_pts(:), Z_pts(:));
    if xp_at_ref
        % reference X-point at this time, interpolated as the From Workspace blocks Rxp_ref/Zxp_ref do
        % (linear, extrapolated past the last slice; constant with a single slice)
        [t_u, iu] = unique(evalin('base', 't_ref'));
        rz_ref = [reshape(evalin('base', 'rxp_ref'), [], 1), reshape(evalin('base', 'zxp_ref'), [], 1)];
        rz_ref = rz_ref(iu, :);
        if numel(t_u) > 1
            rz = interp1(t_u, rz_ref, t_out - evalin('base', 't_start'), 'linear', 'extrap');
        else
            rz = rz_ref(1, :);
        end
        psi_xp = ip_sign * F_flux(rz(1), rz(2));
        if ~xp_found        % no saddle in this step: the X-point position rows see the reference (zero error)
            r_xp = rz(1);
            z_xp = rz(2);
        end
    end
    if evalin('base', 'ShapeIpControlMode') == 2 && ~evalin('base', 'lim_point')
        PsiLimiter = ip_sign * ts.global_quantities.psi_boundary;
    else
        PsiLimiter = ip_sign * F_flux(evalin('base', 'R_lim'), evalin('base', 'Z_lim'));
    end
end
if xp_at_ref && isnan(psi_xp)
    error('muscle_nice2control: xp_flux_at_ref needs the flux map of NICE (ggd), missing at t=%g', t_out);
end

f = fopen(fullfile(getenv('RUN_DIR'), 'isoflux_log.csv'), 'a');
if f >= 0
    fprintf(f, '%.6f', t_out);  fprintf(f, ',%.9g', flux_contr(1:8), psi_xp, r_xp, z_xp, zpl, ip, PsiLimiter);  fprintf(f, '\n');  fclose(f);
end
last = {coils_current,rpl,zpl,zdot,ip,r_xp,z_xp,psi_xp,PsiLimiter,flux_contr,gaps,li3,schedulingVar_RZIp};
end
