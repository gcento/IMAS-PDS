%% plot_run: closed-loop runs of the CCS controller against NICE, several runs overlaid (one colour per run).
% Figures: plasma (Ip, R, Z), vertical stabilisation, gaps (poloidal section of the first run and traces),
% isoflux (modes 2 and 3), coil currents with the KCURR reference, KCURR error, supply voltages.
% Measurements come from the plant (out_nice); references, commands and errors from the controller log
% (simulink_log.csv). Set run_dir and run plot_run.
%
% Signs: everything is shown in the IMAS frame of NICE (Ip < 0 for 105073).
% - The controller works in the frame with Ip > 0: its log has currents, Ip and voltage commands times ip_sign
%   and psi times the factor of isoflux_refs.mat. controller_frame and to_imas_frame (end of the file) undo
%   both. Runs older than the ip_sign frame logged NICE's current signs: their current factor is 1.
% - Voltages are the ones NICE applies, so that a positive voltage drives a positive current: with
%   nice_signs = 1, the opposite of the ones sent (out_controller, Output_total; controller_output_preprocess.m).

close all
clc
clear

run_dir = {'/home/ITER/centong/Desktop/PDSs/IMAS-PDS-torax/cases/runs/evolutive_controller_105073_20261124_084943'
           '/home/ITER/centong/Desktop/PDSs/IMAS-PDS-torax/cases/runs/evolutive_controller_105073_20261124_103707'};

t_section = [];                    % time of the poloidal section [s], [] = last slice of the first run
GAP_LIST  = [3 8 15 21 25 48];     % as in muscle_controller_NICE_IMAS_iter_init.m: gaps_1..6 of the csv

ts = '/equilibrium/time_slice[]&';
eq_signals = {
    'ip',    [ts 'global_quantities&ip']
    'r_cc',  [ts 'global_quantities&current_centre&r']
    'z_cc',  [ts 'global_quantities&current_centre&z']
    'psi_b', [ts 'boundary&psi']
    'btype', [ts 'boundary&type']                    % 0 limited, 1 with an X-point
};
pfa_signals = {
    'coil_current', '/pf_active/coil[]&current&data'
    'coil_names',   '/pf_active/coil[]&name'
};

% ===== Read every run, the controller log taken to the IMAS frame =====
run_dirs = cellstr(run_dir);
runs = cell(1, numel(run_dirs));
for i = 1:numel(run_dirs)
    r = struct();
    r.run_dir = strip(run_dirs{i}, 'right', '/');
    [runs_dir, run_name] = fileparts(r.run_dir);
    case_name = regexprep(run_name, '_\d{8}_\d{6}.*$', '');
    r.case_dir = fullfile(fileparts(runs_dir), case_name);
    r.label = strrep(erase(run_name, [case_name '_']), '_', '\_');
    r.frame = controller_frame(r.run_dir);
    r.csv = to_imas_frame(read_csv_signals(fullfile(r.run_dir, 'simulink_log.csv')), r.frame);
    plant_dir = first_folder(r.run_dir, {'out_nice', 'out_plant'});
    r.eq_file = fullfile(plant_dir, 'equilibrium.h5');
    r.eq  = read_h5_signals(r.eq_file, '/equilibrium/time', eq_signals);
    r.pfa = read_h5_signals(fullfile(plant_dir, 'pf_active.h5'), '/pf_active/time', pfa_signals);
    r.ctrl = read_h5_signals(fullfile(r.run_dir, 'out_controller', 'pf_active.h5'), '/pf_active/time', ...
        {'supply_voltage', '/pf_active/supply[]&voltage&data'});        % voltages sent by the controller
    r.ctrl.supply_voltage = r.frame.voltage * r.ctrl.supply_voltage;     % the ones NICE applies
    fprintf('%s: controller log to the IMAS frame: currents x%g, fluxes x%g; sent voltages x%g\n', ...
        run_name, r.frame.current, r.frame.flux, r.frame.voltage);
    if abs(r.csv.time(end) - r.pfa.time(end)) > 0.01
        warning('%s: simulink_log.csv ends at t = %.3f s but the run at %.3f s: it is from another run', ...
            run_name, r.csv.time(end), r.pfa.time(end));
    end
    runs{i} = r;
end
R = [runs{:}];
colors = lines(numel(R));
t_lim = [min(arrayfun(@(r) r.csv.time(1), R)), max(arrayfun(@(r) r.csv.time(end), R))];

% ===== Constants =====
% the 11 KCURR circuits (columns of KCURR_actual_ref and KCURR_CSPF_ctrl_e)
circ_names = {'CS3U','CS2U','CS1','CS2L','CS3L','PF1','PF2','PF3','PF4','PF5','PF6'};
% the 14 coils of pf_active and the KCURR circuit of each (0 = none)
coil_names = {'CS3U','CS2U','CS1U','CS1L','CS2L','CS3L','PF1','PF2','PF3','PF4','PF5','PF6','VS3U','VS3L'};
coil_circ  = [1 2 3 3 4 5 6 7 8 9 10 11 0 0];
% the 14 power supplies (pf_active.supply of out_controller, columns of Output_total) and their voltage limits
supply_names = {'CS3U','CS2U','CS1U','CS1L','CS2L','CS3L','PF1','PF2','PF3','PF4','PF5','VS1','PF6','VS3'};
vlim_supply  = [2100 2100 2100 2100 2100 2100 2100 3150 3150 3150 3150 6300 2100 2400];
coil = @(s, name) s.coil_current(:, strcmp(cellstr(string(s.coil_names)), name));

% ===== Plasma: Ip and current centre of the plant (solid), controller references (dashed) =====
figure('Name', 'Plasma');
tiledlayout(3, 1, 'TileSpacing', 'compact');
ax = [nexttile, nexttile, nexttile];
hold(ax, 'on');
h_run = gobjects(1, numel(R));
for i = 1:numel(R)
    e = R(i).eq;
    c = R(i).csv;
    h_run(i) = plot(ax(1), e.time, per_time(e.ip, e.time) / 1e6, 'Color', colors(i,:));
    plot(ax(1), c.time, c.Ip_ref / 1e6, '--', 'Color', colors(i,:));
    if isfield(e, 'z_cc')
        plot(ax(2), e.time, per_time(e.r_cc, e.time), 'Color', colors(i,:));
        plot(ax(3), e.time, 1e3 * per_time(e.z_cc, e.time), 'Color', colors(i,:));
    end
    plot(ax(2), c.time, c.R_ref, '--', 'Color', colors(i,:));
    plot(ax(3), c.time, 1e3 * c.Z_ref, '--', 'Color', colors(i,:));
end
xlim(ax, t_lim); grid(ax, 'on'); linkaxes(ax, 'x');
title(ax(1), 'I_p: plant solid, reference dashed'); ylabel(ax(1), 'I_p [MA]');
title(ax(2), 'R of the current centre');            ylabel(ax(2), 'R [m]');
title(ax(3), 'Z of the current centre');            ylabel(ax(3), 'Z [mm]');
xlabel(ax(3), 't [s]');
lg = legend(ax(1), h_run, {R.label});
lg.Layout.Tile = 'east';

% ===== Vertical stabilisation: Zdot seen by VS, VS3 current, VS3 and VS1 voltages =====
% I_VS3 = (I_VS3L - I_VS3U)/2 is the current of the VS3 supply in the IMAS frame, driven up by a positive VS3
% voltage. Voltages: applied by NICE (solid), VS module command (dashed, before the power supply model).
figure('Name', 'Vertical stabilisation');
tiledlayout(4, 1, 'TileSpacing', 'compact');
ax = [nexttile, nexttile, nexttile, nexttile];
hold(ax, 'on');
for i = 1:numel(R)
    c = R(i).csv;
    p = R(i).pfa;
    v = R(i).ctrl;
    h_run(i) = plot(ax(1), c.time, 1e3 * c.Zdot, 'Color', colors(i,:));
    plot(ax(2), p.time, (coil(p, 'VS3L') - coil(p, 'VS3U')) / 2e3, 'Color', colors(i,:));
    plot(ax(3), v.time, v.supply_voltage(:, strcmp(supply_names, 'VS3')) / 1e3, 'Color', colors(i,:));
    plot(ax(4), v.time, v.supply_voltage(:, strcmp(supply_names, 'VS1')) / 1e3, 'Color', colors(i,:));
    if isfield(c, 'VS_VS3_cmd')
        plot(ax(3), c.time, c.VS_VS3_cmd / 1e3, '--', 'Color', colors(i,:));
        plot(ax(4), c.time, c.VS_VS1_cmd / 1e3, '--', 'Color', colors(i,:));
    end
end
xlim(ax, t_lim); grid(ax, 'on'); linkaxes(ax, 'x');
title(ax(1), 'Zdot seen by VS');                                 ylabel(ax(1), 'dZ/dt [mm/s]');
title(ax(2), 'VS3 current (I_{VS3L} - I_{VS3U})/2');             ylabel(ax(2), 'I [kA]');
title(ax(3), 'VS3 voltage: applied solid, VS command dashed');   ylabel(ax(3), 'V [kV]');
title(ax(4), 'VS1 voltage: applied solid, VS command dashed');   ylabel(ax(4), 'V [kV]');
xlabel(ax(4), 't [s]');
lg = legend(ax(1), h_run, {R.label});
lg.Layout.Tile = 'east';

% ===== Gaps: poloidal section of the first run, and each controlled gap against its reference =====
% The section shows the plant at t_section: boundary.gap written by the diagnostics actor,
% from the wall point of each line to the separatrix, a tick where the separatrix should cross
% (the reference the controller tracked) and both values in the labels. At this scale a few mm
% are not visible: read the labels. The traces show what the controller saw (gaps) and what it
% tracked (gaps_ref), both from simulink_log.csv.
e1 = R(1).eq;
if ~has_dataset(R(1).eq_file, '/equilibrium', 'time_slice[]&boundary&gap[]&value')
    warning('%s: no boundary.gap in the plant equilibrium (run without the diagnostics actor): no gap figure', ...
        R(1).label);
else
    if isempty(t_section)
        t_section = e1.time(end);
    end
    [~, k] = min(abs(e1.time - t_section));
    t_k = e1.time(k);

    % plasma boundary and gap lines at slice k; the DD angle is clockwise, theta counter-clockwise from +R
    out_r = slice_of(R(1).eq_file, [ts 'boundary&outline&r'], k);
    out_z = slice_of(R(1).eq_file, [ts 'boundary&outline&z'], k);
    ok = isfinite(out_r) & isfinite(out_z);
    out_r = out_r(ok);
    out_z = out_z(ok);
    g_r   = slice_of(R(1).eq_file, [ts 'boundary&gap[]&r'], k);
    g_z   = slice_of(R(1).eq_file, [ts 'boundary&gap[]&z'], k);
    theta = -slice_of(R(1).eq_file, [ts 'boundary&gap[]&angle'], k);
    gap   = slice_of(R(1).eq_file, [ts 'boundary&gap[]&value'], k);   % NaN where not found
    sep_r = g_r + gap .* cos(theta);                                  % where each line meets the separatrix
    sep_z = g_z + gap .* sin(theta);

    % reference of the controlled gaps at t_k: what the controller tracked (NaN if not logged)
    g_ref = nan(size(gap));
    if isfield(R(1).csv, 'gaps_ref')
        [~, kc] = min(abs(R(1).csv.time - t_k));
        g_ref(GAP_LIST) = R(1).csv.gaps_ref(kc, :);
    end
    ref_r = g_r + g_ref .* cos(theta);                                % where the separatrix should be
    ref_z = g_z + g_ref .* sin(theta);

    % first wall: first wall + divertor outlines, which together close the vessel
    wall_file = fullfile(first_folder(fullfile(R(1).case_dir, 'preprocess'), {'in', 'initial_state', 'dina_in'}), 'wall.h5');
    wr = h5read(wall_file, '/wall/description_2d[]&limiter&unit[]&outline&r');
    wz = h5read(wall_file, '/wall/description_2d[]&limiter&unit[]&outline&z');
    valid = abs(wr) < 1e30;                                           % -9e40 = IMAS empty
    wall_r = [wr(valid(:, 1), 1); wr(valid(:, 2), 2); wr(1, 1)];
    wall_z = [wz(valid(:, 1), 1); wz(valid(:, 2), 2); wz(1, 1)];

    % separatrix with its legs (where the divertor gaps end): psi = psi_b on NICE's mesh, linear
    % on the triangles of the nodes as in the diagnostics actor, on a 2 cm grid inside the first
    % wall (the mesh is much larger than the vessel)
    m_r   = slice_of(R(1).eq_file, [ts 'ggd[]&r[]&values'], k);
    m_z   = slice_of(R(1).eq_file, [ts 'ggd[]&z[]&values'], k);
    m_psi = slice_of(R(1).eq_file, [ts 'ggd[]&psi[]&values'], k);
    ok = isfinite(m_r) & isfinite(m_z) & isfinite(m_psi);
    F = scatteredInterpolant(m_r(ok), m_z(ok), m_psi(ok), 'linear', 'none');
    [RR, ZZ] = meshgrid(min(wall_r):0.02:max(wall_r), min(wall_z):0.02:max(wall_z));
    PSI = F(RR, ZZ);
    PSI(~inpolygon(RR, ZZ, wall_r, wall_z)) = NaN;

    figure('Name', 'Gaps', 'Units', 'normalized', 'Position', [0.1 0.05 0.6 0.85]);
    n_gaps = numel(GAP_LIST);
    tiledlayout(n_gaps, 3, 'TileSpacing', 'compact');

    ax_s = nexttile(1, [n_gaps 2]);                                   % section: first two columns
    hold(ax_s, 'on'); box(ax_s, 'on');
    h = gobjects(1, 7);
    h(1) = plot(ax_s, wall_r, wall_z, 'Color', '#52514e', 'LineWidth', 1.5);
    h(2) = patch(ax_s, 'XData', out_r, 'YData', out_z, 'FaceColor', '#2a78d6', 'FaceAlpha', 0.12, ...
                 'EdgeColor', 'none');
    [~, h(3)] = contour(ax_s, RR, ZZ, PSI, [1 1] * e1.psi_b(k), 'LineColor', '#2a78d6', 'LineWidth', 1);
    others = setdiff(find(isfinite(gap)), GAP_LIST);
    if ~isempty(others)
        h_other = plot(ax_s, [g_r(others) sep_r(others)]', [g_z(others) sep_z(others)]', ...
                       'Color', '#c3c2b7', 'LineWidth', 0.75);
        h(6) = h_other(1);
    end
    for j = GAP_LIST
        h(4) = plot(ax_s, [g_r(j) sep_r(j)], [g_z(j) sep_z(j)], 'Color', '#eb6834', 'LineWidth', 2.5);
        plot(ax_s, g_r(j), g_z(j), 'o', 'MarkerSize', 5, 'MarkerEdgeColor', '#eb6834', 'MarkerFaceColor', '#eb6834');
        if isfinite(g_ref(j))                                         % tick across the line at the reference
            p = 0.1 * [-sin(theta(j)) cos(theta(j))];
            h(5) = plot(ax_s, ref_r(j) + [-1 1]*p(1), ref_z(j) + [-1 1]*p(2), 'Color', '#0b0b0b', 'LineWidth', 1.5);
        end
        % label outside the wall, behind the start point, on the side away from the gap line
        d = [cos(theta(j)) sin(theta(j))];
        ha = 'center';  va = 'middle';
        if d(1) > 0.5,  ha = 'right';  elseif d(1) < -0.5,  ha = 'left';    end
        if d(2) > 0.5,  va = 'top';    elseif d(2) < -0.5,  va = 'bottom';  end
        label = sprintf('GAP%02d\n%.1f cm (ref %.1f)', j, 100*gap(j), 100*g_ref(j));
        text(ax_s, g_r(j) - 0.15*d(1), g_z(j) - 0.15*d(2), label, 'Color', '#0b0b0b', 'FontSize', 9, ...
             'HorizontalAlignment', ha, 'VerticalAlignment', va, 'BackgroundColor', 'w', 'Margin', 1);
    end
    h(7) = plot(ax_s, e1.r_cc(k), e1.z_cc(k), '+', 'Color', '#0b0b0b', 'MarkerSize', 10, 'LineWidth', 1.5);
    axis(ax_s, 'equal');
    axis(ax_s, [2.2 10.5 -5.0 5.4]);                                 % room for the labels outside the wall
    xlabel(ax_s, 'R [m]'); ylabel(ax_s, 'Z [m]');
    title(ax_s, sprintf('%s, t = %.3f s: gap now, black tick at its reference', R(1).label, t_k));
    names = {'first wall', 'plasma', 'separatrix', 'controlled gaps', 'gap reference', 'other gaps', 'current centre'};
    shown = isgraphics(h);
    lg = legend(ax_s, h(shown), names(shown), 'NumColumns', 4);
    lg.Layout.Tile = 'south';

    ax_g = gobjects(1, n_gaps);                                       % traces: third column
    h_gap = gobjects(1, numel(R));
    for n = 1:n_gaps
        ax_g(n) = nexttile(3*n);
        hold(ax_g(n), 'on');
        for i = 1:numel(R)
            c = R(i).csv;
            if isfield(c, 'gaps')
                h_gap(i) = plot(ax_g(n), c.time, 100 * c.gaps(:, n), 'Color', colors(i,:));
            end
            if isfield(c, 'gaps_ref')
                plot(ax_g(n), c.time, 100 * c.gaps_ref(:, n), '--', 'Color', colors(i,:));
            end
        end
        xline(ax_g(n), t_k, ':', 'Color', '#898781');                % time of the section
        grid(ax_g(n), 'on');
        ylabel(ax_g(n), sprintf('GAP%02d [cm]', GAP_LIST(n)));
        if n == 1 && numel(R) > 1
            shown = isgraphics(h_gap);
            lg = legend(ax_g(1), h_gap(shown), {R(shown).label});
            lg.Layout.Tile = 'east';
        end
    end
    xlim(ax_g, t_lim); linkaxes(ax_g, 'x');
    title(ax_g(1), 'seen by the controller, reference dashed');
    xlabel(ax_g(end), 't [s]');
end

% ===== Isoflux: what shapeIp_CCS saw, flux at each control point minus the reference flux =====
% Mode 2 drives FluxContr - FluxLimiterPoint to dPsi_ref = 0, mode 3 FluxContr - FluxXp (and Rxp, Zxp
% to their references: last two panels, runs in mode 3 only; FluxXp is NaN while limited).
% From simulink_log.csv (flux_ctrl_1..8, flux_lim, flux_xp, in the IMAS frame); mode and points from the
% run's settings (magnetic_controller.shape_control, isoflux_r/z). Runs without them are left out.
% Dotted: the X-point appears.
for i = 1:numel(R)
    cfg = fullfile(R(i).run_dir, 'configuration.ymmsl');
    R(i).mode  = ymmsl_setting(cfg, 'magnetic_controller.shape_control');
    R(i).iso_r = ymmsl_setting(cfg, 'magnetic_controller.isoflux_r');
    R(i).iso_z = ymmsl_setting(cfg, 'magnetic_controller.isoflux_z');
    R(i).iso_e = [];
    c = R(i).csv;
    if isfield(c, 'flux_ctrl') && isfield(c, 'flux_lim') && isfield(c, 'flux_xp') && ~isempty(R(i).iso_r)
        n = min(numel(R(i).iso_r), size(c.flux_ctrl, 2));             % 6 points in mode 2, 8 in mode 3
        if isequal(R(i).mode, 3)
            R(i).iso_e = c.flux_ctrl(:, 1:n) - c.flux_xp;             % time x points
        else
            R(i).iso_e = c.flux_ctrl(:, 1:n) - c.flux_lim;
        end
    end
end

has_iso = arrayfun(@(r) ~isempty(r.iso_e), R);
if ~any(has_iso)
    fprintf('\nNo flux_ctrl in simulink_log.csv or no isoflux points in the settings: no isoflux figure\n');
else
    n_pts = max(arrayfun(@(r) size(r.iso_e, 2), R(has_iso)));
    in_xp = has_iso & arrayfun(@(r) isequal(r.mode, 3), R);           % runs in mode 3
    figure('Name', 'Isoflux', 'Units', 'normalized', 'Position', [0.1 0.1 0.7 0.75]);
    tl = tiledlayout('flow', 'TileSpacing', 'compact');
    title(tl, {'seen by shapeIp\_CCS: flux at each point minus the reference flux, driven to 0', ...
               'reference: limiter point in mode 2, X-point in mode 3'});
    h_iso = gobjects(1, numel(R));
    for p = 1:n_pts
        i1 = find(has_iso & arrayfun(@(r) size(r.iso_e, 2) >= p, R), 1);   % a run with this point
        ax = nexttile;
        hold(ax, 'on');  grid(ax, 'on');
        for i = find(has_iso)
            if p <= size(R(i).iso_e, 2)
                h_iso(i) = plot(ax, R(i).csv.time, R(i).iso_e(:, p), 'Color', colors(i,:));
            end
            if isfield(R(i).eq, 'btype') && R(i).eq.btype(1) == 0 && any(R(i).eq.btype == 1)
                xline(ax, R(i).eq.time(find(R(i).eq.btype == 1, 1)), ':', 'Color', colors(i,:));
            end
        end
        yline(ax, 0, 'k');
        xlim(ax, t_lim);
        title(ax, sprintf('P%d (%.2f, %.2f) m', p, R(i1).iso_r(p), R(i1).iso_z(p)));
        ylabel(ax, '\Delta\psi [Wb]');
        xlabel(ax, 't [s]');
    end
    if any(in_xp)
        for v = {'rxp', 'zxp'}                                        % X-point position, mode 3
            ax = nexttile;
            hold(ax, 'on');  grid(ax, 'on');
            for i = find(in_xp)
                plot(ax, R(i).csv.time, R(i).csv.(v{1}), 'Color', colors(i,:));
                plot(ax, R(i).csv.time, R(i).csv.([v{1} '_ref']), '--', 'Color', colors(i,:));
            end
            xlim(ax, t_lim);
            title(ax, sprintf('%s, reference dashed', v{1}));
            ylabel(ax, [v{1} ' [m]']);
            xlabel(ax, 't [s]');
        end
    end
    labels = arrayfun(@(r) sprintf('%s, mode %s', r.label, num2str(r.mode)), R(has_iso), 'UniformOutput', false);
    lg = legend(h_iso(has_iso), labels);
    lg.Layout.Tile = 'east';
end

% ===== Coil currents of the plant (solid) and the reference KCURR tracks (dashed) =====
% CS1U and CS1L share the CS1 reference, VS3U and VS3L have none; the reference moves with the Ip control
figure('Name', 'Coil currents');
tiledlayout(4, 4, 'TileSpacing', 'compact');
h_ref = gobjects(1, numel(R));
for k = 1:numel(coil_names)
    nexttile
    hold on
    for i = 1:numel(R)
        h_run(i) = plot(R(i).pfa.time, coil(R(i).pfa, coil_names{k}) / 1e3, 'Color', colors(i,:));
        if coil_circ(k) > 0
            h_ref(i) = plot(R(i).csv.time, R(i).csv.KCURR_actual_ref(:, coil_circ(k)) / 1e3, '--', ...
                'Color', colors(i,:));
        end
    end
    xlim(t_lim); grid on
    title(coil_names{k}); ylabel('I [kA]');
    if k == 1
        lg = legend([h_run, h_ref], [strcat({'plant '}, {R.label}), strcat({'KCURR reference '}, {R.label})]);
        lg.Layout.Tile = 'east';
    end
end
xlabel('t [s]');

% ===== KCURR control error of each circuit (reference - plant) =====
figure('Name', 'KCURR error');
tiledlayout(3, 4, 'TileSpacing', 'compact');
for k = 1:numel(circ_names)
    nexttile
    hold on
    for i = 1:numel(R)
        h_run(i) = plot(R(i).csv.time, R(i).csv.KCURR_CSPF_ctrl_e(:, k), 'Color', colors(i,:));
    end
    xlim(t_lim); grid on
    title(circ_names{k}); ylabel('e [A]');
    if k == 1
        lg = legend(h_run, {R.label});
        lg.Layout.Tile = 'east';
    end
end
xlabel('t [s]');

% ===== Power-supply voltages applied by NICE, with their limits =====
figure('Name', 'Supply voltages');
tiledlayout(4, 4, 'TileSpacing', 'compact');
for k = 1:numel(supply_names)
    nexttile
    hold on
    for i = 1:numel(R)
        h_run(i) = plot(R(i).ctrl.time, R(i).ctrl.supply_voltage(:, k), 'Color', colors(i,:));
    end
    h_lim = yline([-1 1] * vlim_supply(k), 'r--');
    xlim(t_lim); grid on
    title(supply_names{k}); ylabel('V [V]');
    if k == 1
        lg = legend([h_run, h_lim(1)], [{R.label}, {'limit'}]);
        lg.Layout.Tile = 'east';
    end
end
xlabel('t [s]');

% ===== Local functions (must stay at the end of the script) =====

function f = controller_frame(run_dir)
    % Factors taking the controller log of a run to the IMAS frame of NICE:
    % - current: currents, Ip and voltage commands, logged times ip_sign (isoflux_refs.mat); 1 in older runs,
    %   whose controller saw NICE's current signs;
    % - flux: fluxes, logged times ip_sign, or times flux_scale in older runs;
    % - voltage: the sent voltages -> the applied ones, -1 with nice_signs = 1 (its default).
    f = struct('current', 1, 'flux', 1, 'voltage', 1);
    file = fullfile(run_dir, 'isoflux_refs.mat');
    if isfile(file)
        s = load(file);
        if isfield(s, 'ip_sign')
            f.current = s.ip_sign;
            f.flux = s.ip_sign;
        elseif isfield(s, 'flux_scale')
            f.flux = 1 / s.flux_scale;
        end
    else
        warning('%s: no isoflux_refs.mat, the controller log is taken as already in the IMAS frame', run_dir);
    end
    nice_signs = ymmsl_setting(fullfile(run_dir, 'configuration.ymmsl'), 'magnetic_controller.nice_signs');
    if isempty(nice_signs) || nice_signs == 1
        f.voltage = -1;
    end
end

function c = to_imas_frame(c, f)
    % The signals of simulink_log.csv in the IMAS frame (controller_frame). Positions, gaps, Zdot and the PSC
    % flags are the same in both frames; Output_total holds the sent voltages.
    for n = {'KCURR_CSPF_ctrl_e', 'KCURR_Ip_ctrl_e', 'KCURR_residual_ctrl_e', 'KCURR_actual_ref', ...
             'RZIp_CSPF_curr_cmd', 'RZIp_VS1_cmd', 'VS_VS1_cmd', 'VS_VS3_cmd', 'Shape_CSPF_curr_cmd', ...
             'Shape_VS_decoupling', 'CSPF_curr_ref', 'Ip_ref', 'ip'}
        if isfield(c, n{1}), c.(n{1}) = f.current * c.(n{1}); end
    end
    for n = {'flux_ctrl', 'flux_lim', 'flux_xp'}
        if isfield(c, n{1}), c.(n{1}) = f.flux * c.(n{1}); end
    end
    if isfield(c, 'Output_total'), c.Output_total = f.voltage * c.Output_total; end
end

function data = read_h5_signals(file, time_path, signals)
    % Reads a list of signals from an IMAS HDF5 file into a struct.
    fprintf('\nReading %s\n', file);

    data = struct();
    data.time = double(h5read(file, time_path));
    data.time = data.time(:);
    nt = numel(data.time);
    fprintf('  %-15s %s\n', 'time', mat2str(size(data.time)));

    for k = 1:size(signals, 1)
        name = signals{k, 1};
        path = signals{k, 2};

        try
            x = h5read(file, path);
        catch err
            warning('Skipped "%s": cannot read %s\n  %s', name, path, err.message);
            continue
        end

        data.(name) = clean_imas(x, nt);
        fprintf('  %-15s %s\n', name, mat2str(size(data.(name))));
    end
end

function v = slice_of(file, path, k)
    % Time slice k of an IMAS dataset (time is its last dimension in MATLAB) as a column,
    % without reading the other slices.
    info = h5info(file, path);
    n = info.Dataspace.Size;       % in the order of h5read: [values, (1, 1,) time]
    v = double(h5read(file, path, [ones(1, numel(n) - 1) k], [n(1) ones(1, numel(n) - 1)]));
    v = v(:);
    v(v < -1e30) = NaN;            % IMAS "empty" value is -9e40
end

function tf = has_dataset(file, group, name)
    info = h5info(file, group);
    tf = any(strcmp({info.Datasets.Name}, name));
end

function data = read_csv_signals(file)
    % Reads the CSV into a struct. Columns like name_1, name_2, ... are grouped into one matrix.
    fid = fopen(file, 'r');
    headers = strsplit(strtrim(fgetl(fid)), ',');
    fclose(fid);

    M = readmatrix(file, 'NumHeaderLines', 1);
    base = regexprep(headers, '_\d+$', '');
    [names, ~, idx] = unique(base, 'stable');

    data = struct();
    data.time = M(:, 1);
    for k = 1:numel(names)
        field = matlab.lang.makeValidName(names{k});
        data.(field) = M(:, idx == k);
    end
end

function x = clean_imas(x, nt)
    % Numbers as double with the IMAS empty value as NaN, time on rows.
    if isnumeric(x)
        x = double(x);
        x(x < -1e30) = NaN;        % IMAS "empty" value is -9e40
    end
    if isvector(x)
        x = x(:);                  % always a column vector
    elseif size(x, 2) == nt && size(x, 1) ~= nt
        x = x.';                   % MATLAB reads HDF5 with dimensions reversed
    end
end

function v = ymmsl_setting(file, name)
    % Numbers of a setting in a run's configuration.ymmsl, [] if absent: a value, a list [a, b]
    % or a list written one "- a" per line.
    v = [];
    if ~isfile(file), return, end
    tok = regexp(fileread(file), ['\n\s*' regexptranslate('escape', name) ':([^\n]*)((?:\n\s+- [^\n]*)*)'], ...
                 'tokens', 'once');
    if ~isempty(tok)
        v = str2double(regexp([regexprep(tok{1}, '#.*', '') tok{2}], '[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?', 'match'));
    end
end

function d = first_folder(root, names)
    % the first of root/names{k} that exists (the first name when none does)
    d = fullfile(root, names{1});
    for k = 1:numel(names)
        if isfolder(fullfile(root, names{k}))
            d = fullfile(root, names{k});
            return
        end
    end
end

function v = per_time(v, t)
    % one value per time: a file read while the run writes it can have one more or one less
    % sample than time (cut, or padded with NaN)
    v = v(:);
    n = numel(t);
    if numel(v) >= n
        v = v(1:n);
    else
        v(end+1:n) = NaN;
    end
end
