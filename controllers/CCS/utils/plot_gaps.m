% Poloidal section of ITER with the controlled gaps: those used by the gap controller (Proj) and
% those of the IDM document. Gap lines, gap values and separatrix come from the CREATE model in
% PCSSP (MScen01, 15 MA flattop, t = 80 s); the first wall from the wall IDS of the case.
% Each gap is drawn from its point on the wall to the separatrix. No IMAS needed: run  plot_gaps

%% parameters
gaps_controller = [3 8 15 21 25 48];       % nonzero columns of Proj in Gains_GapControl
gaps_document   = [21 25 15 12 10 3];      % cg1 ... cg6 of the IDM document, in this order
show_all        = true;                    % the other CREATE gap lines too, in light grey
case_name       = 'evolutive_controller_105073';   % only for the first wall (preprocess/in/wall.h5)

%% files
repo = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
S = load(fullfile(repo, 'local_installs', 'pcs', 'modules', 'PlasmaTokamakLinearModel', 'modelFiles', ...
                  'Config_EQ_MScen01_t80_Ip15MA_Div_Flattop.mat'), 'Config');
c = S.Config;
names = strtrim(string(c.Ytype));
y = c.Yequil(:);                           % model outputs at the equilibrium
wall_file = fullfile(repo, 'cases', case_name, 'preprocess', 'in', 'wall.h5');

%% gap lines: from (rgap, zgap) along thetagap (counter-clockwise from +R), as long as the gap
gap = y(startsWith(names, "GAP"));         % GAP01 ... GAP62 [m]
th = c.thetagap(:);
r1 = c.rgap(:);              z1 = c.zgap(:);
r2 = r1 + gap.*cos(th);      z2 = z1 + gap.*sin(th);              % on the separatrix

%% first wall: first wall + divertor outlines, which together close the vessel
wr = h5read(wall_file, '/wall/description_2d[]&limiter&unit[]&outline&r');
wz = h5read(wall_file, '/wall/description_2d[]&limiter&unit[]&outline&z');
valid = abs(wr) < 1e30;                    % -9e40 = IMAS empty
wall_r = [wr(valid(:, 1), 1); wr(valid(:, 2), 2); wr(1, 1)];
wall_z = [wz(valid(:, 1), 1); wz(valid(:, 2), 2); wz(1, 1)];

%% separatrix: psi on the CREATE grid (last block of outputs), contour at the X-point flux
R = c.RgridFWHLC(:);  Z = c.ZgridFWHLC(:);
[Rg, ~, iR] = unique(round(R, 6));  [Zg, ~, iZ] = unique(round(Z, 6));
psi = accumarray([iZ iR], y(end-numel(R)+1:end), [numel(Zg) numel(Rg)]);
psi_x = interp2(Rg, Zg, psi, y(names == "CV-RX"), y(names == "CV-ZX"));
[RR, ZZ] = meshgrid(Rg, Zg);
psi(~inpolygon(RR, ZZ, wall_r, wall_z)) = NaN;      % only inside the vessel

%% plot
figure('Color', '#fcfcfb');
axes('Color', '#fcfcfb', 'XColor', '#52514e', 'YColor', '#52514e'); hold on; box on
h_wall = plot(wall_r, wall_z, 'Color', '#52514e', 'LineWidth', 1.5);
[~, h_sep] = contour(Rg, Zg, psi, [psi_x psi_x], 'LineColor', '#898781', 'LineWidth', 1);

if show_all
    k = setdiff(1:numel(gap), [gaps_controller gaps_document]);
    plot([r1(k) r2(k)]', [z1(k) z2(k)]', 'Color', '#c3c2b7', 'LineWidth', 0.75);
end

both = intersect(gaps_controller, gaps_document);
groups = {both,                         'used by both',           '#2a78d6', '-';
          setdiff(gaps_controller, both), 'controller only (Proj)', '#eb6834', '-';
          setdiff(gaps_document, both),   'IDM document only',      '#1baf7a', '--'};
h = gobjects(1, size(groups, 1));
for i = 1:size(groups, 1)
    for k = groups{i, 1}
        h(i) = plot([r1(k) r2(k)], [z1(k) z2(k)], groups{i, 4}, 'Color', groups{i, 3}, 'LineWidth', 2.5);
        plot(r1(k), z1(k), 'o', 'MarkerSize', 6, 'MarkerEdgeColor', groups{i, 3}, 'MarkerFaceColor', groups{i, 3});
        label = sprintf('GAP%02d', k);
        j = find(gaps_document == k);
        if ~isempty(j)
            label = sprintf('%s (cg%d)', label, j);
        end
        % label outside the wall, behind the start point, on the side away from the gap line
        d = [cos(th(k)) sin(th(k))];
        ha = 'center';  va = 'middle';
        if d(1) > 0.5,  ha = 'right';  elseif d(1) < -0.5,  ha = 'left';    end
        if d(2) > 0.5,  va = 'top';    elseif d(2) < -0.5,  va = 'bottom';  end
        text(r1(k) - 0.15*d(1), z1(k) - 0.15*d(2), label, 'Color', '#0b0b0b', 'FontSize', 9, ...
             'HorizontalAlignment', ha, 'VerticalAlignment', va, 'BackgroundColor', '#fcfcfb', 'Margin', 1);
    end
end
ok = isgraphics(h);
legend([h(ok), h_wall, h_sep], [groups(ok, 2)', {'first wall', 'separatrix (CREATE)'}], ...
       'Location', 'eastoutside', 'TextColor', '#0b0b0b');
axis equal
axis([2.9 9.6 -5.0 5.4])
xlabel('R [m]'); ylabel('Z [m]');
title('Controlled gaps: controller (Proj) vs IDM document', 'Color', '#0b0b0b');
