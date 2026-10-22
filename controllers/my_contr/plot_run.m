close all
clc
clear

run_dir  = '/home/ITER/centong/Desktop/PDSs/IMAS-PDS/cases/runs/linear_plant_controller_105073_20261018_131731_Vert_VS1';
csv_file = fullfile(run_dir, 'simulink_log.csv');

[runs_dir, run_name] = fileparts(strip(run_dir, 'right', '/'));
case_dir = fullfile(fileparts(runs_dir), regexprep(run_name, '_\d{8}_\d{6}.*$', ''));
plant_dir = fullfile(run_dir, 'out_plant');               % plant outputs

eq_file   = fullfile(plant_dir, 'equilibrium.h5');
pfa_file  = fullfile(plant_dir, 'pf_active.h5');      % plant: measured coil currents
ppas_file = fullfile(plant_dir, 'pf_passive.h5');     % plant: vessel loop currents
ctrl_file = fullfile(run_dir, 'out_controller', 'pf_active.h5');   % voltages sent by the controller

ref_file  = fullfile(case_dir, 'preprocess', 'initial_state', 'pf_active.h5');

ts = '/equilibrium/time_slice[]&';
eq_signals = {
    'ip', [ts 'global_quantities&ip']
    'zgeo', [ts 'boundary&geometric_axis&z']
    'rgeo', [ts 'boundary&geometric_axis&r']
    'rmag', [ts 'global_quantities&magnetic_axis&r']
    'zmag', [ts 'global_quantities&magnetic_axis&z']
    %'li_3',[ts 'global_quantities&li_3']
    %'psi_boundary',[ts 'boundary&psi']
};

pfa_signals = {
    'coil_current', '/pf_active/coil[]&current&data'
    'coil_names', '/pf_active/coil[]&name'
};

ref_signals = {
    'coil_current', '/pf_active/coil[]&current&data'
    'coil_voltage', '/pf_active/coil[]&voltage&data'
    'coil_names', '/pf_active/coil[]&name'
};

ctrl_signals = {
    'coil_voltage', '/pf_active/coil[]&voltage&data'
};

ppas_signals = {
    'loop_current', '/pf_passive/loop[]&current'
};

csv_ctrl = read_csv_signals(csv_file);
eq   = read_h5_signals(eq_file,   '/equilibrium/time', eq_signals);
pfa  = read_h5_signals(pfa_file,  '/pf_active/time',   pfa_signals);
ref  = read_h5_signals(ref_file,  '/pf_active/time',   ref_signals);
ctrl = read_h5_signals(ctrl_file, '/pf_active/time',   ctrl_signals);
ppas = read_h5_signals(ppas_file, '/pf_passive/time',  ppas_signals);

if abs(csv_ctrl.time(end) - pfa.time(end)) > 0.01
    warning('simulink_log.csv ends at t = %.3f s but the run at %.3f s: it is from another run', ...
        csv_ctrl.time(end), pfa.time(end));
end

% ===== Constants =====
circ_names = {'CS3U','CS2U','CS1','CS2L','CS3L','PF1','PF2','PF3','PF4','PF5','PF6'};
idx11 = [1 2 3 5 6 7 8 9 10 11 12];   
% the 14 voltages sent by the controller are power-supply voltages, coil{k}.voltage = supply k
% (controller_output_preprocess); limits: voltage_limit_max of pf_active.supply in the F_INIT entry
supply_names = {'CS3U','CS2U','CS1U','CS1L','CS2L','CS3L','PF1','PF2','PF3','PF4','PF5','VS1','PF6','VS3'};
vlim_supply = [2100 2100 2100 2100 2100 2100 2100 3150 3150 3150 3150 6300 2100 2400];
tau = 0.1;

t_csv = csv_ctrl.time;
[~, dIref_dt] = gradient(csv_ctrl.CSPF_curr_ref, 1, t_csv);   % derivative along time
slices = ref.time(ref.time >= t_csv(1) & ref.time <= t_csv(end));

figure('Name','Coil error')
tiledlayout(2, 1, 'TileSpacing', 'compact');
nexttile
hh = plot(csv_ctrl.t,csv_ctrl.KCURR_CSPF_ctrl_e);
title('KCURR control error')
ylabel('I [A]')
nexttile
% normalised to each circuit's max |I_ref| over the run: dividing by I_ref itself blows up
% where a reference crosses zero (CS2L at 70.49 s in linear_plant_controller_105073)
I_norm = max(abs(csv_ctrl.CSPF_curr_ref), [], 1);
h = plot(csv_ctrl.t, 100*csv_ctrl.KCURR_CSPF_ctrl_e ./ I_norm);
title('KCURR control error [% of max |I_{ref}|]')
%h =plot(csv_ctrl.t,(csv_ctrl.KCURR_CSPF_ctrl_e));
legend(circ_names)
xlabel('t [s]')
%ylabel('[A]')
set(h([1 2 3 4 5]), 'LineWidth', 1.7)
set(hh([1 2 3 4 5]), 'LineWidth', 1.7)
% ===== 1) Coil currents v*s reference =====
figure('Name', 'Coil currents vs reference');
tiledlayout(3, 4, 'TileSpacing', 'compact');
for k = 1:11
    nexttile
    plot(pfa.time, pfa.coil_current(:, idx11(k)) / 1e3, 'b'); hold on
    plot(t_csv, csv_ctrl.CSPF_curr_ref(:, k) / 1e3, 'r--');
    plot(ref.time, ref.coil_current(:, idx11(k)) / 1e3, 'ko', 'MarkerSize', 4);
    xlim([t_csv(1) t_csv(end)]); grid on
    title(circ_names{k}); xlabel('t [s]'); ylabel('I [kA]');
end
legend('measured (plant)', 'reference seen by KCURR', 'F\_INIT slices', 'Location', 'best');
%VS3 coils (not in KCURR): VS3U, -VS3L and the circuit current (I_U - I_L)/2 seen by VS
nexttile
Iu = pfa.coil_current(:, 13); Il = pfa.coil_current(:, 14);
plot(pfa.time, Iu / 1e3, 'b', pfa.time, -Il / 1e3, 'm', pfa.time, (Iu - Il) / 2e3, '--'); hold on
plot(ref.time, ref.coil_current(:, 13) / 1e3, 'ko', 'MarkerSize', 4);
xlim([t_csv(1) t_csv(end)]); grid on
title('VS3'); xlabel('t [s]'); ylabel('I [kA]');
legend('VS3U', '-VS3L', '(I_U - I_L)/2', 'F\_INIT VS3U', 'Location', 'best');

% % ===== 2) KCURR current error vs the P-only prediction =====
% figure('Name', 'KCURR current error');
% tiledlayout(3, 4, 'TileSpacing', 'compact');
% for k = 1:11
%     nexttile
%     plot(t_csv, csv_ctrl.KCURR_CSPF_ctrl_e(:, k), 'b'); hold on
%     plot(t_csv, tau * dIref_dt(:, k), 'r--');
%     xline(slices, ':', 'Color', [0.6 0.6 0.6]);
%     grid on
%     title(circ_names{k}); xlabel('t [s]'); ylabel('I_{ref} - I [A]');
% end
% legend('error (KCURR)', '\tau dI_{ref}/dt (P-only prediction)', 'F\_INIT slices', 'Location', 'best');
% nexttile
% plot(t_csv, sqrt(mean(csv_ctrl.KCURR_CSPF_ctrl_e .^ 2, 2))); grid on
% title('RMS over the 11 circuits'); xlabel('t [s]'); ylabel('[A]');
% 
% ===== 3) Power-supply voltages sent by the controller vs limits and F_INIT reference =====
% coil{k}.voltage of out_controller is the voltage of power supply k. The F_INIT entry has
% coil voltages only: the reference of a supply is the coil of the same name, V(VS3L) -
% V(VS3U) for VS3 (anti-series pair) and none for VS1 (in series with PF2..PF5). They are
% instantaneous values at the slices and the PF2..PF5 coil voltages include the VS1 action
% (about +-1 kV) and the plasma: compare orders of magnitude and signs, not values.
ref_coil_names = cellstr(string(ref.coil_names));
ref_supply_voltage = nan(numel(ref.time), numel(supply_names));
for k = 1:numel(supply_names)
    switch supply_names{k}
        case 'VS1'
            % no coil of its own
        case 'VS3'
            ref_supply_voltage(:, k) = ref.coil_voltage(:, strcmp(ref_coil_names, 'VS3L')) ...
                - ref.coil_voltage(:, strcmp(ref_coil_names, 'VS3U'));
        otherwise
            ref_supply_voltage(:, k) = ref.coil_voltage(:, strcmp(ref_coil_names, supply_names{k}));
    end
end
figure('Name', 'Power-supply voltages');
tiledlayout(4, 4, 'TileSpacing', 'compact');
for k = 1:numel(supply_names)
    nexttile
    plot(ctrl.time, ctrl.coil_voltage(:, k), 'b'); hold on
    plot(ref.time, ref_supply_voltage(:, k), 'ko', 'MarkerSize', 4);
    yline([-1 1] * vlim_supply(k), 'r--');
    xlim([ctrl.time(1) ctrl.time(end)]); grid on
    title(supply_names{k}); xlabel('t [s]'); ylabel('V [V]');
end
legend('sent by controller', 'F\_INIT reference (from coil voltages)', 'limit', 'Location', 'best');
% 
% % ===== 4) Series (CS1) and anti-series (VS3) pairs in the plant =====
% % In ITER CS1U+CS1L are one series circuit and VS3U/VS3L one anti-series circuit; the
% % plasmaless model drives each coil with its own voltage, so these pairs are not forced.
% Iu = pfa.coil_current(:, 13); Il = pfa.coil_current(:, 14);
% figure('Name', 'CS1 and VS3 coil pairs');
% tiledlayout(2, 2, 'TileSpacing', 'compact');
% nexttile
% plot(pfa.time, pfa.coil_current(:, 3) / 1e3, pfa.time, pfa.coil_current(:, 4) / 1e3, '--'); grid on
% legend('CS1U', 'CS1L'); title('CS1 coils (series in ITER)'); xlabel('t [s]'); ylabel('I [kA]');
% nexttile
% plot(pfa.time, pfa.coil_current(:, 3) - pfa.coil_current(:, 4)); grid on
% title('I_{CS1U} - I_{CS1L} (0 in a series circuit)'); xlabel('t [s]'); ylabel('[A]');
% nexttile
% plot(pfa.time, Iu, pfa.time, -Il, '--'); grid on
% legend('VS3U', '-VS3L'); title('VS3 coils (anti-series in ITER)'); xlabel('t [s]'); ylabel('I [A]');
% nexttile
% plot(pfa.time, (Iu - Il) / 2, pfa.time, (Iu + Il) / 2); grid on
% legend('(I_U - I_L)/2 = I_{VS3} seen by VS', '(I_U + I_L)/2 common mode (0 in ITER)');
% xlabel('t [s]'); ylabel('[A]');
% 
% % ===== 5) VS controller, plasma pass-through, vessel =====
% figure('Name', 'VS, plasma and vessel');
% tiledlayout(2, 2, 'TileSpacing', 'compact');
% nexttile
% plot(t_csv, csv_ctrl.VS_VS1_cmd, t_csv, csv_ctrl.VS_VS3_cmd); grid on
% legend('VS1 cmd', 'VS3 cmd'); title('VS controller output (0 if enable\_VS = 0)');
% xlabel('t [s]'); ylabel('[V]');
% nexttile
% plot(t_csv, csv_ctrl.Zdot); grid on
% title('dZ/dt fed to VS (from the pass-through reference)'); xlabel('t [s]'); ylabel('[m/s]');
% nexttile
% plot(eq.time, eq.ip / 1e6, t_csv, csv_ctrl.Ip_ref / 1e6, '--'); grid on
% legend('I_p plant (pass-through)', 'I_p ref'); xlabel('t [s]'); ylabel('[MA]');
% nexttile
% plot(ppas.time, sum(ppas.loop_current, 2) / 1e3); grid on
% title('Sum of vessel loop currents (pf\_passive)'); xlabel('t [s]'); ylabel('[kA]');


% ===== Local functions (must stay at the end of the script) =====

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
    % Converts to double, replaces IMAS empty values with NaN, puts time on rows.
    %x = double(x);
    %x(x < -1e40) = NaN;            % IMAS "empty" value is -9e40

    if isvector(x)
        x = x(:);                  % always a column vector
    elseif size(x, 2) == nt && size(x, 1) ~= nt
        x = x.';                   % MATLAB reads HDF5 with dimensions reversed
    end
end