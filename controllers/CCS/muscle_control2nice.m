function muscle_control2nice(t_cur, voltage)
% Sends the 14 supply voltages to NICE through MUSCLE3, as a pf_active IDS.
% Called by the CODAC output block at every step.

instance = evalin('base','instance');

if numel(voltage) ~= 14
    error('muscle_control2nice: expected 14 voltages, got %s', mat2str(size(voltage)));
end

% supplies and coils in NICE's order
supply_names = {'CS3U','CS2U','CS1U','CS1L','CS2L','CS3L','PF1','PF2','PF3','PF4','PF5','VS1','PF6','VS3'};
coil_names = {'CS3U','CS2U','CS1U','CS1L','CS2L','CS3L','PF1','PF2','PF3','PF4','PF5','PF6','VS3U','VS3L'};
pfa = ids_init('pf_active');
pfa.ids_properties.homogeneous_time = 1;
pfa.time = t_cur;
pfa.supply = ids_allocate('pf_active', 'supply', 14);
pfa.coil = ids_allocate('pf_active', 'coil', 14);
for i = 1:14
    pfa.supply{i}.name = supply_names{i};
    pfa.supply{i}.voltage.data = voltage(i);
    pfa.coil{i}.name = coil_names{i};
end

% MATLAB uint8 -> Python bytes, through a text of numbers
s = sprintf('%d ', uint8(imas_serialize(pfa, 'pf_active')));
data = py.bytes(py.numpy.fromstring(s, py.numpy.int8, int8(-1), char(' ')));
instance.send("pf_active_out_i", py.libmuscle.Message(t_cur, py.None, data));
end
