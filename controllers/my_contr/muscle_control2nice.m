function muscle_control2nice(t_cur,voltage)
% This function receives the current simulation time and
% a voltage from simulink.

% retrieve instance object from base WS
instance = evalin('base','instance');
logger = evalin('base','logger');

if numel(voltage) ~= 14
    error('muscle_control2nice:unexpectedVoltageWidth', ...
        'expected the 14-wide voltage command, got %s', mat2str(size(voltage)));
end

%kick = evalin('base','kick');
%if t_cur >= kick.t_on && t_cur < kick.t_on + kick.width
 %   voltage(kick.supply) = voltage(kick.supply) + kick.voltage;
%end

pfa = ids_init('pf_active');
pfa.ids_properties.homogeneous_time = 1;
pfa.time = t_cur;
pfa.coil=ids_allocate('pf_active', 'coil', 14);
for i = 1:14
    pfa.coil{i}.voltage.data=voltage(i);
end

pfa_serialized=imas_serialize(pfa, 'pf_active');
s = sprintf('%d ' ,uint8(pfa_serialized));
x = py.numpy.fromstring(s, py.numpy.int8, int8(-1), char(' '));
pfa_serialized=py.bytes(x);
msg =  py.libmuscle.Message(t_cur,py.None,pfa_serialized);
instance.send("pf_active_out_i", msg);
