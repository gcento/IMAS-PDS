function [Icspf,Ivs3,rpl,zpl,zdot,ip,gaps,li3,schedulingVar_RZIp] = controller_input_preprocess(coil_currents,rpl,zpl,ip,gaps,li3,schedulingVar_RZIp,t)
% Turns the plant measurements into the controller inputs:
% - Icspf: the 11 CS/PF circuit currents, CS1 = mean of CS1U and CS1L;
% - Ivs3: the VS3 circuit current, (I_VS3U - I_VS3L)/2 (anti-series pair);
% - zdot: vertical velocity, filtered derivative of zpl, s/(tau*s + 1) discretised with
%   Tustin at the actual time step (tau = 7.5 ms). The first call returns 0; a call with
%   the same or an earlier time holds the previous value.

persistent last_t last_zpl last_zdot
tau = 0.0075;

coil_currents = coil_currents(:);
cs1mean = (coil_currents(3) + coil_currents(4))/2;
Icspf = [coil_currents(1:2); cs1mean; coil_currents(5:12)];
Ivs3 = 0.5*(coil_currents(13) - coil_currents(14));

if isempty(last_t)
    zdot = 0;
else
    dt = t - last_t;
    if dt > 0
        a = (2*tau - dt) / (2*tau + dt);
        b = 2 / (2*tau + dt);
        zdot = a*last_zdot + b*(zpl - last_zpl);
    else
        zdot = last_zdot;
    end
end

last_t = t;
last_zpl = zpl;
last_zdot = zdot;

if abs(ip) <= 100
    error("plasma current is almost 0");
end
