function total_cmd = controller_output_preprocess(cspfout, vs1out, vs3out)
% Controller outputs -> the 14 supply voltages in NICE's order:
% [CS3U CS2U CS1U CS1L CS2L CS3L PF1 PF2 PF3 PF4 PF5 VS1 PF6 VS3].
% The CS1 circuit (two coils in series) gets half of its voltage on CS1U and half on CS1L.

cspf = [cspfout(1:2); cspfout(3)/2; cspfout(3)/2; cspfout(4:11)];   % 11 circuits -> 12 coils
total_cmd = [cspf(1:11); vs1out; cspf(12); vs3out];
if evalin('base','nice_signs')
    total_cmd = -total_cmd;     % NICE evo_rd's voltage sign convention
end
end
