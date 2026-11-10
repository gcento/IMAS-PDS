
clear
clc

VS_control_mode = 1;      % 0 off, 1 VS1, 2 VS3, 3 VS1+VS3
% KCURR on/off, control_Ip and dt are at the top of build_cl.m

prm = fullfile(fileparts(mfilename('fullpath')), '..', 'params');
K = load(fullfile(prm, 'controller.mat'));
F = load(fullfile(prm, 'FC016_Config_VS.mat'), 'Config_VS');

%% gains
Gold = K.VS;    % old gains
% FC016 (continuous A/B/C/D, one set), the ones the CCS controller loads; Ip_ref stays the one of controller.mat
Gnew = Gold;
for f = fieldnames(F.Config_VS.Gains)'
    Gnew.(f{1}) = F.Config_VS.Gains.(f{1});
end

%% closed loops
[Acl_old, sys_old] = build_cl(Gold, VS_control_mode);
[Acl_new, sys_new] = build_cl(Gnew, VS_control_mode);

[~,~, Lold] = build_cl(Gold, VS_control_mode,'VS1');
[~, ~,Lnew] = build_cl(Gnew, VS_control_mode,'VS1');

%% test 1: poles before and after
eo = eig(Acl_old);  
en = eig(Acl_new);
figure; plot(real(eo), imag(eo), 'bx', real(en), imag(en), 'rx'); grid on
xlim([-200 20]); xline(0); legend('old', 'FC016'); xlabel('Re \lambda [1/s]'); ylabel('Im \lambda [rad/s]')
fprintf('max Re: old %.3f   FC016 %.3f   [1/s]\n', max(real(eo)), max(real(en)));