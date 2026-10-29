
close all
clear
clc

model           = 'vertical_circuit_105073_t84.3176';
use_KCURR       = true;
VS_control_mode = 3;       
use_allocator   = true;       
use_PS          = true;       % power-supply model, as magnetic_controller.PSModel 1
gain_set        = 'noIp';     % KCURR gains: 'noIp' (control_Ip 0) or 'Ip' (control_Ip 1)
Ip_ref          = -15.063e6;  % [A] Ipl_ref of the VS laws (F_INIT of 105073, IMAS sign)
tau             = 7.5e-3;     % [s] Zdot filter 1/(tau s + 1): 7.5e-3 as zdot_filter 1,
                              % about 0.5e-3 for zdot_filter 0 (backward difference at 1 ms)
md_uri          = '';         % machine description ('' = data/md_dd4 of the fork)

%% plant, controllers and allocator
pds_root = fileparts(mfilename('fullpath'));            % this file is in the IMAS-PDS root
fork = fullfile(pds_root,'local_installs','plasmaless-tokamak-circuits');
pcs = fullfile(pds_root,'local_installs','pcs','modules');   % PCSSP of my_contr (PCS_PATH)
addpath(fullfile(fork,'plant'), fullfile(fork,'utils'), fullfile(pcs,'KCURR_PFPO1'), ...
    fullfile(pcs,'VS'), fullfile(pcs,'VS1allocator'));
if ~startsWith(model, filesep)
    model = fullfile(fork,'linear_models',model);
end

% plant: the continuous matrices of the loader (its dt only sets the actor's discrete
% model, not used here)
[~, pf_active_md] = load_machine_description(md_uri);
warn_state = warning('off','linear_model:unstable');      % vertical models: unstable by design
plant = configure_linear_model(model, 1e-3, pf_active_md);
warning(warn_state);
sys_plant = ss(plant.A_c, plant.B_c, plant.C_c, 0);

% measurements as my_contr takes them: I = S*x_sys the 11 KCURR currents (CS1 = mean of
% CS1U and CS1L), I_VS3 = c_VS3*x_sys = (I_VS3U - I_VS3L)/2, Zpl = c_z*x_sys,
% Ip = c_Ip*x_sys (zero if the model has no Ip)
[~, X2I] = circuit_currents(plant.x_to_coil', plant.coil_names);   % one row per state
S = X2I(:,1:11)';
c_VS3 = X2I(:,12)';
c_z = plant.C_c(strcmp(plant.y_names,'zmag'),:);
c_Ip = zeros(1,plant.n_states);
if any(strcmp(plant.y_names,'Ip'))
    c_Ip = plant.C_c(strcmp(plant.y_names,'Ip'),:);
end

% KCURR: K_e is 11 x 12, the last column on the Ip error (zero in the noIp set)
G = KCURR_PFPO1_load_gains_tp();
K_e = G.(['Ke_' gain_set]);
I_ref = zeros(11,1);
dIp_ref = 0;

% VS: the continuous filters and gains of Config_VS_structure.m (pcssp_VS_loadfp.m
% discretises the filters with Tustin at 1 ms for the Simulink blocks and sets KZ_VS1).
% Mode 3 ("VS1+VS3 controller") has its own _VS13 gains for VS3, and VS1 = KIVS1_VS13*I_VS3.
VS = Config_VS_structure();
g = VS.Gains(1);                                          % the 3 scheduling levels are equal
KZ_VS1 = 15000/15e6;
sys_VS1 = ss(g.AVS1, g.BVS1, g.CVS1, g.DVS1);
if VS_control_mode == 3
    KZ_VS3 = g.KZ_VS13;
    KI_VS3 = g.KIVS3_VS13;
    sys_VS3 = ss(g.AVS13, g.BVS13, g.CVS13, g.DVS13);
else
    KZ_VS3 = g.KZ_VS3;
    KI_VS3 = g.KI_VS3;
    sys_VS3 = ss(g.AVS3, g.BVS3, g.CVS3, g.DVS3);
end
KIVS1_VS13 = g.KIVS1_VS13;

% allocator (VS1_allocator.slx): projVS1, Fovs1, KV2toV5 and the low pass 1/(tauF s + 1)
alloc = VS1_allocator_loadtp();

% power supplies (pcssp powerSupplies, PSModel 1), linear part only: a delay of 3 ms on the
% CS/PF and VS1 supplies, 2.5 ms plus a 7.5 ms lag on VS3; delays as 2nd-order Pade
if use_PS
    ps_delay = ss(pade(tf(1,1,'IODelay',3e-3), 2));
    ps_VS3 = ss(pade(tf(1,1,'IODelay',2.5e-3), 2)*tf(1,[7.5e-3 1]));
    blocks = [repmat({ps_delay},1,13), {ps_VS3}];        % VS3 is the 14th supply
    sys_PS = append(blocks{:});
else
    sys_PS = ss(eye(14));                                 % no states, voltages pass through
end

%% the loop
%            ┌────┐      ┌───────────────────────┐ Zdot        ┌──────────┐
%       ┌───►│ VS ├─────►│ PLANT (with KCURR, PS,├───────►(+)─►│ 1/(1+τs) ├──┬──► y
%       │    └────┘      │ allocator) + w        │         ▲   └──────────┘  │
%       │                └───────────────────────┘         d                 │
%       └────────────────────────────────────────────────────────────────────┘
% d: noise on the measured Zdot (controller_input_preprocess.m); w: voltages added to the
% 14 supplies after the power supplies (the kick of muscle_control2nice.m).

%% PLANT block: linear model, power supplies, KCURR and allocator, already closed
% Inputs: the two VS commands and w. Outputs: Zdot (exact derivative of zmag), I_VS3 (as
% the VS law measures it), zmag, the 14 coil currents and V_VS3, the voltage applied to the
% VS3 supply.
% [X; inputs] is set to the identity, one unit vector per column, with X = [x_sys; x_alloc;
% x_PS]: every signal below is then the matrix that gives it from [X; inputs], and dX/dt
% and the outputs computed this way are [A B] and [C D] of the block.
sizes = [plant.n_states, use_allocator, size(sys_PS.A,1)];
N = sum(sizes);
nc = N + 16;                                              % columns: states, 2 VS commands, 14 w
parts = mat2cell(eye(nc), [sizes 1 1 14], nc);
[x_sys, x_alloc, x_PS, VS1_cmd, VS3_cmd, w] = parts{:};   % x_alloc, x_PS empty if off

% KCURR (P = 1, I = D = 0), Ip error as a deviation from the reference (Ip gains only)
I = S*x_sys;
Ip = c_Ip*x_sys;
V_c = use_KCURR*K_e*[I_ref - I; dIp_ref - Ip];

% allocator: the KCURR voltage along VS1 (p) goes to the VS1 supply, its low-pass part (f)
% back to the PF supplies
V_VS1 = VS1_cmd;
if use_allocator
    p = alloc.projVS1*V_c;
    f = x_alloc;
    dx_alloc = (V_VS1 - p - x_alloc)/alloc.tauF;
    V_c = alloc.Fovs1*alloc.Fovs1'*V_c - alloc.KV2toV5*f;
    V_VS1 = V_VS1 - p - f;
else
    dx_alloc = zeros(0,nc);
end

% the 14 supply voltages (controller_output_preprocess.m), CS1 half on CS1U and CS1L
V_Cspf = [V_c(1:2,:); V_c(3,:)/2; V_c(3,:)/2; V_c(4:11,:)];
V_cmd = [V_Cspf(1:11,:); V_VS1; V_Cspf(12,:); VS3_cmd];

% power supplies, then w added to the voltages they apply
dx_PS = sys_PS.A*x_PS + sys_PS.B*V_cmd;
V = sys_PS.C*x_PS + sys_PS.D*V_cmd + w;

% linear model
dx_sys = sys_plant.A*x_sys + sys_plant.B*V;

dX = [dx_sys; dx_alloc; dx_PS];
out = [c_z*dx_sys; c_VS3*x_sys; c_z*x_sys; plant.x_to_coil*x_sys; V(14,:)];
P = ss(dX(:,1:N), dX(:,N+1:end), out(:,1:N), out(:,N+1:end), ...
    'InputName', [{'VS1_cmd', 'VS3_cmd'}, strcat('w_', plant.input_names)], ...
    'OutputName', [{'Zdot', 'I_VS3', 'zmag'}, strcat('I_', plant.coil_names), {'V_VS3'}]);

%% VS block (pcssp_VS.slx): inputs y (measured, filtered Zdot) and I_VS3
% "VS1 controller" (mode 1), "VS3 controller" (mode 2), "VS1+VS3 controller" (mode 3)
VS1_law = sys_VS1*[KZ_VS1*Ip_ref, 0];                     % F1(KZ_VS1*Ip_ref*y)
VS3_law = sys_VS3*[-KZ_VS3*Ip_ref, -KI_VS3];              % F3(-(KZ_VS3*Ip_ref*y + KI_VS3*I_VS3))
switch VS_control_mode
    case 1
        VS_blk = [VS1_law; ss([0 0])];
    case 2
        VS_blk = [ss([0 0]); VS3_law];
    case 3
        VS_blk = [ss([0 KIVS1_VS13]); VS3_law];           % VS1 brings the VS3 current to zero
    otherwise
        VS_blk = ss(zeros(2));
end
VS_blk.InputName = {'y', 'I_VS3'};
VS_blk.OutputName = {'VS1_cmd', 'VS3_cmd'};

%% Zdot measurement: noise d added to Zdot, then the filter
sum_d = sumblk('Zdot_m = Zdot + d');
F = tf(1, [tau 1]);
F.InputName = 'Zdot_m';
F.OutputName = 'y';

%% closed loop, with an analysis point on y
% sys_cl from d to y is 1/(1 + tau s - C(s) P(s)), C the VS law on Zdot and P the PLANT
% block from the VS command to Zdot.
T = connect(P, VS_blk, sum_d, F, [{'d'}, strcat('w_', plant.input_names)], ...
    [{'y', 'zmag'}, strcat('I_', plant.coil_names), {'V_VS3'}], 'y');
L_y = getLoopTransfer(T, 'y', -1);    % loop opened at y (convention: closed loop = 1/(1+L_y))
S_y = getSensitivity(T, 'y');         % sensitivity at y, 1/(1 + L_y)
sys_cl = ss(T);                       % plain closed loop for eig, damp, lsim, covar
