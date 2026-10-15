function export_plasmaless_model(plasmaless_repo, out_dir, md_uri)
% Export the plasmaless coil+vessel circuit model (no plasma) as a linear model folder
% (linear_model.json + A, B, C .npy) for the actor setting linear_model_dir, so that it is
% selected like the NICE models of linear_models/.
%   export_plasmaless_model(plasmaless_repo, out_dir)          % md: <repo>/data/md_dd4
%   export_plasmaless_model(plasmaless_repo, out_dir, md_uri)
% Needs IMAS-MATLAB (reads the machine description) and the Control System Toolbox; run
% it on a compute node.
%
% The circuit model is absolute (states = coil and passive currents), so the descriptor
% has coordinates 'absolute' and no operating point: the actor starts it from the F_INIT
% coil currents at t0, passive currents 0, as with linear_model_dir = '' (needs the
% plasmaless actor branch IMAS_Muscle3_absolute).
% After writing, the folder is read back with configure_linear_model and checked against
% configure_plasmaless_model over 1000 steps with random voltages, from the same initial
% coil currents.

plasmaless_repo = char(plasmaless_repo);
out_dir = char(out_dir);
if nargin < 3 || isempty(md_uri)
    md_uri = ['imas:hdf5?path=' fullfile(plasmaless_repo,'data','md_dd4')];
end
md_uri = char(md_uri);   % imas_open needs a char URI
addpath(fullfile(plasmaless_repo,'imas_model'));

dt = 0.002;   % only used to build the circuit model's M and R and for the check
[em_coupling, pf_active, pf_passive] = load_plasmaless_md(md_uri);
native = configure_plasmaless_model(em_coupling, pf_active, pf_passive, dt);
n = native.n_states;
n_coils = native.n_coils;
% continuous matrices exactly as configure_plasmaless_model builds them
Ac = -native.M\native.R;
Bc = native.M\eye(n,n_coils);
r_coil = diag(native.R(1:n_coils,1:n_coils));

if ~isfolder(out_dir)
    mkdir(out_dir);
end
write_npy(fullfile(out_dir,'A.npy'), Ac);
write_npy(fullfile(out_dir,'B.npy'), Bc);
write_npy(fullfile(out_dir,'C.npy'), zeros(0,n));

d.description = ['plasmaless coil+vessel circuit model (no plasma), ', num2str(n_coils), ...
    ' coils + ', num2str(native.n_loops), ' passive loops: dI/dt = A I + B U, ', ...
    'absolute currents and voltages'];
d.shot = 0;
d.time = 0;
d.files = struct('A','A.npy','B','B.npy','C','C.npy');
d.states = native.state_names;
d.inputs = native.coil_names;
d.coil_states = struct();
for ii = 1:n_coils
    d.coil_states.(native.coil_names{ii}) = native.coil_names{ii};
end
d.passive_states = native.loop_names;
d.outputs = {};
d.coordinates = 'absolute';
d.units = struct('currents','A','voltages','V');
d.provenance = sprintf(['machine description %s; A = -M\\R, B = M\\eye(n,n_coils) as ', ...
    'configure_plasmaless_model.m; written by linear_models/export_plasmaless_model.m on %s'], ...
    md_uri, char(datetime('today','Format','yyyy-MM-dd')));
d.notes = {['Not tied to a pulse: shot and time are 0. Absolute model: the actor starts ', ...
    'it from the F_INIT coil currents at t0 (passive currents 0) and ignores the setting ', ...
    'operating_point, exactly as the circuit model of linear_model_dir = ''''.']};
fid = fopen(fullfile(out_dir,'linear_model.json'),'w');
fprintf(fid,'%s\n',jsonencode(d,'PrettyPrint',true));
fclose(fid);

% check: the folder read by the actor's own loader against the circuit model
lin = configure_linear_model(out_dir, dt, pf_active);
if ~isequal(lin.A, native.A) || ~isequal(lin.B, native.B)
    error('export_plasmaless_model:check', ['discrete matrices differ: max |dA| %g, ', ...
        'max |dB| %g.'], max(abs(lin.A(:)-native.A(:))), max(abs(lin.B(:)-native.B(:))));
end
if ~strcmp(lin.coordinates,'absolute') || any(lin.op.coil_current) || any(lin.op.voltage)
    error('export_plasmaless_model:check', ['the actor does not read the folder as an ', ...
        'absolute model (plasmaless repository without coordinates support?).']);
end
rng(0);
i0 = 1e4*randn(n_coils,1);
i0(r_coil > 0) = 200*randn(nnz(r_coil > 0),1);
x_native = [i0; zeros(native.n_loops,1)];
x_lin = lin.x_to_coil'*i0;          % initial state as muscle_plasmaless_actor builds it
err = 0;
for k = 1:1000
    u = 100*randn(n_coils,1);
    x_native = native.A*x_native + native.B*u;
    x_lin = lin.A*x_lin + lin.B*u;
    i_lin = [lin.x_to_coil*x_lin; lin.x_to_loop*x_lin];
    err = max(err, max(abs(i_lin - x_native))/max(abs(x_native)));
end
if err > 1e-9
    error('export_plasmaless_model:check','currents differ from the circuit model (relative %g).', err);
end
fprintf(['export_plasmaless_model: wrote %s (%d states, %d inputs); same discrete A and B ', ...
    'as the circuit model, currents equal over 1000 steps (max relative difference %.1e)\n'], ...
    out_dir, n, n_coils, err);
end

function write_npy(filename, M)
% minimal .npy writer: float64, 2-D, Fortran order (MATLAB's column-major storage)
header = sprintf('{''descr'': ''<f8'', ''fortran_order'': True, ''shape'': (%d, %d), }', ...
    size(M,1), size(M,2));
pad = mod(-(10 + numel(header) + 1), 64);    % total header length a multiple of 64
header = [header repmat(' ',1,pad) newline];
fid = fopen(filename,'w','ieee-le');
fwrite(fid,[147 double('NUMPY') 1 0],'uint8');
fwrite(fid,numel(header),'uint16');
fwrite(fid,header,'char');
fwrite(fid,M(:),'float64');
fclose(fid);
end
