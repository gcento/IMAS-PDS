function log_row(t, u)
% Appends one line [t, u] to simulink_log.csv in the run directory (the header is written by the init).
f = fopen(fullfile(getenv('RUN_DIR'), 'simulink_log.csv'), 'a');
if f < 0, return; end
fprintf(f, '%.6f', t);
fprintf(f, ',%.9g', u(:));
fprintf(f, '\n');
fclose(f);
end
