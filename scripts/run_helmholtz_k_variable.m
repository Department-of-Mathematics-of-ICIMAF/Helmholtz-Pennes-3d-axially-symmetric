
%RUN_HELMHOLTZ_K_VARIABLE Ejecuta en secuencia los casos con k variable.
%   RESULTS = RUN_HELMHOLTZ_K_VARIABLE() lee inputs/helmholtz_k_variable.json.
%   RESULTS = RUN_HELMHOLTZ_K_VARIABLE(INPUT_FILE) usa otro fichero JSON.
%
%   Cada elemento de RESULTS contiene last_max, preim, XYQ, params,
%   metadata y output_dir. Las corridas no se ejecutan en paralelo.
input_file = 'inputs/helmholtz_k_variable.json';
root_dir = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root_dir, 'utils'));
results = run_helmholtz_cases('variable', input_file);
