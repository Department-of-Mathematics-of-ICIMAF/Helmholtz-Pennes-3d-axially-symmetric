
%RUN_HELMHOLTZ_K_CONSTANT Ejecuta en secuencia los casos con k constante.
%   RESULTS = RUN_HELMHOLTZ_K_CONSTANT() lee inputs/helmholtz_k_constant.json.
%   RESULTS = RUN_HELMHOLTZ_K_CONSTANT(INPUT_FILE) usa otro fichero JSON.
%
%   Cada elemento de RESULTS contiene last_max, preim, XYQ, params,
%   metadata y output_dir. Las corridas no se ejecutan en paralelo.
input_file = 'inputs/helmholtz_k_constant.json';
root_dir = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root_dir, 'utils'));
results = run_helmholtz_cases('constant', input_file);
