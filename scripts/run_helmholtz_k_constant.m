function results = run_helmholtz_k_constant(input_file)
%RUN_HELMHOLTZ_K_CONSTANT Ejecuta en secuencia los casos con k constante.
%   RESULTS = RUN_HELMHOLTZ_K_CONSTANT() lee inputs/helmholtz_k_constant.json.
%   RESULTS = RUN_HELMHOLTZ_K_CONSTANT(INPUT_FILE) usa otro fichero JSON.
%
%   Cada elemento de RESULTS contiene last_max, preim, XYQ, params,
%   metadata y output_dir. Las corridas no se ejecutan en paralelo.
    root_dir = fileparts(fileparts(mfilename('fullpath')));
    addpath(fullfile(root_dir, 'utils'));
    if nargin < 1
        input_file = '';
    end
    results = run_helmholtz_cases('constant', input_file);
end
