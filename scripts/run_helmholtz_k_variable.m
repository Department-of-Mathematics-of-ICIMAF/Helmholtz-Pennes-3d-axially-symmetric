
%RUN_HELMHOLTZ_K_VARIABLE Ejecuta en secuencia los casos con k variable.
%   RESULTS = RUN_HELMHOLTZ_K_VARIABLE() lee inputs/helmholtz_k_variable.json.
%   RESULTS = RUN_HELMHOLTZ_K_VARIABLE(INPUT_FILE) usa otro fichero JSON.
%
%   Cada elemento de RESULTS contiene last_max, preim, XYQ, params,
%   metadata y output_dir. Las corridas no se ejecutan en paralelo.

% ---------------------------------------------------------------------
% ETAPA 1: Helmholtz
% ---------------------------------------------------------------------
input_file = 'inputs/helmholtz_k_variable.json';
root_dir = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root_dir, 'utils'));
addpath(fullfile(root_dir, 'inputs'));
results = run_helmholtz_cases('variable', input_file);

% ---------------------------------------------------------------------
% ETAPA 2: Pennes
% ---------------------------------------------------------------------
n = numel(results);
fprintf('\n=== ETAPA 2: Pennes ===\n');
P = cell(n,1);

case_name = 'Pennes_3_layers';
% Raiz de salida: <pennes_dir>/outputs
pennes_outputs_root = fullfile(root_dir, 'Pennes', case_name, 'outputs');
if ~exist(pennes_outputs_root, 'dir')
    mkdir(pennes_outputs_root);
end

addpath(fullfile(root_dir, 'Pennes', case_name));
addpath(fullfile(root_dir, 'Pennes', case_name, 'utils'));

for i = 1:n
    fmat = fullfile(results{i}.output_dir, 'helmholtz_results.mat');
    if exist(fmat, 'file') ~= 2
        error('No existe: %s', fmat);
    end

    S = load(fmat, 'params', 'XYQ');
    p = S.params;
    R  = p.r;
    a  = p.a;
    if ~isfield(p,'rl')
        % Intentar z_c con los nombres habituales
        if ~isempty(p.zc)
            zc = p.zc;
            rl = sqrt(zc^2 + a^2);
            fprintf('  [rl] calculado desde z_c=%.6g -> rl=%.6g\n', zc, rl);
        else
            error('  [rl] no encontrado (ni rl ni z_c)');
        end
    else
        rl = p.rl;
    end
    
    % ---- Construccion del nombre de carpeta unico ----
    base_folder_name = sprintf('%g_per_layers_%g', a, rl);
    outdir_pennes = fullfile(pennes_outputs_root, base_folder_name);
    
    if exist(outdir_pennes, 'dir')
        % Buscar el primer sufijo _01, _02, ... disponible
        idx = 1;
        while true
            candidate = fullfile(pennes_outputs_root, ...
                                 sprintf('%s_%02d', base_folder_name, idx));
            if ~exist(candidate, 'dir')
                outdir_pennes = candidate;
                break;
            end
            idx = idx + 1;
        end
    end
    % mainPennesODE crea la carpeta si no existe
    % (aqui aun no existe, pero la pasamos como outdir)
    
    
    fprintf('\n[%d/%d] Pennes: %s\n', i, n, outdir_pennes);
    fprintf('  R=%.6g  a=%.6g  rl=%.6g\n', ...
            R, a, rl);

    % ---- LLAMADA A PENNES (nueva firma de 9 argumentos) ----
    pennes_results = mainPennesODE(R, a, rl, fmat, ...
                                   48, 30, ...
                                   outdir_pennes, ...
                                   1, ...
                                   3*pi/2);

    % Metadatos del lado Helmholtz para trazabilidad
    pennes_meta = struct();
    pennes_meta.case_name          = case_name;
    pennes_meta.helmholtz_outdir   = results{i}.output_dir;
    pennes_meta.helmholtz_metadata = results{i}.metadata;
    pennes_meta.helmholtz_params   = p;
    pennes_meta.pennes_filenameQ   = fmat;
    save(fullfile(outdir_pennes,'pennes_meta.mat'), 'pennes_meta');

    P{i} = struct('pennes_results',   pennes_results, ...
                  'outdir',           outdir_pennes, ...
                  'case_name',        case_name, ...
                  'helmholtz_outdir', results{i}.output_dir);
end


fprintf('\nIntegracion Helmholtz -> Pennes finalizada (%d casos).\n', n);
fprintf('Resultados en: %s\n', pennes_outputs_root);