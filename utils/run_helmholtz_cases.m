function results = run_helmholtz_cases(model_type, input_file)
%RUN_HELMHOLTZ_CASES Ejecuta secuencialmente un conjunto de casos Helmholtz.
    root_dir = get_project_root();
    utils_dir = fullfile(root_dir, 'utils');
    original_path = path;
    restore_path = onCleanup(@() path(original_path));
    addpath(utils_dir);

    switch model_type
        case 'constant'
            case_dir = fullfile(root_dir, 'Helmholtz', 'Homogeneous Helmholtz');
            default_input = fullfile(root_dir, 'inputs', 'helmholtz_k_constant.json');
        case 'variable'
            case_dir = fullfile(root_dir, 'Helmholtz', 'Variable k Helmholtz');
            default_input = fullfile(root_dir, 'inputs', 'helmholtz_k_variable.json');
        otherwise
            error('run_helmholtz_cases:UnknownModel', ...
                  'Modelo desconocido: %s', model_type);
    end
    if nargin < 2 || isempty(input_file)
        input_file = default_input;
    end

    cases = load_parameter_sets(input_file);
    results = cell(numel(cases), 1);
    addpath(case_dir);
    % Los dos modelos usan los mismos nombres de funcion: limpiar versiones
    % cacheadas para garantizar que se ejecute la carpeta seleccionada.
    clear mainHelmholtzMultipatch_3Daxial ...
          radiationmatricesHelmholtz_3Daxial ParametMultipatch_mediodominio;

    for i = 1:numel(cases)
        [params, metadata] = prepare_parameter_set(cases{i});
        output_dir = create_run_output_dir(fullfile(case_dir, 'outputs'), ...
                                           params, metadata, i, model_type);
        write_parameter_summary(output_dir, params, metadata, model_type);
        fprintf('\n[%d/%d] Ejecutando Helmholtz %s: %s\n', ...
                i, numel(cases), model_type, output_dir);

        try
            [last_max, preim, XYQ, quarter_points] = ...
                mainHelmholtzMultipatch_3Daxial(params, output_dir);
        catch exception
            close(findall(0, 'Type', 'figure'));
            rethrow(exception);
        end

        % Valores derivados disponibles para Pennes en params y en results.
        params.quarter_points = quarter_points;
        horizontal_semi_axis = abs(quarter_points(1,1));
        vertical_semi_axis = mean(abs(quarter_points(2:3,2)-last_max(1)));
        params.quarter_semi_axes = [horizontal_semi_axis, vertical_semi_axis];
        write_parameter_summary(output_dir, params, metadata, model_type);

        save(fullfile(output_dir, 'helmholtz_results.mat'), ...
             'last_max', 'preim', 'XYQ', 'quarter_points', 'params');
        result = struct('last_max', last_max, 'preim', preim, 'XYQ', XYQ, ...
                        'quarter_points', quarter_points, 'params', params, ...
                        'metadata', metadata, 'output_dir', output_dir);
        results{i} = result;

        % Punto de integracion para Pennes: consumir aqui last_max, preim,
        % XYQ y quarter_points (o results{i}) al implementar la segunda etapa.
        fprintf('Corrida finalizada. Resultados: %s\n', ...
                fullfile(output_dir, 'helmholtz_results.mat'));
    end
end
