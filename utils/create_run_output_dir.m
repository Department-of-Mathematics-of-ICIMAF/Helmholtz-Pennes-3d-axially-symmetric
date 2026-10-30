function output_dir = create_run_output_dir(base_output_dir, params, metadata, run_index, model_type)
%CREATE_RUN_OUTPUT_DIR Crea una subcarpeta breve y unica para una corrida.
    if ~isempty(metadata.name)
        label = metadata.name;
    elseif ~isempty(metadata.tissue)
        label = metadata.tissue;
    elseif strcmp(model_type, 'variable')
        label = 'capas';
    else
        label = 'personalizado';
    end
    label = regexprep(strtrim(label), '[^a-zA-Z0-9_-]+', '_');
    label = regexprep(label, '^_+|_+$', '');
    if isempty(label)
        label = 'corrida';
    end

    folder_name = sprintf('%03d_%s_n%d', run_index, label, params.n1);
    output_dir = fullfile(base_output_dir, folder_name);
    suffix = 2;
    while exist(output_dir, 'dir')
        output_dir = fullfile(base_output_dir, sprintf('%s_%02d', folder_name, suffix));
        suffix = suffix + 1;
    end
    mkdir(output_dir);
end
