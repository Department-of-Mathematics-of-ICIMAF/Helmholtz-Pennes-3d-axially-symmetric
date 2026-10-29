function write_parameter_summary(output_dir, params, metadata, model_type)
%WRITE_PARAMETER_SUMMARY Guarda todos los parametros usados en la corrida.
    summary_file = fullfile(output_dir, 'parametros.txt');
    fid = fopen(summary_file, 'w');
    if fid < 0
        error('write_parameter_summary:OpenFailed', ...
              'No se pudo crear el fichero: %s', summary_file);
    end
    cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>

    fprintf(fid, 'Resumen de parametros de Helmholtz\n');
    fprintf(fid, 'Modelo: %s\n', model_type);
    fprintf(fid, 'Fecha: %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
    if ~isempty(metadata.name)
        fprintf(fid, 'Nombre: %s\n', metadata.name);
    end
    if ~isempty(metadata.tissue)
        fprintf(fid, 'Tejido: %s\n', metadata.tissue);
    end
    fprintf(fid, '\n');

    fields = sort(fieldnames(params));
    for i = 1:numel(fields)
        field = fields{i};
        value = params.(field);
        if isnumeric(value)
            text = mat2str(value, 16);
        elseif islogical(value)
            text = mat2str(value);
        elseif ischar(value)
            text = value;
        elseif isstring(value) && isscalar(value)
            text = char(value);
        else
            text = sprintf('<%s>', class(value));
        end
        fprintf(fid, '%s = %s\n', field, text);
    end
end
