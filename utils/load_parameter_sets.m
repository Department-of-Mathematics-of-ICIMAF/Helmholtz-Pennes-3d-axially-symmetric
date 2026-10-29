function cases = load_parameter_sets(input_file)
%LOAD_PARAMETER_SETS Lee uno o varios juegos de parametros desde JSON.
    if ~exist(input_file, 'file')
        error('load_parameter_sets:FileNotFound', ...
              'No existe el fichero de parametros: %s', input_file);
    end

    decoded = jsondecode(fileread(input_file));
    if iscell(decoded)
        cases = decoded;
    elseif isstruct(decoded)
        cases = cell(1, numel(decoded));
        for i = 1:numel(decoded)
            cases{i} = decoded(i);
        end
    else
        error('load_parameter_sets:InvalidFormat', ...
              'El JSON debe contener un objeto o un arreglo de objetos.');
    end

    if isempty(cases)
        error('load_parameter_sets:EmptyInput', ...
              'El fichero no contiene juegos de parametros.');
    end
    for i = 1:numel(cases)
        if ~isstruct(cases{i})
            error('load_parameter_sets:InvalidCase', ...
                  'El juego de parametros %d no es un objeto JSON.', i);
        end
    end
end
