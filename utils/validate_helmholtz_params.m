function validate_helmholtz_params(params, required_fields)
%VALIDATE_HELMHOLTZ_PARAMS Comprueba los parametros requeridos por un main.
    if ~isstruct(params)
        error('validate_helmholtz_params:InvalidParameters', ...
              'Los parametros deben estar en una estructura.');
    end
    for i = 1:numel(required_fields)
        field = required_fields{i};
        if ~isfield(params, field) || isempty(params.(field))
            error('validate_helmholtz_params:MissingParameter', ...
                  'Falta el parametro requerido: %s', field);
        end
    end
end
