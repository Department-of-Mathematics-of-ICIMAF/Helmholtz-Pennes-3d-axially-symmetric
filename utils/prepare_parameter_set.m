function [params, metadata] = prepare_parameter_set(case_def)
%PREPARE_PARAMETER_SET Separa metadatos y expande el preset de tejido.
    metadata = struct('name', '', 'tissue', '');
    if isfield(case_def, 'name')
        metadata.name = case_def.name;
    end
    if isfield(case_def, 'tissue')
        metadata.tissue = case_def.tissue;
    end

    params = case_def;
    fields_to_remove = {};
    if isfield(params, 'name')
        fields_to_remove{end+1} = 'name'; %#ok<AGROW>
    end
    if isfield(params, 'tissue')
        fields_to_remove{end+1} = 'tissue'; %#ok<AGROW>
    end
    if ~isempty(fields_to_remove)
        params = rmfield(params, fields_to_remove);
    end

    if isempty(metadata.tissue)
        return
    end

    tissue = lower(strtrim(metadata.tissue));
    presets = struct( ...
        'phantom', struct('c', 1568, 'aten', 8.55, 'rho', 1044), ...
        'agua',    struct('c', 1500, 'aten', 0,    'rho', 1000), ...
        'piel',    struct('c', 1558, 'aten', 24,   'rho', 1200), ...
        'grasa',   struct('c', 1478, 'aten', 5.58, 'rho', 950), ...
        'musculo', struct('c', 1547, 'aten', 12.7, 'rho', 1050));

    if isfield(presets, tissue)
        preset = presets.(tissue);
        overrides = params;
        params = preset;
        override_fields = fieldnames(overrides);
        for i = 1:numel(override_fields)
            field = override_fields{i};
            params.(field) = overrides.(field);
        end
    elseif ~strcmp(tissue, 'personalizado')
        error('prepare_parameter_set:UnknownTissue', ...
              ['Tejido desconocido: %s. Use phantom, agua, piel, grasa, ' ...
               'musculo o personalizado (con c, aten y rho explicitos).'], ...
              metadata.tissue);
    end
end
