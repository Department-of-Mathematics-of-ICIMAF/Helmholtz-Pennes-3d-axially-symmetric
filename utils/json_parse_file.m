function data = json_parse_file(filename)
%JSON_PARSE_FILE Parser JSON minimalista en MATLAB puro (sin jsondecode).
%   Soporta: objetos {}, arrays [], strings "...", numeros, true/false/null.
%   No soporta: \uXXXX, notacion cientifica rara, comentarios (los quita).

    fid = fopen(filename, 'r');
    if fid < 0
        error('json_parse_file:open', 'No se pudo abrir: %s', filename);
    end
    lines = {};
    while true
        ln = fgetl(fid);
        if ~ischar(ln), break; end
        lines{end+1} = ln; %#ok<AGROW>
    end
    fclose(fid);
    txt = strjoin(lines, sprintf('\n'));

    % Quitar comentarios // ... y /* ... */
    txt = regexprep(txt, '/\*.*?\*/', '', 'dotall');
    txt = regexprep(txt, '//[^\n]*', '');

    % Tokenizar y parsear
    state.txt = txt;
    state.pos = 1;
    state.len = numel(txt);
    [data, state] = parse_value(state);

    skip_ws(state);
    if state.pos <= state.len
        error('json_parse_file:trailing', ...
              'Texto sobrante en posicion %d', state.pos);
    end
end

% -------------------------------------------------------------------------
function [val, s] = parse_value(s)
    s = skip_ws(s);
    if s.pos > s.len
        error('json_parse_file:eof', 'Fin de fichero inesperado');
    end
    c = s.txt(s.pos);
    switch c
        case '{'
            [val, s] = parse_object(s);
        case '['
            [val, s] = parse_array(s);
        case '"'
            [val, s] = parse_string(s);
        case 't'
            expect(s, 'true');  val = true;  s.pos = s.pos + 4;
        case 'f'
            expect(s, 'false'); val = false; s.pos = s.pos + 5;
        case 'n'
            expect(s, 'null');  val = [];    s.pos = s.pos + 4;
        otherwise
            [val, s] = parse_number(s);
    end
end

% -------------------------------------------------------------------------
function [obj, s] = parse_object(s)
    s.pos = s.pos + 1;              % consume '{'
    obj = struct();
    s = skip_ws(s);
    if s.pos <= s.len && s.txt(s.pos) == '}'
        s.pos = s.pos + 1; return
    end
    while true
        [key, s] = parse_string(s);         % clave
        s = skip_ws(s);
        if s.txt(s.pos) ~= ':'
            error('json_parse_file:colon', ...
                  'Se esperaba ":" en posicion %d', s.pos);
        end
        s.pos = s.pos + 1;                  % consume ':'
        [val, s] = parse_value(s);
        obj.(key) = val;
        s = skip_ws(s);
        if s.pos > s.len, error('json_parse_file:eofobj','Objeto sin cerrar'); end
        c = s.txt(s.pos);
        if c == ','
            s.pos = s.pos + 1; s = skip_ws(s);
        elseif c == '}'
            s.pos = s.pos + 1; return
        else
            error('json_parse_file:objsep', ...
                  'Se esperaba "," o "}" en posicion %d', s.pos);
        end
    end
end

% -------------------------------------------------------------------------
function [arr, s] = parse_array(s)
    s.pos = s.pos + 1;              % consume '['
    arr = {};
    s = skip_ws(s);
    if s.pos <= s.len && s.txt(s.pos) == ']'
        s.pos = s.pos + 1; arr = []; return
    end
    while true
        [val, s] = parse_value(s);
        arr{end+1} = val; %#ok<AGROW>
        s = skip_ws(s);
        if s.pos > s.len, error('json_parse_file:eofarr','Array sin cerrar'); end
        c = s.txt(s.pos);
        if c == ','
            s.pos = s.pos + 1; s = skip_ws(s);
        elseif c == ']'
            s.pos = s.pos + 1; break
        else
            error('json_parse_file:arrsep', ...
                  'Se esperaba "," o "]" en posicion %d', s.pos);
        end
    end
    % Si todos son numeros -> vector double
    if ~isempty(arr) && all(cellfun(@(x) isnumeric(x) && isscalar(x), arr))
        arr = cell2mat(arr);
    end
end

% -------------------------------------------------------------------------
function [str, s] = parse_string(s)
    if s.txt(s.pos) ~= '"'
        error('json_parse_file:str', 'Se esperaba string en %d', s.pos);
    end
    s.pos = s.pos + 1;
    out = '';
    while s.pos <= s.len
        c = s.txt(s.pos);
        if c == '"'
            s.pos = s.pos + 1; str = out; return
        elseif c == '\'
            s.pos = s.pos + 1;
            e = s.txt(s.pos);
            switch e
                case '"',  out(end+1) = '"';
                case '\',  out(end+1) = '\';
                case '/',  out(end+1) = '/';
                case 'n',  out(end+1) = sprintf('\n');
                case 't',  out(end+1) = sprintf('\t');
                case 'r',  out(end+1) = sprintf('\r');
                case 'b',  out(end+1) = sprintf('\b');
                case 'f',  out(end+1) = sprintf('\f');
                otherwise
                    error('json_parse_file:esc', ...
                          'Escape no soportado: \\%s', e);
            end
            s.pos = s.pos + 1;
        else
            out(end+1) = c; %#ok<AGROW>
            s.pos = s.pos + 1;
        end
    end
    error('json_parse_file:stropen', 'String sin cerrar');
end

% -------------------------------------------------------------------------
function [num, s] = parse_number(s)
    start = s.pos;
    while s.pos <= s.len && any(s.txt(s.pos) == '+-0123456789.eE')
        s.pos = s.pos + 1;
    end
    tok = s.txt(start:s.pos-1);
    num = str2double(tok);
    if isnan(num)
        error('json_parse_file:num', 'Numero invalido: %s', tok);
    end
end

% -------------------------------------------------------------------------
function expect(s, word)
    if s.pos + numel(word) - 1 > s.len || ...
       ~strcmp(s.txt(s.pos:s.pos+numel(word)-1), word)
        error('json_parse_file:kw', 'Se esperaba "%s"', word);
    end
end

% -------------------------------------------------------------------------
function s = skip_ws(s)
    while s.pos <= s.len && any(s.txt(s.pos) == sprintf(' \t\r\n'))
        s.pos = s.pos + 1;
    end
end