function [horizontal, vertical_up, vertical_down] = ...
    quarter_pressure_points(z_f, preim, u, space, geometry, npoints_axis)
%QUARTER_PRESSURE_POINTS Busca los primeros puntos cercanos a |u_f|/4.
%   [HORIZONTAL, VERTICAL_UP, VERTICAL_DOWN] =
%   QUARTER_PRESSURE_POINTS(Z_F, PREIM, U, SPACE, GEOMETRY, NPOINTS_AXIS)
%   evalua la solucion sobre la linea focal y sobre el eje de simetria.
%
%   En cada direccion se parte del foco y se toma el primer par de muestras
%   consecutivas cuyos valores de |u| encierran |u(0,z_f)|/4. De esos dos
%   extremos se devuelve el mas cercano al valor objetivo. La busqueda
%   horizontal recorre PREIM hacia rho creciente; las busquedas verticales
%   usan puntos parametricos (0,eta) y recorren z creciente o decreciente.
%
%   Cada salida es un vector fila con columnas:
%       [rho, z, |u|, numero_de_parche, xi, eta].

    if nargin < 6 || isempty(npoints_axis)
        npoints_axis = 1000;
    end
    if npoints_axis < 2
        error('quarter_pressure_points:InvalidSampling', ...
              'npoints_axis debe ser al menos 2.');
    end

    % Evaluacion sobre las preimagenes de la linea horizontal z = z_f.
    h_rho = [];
    h_z = [];
    h_abs = [];
    h_patch = [];
    h_xi = [];
    h_eta = [];
    for iptc = 1:numel(geometry)
        if preim(iptc).n == 0
            continue
        end
        values = zeros(preim(iptc).n, 1);
        for k = 1:preim(iptc).n
            uv = sp_eval(u(space.gnum{iptc}), space.sp_patch{iptc}, ...
                         geometry(iptc), ...
                         {preim(iptc).xi(k), preim(iptc).eta(k)});
            values(k) = abs(uv);
        end
        n = preim(iptc).n;
        h_rho = [h_rho; preim(iptc).rho(:)]; %#ok<AGROW>
        h_z = [h_z; preim(iptc).z(:)]; %#ok<AGROW>
        h_abs = [h_abs; values]; %#ok<AGROW>
        h_patch = [h_patch; repmat(iptc, n, 1)]; %#ok<AGROW>
        h_xi = [h_xi; preim(iptc).xi(:)]; %#ok<AGROW>
        h_eta = [h_eta; preim(iptc).eta(:)]; %#ok<AGROW>
    end
    if isempty(h_rho)
        error('quarter_pressure_points:EmptyHorizontalLine', ...
              'No hay preimagenes para evaluar la linea focal.');
    end

    [h_rho, order] = sort(h_rho);
    h_z = h_z(order);
    h_abs = h_abs(order);
    h_patch = h_patch(order);
    h_xi = h_xi(order);
    h_eta = h_eta(order);

    rho_scale = max(1, max(abs(h_rho)));
    z_scale = max(1, abs(z_f));
    axis_tol = 1e-9*rho_scale;
    line_tol = 1e-8*z_scale;
    focal_candidates = find(abs(h_rho) <= axis_tol & ...
                            abs(h_z-z_f) <= line_tol);
    if isempty(focal_candidates)
        error('quarter_pressure_points:MissingFocus', ...
              'No se encontro la preimagen del foco (0,z_f).');
    end
    [~, focal_local] = min(abs(h_rho(focal_candidates)));
    focal_index = focal_candidates(focal_local);
    focal_abs_u = h_abs(focal_index);
    if ~isfinite(focal_abs_u) || focal_abs_u <= 0
        error('quarter_pressure_points:InvalidFocusValue', ...
              'El valor absoluto de la presion focal debe ser positivo.');
    end
    target = focal_abs_u/4;

    % El foco es la primera muestra; despues se avanza hacia rho creciente.
    outward = find(h_rho > h_rho(focal_index) + axis_tol);
    horizontal_indices = [focal_index; outward];
    selected = first_bracket(h_abs(horizontal_indices), target, 'horizontal');
    index = horizontal_indices(selected);
    horizontal = [h_rho(index), h_z(index), h_abs(index), ...
                  h_patch(index), h_xi(index), h_eta(index)];

    % Evaluacion vertical sobre (0,eta) en todos los parches.
    eta_samples = linspace(0, 1, npoints_axis);
    v_rho = [];
    v_z = [];
    v_abs = [];
    v_patch = [];
    v_xi = [];
    v_eta = [];
    for iptc = 1:numel(geometry)
        [uv, Fv] = sp_eval(u(space.gnum{iptc}), space.sp_patch{iptc}, ...
                           geometry(iptc), {0, eta_samples});
        rho_values = squeeze(Fv(1,:,:));
        z_values = squeeze(Fv(2,:,:));
        abs_values = abs(squeeze(uv));
        n = numel(z_values);
        v_rho = [v_rho; rho_values(:)]; %#ok<AGROW>
        v_z = [v_z; z_values(:)]; %#ok<AGROW>
        v_abs = [v_abs; abs_values(:)]; %#ok<AGROW>
        v_patch = [v_patch; repmat(iptc, n, 1)]; %#ok<AGROW>
        v_xi = [v_xi; zeros(n, 1)]; %#ok<AGROW>
        v_eta = [v_eta; eta_samples(:)]; %#ok<AGROW>
    end

    % Conservar unicamente el eje fisico rho = 0 y ordenarlo por z.
    vertical_axis_tol = 1e-9*max(1, max(abs(v_rho)));
    on_axis = find(abs(v_rho) <= vertical_axis_tol);
    if isempty(on_axis)
        error('quarter_pressure_points:MissingVerticalAxis', ...
              'No se encontraron puntos (0,eta) sobre el eje fisico.');
    end
    v_rho = v_rho(on_axis);
    v_z = v_z(on_axis);
    v_abs = v_abs(on_axis);
    v_patch = v_patch(on_axis);
    v_xi = v_xi(on_axis);
    v_eta = v_eta(on_axis);
    [v_z, order] = sort(v_z);
    v_rho = v_rho(order);
    v_abs = v_abs(order);
    v_patch = v_patch(order);
    v_xi = v_xi(order);
    v_eta = v_eta(order);

    % Hacia arriba: foco seguido de muestras con z creciente.
    upward = find(v_z > z_f + line_tol);
    selected = first_bracket([focal_abs_u; v_abs(upward)], target, ...
                             'vertical superior');
    if selected == 1
        error('quarter_pressure_points:UnexpectedFocusBracket', ...
              'El foco no puede ser el punto vertical superior.');
    end
    index = upward(selected-1);
    vertical_up = [v_rho(index), v_z(index), v_abs(index), ...
                   v_patch(index), v_xi(index), v_eta(index)];

    % Hacia abajo: foco seguido de muestras con z decreciente.
    downward = flipud(find(v_z < z_f - line_tol));
    selected = first_bracket([focal_abs_u; v_abs(downward)], target, ...
                             'vertical inferior');
    if selected == 1
        error('quarter_pressure_points:UnexpectedFocusBracket', ...
              'El foco no puede ser el punto vertical inferior.');
    end
    index = downward(selected-1);
    vertical_down = [v_rho(index), v_z(index), v_abs(index), ...
                     v_patch(index), v_xi(index), v_eta(index)];

    fprintf('Puntos cercanos a |u(0,z_f)|/4 (objetivo %.16e):\n', target);
    fprintf('  Horizontal: rho=%.16e, z=%.16e, |u|=%.16e\n', ...
            horizontal(1), horizontal(2), horizontal(3));
    fprintf('  Vertical superior: rho=%.16e, z=%.16e, |u|=%.16e\n', ...
            vertical_up(1), vertical_up(2), vertical_up(3));
    fprintf('  Vertical inferior: rho=%.16e, z=%.16e, |u|=%.16e\n', ...
            vertical_down(1), vertical_down(2), vertical_down(3));
end

function selected = first_bracket(values, target, direction)
%FIRST_BRACKET Primer par consecutivo que encierra target; extremo mas cerca.
    for k = 1:numel(values)-1
        first_distance = values(k)-target;
        second_distance = values(k+1)-target;
        if first_distance == 0
            selected = k;
            return
        end
        if first_distance*second_distance < 0 || second_distance == 0
            if abs(first_distance) <= abs(second_distance)
                selected = k;
            else
                selected = k+1;
            end
            return
        end
    end
    error('quarter_pressure_points:NoBracket', ...
          ['No se encontro un par consecutivo que encierre el objetivo ' ...
           'en la direccion %s. Aumente el muestreo o revise el dominio.'], ...
          direction);
end
