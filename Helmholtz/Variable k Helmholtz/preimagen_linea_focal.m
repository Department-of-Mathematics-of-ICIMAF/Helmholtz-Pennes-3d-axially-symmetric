function preim = preimagen_linea_focal(geometry, z_f, npoints)
% PREIMAGEN_LINEA_FOCAL  Calcula las preimagenes, en el dominio parametrico
% [0,1]^2 de cada parche, de una muestra densa de puntos sobre la linea
% horizontal z = z_f del dominio fisico (plano meridiano (rho, z)).
%
% USO:
%   preim = preimagen_linea_focal(geometry, z_f)
%   preim = preimagen_linea_focal(geometry, z_f, npoints)
%
% ENTRADA:
%   geometry : arreglo de estructuras con la parametrizacion de cada parche
%              (salida de mp_geo_load, tal como la devuelve
%              ParametMultipatch_mediodominio).
%   z_f      : altura z de la linea horizontal. En el problema corresponde
%              a la altura del foco, z_f = Z(ind).
%   npoints  : numero aproximado de puntos de la muestra (por defecto 1000).
%
% SALIDA:
%   preim : arreglo de estructuras con un elemento por parche. Para el
%           parche iptc:
%     preim(iptc).xi  : coordenada xi de las preimagenes (vector columna)
%     preim(iptc).eta : coordenada eta de las preimagenes
%     preim(iptc).rho : coordenada rho fisica de cada punto, F(xi,eta)(1)
%     preim(iptc).z   : coordenada z fisica de cada punto (debe ser ~ z_f)
%     preim(iptc).n   : cantidad de puntos de la muestra que cayeron en
%                       ese parche
%
% METODO:
%   Se recorre la linea con un paso pequeno en rho desde rho = 0 hasta
%   salir del dominio y, para cada punto (rho, z_f), se resuelve el sistema
%   no lineal
%        F_i(xi, eta) = (rho, z_f)
%   con el metodo de Newton-Raphson (con proyeccion al cuadrado [0,1]^2)
%   usando como semilla la solucion convergida del punto anterior en ese
%   parche (estrategia de continuacion, como en los codigos new_raph de la
%   tesis). Un punto se acepta en el parche i solo si el residuo
%   ||F_i(xi,eta) - (rho, z_f)|| es menor que una tolerancia, lo que
%   garantiza que la preimagen pertenece a ese parche. Cada punto fisico
%   se asigna solamente al primer parche que lo reclama, de modo que los
%   perfiles no se duplican sobre las interfaces entre parches.
%   El barrido se hace en dos pasadas: una gruesa que busca, para cada
%   parche activo, una semilla que converja (solo en una rejilla de rhos),
%   y una fina con continuacion que genera la muestra densa. Los parches
%   cuyo rango de z no contiene z_f se descartan de antemano.
%
% Nota: space no se necesita como argumento; la inversion solo usa el
%       mapeo geometrico F de cada parche.

% VHM octubre 2026

if nargin < 3 || isempty(npoints)
    npoints = 1000;
end

npatch   = numel (geometry);
tol_nr   = 1e-12;    % tolerancia del paso de Newton
maxit    = 60;       % maximo de iteraciones de Newton por punto
tol_res  = 1e-9;     % tolerancia del residuo para aceptar la preimagen
n_fallo  = 12;       % fallos consecutivos para dar por terminada la linea
n_arista = 257;      % puntos por arista para estimar rangos de rho y z
n_grueso = 80;       % puntos de la pasada gruesa (busqueda de semillas)
tol_z    = 1e-6;     % holgura en el filtro del rango de z

%--------------------------------------------------------------------------
% 1. Rangos de rho y z de cada parche: se evaluan las aristas parametricas
%    (xi=0, xi=1, eta=0, eta=1). De ahi salen rho_max del dominio y el
%    filtro de parches activos (z_f dentro del rango de z del parche).
tt = linspace (0, 1, n_arista);
rho_max = 0;
zmin_p = inf (1, npatch);
zmax_p = -inf (1, npatch);
for iptc = 1:npatch
    for arista = 1:4
        switch arista
            case 1, pts = [tt; zeros(size(tt))];   % eta = 0
            case 2, pts = [tt; ones(size(tt))];    % eta = 1
            case 3, pts = [zeros(size(tt)); tt];   % xi = 0
            case 4, pts = [ones(size(tt)); tt];    % xi = 1
        end
        e = geometry(iptc).map (pts);
        rho_max = max (rho_max, max (e(1,:)));
        zmin_p(iptc) = min (zmin_p(iptc), min (e(2,:)));
        zmax_p(iptc) = max (zmax_p(iptc), max (e(2,:)));
    end
end
if rho_max <= 0
    error ('preimagen_linea_focal: no se pudo estimar la extension de la linea');
end

activo = (z_f >= zmin_p - tol_z) & (z_f <= zmax_p + tol_z);

%--------------------------------------------------------------------------
% 2. Semillas iniciales y pasada gruesa: para cada parche activo se busca
%    una semilla que converja en algun punto de una rejilla gruesa de rhos.
%    A partir de ahi la continuacion hace el resto.
seeds = [ 0.5  0.10 0.90 0.5  0.5  0.0  1.0  0.0  1.0  0.02 0.98; ...
          0.5  0.5  0.5  0.10 0.90 0.0  0.0  1.0  1.0  0.5  0.5 ];

ders = cell (1, npatch);
for iptc = 1:npatch
    ders{iptc} = nrbderiv (geometry(iptc).nurbs);
end

rhos_grueso = linspace (0, rho_max, n_grueso);
seed_prev = cell (1, npatch);   % ultima preimagen convergida de cada parche
for iptc = 1:npatch
    if ~activo(iptc)
        continue
    end
    for rg = rhos_grueso
        target = [rg; z_f];
        for s = 1:size (seeds, 2)
            [x_e, res] = newton_invierte (geometry(iptc), ders{iptc}, ...
                               target, seeds(:,s), tol_nr, maxit);
            if res < tol_res
                seed_prev{iptc} = x_e;
                break
            end
        end
        if ~isempty (seed_prev{iptc})
            break
        end
    end
end

%--------------------------------------------------------------------------
% 3. Pasada fina: barrido global en rho con continuacion. Cada punto se
%    asigna solo al primer parche que lo reclama.
drho = rho_max / (3*npoints);   % sobremuestreo ~3x, se recortan los fallos

xi   = cell (1, npatch);
eta  = cell (1, npatch);
rhov = cell (1, npatch);
zv   = cell (1, npatch);

rho   = 0;
nfail = 0;
while rho <= rho_max && nfail < n_fallo
    target = [rho; z_f];
    reclamado = false;

    for iptc = 1:npatch
        if ~activo(iptc) || isempty (seed_prev{iptc})
            continue
        end
        [x_e, res] = newton_invierte (geometry(iptc), ders{iptc}, ...
                           target, seed_prev{iptc}, tol_nr, maxit);
        if res < tol_res
            seed_prev{iptc} = x_e;
            if ~reclamado
                xi{iptc}  (end+1,1) = x_e(1);
                eta{iptc} (end+1,1) = x_e(2);
                pto = geometry(iptc).map (x_e);
                rhov{iptc}(end+1,1) = pto(1);
                zv{iptc}  (end+1,1) = pto(2);
                reclamado = true;
            end
        end
    end

    if reclamado
        nfail = 0;
    else
        nfail = nfail + 1;
    end
    rho = rho + drho;
end

%--------------------------------------------------------------------------
% 4. Salida
for iptc = 1:npatch
    preim(iptc).xi  = xi{iptc};
    preim(iptc).eta = eta{iptc};
    preim(iptc).rho = rhov{iptc};
    preim(iptc).z   = zv{iptc};
    preim(iptc).n   = numel (xi{iptc});
end

fprintf ('-----------------------------------------------------------------\n')
fprintf ('Preimagenes sobre la linea z = z_f = %e:\n', z_f)
for iptc = 1:npatch
    if preim(iptc).n > 0
        fprintf ('  parche %i: %i puntos, rho en [%e, %e]\n', iptc, ...
                 preim(iptc).n, min (preim(iptc).rho), max (preim(iptc).rho));
    else
        fprintf ('  parche %i: sin puntos en la linea\n', iptc);
    end
end

end

%==========================================================================
function [x_e, res] = newton_invierte (geom, der, target, seed, tol_nr, maxit)
% NEWTON_INVIERTE  Invierte el mapeo del parche: resuelve F(xi,eta) = target
% por Newton-Raphson con proyeccion al cuadrado [0,1]^2. Devuelve la
% aproximacion final y el residuo ||F(xi,eta) - target||.

    x_e = min (max (seed(:), 0), 1);
    for it = 1:maxit
        [Fv, Jm] = geo_nurbs (geom.nurbs, der, [], x_e, 1, 2);
        Fv = Fv - target;
        Jm = reshape (Jm, 2, 2);
        if rcond (Jm) < 1e-14
            break
        end
        delta = -Jm \ Fv;
        x_e = min (max (x_e + delta, 0), 1);
        if norm (delta) < tol_nr
            break
        end
    end
    res = norm (geom.map (x_e) - target);

end
