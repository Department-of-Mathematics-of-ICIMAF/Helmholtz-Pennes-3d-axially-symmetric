function [last_max, preim, XYQ] = mainHelmholtzMultipatch_3Daxial(params, output_dir)
%MAINHELMHOLTZMULTIPATCH_3DAXIAL Resuelve Helmholtz 3D axial con k constante.
%   [LAST_MAX, PREIM, XYQ] = MAINHELMHOLTZMULTIPATCH_3DAXIAL(PARAMS,
%   OUTPUT_DIR) ejecuta el modelo con los parametros de PARAMS y guarda
%   figuras y resultados intermedios en OUTPUT_DIR.
%
%   LAST_MAX contiene [z_foco, |u|_foco], PREIM las preimagenes de la linea
%   focal y XYQ la salida de evaluaQ: columnas [X, Y, Q].

required_fields = {'r', 'a', 'c', 'fr', 'zc', 'aten', 'C', 'rho', ...
                   'alpha', 'Q', 'n1', 'm1', 'n3', 'npoints_axis', ...
                   'npoints_focal', 'npoints_surface', 'ang', 'q_npoints'};
validate_helmholtz_params(params, required_fields);
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

fprintf('------------------------------------------------------------------------\n');
fprintf('Solucion IgA de la ecuacion de Helmholtz 3D (simetria axial)\n');
fprintf('Lente curva y k constante\n');
fprintf('------------------------------------------------------------------------\n');

% Datos del problema y magnitudes constantes en todo el dominio.
r = params.r;
a = params.a;
c = params.c;
fr = params.fr;
zc = params.zc;
rl = sqrt(a^2 + zc^2);
aten = params.aten;
rho = params.rho;
C = params.C;
kval = 2*pi*fr/c - 1i*aten;

fprintf('Wave number k=%e\n', kval);
fprintf('r=%e\n', r);
fprintf('Velocidad c=%g, Atenuacion alpha=%g, Densidad rho=%g\n', c, aten, rho);

% Parametrizacion multiparche del dominio.
alpha = params.alpha;
Q = params.Q;
[geometry, boundaries, interfaces, boundary_interfaces] = ...
    ParametMultipatch_mediodominio(r, a, rl, alpha, Q, zc, output_dir);

% Espacio spline de aproximacion.
n1 = params.n1;
m1 = params.m1;
n3 = params.n3;
fprintf('-----------------------------------------------------------------\n');
fprintf('Dimension del espacio: n1=%i, m1=%i, n3=%i\n', n1, m1, n3);
[space, msh] = approxspace_multipatch(geometry, boundaries, interfaces, ...
                                      boundary_interfaces, n1, m1, n3);
fprintf('Cantidad de grados de libertad: %d\n', space.ndof);

% Matrices IgA, solucion y evaluacion de Q.
[u, XYQ] = radiationmatricesHelmholtz_3Daxial( ...
    geometry, space, msh, kval, C, r, alpha, c, aten, rho, ...
    output_dir, params.q_npoints);

% Postprocesamiento y salidas para la posterior integracion con Pennes.
[last_max, preim] = postprocess_helmholtz(u, space, geometry, params, output_dir);
end
