function [last_max, preim, XYQ] = mainHelmholtzMultipatch_3Daxial(params, output_dir)
%MAINHELMHOLTZMULTIPATCH_3DAXIAL Resuelve Helmholtz 3D axial con k variable.
%   [LAST_MAX, PREIM, XYQ] = MAINHELMHOLTZMULTIPATCH_3DAXIAL(PARAMS,
%   OUTPUT_DIR) ejecuta el modelo por capas y guarda figuras y resultados
%   intermedios en OUTPUT_DIR.
%
%   LAST_MAX contiene [z_foco, |u|_foco], PREIM las preimagenes de la linea
%   focal y XYQ la salida de evaluaQ: columnas [X, Y, Q].

required_fields = {'r', 'a', 'c', 'fr', 'zc', 'C', 'z1', 'z2', 'z3', ...
                   'alpha', 'Q', 'n1', 'm1', 'n3', 'tissue_n', 'tissue_m', ...
                   'spline_n', 'spline_order', 'npoints_axis', ...
                   'npoints_focal', 'npoints_surface', 'ang', 'q_npoints'};
validate_helmholtz_params(params, required_fields);
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

fprintf('------------------------------------------------------------------------\n');
fprintf('Solucion IgA de la ecuacion de Helmholtz 3D (simetria axial)\n');
fprintf('Lente curva y k variable por capas\n');
fprintf('------------------------------------------------------------------------\n');

% Datos geometricos y del modelo.
r = params.r;
a = params.a;
fr = params.fr;
zc = params.zc;
rl = sqrt(a^2 + zc^2);
zmin = zc - rl;
z1 = params.z1;
z2 = params.z2;
z3 = params.z3;
C = params.C;

% k en la frontera de Robin y aproximacion suave en el dominio.
kfunR = @(rho, z) kfuntramos(rho, z, r, fr, z1, z2, z3);
graficatejido_lentecurva(r, z1, z2, params.tissue_n, params.tissue_m, ...
                         fr, zmin, rl, output_dir);
[ccoef, alphacoef, rhocoef, knots] = ...
    ConstruyeksmoothShoenberg_lentecurva(params.spline_n, params.spline_order, ...
                                         z1, z2, z3, zmin, r, params.c);
kfun2 = @(rho, z) ksmoothShoenbergPto(rho, z, fr, params.spline_order, ...
                                      knots, ccoef, alphacoef);

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
    geometry, space, msh, C, r, alpha, kfunR, kfun2, ...
    ccoef, alphacoef, rhocoef, knots, output_dir, params.q_npoints);

% Postprocesamiento y salidas para la posterior integracion con Pennes.
[last_max, preim] = postprocess_helmholtz(u, space, geometry, params, output_dir);
end
