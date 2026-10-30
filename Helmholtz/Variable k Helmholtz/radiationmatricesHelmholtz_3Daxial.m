function [u, XYQ] = radiationmatricesHelmholtz_3Daxial(geometry, space, msh, C, r, alpha, KfunR, Kfun2, ccoef, alphacoef, rhocoef, knots, output_dir, q_npoints)
if nargin < 13 || isempty(output_dir)
    output_dir = pwd;
end
if nargin < 14 || isempty(q_npoints)
    q_npoints = 250;
end
% Esta funcion calcula las matrices de masa y rigidez para el problema de 
% Helmholtz en un dominio 3D con simetria axial. Ademas calcula la matriz
% A y el lado derecho b del sistema lineal, cuya solucion son los
% coeficientes en la base B-spline (multipatch y producto tensorial bicubica)
% de la solucion aproximada de la ecuacion de Helmholtz.
%
% Entrada
% space - objecto que representa el espacio multipatch al que petenece la solucion aproximada 
% msh- malla multipatch
% kval- numero de onda
% C - constante de Dirichlet
% r- radio del dominio 
% alpha- algulo que define la parametrizacion. En particular (r*cos(alpha),r*sin(alpha))
%        es el punto sobre la frontera que divide el segmento de arco de
%        circulo del parche y del parche 3
% Salida
% u-coeficientes en la base B-spline (multipatch y producto tensorial bicubica)
% de la solucion aproximada de la ecuacion de Helmholtz.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Chequeo 1: usamos la funcion op_f_v_tp (usualmente pensada para calcular una 
% integral doble) para calcular una integral de linea, que es la longitud de 
% un arco de circunferencia. El resultado se basa en que la suma de las
% funciones basicas es 1.
% gfun=@(rho,z)ones(length(z),1);
% v4= op_f_v_tp (space.sp_patch{2}.boundary(4), msh.msh_patch{2}.boundary(4),gfun);
%                                                  
% fprintf('long arco exacta frontera 4=%e\n',r*(pi/2-alpha))                                                    
% fprintf('long arco aprox. frontera 4=%e\n',sum(v4))   
% fprintf('----------------------------------------------------\n')
% v5= op_f_v_tp (space.sp_patch{3}.boundary(2), msh.msh_patch{3}.boundary(2),gfun);
%                                                  
% fprintf('long arco exacta frontera 5=%e\n',r*alpha)                                                    
% fprintf('long arco aprox. frontera 5=%e\n',sum(v5))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('----------------------------------------------------------------\n')
fprintf('Calculando las matrices IgA \n')
tic
% Calculo de la matriz de rigidez de orden Ndof x Ndof. Se almacena de forma sparse
frho = @(rho,z) rho;
stiff_mat = op_gradu_gradv_mp (space, space, msh, frho);


% Calculo de la matriz de masa de orden Ndof x Ndof. Se almacena de forma sparse
k2frho = @(rho,z) Kfun2(rho,z).*rho;
mass_mat = op_u_v_mp (space,space, msh, k2frho);

% Matriz del sist. lineal antes de considerar las condiciones de frontera 
mat = stiff_mat -mass_mat;  % mass_mat ya incluye el signo negativo y k^2

Ndof = size(mat,1);
% Las fronteras de cada parche se numeran localmente en el orden:
% frontera 1 = Oeste
% frontera 2 = Este
% frontera 3 = Sur
% frontera 4 = Norte
% 
% El orden de la numeracion global de las fronteras es por parches y dentro
% de cada parche siguiendo el orden de numeracion local. Por eso en la
% numeracion global tenemos que
%   Frontera 1 global: frontera Oeste del parche 1- condicion de Neuman
%   Frontera 2 global: frontera Sur del parche 1- condicion de Dirichlet (arco inferior, transductor)
%   Frontera 3 global: frontera Oeste del parche 2- condicion de Neuman
%   Frontera 4 global: frontera Norte del parche 2- condicion de Robin (arco superior)
%   Frontera 5 global: frontera Este del parche 3- condicion de Robin (arco derecho)
%   Frontera 6 global: frontera Sur del parche 3- condicion de Dirichlet (arco inferior)

dofs_dirich2 = space.boundary.dofs(space.boundary.gnum{2}); % Indices de las funciones basicas que no se anulan sobre 
                                                            % la frontera de Dirichlet (arco inferior)
dofs_robin4 = space.boundary.dofs(space.boundary.gnum{4}); % Indices de las funciones basicas que no se anulan sobre la 
                                                           % frontera de Robin superior
dofs_robin5 = space.boundary.dofs(space.boundary.gnum{5}); % Indices de las funciones basicas que no se anulan sobre la 
                                                           % frontera de Robin derecha


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                           
% Matrices de masa en fronteras de Robin 
frhoR = @(rho,z) KfunR(rho,z).*rho;
%robin4_mat = op_u_v_tp (space.boundary.sp_patch{4}, space.boundary.sp_patch{4}, msh.boundary.msh_patch{4},frho);
% El resultado anterior es el mismo que el siguiente
 robin4_mat= op_u_v_tp(space.sp_patch{2}.boundary(4), space.sp_patch{2}.boundary(4), msh.msh_patch{2}.boundary(4),frhoR);
% o sea integramos sobre la frontera 4 del parche 2

%robin5_mat = op_u_v_tp (space.boundary.sp_patch{5}, space.boundary.sp_patch{5}, msh.boundary.msh_patch{5},frho);
% El resultado anterior es el mismo que el siguiente
robin5_mat= op_u_v_tp(space.sp_patch{3}.boundary(2), space.sp_patch{3}.boundary(2), msh.msh_patch{3}.boundary(2),frhoR);
% o sea integramos sobre la frontera 2 del parche 3

% Agregar contribucion de Robin: i * k * (matriz de masa en frontera)
mat(dofs_robin4, dofs_robin4) = mat(dofs_robin4, dofs_robin4) + 1i  * robin4_mat;
mat(dofs_robin5, dofs_robin5) = mat(dofs_robin5, dofs_robin5) + 1i  * robin5_mat;
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Condicion de Dirichlet
u(dofs_dirich2) = C;

% Lado derecho (cero en este problema)
rhs = zeros(Ndof, 1);

% Ajustar lado derecho por condición de Dirichlet
int_dofs = setdiff (1:space.ndof, dofs_dirich2);
rhs(int_dofs) = rhs(int_dofs) - mat(int_dofs, dofs_dirich2) * u(dofs_dirich2)';

fprintf('-----------------------------------------------------------------\n')
fprintf('Resolviendo el sistema de ecuaciones lineales\n')
% Resolver sistema para grados de libertad interiores
u(int_dofs) = mat(int_dofs, int_dofs) \ rhs(int_dofs);
fprintf('-----------------------------------------------------------------\n')

XYQ = evaluaQ(knots, ccoef, alphacoef, rhocoef, u, space, geometry, ...
              output_dir, q_npoints);
save(fullfile(output_dir, 'Q_values.mat'), 'XYQ');

% Q queda guardado automaticamente en la carpeta de esta corrida.

