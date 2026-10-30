function [u, XYQ] = radiationmatricesHelmholtz_3Daxial(geometry, space, msh, kval, C, r, alpha, c_const, alpha_const, rho_const, output_dir, q_npoints)
if nargin < 11 || isempty(output_dir)
    output_dir = pwd;
end
if nargin < 12 || isempty(q_npoints)
    q_npoints = 250;
end
% Esta funcion calcula las matrices de masa y rigidez para el problema de 
% Helmholtz en un dominio 3D con simetria axial. Ademas calcula la matriz
% A y el lado derecho b del sistema lineal, cuya solucion son los
% coeficientes en la base B-spline (multipatch y producto tensorial bicubica)
% de la solucion aproximada de la ecuacion de Helmholtz.
%
% Entrada
% geometry - estructura que describe el mapeo multipatch, parametro de
%            salida de mp_geo_load
% space    - objeto que representa el espacio multipatch al que pertenece
%            la solucion aproximada
% msh      - malla multipatch
% kval     - numero de onda (constante en todo el dominio)
% C        - constante de Dirichlet
% r        - radio del dominio
% alpha    - angulo que define la parametrizacion. En particular
%            (r*cos(alpha), r*sin(alpha)) es el punto sobre la frontera que
%            divide el segmento de arco de circulo del parche 2 y del parche 3
% c_const     - (OPCIONAL) valor constante de la velocidad de propagacion c
%               del ultrasonido en todo el dominio (m/s). Necesario para
%               evaluar Q.
% alpha_const - (OPCIONAL) valor constante del coeficiente de atenuacion
%               alpha en todo el dominio (Np/m). Necesario para evaluar Q.
% rho_const   - (OPCIONAL) valor constante de la densidad rho del tejido en
%               todo el dominio (kg/m^3). Necesario para evaluar Q.
%
% Salida
% u - coeficientes en la base B-spline (multipatch y producto tensorial
%     bicubica) de la solucion aproximada de la ecuacion de Helmholtz.


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
frho = @(rho,z)rho;
stiff_mat = op_gradu_gradv_mp (space, space, msh, frho);

% Calculo de la matriz de masa de orden Ndof x Ndof. Se almacena de forma sparse
mass_mat = op_u_v_mp (space,space, msh, frho);

% Matriz del sist. lineal antes de considerar las condiciones de frontera 
mat = stiff_mat -kval^2* mass_mat;  % mass_mat ya incluye el signo negativo y k^2

Ndof = size(mat,1);
fprintf('Cantidad de grados de libertad: %.3f \n', Ndof)
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
%robin4_mat = op_u_v_tp (space.boundary.sp_patch{4}, space.boundary.sp_patch{4}, msh.boundary.msh_patch{4},frho);
% El resultado anterior es el mismo que el siguiente
 robin4_mat= op_u_v_tp(space.sp_patch{2}.boundary(4), space.sp_patch{2}.boundary(4), msh.msh_patch{2}.boundary(4),frho);
% o sea integramos sobre la frontera 4 del parche 2

%robin5_mat = op_u_v_tp (space.boundary.sp_patch{5}, space.boundary.sp_patch{5}, msh.boundary.msh_patch{5},frho);
% El resultado anterior es el mismo que el siguiente
robin5_mat= op_u_v_tp(space.sp_patch{3}.boundary(2), space.sp_patch{3}.boundary(2), msh.msh_patch{3}.boundary(2),frho);
% o sea integramos sobre la frontera 2 del parche 3

% Agregar contribucion de Robin: i * k * (matriz de masa en frontera)
mat(dofs_robin4, dofs_robin4) = mat(dofs_robin4, dofs_robin4) + 1i * kval * robin4_mat;
mat(dofs_robin5, dofs_robin5) = mat(dofs_robin5, dofs_robin5) + 1i * kval * robin5_mat;
time = toc;
fprintf('Tiempo para armar las matrices: %.3f segundos\n', time)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Condicion de Dirichlet
u(dofs_dirich2) = C;

% Lado derecho (cero en este problema)
rhs = zeros(Ndof, 1);

% Ajustar lado derecho por condicion de Dirichlet
int_dofs = setdiff (1:space.ndof, dofs_dirich2);
rhs(int_dofs) = rhs(int_dofs) - mat(int_dofs, dofs_dirich2) * u(dofs_dirich2)';

fprintf('-----------------------------------------------------------------\n')
fprintf('Resolviendo el sistema de ecuaciones lineales\n')
% Resolver sistema para grados de libertad interiores
tic
u(int_dofs) = mat(int_dofs, int_dofs) \ rhs(int_dofs);
time = toc;
fprintf('Tiempo para resolver el sistema: %.3f segundos\n', time)
fprintf('-----------------------------------------------------------------\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculo de Q = alpha * |u|^2 / (rho * c)
% Solo se realiza si se han proporcionado los parametros constantes
% c_const, alpha_const y rho_const. Si no se proporcionan, se omite el
% calculo de Q (la funcion u ya esta calculada y se devuelve normalmente).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin >= 8 && ~isempty(c_const) && ~isempty(alpha_const) && ~isempty(rho_const)
    XYQ = evaluaQ(u, space, geometry, c_const, alpha_const, rho_const, ...
                  output_dir, q_npoints);
    save(fullfile(output_dir, 'Q_values.mat'), 'XYQ');
else
    XYQ = [];
    fprintf('No se calcularon valores de Q: faltan c_const, alpha_const o rho_const.\n');
end
