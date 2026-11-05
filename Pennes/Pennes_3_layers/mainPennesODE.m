function alpha = mainPennesODE(varargin)
% MainPennesODE
% Octubre 2026
%
% Programa para resolver la ecuacion de Pennes en un tejido de 3 capas.
%
% rho_d(x,y,z) c_d(x,y,z) (dT/dt) - k_d(x,y,z) Delta T = Q_e(x,y,z,t)
% (x,y,z) en Omega, t en [0,t_f]
% Q_e(x,y,z,t) = Q(x,y,z) I_s(t)
% Q(x,y,z) = mu_d(x,y,z) |u(x,y,z)|^2 / (rho_d(x,y,z) c_d(x,y,z))
%
% Condicion inicial: T(x,y,z,0) = T_b
% Frontera Dirichlet: T = T_b en Abajo, Derecha, Arriba.
%
% -------------------------------------------------------------------------
% SINTAXIS
% -------------------------------------------------------------------------
%   alpha = mainPennesODE()
%   alpha = mainPennesODE(R)
%   alpha = mainPennesODE(R, a)
%   alpha = mainPennesODE(R, a, rl)
%   alpha = mainPennesODE(R, a, rl, filenameQ)
%   alpha = mainPennesODE(R, a, rl, filenameQ, tf)
%   alpha = mainPennesODE(R, a, rl, filenameQ, tf, tstop)
%   alpha = mainPennesODE(R, a, rl, filenameQ, tf, tstop, show_figure)
%   alpha = mainPennesODE(R, a, rl, filenameQ, tf, tstop, show_figure, theta)
%
% Todos los argumentos son opcionales. Si un argumento no se pasa (o se
% pasa como []), se usa el valor por defecto o se solicita por teclado,
% segun se indica mas abajo.
%
% -------------------------------------------------------------------------
% PARAMETROS DE ENTRADA (en orden posicional, todos opcionales)
% -------------------------------------------------------------------------
%   R          : (double, escalar, > 0)
%                Radio de la semiesfera que define el dominio fisico Omega.
%                Por defecto: R = 0.133 [m].
%                Si no se pasa, se usa 0.133.
%
%   a          : (double, escalar, > 0)
%                Semiapertura de la lente acustica.
%                Por defecto: a = 0.01 [m].
%                Si no se pasa, se usa 0.01.
%
%   rl         : (double, escalar, > a)
%                Radio de curvatura de la lente.
%                Por defecto: rl = 0.0672 [m].
%                Si no se pasa, se usa 0.0672.
%
%   filenameQ  : (char | string)
%                Nombre del fichero .mat que contiene la variable XYQ con
%                los valores del campo de presion acustica Q(x,y). Se busca
%                en la carpeta ./inputs/. XYQ debe ser una matriz N x 3:
%                [x, y, Q(x,y)].
%                Por defecto: se solicita por teclado con input().
%
%   tf         : (double, escalar, > 0)
%                Tiempo final de simulacion [s].
%                Por defecto: se solicita por teclado con input().
%
%   tstop      : (double, escalar, 0 < tstop < tf)
%                Instante en el que se detiene la intensidad I_s(t) [s].
%                Por defecto: se solicita por teclado con input().
%
%   show_figure: (logical | 0 | 1)
%                Bandera para mostrar (1/true) u omitir (0/false) las
%                figuras de resultados:
%                  - Figura 500: T(z) a lo largo del eje Z para t = 5, 10,
%                    15, 25 y 30 s.
%                  - Figura 602: superficie 3D rotada de T en t = 30 s.
%                  - Figura 603: seccion 2D (plano XZ) de T en t = 30 s.
%                Por defecto: 0 (no mostrar figuras).
%
%   theta      : (double, escalar, en radianes)
%                Angulo de rotacion (alrededor del eje Z) usado para
%                reconstruir el solido completo a partir del cuarto de
%                circulo en las figuras 3D y 2D.
%                Por defecto: theta = 3*pi/2.
%
% -------------------------------------------------------------------------
% PARAMETRO DE SALIDA
% -------------------------------------------------------------------------
%   alpha      : (double, matriz Nt x Ndof)
%                Matriz de coeficientes de la aproximacion IgA de la
%                temperatura en cada instante de tiempo devuelto por el
%                solver de ODEs.
%
%                  - Filas    : Nt = length(tiempo), numero de instantes
%                               de tiempo resueltos por ode23t.
%                  - Columnas : Ndof = space.ndof, numero total de grados
%                               de libertad del espacio spline multipatch.
%
%                El elemento alpha(j, i) es el coeficiente alpha_i(t_j) en
%                la representacion
%
%                    Th(x,y,z,t_j) = sum_{i=1..Ndof} alpha(j,i) * Phi_i(x,y,z)
%
%                donde Phi_i son las B-splines del espacio approxspace.
%                Las columnas correspondientes a grados de libertad de
%                frontera de Dirichlet (drchlt_dofs) contienen el valor
%                constante T_b = 35 (condicion de frontera).
%
%                Para evaluar la temperatura en un patch concreto se usa
%                sp_eval(alpha(j, space.gnum{patch}), space.sp_patch{patch},
%                        geometry(patch), pts)
%
% -------------------------------------------------------------------------
% OTROS DATOS INTERNOS RELEVANTES
% -------------------------------------------------------------------------
%   T_b        : 35  (temperatura basal del tejido, condicion inicial y de
%                     frontera de Dirichlet).
%   z1, z2, z3 : 0.010, 0.032, 0.0675  (grosores de las tres capas).
%   lambda2    : 1   (factor de escala entre la C de Helmholtz y la de
%                     Pennes; ver bloque comentado en el codigo).
%
% -------------------------------------------------------------------------
% EJEMPLO DE USO
% -------------------------------------------------------------------------
%   % Usar todos los valores por defecto (pide filenameQ, tf y tstop)
%   alpha = mainPennesODE();
%
%   % Especificar todo salvo theta, sin figuras
%   alpha = mainPennesODE(0.133, 0.01, 0.0672, 'Q_data.mat', 30, 25, 0);
%
%   % Especificar todo y mostrar figuras con rotacion por defecto
%   alpha = mainPennesODE(0.133, 0.01, 0.0672, 'Q_data.mat', 30, 25, 1);
%
%   % Especificar todo, mostrar figuras y rotar 90 grados
%   alpha = mainPennesODE(0.133, 0.01, 0.0672, 'Q_data.mat', 30, 25, 1, pi/2);
%
% -------------------------------------------------------------------------
% DEPENDENCIAS (deben estar en ./utils o en el path de MATLAB)
% -------------------------------------------------------------------------
%   paramCuartCirculoLente, approxspace_multipatch,
%   ShoenbergparaPennes_lentecurva, k_dSmoothShoenbergPto,
%   rho_dxc_dSmoothShoenbergPto, matricesIgAbasicas3D,
%   select_drchlt_dofs, eval_q, eval_lado_derecho, sp_eval,
%   plot_3D_cylinder_1
%
% -------------------------------------------------------------------------

addpath(genpath('./utils'));

% -------------------------------------------------------------------------
% Parseo de argumentos
% Orden esperado (tras el posible axes):
%   args{1}=R, args{2}=a, args{3}=rl, args{4}=filenameQ,
%   args{5}=tf, args{6}=tstop, args{7}=show_figure, args{8}=theta
% -------------------------------------------------------------------------
narginchk(0,8);
args  = varargin;
nargs = numel(args);

R = 0.133;
if nargs > 0 && ~isempty(args{1}), R = args{1}; end

a = 0.01;
if nargs > 1 && ~isempty(args{2}), a = args{2}; end

rl = 0.0672;
if nargs > 2 && ~isempty(args{3}), rl = args{3}; end

if nargs > 3 && ~isempty(args{4})
    filenameQ = args{4};
    filename  = sprintf('inputs//%s', filenameQ);
else
    filenameQ = input('diga el nombre del fichero de valores de Q\n','s');
    filename  = sprintf('inputs//%s', filenameQ);
end
S = load(filename,'XYQ');
XYQ = S.XYQ;

% -------------------------------------------------------------------------
% Parametros de tiempo y visualizacion (leer ANTES de usar show_figure)
% -------------------------------------------------------------------------
if nargs > 4 && ~isempty(args{5})
    tf = args{5};
    fprintf('El problema se resuelve en [0,%f]\n', tf)
else
    tf = input('tf=');
end

if nargs > 5 && ~isempty(args{6})
    tstop = args{6};
    fprintf('La intensidad se detiene a los: %f sec.\n', tstop)
else
    fprintf('Diga en que segundo se detiene la intensidad\n')
    tstop = input('tstop=');
end

show_figure = 0;
if nargs > 6 && ~isempty(args{7})
    show_figure = args{7};
end

% -------------------------------------------------------------------------
% Geometria
% -------------------------------------------------------------------------
fprintf('------------------------------------------------------------------\n')
fprintf('Solucion IgA de la ecuacion de Pennes en un problema\n')
fprintf('donde se conoce la solucion exacta\n')
fprintf('------------------------------------------------------------------\n')

z_c = sqrt(rl^2 - a^2);

z1 = 0.010; z2 = 0.032; z3 = 0.0675;
zmin = z_c - rl;

[geometry,boundaries,interfaces,boundary_interfaces] = paramCuartCirculoLente(R, a, z_c, rl, show_figure);

% -------------------------------------------------------------------------
% Espacio spline
% -------------------------------------------------------------------------
fprintf('------------------------------------------------------------------\n')

n1 = 50;
m1=3*n1; n3=4*n1;

[space,msh]=approxspace_multipatch(geometry,boundaries,interfaces,boundary_interfaces,n1,m1,n3);
Ndof = space.ndof;

% -------------------------------------------------------------------------
% Matrices IgA
% -------------------------------------------------------------------------
fprintf('----------------------------------------------------\n')
fprintf('Calculo de las matrices IgA basicas\n')
fprintf('orden de las matrices=%i\n',Ndof)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construir el aproximante de Shoenberg para las funciones constantes por
% capas que definen la densidad rho, el calor especifico c_s y la
% conductividad termica k_s
[rho_dcoef,c_dcoef,k_dcoef,knots]=ShoenbergparaPennes_lentecurva(50,3,z1,z2,z3,zmin,R,show_figure);

%
% Calcular las matrices IgA de masa y rigidez del espacio de aproximacion multipatch
fprintf('----------------------------------------------------------------\n')
fprintf('Calculo de las matrices IgA multipatch de masa y rigidez\n')
fprintf('----------------------------------------------------------------\n')
tic
k_dfun       = @(rho,z) k_dSmoothShoenbergPto(rho,z,3,knots,k_dcoef);
rho_dxc_dfun = @(rho,z) rho_dxc_dSmoothShoenbergPto(rho,z,3,knots,rho_dcoef,c_dcoef);

[M,K] = matricesIgAbasicas3D(space,msh,rho_dxc_dfun,k_dfun);
toc

%% Reducir las matrices M y K a los bloques correspondientes a los nodos que no estan
%% relacionados con la frontera de Dirichlet, que son los coeficientes que se tienen 
%% que calcular.
[drchlt_dofs, free_dofs] = select_drchlt_dofs(space);

M_0 = M(free_dofs, free_dofs);
K_0 = K(free_dofs, free_dofs);
K_d = K(free_dofs, drchlt_dofs);

%% Calculo de la condicion inicial
T_b = 35;
% T_0 = @(rho,z) T_b;
alpha_0_inicial = T_b*ones(length(free_dofs),1);
n_d=length(drchlt_dofs);
alpha_d = T_b*ones(n_d,1);

%% Calculo del vector del lado derecho

%%% La escala lambda2 es el cociente entre
%%% la C que viene de Helmotz dividido por el nuevo C que se quiere para
%%% Pennes
% C_Helmotz = 1e6;
% C_Pennes = 4e5;
% lambda2 = (C_Pennes/C_Helmotz)^2;

lambda2 = 1;
Qval = scatteredInterpolant(XYQ(:,1),XYQ(:,2),lambda2*XYQ(:,3),'nearest');

q = eval_q(space, msh, free_dofs, Qval);  % no depende del tiempo

%% Calculo del sistema de ODEs
fprintf('-------------------------------------------------------------\n')
fprintf('Calculando la solucion del sistema ODE\n')
ATol = 1.0e-12; 
RTol = 1.0e-10; 
opts = odeset('RelTol',RTol,'AbsTol',ATol,'Mass',M_0,'Jacobian',-K_0);

t_in = 0;
N= 100;
step = (tf- t_in)/N;

tic;
solver_name = 'ode23t';
[tiempo,alphaI0,S]= ode23t(@(t, alpha_0) eval_lado_derecho(t, alpha_0, tstop, K_0, K_d, alpha_d, q),0:step:tf, alpha_0_inicial, opts);

% solver_name = 'ode23s';
% [tiempo,alphaI0,S]= ode23s(@(t, alpha_0) eval_lado_derecho(t, alpha_0, tstop, K_0, K_d, alpha_d, q),0:step:tf, alpha_0_inicial, opts);

% solver_name = 'ode45';
% [tiempo,alphaI0,S]= ode45(@(t, alpha_0) eval_lado_derecho(t, alpha_0, tstop, K_0, K_d, alpha_d, q),0:step:tf, alpha_0_inicial, opts);

% solver_name = 'ode15s';
% [tiempo,alphaI0,S]= ode15s(@(t, alpha_0) eval_lado_derecho(t, alpha_0, tstop, K_0, K_d, alpha_d, q),0:step:tf, alpha_0_inicial, opts);

% solver_name = 'ode23';
% [tiempo,alphaI0,S]= ode23(@(t, alpha_0) eval_lado_derecho(t, alpha_0, tstop, K_0, K_d, alpha_d, q),0:step:tf, alpha_0_inicial, opts);

% solver_name = 'ode113';
% [tiempo,alphaI0,S]= ode113(@(t, alpha_0) eval_lado_derecho(t, alpha_0, tstop, K_0, K_d, alpha_d, q),0:step:tf, alpha_0_inicial, opts);

% solver_name = 'ode23tb';
% [tiempo,alphaI0,S]= ode23tb(@(t, alpha_0) eval_lado_derecho(t, alpha_0, tstop, K_0, K_d, alpha_d, q),0:step:tf, alpha_0_inicial, opts);
solver_time = toc;

fprintf('Metodo de solucion %s\n',solver_name)
fprintf('-------------------------------------------------------------\n')
fprintf('Results of ODE solver\n')
fprintf('Number of successful steps = %i\n',S(1))
fprintf('Number of failed attemps = %i\n',S(2))
fprintf('Number of evaluations of the right hand side F(t,y)= %i\n',S(3))
fprintf('Number of evaluations of the Jacobian dF(t,y)/dy =%i\n',S(4))
fprintf('Number of LU decompositions = %i\n',S(5))
fprintf('Number of solutions of linear systems = %i\n',S(6))
fprintf('time to solve de ODE systems = %e\n',solver_time)
fprintf('-------------------------------------------------------------\n')

% -------------------------------------------------------------------------
% Reconstruir el vector completo alpha(t)
% -------------------------------------------------------------------------
Nt=length(tiempo);
alpha = zeros(Nt, Ndof);
alpha(:, free_dofs) = alphaI0;
alpha(:, drchlt_dofs) = T_b;

if show_figure == 1
    % ---------------------------------------------------------------------
    % Carpeta de salida: outputs/a_<a>_rl_<rl>/
    % Nombres de figura fijos dentro de la carpeta (Tz.png, 3D.png, 2D.png)
    % ---------------------------------------------------------------------
    outdir = fullfile('outputs', sprintf('a_%g_rl_%g', a, rl));
    if ~exist(outdir,'dir')
        mkdir(outdir);
    end
    fprintf('Figuras se guardaran en: %s\n', outdir);
    
    %% ---- Solucion a lo largo del eje Z ----
    figure(500); clf; grid on; hold on
    npoints = 800;
    z_p     = linspace(0,1,npoints);
    ptsejez = {0, z_p};
    max_value_g = 0;

    tiempos_plot = [5 10 15 25 30];
    colores      = {'-m','-k','-r','-g','-b'};
    leyendas     = cell(1,numel(tiempos_plot));

    for k = 1:numel(tiempos_plot)
        tk = tiempos_plot(k);
        [~,jk] = min(abs(tiempo - tk));

        [T1, F1] = sp_eval(alpha(jk,space.gnum{1}), space.sp_patch{1},geometry(1), ptsejez);
        Z1 = squeeze(F1(2,:,:))';

        [T2, F2] = sp_eval(alpha(jk,space.gnum{2}), space.sp_patch{2},geometry(2), ptsejez);
        Z2 = squeeze(F2(2,:,:))';

        Zplot  = [Z1, Z2];
        Tplot  = [T1, T2];
        plot(Zplot, Tplot, colores{k}, 'LineWidth', 2, 'MarkerSize', 6);
        max_value_g = max(max(Tplot), max_value_g);
        leyendas{k} = sprintf('t=%d sec', tk);
    end

    plot([z1,z1],[35,max_value_g],'--c','LineWidth',2,'HandleVisibility','off')
    plot([z2,z2],[35,max_value_g],'--c','LineWidth',2,'HandleVisibility','off')

    legend(leyendas,'Location','best')
    xlabel('Temp.'); ylabel('Z');
    ax = gca; ax.FontSize = 12; ax.FontWeight = 'bold';
    hold off
    
    % Guardar figura 500 -> outputs/a_<a>_rl_<rl>/Tz.png
    fname500 = fullfile(outdir, 'T_vs_z.fig');
    exportgraphics(figure(500), fname500, 'Resolution', 300);
    fprintf('Guardada: %s\n', fname500);
    
    %% ---- Solucion en todo Omega ----
    theta = 3*pi/2;
    if nargs > 7 && ~isempty(args{8})
        theta = args{8};
    end
    
    figure(602); clf; hold on; grid on
    figure(603); clf; hold on; grid on
    
    npoints=100;
    pts = {linspace(0, 1, npoints), linspace(0, 1, npoints)};
    
    n_patch = space.npatch;
    
    min_value_g = inf;
    max_value_g = 0;
    
    t_eval = 30;
    [~,j]  = min(abs(tiempo - t_eval));  
    for patch = 1:n_patch
        [Taprox, F] = sp_eval (alpha(j,space.gnum{patch}), space.sp_patch{patch}, geometry(patch), pts);
        [Rho, Z]  = deal(squeeze(F(1,:,:)), squeeze(F(2,:,:)));
        
        %% Figura 3D con la solucion calculada en un tiempo fijo
        min_value_g = min(min(Taprox(:)), min_value_g);
        max_value_g = max(max(Taprox(:)), max_value_g);

        figure(602), hold on,grid on
        plot_3D_cylinder_1(Rho, Z, Taprox, theta, min_value_g, max_value_g, ...
            ['Temperatura para t = ', num2str(tiempo(j)), 'sec'],...
            z1,z2,z3);
        title(['Superficie rotada para el timepo t=', num2str(tiempo(j)), 'sec']);
        
        %%% Grafico las curvas de los bordes
        hold on
        x_rho = R*sin(0:0.001:pi/2);
        x_z = R*cos(0:0.001:pi/2);
        x_zero = zeros(1,length(x_rho));
        plot3(x_rho,x_zero,x_z, 'LineStyle', '--', 'LineWidth', 2, 'Color', 'w'); % 
        plot3(x_zero, -x_rho,x_z, 'LineStyle', '--', 'LineWidth', 2, 'Color', 'w'); % 
        hold off
        ax = gca;ax.FontSize = 12;ax.FontWeight = 'bold';
        
        %% Figura 2D sobre el Plano XZ
        figure(603), hold on,grid on
        surf(Rho, Z, Taprox,'EdgeColor', 'flat')
        surf(-Rho, Z, Taprox,'EdgeColor', 'flat')
        title(['Plano 2D para el timepo t=', num2str(tiempo(j)), 'sec']);
        ax = gca; ax.FontSize = 12; ax.FontWeight = 'bold';
        cb = colorbar;
        cb.Limits = [min_value_g, max_value_g];
        xlabel('X'); ylabel('Z');
        hold off
        
    end
end
end