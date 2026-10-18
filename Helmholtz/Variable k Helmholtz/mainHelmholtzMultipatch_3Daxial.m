%mainHelmholtzMultipatch_3Daxial
% VHM marzo 2026
% programa para resolver la ecuacion de Helmholtz en un dominio 3D con 
% simetria axial, lente curva y k constante.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('------------------------------------------------------------------------\n')
fprintf('Solucion IgA de la ecuacion de Helmholtz 3D (simetria axial) en un problema\n')
fprintf('con lente curva \n')
fprintf('-------------------------------------------------------------------------\n')

% Datos del problema
r=0.133;
a=0.01;
c=1500;
fr=1e6;
%zc=input('zc=');
zc=r/2;
rl=sqrt(a^2+zc^2);
zmin = zc - rl;

%%%% k variable %%%%%%%%%%%%%%%%%%%%%%%%
fprintf('k es una funcion por capas\n') 
kfunR = @(rho,z)kfuntramos(rho,z,r,fr); % Esta funcion se utiliza para evaluar
%       exclusivamente la funcion k(x,y) restringida a la frontera del
%       semicirculo donde se impone la condicion de Robin

z1=0.010; z2=0.032; z3=0.0675;
% Grafico del tejid
graficatejido_lentecurva(r,z1,z2,1000,500,fr,zmin, rl);  
[ccoef,alphacoef,rhocoef,knots]=ConstruyeksmoothShoenberg_lentecurva(50,3,z1,z2,z3,zmin,r);

kfun2=@(rho,z) ksmoothShoenbergPto(rho,z,fr,3,knots,ccoef,alphacoef);
%kfunR = @(rho,z) sqrt(kfun2(rho,z));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C=5.0e5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calcular la parametrizacion Multiparche del dominio con 3 parches 

alpha = 0.7698;
Q = [0.0217,0.0385];
% Parametrizacion Multiparche VHM con nodos simples
[geometry,boundaries,interfaces,boundary_interfaces]=ParametMultipatch_mediodominio(r,a,rl,alpha,Q,zc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construir el espacio spline al que pertenece la solucion 
% 
%n1=40;m1=24;n3=20;
%n1=80;m1=48;n3=40;
%n1=100; m1=300;  n3=400; % N=189201
n1 = 15; 
m1 = 3*n1;
n3 = 4*n1;
%n1=200; m1=600;  n3=800; % N=758401
%n1=50; m1=500;  n3=500; 
fprintf('-----------------------------------------------------------------\n')
fprintf('La dimension del espacio de aproximacion se calcula en terminos de\n')
fprintf('n1=%i,m1=%i,n3=%i\n',n1,m1,n3)

[space,msh]=approxspace_multipatch(geometry,boundaries,interfaces,boundary_interfaces,n1,m1,n3);
Ndof = space.ndof

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculo de las matrices IgA

% Area exacta de Omega
% beta=atan(a/zc);
% Aexacta=pi*r^2/4+beta*rl^2/2-a*zc/2;
% fprintf('Area exacta de Omega=%5e\n', Aexacta);

u=radiationmatricesHelmholtz_3Daxial(geometry,space,msh,C,r,alpha, kfunR, kfun2, ccoef, alphacoef, rhocoef, knots);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Postprocesamiento grafico\n')
% Post procesamiento grafico
npoints=1000;
zp=linspace(0, 1, npoints);
pts = {0, zp};

% Evaluar los 2 parches de la solucion aproximada a lo largo del eje z
[up1, Fzp1] = sp_eval (u(space.gnum{1}), space.sp_patch{1}, geometry(1), pts);
Zp1=deal (squeeze(Fzp1(2,:,:)));

[up2, Fzp2] = sp_eval (u(space.gnum{2}), space.sp_patch{2}, geometry(2), pts);
Zp2=deal (squeeze(Fzp2(2,:,:)));

figure(400),grid on, hold on
plot(Zp1,real(up1),'LineWidth',2,'MarkerSize',6),grid on
plot(Zp2,real(up2),'LineWidth',2,'MarkerSize',6),grid on
title ('Parte Real de la solucion de Helmholtz sobre el eje z')

figure(500),grid on, hold on
plot(Zp1,imag(up1),'LineWidth',2,'MarkerSize',6),grid on
plot(Zp2,imag(up2),'LineWidth',2,'MarkerSize',6),grid on
title ('Parte imaginaria de la solucion de Helmholtz el eje z')

figure(600),grid on, hold on
plot(Zp1,abs(up1),'LineWidth',2,'MarkerSize',6),grid on
plot(Zp2,abs(up2),'LineWidth',2,'MarkerSize',6),grid on
title ('Valor absoluto de la solucion de Helmholtz sobre el eje z')

%saveas(gcf,'multi.fig')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculo del punto focal por la heuristica del ultimo maximo local
Z = [Zp1', Zp2'];
up = [abs(up1), abs(up2)];

ind = ultimo_maximo_local(up);
last_max = [Z(ind) up(ind)];

figure(600), hold on
plot(last_max(1),last_max(2),'ko','MarkerSize',10),hold on

fprintf('Ultimo maximo local: (%e,  %e)\n',last_max(1), last_max(2))
fprintf('Indice: %e\n',ind)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perfil de la solucion sobre la linea horizontal z = z_f (altura del foco)
% Se calculan las preimagenes de una muestra densa de puntos de esa linea
% en cada parche (dominio parametrico) y se evalua la solucion sobre ellas.
z_f = Z(ind);
fprintf('Altura del foco z_f = %e\n', z_f)

npoints_f = 1000;
preim = preimagen_linea_focal(geometry, z_f, npoints_f);

figure(700),grid on, hold on
figure(800),grid on, hold on
figure(900),grid on, hold on

for iptc = 1:numel(geometry)
    if preim(iptc).n == 0
        continue
    end
    % sp_eval trabaja con rejillas tensoriales, asi que sobre la curva
    % preimagen (xi(t), eta(t)) se evalua punto a punto
    upf = zeros(1, preim(iptc).n);
    for k = 1:preim(iptc).n
        upk = sp_eval (u(space.gnum{iptc}), space.sp_patch{iptc}, ...
                       geometry(iptc), {preim(iptc).xi(k), preim(iptc).eta(k)});
        upf(k) = upk;
    end
    rhof = preim(iptc).rho';   % coordenada x (rho) de los puntos de la linea

    figure(700)
    plot(rhof,real(upf),'LineWidth',2,'MarkerSize',6),grid on
    figure(800)
    plot(rhof,imag(upf),'LineWidth',2,'MarkerSize',6),grid on
    figure(900)
    plot(rhof,abs(upf),'LineWidth',2,'MarkerSize',6),grid on
end

figure(700)
xlabel('x'), title ('Parte Real de la solucion de Helmholtz sobre la linea z = z_f (altura del foco)')
figure(800)
xlabel('x'), title ('Parte imaginaria de la solucion de Helmholtz sobre la linea z = z_f (altura del foco)')
figure(900)
xlabel('x'), title ('Valor absoluto de la solucion de Helmholtz sobre la linea z = z_f (altura del foco)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rotacion
npoints=100;
pts=linspace(0, 1, npoints);
pts = {pts, pts};

[up1, Fzp1] = sp_eval (u(space.gnum{1}), space.sp_patch{1}, geometry(1), pts);
Rhop1=deal (squeeze(Fzp1(1,:,:)));
Zp1=deal (squeeze(Fzp1(2,:,:)));

[up2, Fzp2] = sp_eval (u(space.gnum{2}), space.sp_patch{2}, geometry(2), pts);
Rhop2=deal (squeeze(Fzp2(1,:,:)));
Zp2=deal (squeeze(Fzp2(2,:,:)));

[up3, Fzp3] = sp_eval (u(space.gnum{3}), space.sp_patch{3}, geometry(3), pts);
Rhop3=deal (squeeze(Fzp3(1,:,:)));
Zp3=deal (squeeze(Fzp3(2,:,:)));

sol1 = abs(up1);
sol2 = abs(up2);
sol3 = abs(up3);

% Encontrar mínimo y máximo globales
min_global = min([sol1(:); sol2(:); sol3(:)]);
max_global = max([sol1(:); sol2(:); sol3(:)]);
global_clim = [min_global, max_global];

% Crear figura y ejes
figure();
ax = axes;
view(3);
axis equal;
grid on;
xlabel('X');
ylabel('Y');
zlabel('Z');
title('abs(u_h)');
hold(ax, 'on');

ang = (3/2)*pi;
% Agregar cada sólido con la escala global
crear_solido_revolucion2(ax, Rhop1, Zp1, sol1,ang, global_clim);
crear_solido_revolucion2(ax, Rhop2, Zp2, sol2,ang, global_clim);
crear_solido_revolucion2(ax, Rhop3, Zp3, sol3,ang, global_clim);

% Mejorar iluminación (opcional)
camlight(ax, 'left');
lighting(ax, 'gouraud');
material(ax, 'dull');

% Agregar barra de colores (usará la escala global automáticamente)
colorbar(ax);

















