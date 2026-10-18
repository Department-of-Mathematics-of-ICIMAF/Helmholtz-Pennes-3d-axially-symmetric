%mainHelmholtzMultipatch_3Daxial
% VHM marzo 2026 / octubre 2026 (integracion de evaluaQ)
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
c=1568;
fr=1e6;
%zc=input('zc=');
zc=r/2;
rl=sqrt(a^2+zc^2);
aten = 8.55;
kval=2*pi*fr/c-1i*aten;
C=4.0e5;
% C=1;

% Densidad del tejido (constante en todo el dominio). Valor tipico para
% tejido phantom / musculo. Modifiquelo segun el tejido que se este
% simulando (agua: rho=1000, piel: 1200, grasa: 950, musculo: 1050).
rho = 1044;

% En el modelo con k CONSTANTE las magnitudes c, alpha y rho son valores
% constantes en todo el dominio. Estos valores se pasan a
% radiationmatricesHelmholtz_3Daxial y de ahi a evaluaQ para calcular
%   Q = alpha * |u|^2 / (rho * c)
% a partir de la solucion aproximada u de la ecuacion de Helmholtz.
c_const     = c;       % velocidad de propagacion del ultrasonido (m/s)
alpha_const = aten;    % coeficiente de atenuacion (Np/m)
rho_const   = rho;     % densidad del tejido (kg/m^3)

fprintf('Wave number k=%e\n',kval);
fprintf('r=%e\n',r);
fprintf('Velocidad c=%g, Atenuacion alpha=%g, Densidad rho=%g\n', c_const, alpha_const, rho_const);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calcular la parametrizacion Multiparche del dominio con 3 parches 

alpha = 0.7698;
Q = [0.0217,0.0385];
% Parametrizacion Multipatch VHM con nodos simples
[geometry,boundaries,interfaces,boundary_interfaces]=ParametMultipatch_mediodominio(r,a,rl,alpha,Q,zc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construir el espacio spline al que pertenece la solucion 
% 
%n1=40;m1=24;n3=20;
%n1=80;m1=48;n3=40;
n1=100; % N=308502 
m1 = n1*3;
n3 = n1*4; 
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

u=radiationmatricesHelmholtz_3Daxial(geometry,space,msh,kval,C,r,alpha,c_const,alpha_const,rho_const);

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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Z = [Zp1', Zp2'];
up = [abs(up1), abs(up2)];

ind = ultimo_maximo_local(up);
last_max = [Z(ind) up(ind)];

figure(600), hold on
plot(last_max(1),last_max(2),'ko','MarkerSize',10),hold on

ind_vector = ind-10:ind-1;
plot(Z(ind_vector), up(ind_vector), 'o', 'Color', 'green', 'MarkerFaceColor', 'green'), hold on

ind_vector = ind+1:ind+10;
plot(Z(ind_vector), up(ind_vector), 'o', 'Color', 'green', 'MarkerFaceColor', 'green')

fprintf('Ultimo maximo local: (%e,  %e)\n',last_max(1), last_max(2))
fprintf('Indice: %e\n',ind)
saveas(gcf,'multi.fig')

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

% Encontrar mnimo y mximo globales
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
% Agregar cada slido con la escala global
crear_solido_revolucion2(ax, Rhop1, Zp1, sol1,ang, global_clim);
crear_solido_revolucion2(ax, Rhop2, Zp2, sol2,ang, global_clim);
crear_solido_revolucion2(ax, Rhop3, Zp3, sol3,ang, global_clim);

% Mejorar iluminacin (opcional)
camlight(ax, 'left');
lighting(ax, 'gouraud');
material(ax, 'dull');

% Agregar barra de colores (usar la escala global automticamente)
colorbar(ax);
