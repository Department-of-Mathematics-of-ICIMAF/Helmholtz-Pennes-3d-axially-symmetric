function[ccoef,alphacoef,rhocoef,knots]=ConstruyeksmoothShoenberg_lentecurva(n,k,y1,y2,y3,ymin,r,c_water)

% Por definicion
% k(x,y)=2*pi*fr/c(x,y)-1i*alpha(x,y)
% donde c(x,y) es la velocidad del ultrasonido que es una funcion constante
% por tramos y alpha(x,y) es ela atenuacion que tambien es constante por tramos. 
% En este problema estamos considerando que 
% c(x,y)=1558, alpha(x,y)=24 para 0<= y <= 0.010,  zona de la piel
% c(x,y)=1478, alpha(x,y)=5.58 para 0.010 < y <= 0.022, zona de la grasa
% c(x,y)=1547, alpha(x,y)=12.7 para 0.022 < y <= 0.0355, zona del musculo
% c(x,y)=3476, alpha(x,y)=80.2 para 0.0355 < y <= r, zona del hueso
% Se supone que el punto (x,y) donde vamos a evaluar esta en el interior
% del semicirculo con centro en el origen y radio r
%
% Las funciones c(x,y) y alpha(x,y) en este problema solamente dependen de
% la variable y. En esta funcion c(x,y) y alpha(x,y) se aproximan cada una
% por un spline de orden k de Shoenberg que pertenece a un espacio de
% dimension n

% Entrada
% n,k- dimension y orden del espacio spline al que pertenece las funciones 
%    aproximantes de c(x,y) y alpha(x,y) (como funciones de y
%    exclusivamente):
%  Vc(y)=sum_{i=1}^n c(tstar_i)B_i^k(y) - spline que aproxima a c(y)
%  Valpha(y)=sum_{i=1}^n alpha(tstar_i)B_i^k(y) - spline que aproxima a alpha(y)
% donde tstar_i se calculan como el promedio de los nodos definidos en [0,r]
% y1,y2,y3 - valores de y en (0,r) que definen los tramos de las funciones
%           c(ax,y) y alpha(x,y)
% ymin- valor de y para el punto sobre la lente curvas con x=0
% r - radio del semicirculo ( 0< y1 < y2 < y3 < r)
% Salida
% ccoef=c(tstar) - vector de los coeficientes en la base de B-splines del aproximante Vc
% alphacoef= alpha(tstar)- - vector de los coeficientes en la base de B-splines del aproximante Valpha
% knots- nodos de los splines Vc y Valpha
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if nargin < 8 || isempty(c_water)
    c_water = 1500;
end

% Calculo del aproximante Vc
gfun=@(y)speed(y,r,ymin);
[ccoef,knots]=ShoenbergSpline(n,k,gfun,ymin,r);

% % % Evaluamos los B-splines en una secuencia fina de puntos en [a,b]
tau= linspace(ymin, r, 1000);
gtau=gfun(tau);
% figure(201), hold on
% plot(tau,gtau,'.-r')
% 
% Evaluamos el aproximante de Shoenberg
Bmat = spcol(knots,k,tau);
Vc=Bmat*ccoef;
% plot(tau,Vc,'-b')
% title('Aproximante de Shoenberg para velocidad c(0,y)'), axis tight

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculo del aproximante Valpha
gfun=@(y)atenuacion(y,r,ymin);

[alphacoef,knots]=ShoenbergSpline(n,k,gfun,ymin,r);

% % Evaluamos los B-splines en una secuencia fina de puntos en [a,b]
gtau=gfun(tau);
% figure(202), hold on
% plot(tau,gtau,'.-r')

% Evaluamos el aproximante de Shoenberg
Bmat = spcol(knots,k,tau);
Valpha=Bmat*alphacoef;
% plot(tau,Valpha,'-b')
% 
% title('Aproximante de Shoenberg para atenuacion \alpha(0,y)'), axis tight

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculo del aproximante Vrho
gfun=@(y)densidad(y,r,ymin);

[rhocoef,knots]=ShoenbergSpline(n,k,gfun,ymin,r);

% % Evaluamos los B-splines en una secuencia fina de puntos en [a,b]
gtau=gfun(tau);

% figure(203), hold on
% plot(tau,gtau,'.-r')

% Evaluamos el aproximante de Shoenberg
Bmat = spcol(knots,k,tau);
Vrho=Bmat*rhocoef;
% plot(tau,Vrho,'-b')
% title('Aproximante de Shoenberg para la densidad \rho(0,y)'), axis tight

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
function[cy]=speed(y,r,ymin)
ny=length(y);
cy=zeros(ny,1);
 for i=1:ny
    if ymin <= y(i) && y(i) < 0
          cy(i)=c_water; % water
       elseif 0 <= y(i) && y(i) <= y1 
          cy(i)=1558; % skin
       elseif y1 < y(i) && y(i) <= y2 
         cy(i)=1478; % fat
       elseif y2 < y(i) && y(i) <= y3 
         cy(i)=1547; % muscle
        elseif y3 < y(i) && y(i) <= r         
% %         cy(i)=3476; % bone
           cy(i)=1547; % muscle
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[alpha]=atenuacion(y,r,ymin)
ny=length(y);
alpha=zeros(ny,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modificacion: enero 2026
% En la version original y0=0. Aqui cambiamos pensando que la atenuacion 
% muy cerca de la superficie de la piel sea 0 y no 24
y0=0;
%y0=0.2*y1;
%y0=0.5*y1;
%y0=0.9*y1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 for i=1:ny
%    if ymin <= y(i) && y(i) < 0
    if ymin <= y(i) && y(i) < y0   
        alpha(i)=0; %water
       %elseif 0 <= y(i) &&  y(i) <= y1 
       elseif y0 <= y(i) &&  y(i) <= y1 
          alpha(i)=24; %skin
          %alpha(i)=17; %skin
       elseif y1 < y(i) && y(i) <= y2 
         alpha(i)=5.58; %fat
       elseif y2 < y(i) && y(i) <= y3 
         alpha(i)=12.7; %muscle
       elseif y3 < y(i) && y(i) <= r 
% %         alpha(i)=80.2; %bone
           alpha(i)=12.7; %muscle
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[rho]=densidad(y,r,ymin)
ny=length(y);
rho=zeros(ny,1);
 for i=1:ny
    if ymin <= y(i) && y(i) < 0
        rho(i)=1000; %water
       elseif 0 <= y(i) && y(i) <= y1 
          rho(i)=1200; %skin
       elseif y1 < y(i) && y(i) <= y2 
         rho(i)=950; %fat
       elseif y2 < y(i) && y(i) <= y3 
         rho(i)=1050; %muscle
        elseif y3 < y(i) && y(i) <= r 
% %         rho(i)=1975; %bone
           rho(i)=1050; %muscle
    end
  end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
