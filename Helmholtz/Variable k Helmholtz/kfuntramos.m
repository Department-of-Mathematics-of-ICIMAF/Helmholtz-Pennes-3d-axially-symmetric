function kxy = kfuntramos(x, y, r, fr, y1, y2, y3)
% Esta funcion evalua en un punto (x,y) la funcion k(x,y) que aparece en 
% la ecuacion de Helmholtz para el problema de radiacion en un dominio
% semicircular de radio r. 
% La kfuntramos se llama en RadiationMatrixGeoPDEKcteaTramos cuand vamos a 
% calcular la matriz de Robin y entra  como function handle de la funcion de 
% GeoDPE op_u_v_tp .
% Por definicion
% k(x,y)=2*pi*fr/c(x,y)
% donde c(x,y) es la velocidad del ultrasonido que es una funcion constante
% por tramos. En este problema estamos considerando que 
% c(x,y)=1558 para 0<= y <= 0.010,  zona de la piel
% c(x,y)=1478 para 0.010 < y <= 0.022, zona de la grasa
% c(x,y)=1547 para 0.022 < y <= 0.0355, zona del musculo
% c(x,y)=3476 para 0.0355 < y <= r, zona del hueso
% Se supone que el punto (x,y) donde vamos a evaluar esta en el interior
% del semicirculo con centro en el origen y radio r

if nargin < 7
    y1 = 0.010;
    y2 = 0.032;
    y3 = 0.0675;
end
cte=2*pi*fr;
[n,m]=size(x);
kxy=zeros(n,m);
for i=1:n
   for j=1:m 
%          cxy=1478; % fat
%          alpha=5.58; %fat      
      if  y(i,j) <= y1 && abs(x(i,j)) <=r
         cxy=1558; % skin
         alpha=24; %skin
      elseif y1 < y(i,j) && y(i,j) <= y2 && abs(x(i,j)) <=r
         cxy=1478; % fat
         alpha=5.58; %fat
      elseif y2 < y(i,j) && y(i,j) <= y3 && abs(x(i,j)) <=r
        cxy=1547; % muscle
        alpha=12.7; %muscle
        elseif y3 < y(i,j) && y(i,j) <= r && abs(x(i,j)) <=r
%          %cxy=3476; % bone
%          %alpha=80.2; %bone
         cxy=1547; % muscle
         alpha=12.7; %muscle
       end
      kxy(i,j)=cte/cxy-1i*alpha;
%      c=1568;   alpha=8.55; constantes tejido phantom
%      OJO PRUEBA de la condicion de Robin con la parte real de k
%      kxy(i,j)=cte/cxy;
   end
end
