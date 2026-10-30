function[geometry,boundaries,interfaces,boundary_interfaces]=ParametMultipatch_mediodominio(r,a,rl,alpha,Q,zc,output_dir)
if nargin < 7 || isempty(output_dir)
    output_dir = pwd;
end
 % VHM, marzo 2026
% Esta funcion construye una parametrizacion F(xi,eta) del dominio del problema 
% 3D con lente CURVA cuando se utiliza la simetria axial. Este dominio es 2D y esta 
% definido por 1/4 de circulo de radio r al cual se le agrega el sector de circulo
% que define la lente. 

r2=norm(Q);
line_right1 = abline([a, 0],Q);
line_right2 = abline(Q,[r*cos(alpha), r*sin(alpha)]);
line_right3 = abline([a, 0],[r, 0]);
line_left1 = abline([0, -(rl-zc)],[0, r2]);
line_left2 = abline([0, r2], [0, r]);

arc_top = nrbcirc(r,[0 0],alpha,pi/2);
aux=arc_top.coefs(:,1);
arc_top.coefs(:,1)=arc_top.coefs(:,3);
arc_top.coefs(:,3)=aux;

beta=acos(Q(1)/r2);
arc_middle = nrbcirc(r2,[0 0],beta,pi/2);
aux=arc_middle.coefs(:,1);
arc_middle.coefs(:,1)=arc_middle.coefs(:,3);
arc_middle.coefs(:,3)=aux;

arc_right = nrbcirc(r,[0 0],0,alpha);

arc_bottom = nrbcirc(rl,[0 zc],3*pi/2,3*pi/2+atan(a/zc));

figure(10), clf
 nrbctrlplot(line_right1), hold on
 nrbctrlplot(line_right2), hold on
 nrbctrlplot(line_right3), hold on
 nrbctrlplot(line_left1), hold on
 nrbctrlplot(line_left2), hold on
 nrbctrlplot(arc_top), hold on
 nrbctrlplot(arc_middle), hold on
 nrbctrlplot(arc_bottom), hold on
 nrbctrlplot(arc_right)


%----------------------------Bottom Patch-------------------------------
%Right
curve_br = line_right1;

%Top
curve_bt = arc_middle;
%Left
curve_bl = line_left1;
%Bottom
curve_bb = arc_bottom;
%Patch
Fmap_b = nrbcoons(curve_bb,curve_bt,curve_bl,curve_br);
%geometry_b = geo_load (Fmap_b);
figure(2), clf, hold on
nrbctrlplot(Fmap_b)
% Calculos para graficar el Jacobiano del mapeo
[ders, ders2] = nrbderiv (Fmap_b);
npts=299;
pts = {linspace(0, 1, npts), linspace(0, 1, npts)};

% Valor de la parametrizacion evaluada en los ptos pts: matriz de orden 2 x npts^2
 
[F,Jac] = geo_nurbs (Fmap_b, ders, ders2, pts, 1, 2); 

%Mean Ratio Jacobian- medida de la calidad de la parametrizacion
Zvector(:)=Jac(1,1,:).*Jac(2,2,:)-Jac(1,2,:).*Jac(2,1,:);
denom(:)=Jac(1,1,:).*Jac(1,1,:)+Jac(2,1,:).*Jac(2,1,:)+ Jac(1,2,:).*Jac(1,2,:)+Jac(2,2,:).*Jac(2,2,:); 
Jrvector=2*Zvector./denom;
ZJr=reshape(Jrvector,npts,npts);

norminfJac_b=max(max(max(Jac)));
MinMrJ_b =min(min(ZJr));

X=reshape(F(1,:),npts,npts);
Y=reshape(F(2,:),npts,npts);
%figure(115), hold on, plot(X,Y,'m.-',X',Y','m.-')

figure(5),clf, hold on, grid on, view(2)
colorbar
surf (X, Y, ZJr,'EdgeColor','none','Facecolor','interp'), axis tight, axis equal
title('Mean ratio Jacobian')
%---------------------------------Top Patch--------------------------------
%Right
curve_tr = line_right2;
%Top
curve_tt = arc_top;
%Left
curve_tl = line_left2;
%Bottom
curve_tb = arc_middle;
%Patch
Fmap_t = nrbcoons(curve_tb,curve_tt,curve_tl,curve_tr);
%geometry_t = geo_load (Fmap_t);

figure(2), hold on
nrbctrlplot(Fmap_t)
% Calculos para graficar el Jacobiano del mapeo
[ders, ders2] = nrbderiv (Fmap_t);
%npts=199;
pts = {linspace(0, 1, npts), linspace(0, 1, npts)};

% Valor de la parametrizacion evaluada en los ptos pts: matriz de orden 2 x npts^2
 
[F,Jac] = geo_nurbs (Fmap_t, ders, ders2, pts, 1, 2); 

%Mean Ratio Jacobian- medida de la calidad de la parametrizacion
Zvector(:)=Jac(1,1,:).*Jac(2,2,:)-Jac(1,2,:).*Jac(2,1,:);
denom(:)=Jac(1,1,:).*Jac(1,1,:)+Jac(2,1,:).*Jac(2,1,:)+ Jac(1,2,:).*Jac(1,2,:)+Jac(2,2,:).*Jac(2,2,:); 
Jrvector=2*Zvector./denom;
ZJr=reshape(Jrvector,npts,npts);

norminfJac_t=max(max(max(Jac)));
MinMrJ_t =min(min(ZJr));

X=reshape(F(1,:),npts,npts);
Y=reshape(F(2,:),npts,npts);
%figure(115), hold on, plot(X,Y,'m.-',X',Y','m.-')

figure(5),hold on, grid on, view(2)
colorbar
surf (X, Y, ZJr,'EdgeColor','none','Facecolor','interp'), axis tight, axis equal
title('Mean ratio Jacobian')


%---------------------------Right Patch-------------

%Right
curve_rr = arc_right;
%Top
curve_rt = line_right2;
%Left
curve_rl = line_right1;
%Bottom
curve_rb = line_right3;
%Patch
Fmap_r = nrbcoons(curve_rb,curve_rt,curve_rl,curve_rr);
%geometry_r = geo_load (Fmap_r);

figure(2), hold on
nrbctrlplot(Fmap_r)
ax=gca;
ax.FontSize=14;
% Calculos para graficar el Jacobiano del mapeo
[ders, ders2] = nrbderiv (Fmap_r);
%npts=199;
pts = {linspace(0, 1, npts), linspace(0, 1, npts)};

% Valor de la parametrizacion evaluada en los ptos pts: matriz de orden 2 x npts^2
 
[F,Jac] = geo_nurbs (Fmap_r, ders, ders2, pts, 1, 2); 

%Mean Ratio Jacobian- medida de la calidad de la parametrizacion
Zvector(:)=Jac(1,1,:).*Jac(2,2,:)-Jac(1,2,:).*Jac(2,1,:);
denom(:)=Jac(1,1,:).*Jac(1,1,:)+Jac(2,1,:).*Jac(2,1,:)+ Jac(1,2,:).*Jac(1,2,:)+Jac(2,2,:).*Jac(2,2,:); 
Jrvector=2*Zvector./denom;
ZJr=reshape(Jrvector,npts,npts);

norminfJac_r=max(max(max(Jac)));
MinMrJ_r =min(min(ZJr));

X=reshape(F(1,:),npts,npts);
Y=reshape(F(2,:),npts,npts);
%figure(115), hold on, plot(X,Y,'m.-',X',Y','m.-')

figure(5),hold on, grid on, view(2)
colorbar
surf (X, Y, ZJr,'EdgeColor','none','Facecolor','interp'), axis tight, axis equal
title('Mean ratio Jacobian')
%close;
%----------------------------------------
MinMrJ=min([MinMrJ_b, MinMrJ_t, MinMrJ_r])
norminfJac = max([norminfJac_t, norminfJac_b, norminfJac_r]);
Fmap = [Fmap_b, Fmap_t, Fmap_r];
%geometry = [geometry_b, geometry_t, geometry_r];
[geometry, boundaries,interfaces,~,boundary_interfaces]= mp_geo_load(Fmap);

save_figure(10, output_dir, 'domain_control_curves');
save_figure(2, output_dir, 'patch_control_points');
save_figure(5, output_dir, 'mean_ratio_jacobian');
