function [] = graficatejido_lentecurva(r, y1, y2, n, m, fr, zmin, rl, output_dir)
if nargin < 9 || isempty(output_dir)
    output_dir = pwd;
end
cp=1558; % velocidad propagacion de la onda en piel
cg=1478; %  "                          "    en grasa
cm=1547; %  "                          "    en musculo
alphag=5.58; %fat      
alphap=24; %skin
alpham=12.7; %muscle

XF=linspace(-r,r,m);
YF=linspace(zmin,r,n);

[MFX,MFY]=meshgrid(XF,YF);

% Y - son n y estan relacionadas con las filas
% X - son m y estan relacionadas con las columnas

C=zeros(n,m);
A=zeros(n,m);

for i=1:n
  if YF(i)<= y1
      C(i,1:m)=cp;
      A(i,1:m)=alphap;
  elseif YF(i)>y1 && YF(i)<= y2
      C(i,1:m)=cg;
      A(i,1:m)=alphag;
  else   
     C(i,1:m)=cm;  
     A(i,1:m)=alpham;
  end
 end

% Valores de la funcion k(x,y)=2*pi*f/c(x,y)-1i*, en una malla rectangular en 
% [-r,r]x[zmin,r]
K=2*pi*fr./C-1i*A;
npt1=size(MFX,1);
npt2=size(MFX,2);

COLOR=ones(npt1,npt2,3);
for ii=1:npt1
    for jj=1:npt2
       if MFX(ii,jj)>= 0
           if C(ii,jj)==cp
             COLOR(ii,jj,1)=0.9; COLOR(ii,jj,2)=0.45; COLOR(ii,jj,3)=0;  
           elseif C(ii,jj)==cg
             COLOR(ii,jj,1)=1; COLOR(ii,jj,2)=1; COLOR(ii,jj,3)=0.3;  
           elseif C(ii,jj)==cm
             COLOR(ii,jj,1)=1; COLOR(ii,jj,2)=0.22; COLOR(ii,jj,3)=0.22; 
           end  
           if MFX(ii,jj)^2+MFY(ii,jj)^2 > r^2
            COLOR(ii,jj,1)=1; COLOR(ii,jj,2)=1; COLOR(ii,jj,3)=1; 
           end
           if MFY(ii,jj) < 0 && MFX(ii,jj)^2+(MFY(ii,jj)-(zmin + rl))^2 > rl^2
               COLOR(ii,jj,1)=1; COLOR(ii,jj,2)=1; COLOR(ii,jj,3)=1; 
           end
       end
    end
end
fig23 = figure(23); clf(fig23); hold on 
%mesh(MFX,MFY,C,COLOR)
surf(MFX,MFY,C,COLOR,'Facecolor','interp','EdgeColor','none','FaceAlpha',0.5);
ax=gca;
ax.FontSize=16;

ax.XLim = [0 0.1400];
view(2)

fig25 = figure(25); clf(fig25); grid on, hold on 
surf(MFX,MFY,abs(K),COLOR,'Facecolor','interp','EdgeColor','none','FaceAlpha',0.5);
ax=gca;
ax.XLim = [0 0.1400];
ax.FontSize=16;

save_figure(fig23, output_dir, 'tissue_speed');
save_figure(fig25, output_dir, 'tissue_wave_number');
end
