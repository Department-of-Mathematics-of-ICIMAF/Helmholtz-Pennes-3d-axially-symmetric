function XYQ = evaluaQ(knots, ccoef, alphacoef, rhocoef, u, space, geometry, output_dir, npoints)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Q se calcula a partir de la solucion 
% u de la ec de Helmholtz como
%  Q= alpha*|u|^2/(rho*c)
% alpha es el coeficiente de atencuacion, rho es la densidad y c la velocidad de propagacion 
% del ultrasonido. Para un tejido compuesto por 3 capas alpha, rho y c son funciones constantes
% sobre cada CAPA, de modo que en el general son constantes por tramos. 
% Qval=alpha*(real(uval).^2+imag(uval).^2)/(rho*c);

fprintf('--------------------------------------------------------------\n')
fprintf(' Calculos para graficar Q ...\n')


if nargin < 8 || isempty(output_dir)
    output_dir = pwd;
end
if nargin < 9 || isempty(npoints)
    npoints = 250;
end
pts = {linspace(0, 1, npoints), linspace(0, 1, npoints)};

% Evaluamos la funcion Q para graficarla parche por parche
fig351 = figure(351); clf(fig351); hold on, grid on
axis([-0.15 0.15 -0.02 0.15])
colorbar
axis tight, axis equal
ax=gca;
ax.FontSize=14;

XYQ=[];
for iptc = 1:space.npatch
    % Ciclo para cada parche
    [Zval, Fval]=sp_eval(u(space.gnum{iptc}),space.sp_patch{iptc},geometry(iptc),pts);
    [Xval, Yval]  = deal (squeeze(Fval(1,:,:)), squeeze(Fval(2,:,:)));
    % Evaluacion y graficacion de un parche de Q
    [cxy,alphaxy,rhoxy]=EvaluaShoenbergPto(Xval, Yval,knots,ccoef,alphacoef,rhocoef);
    deno=rhoxy.*cxy;
    coefQ=alphaxy./deno;
    Qval=coefQ.*(real(Zval).^2+imag(Zval).^2);
    %Qval=alphaxy.*(real(Zval).^2+imag(Zval).^2)./(rhoxy.*cxy); 
    surf (Xval, Yval,Qval,'EdgeColor','none') 
    XYQnew=[Xval(:),Yval(:),Qval(:)];
    XYQ=[XYQ;XYQnew];
end
save_figure(fig351, output_dir, 'Q_distribution');
end
