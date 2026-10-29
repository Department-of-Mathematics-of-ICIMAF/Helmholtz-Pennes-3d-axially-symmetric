function XYQ = evaluaQ(u, space, geometry, c_const, alpha_const, rho_const, output_dir, npoints)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% evaluaQ.m  -  Caso k CONSTANTE (Homogeneous Helmholtz)
% VHM octubre 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Q se calcula a partir de la solucion u de la ecuacion de Helmholtz como
%
%                  alpha * |u|^2
%            Q = ------------------
%                    rho * c
%
% donde:
%   alpha  es el coeficiente de atenuacion del tejido,
%   rho    es la densidad del tejido,
%   c      es la velocidad de propagacion del ultrasonido,
%   |u|^2  = Re(u)^2 + Im(u)^2  es el cuadrado del modulo de la presion
%           acustica u (solucion de Helmholtz).
%
% En el modelo con k CONSTANTE las magnitudes c, alpha y rho son valores
% constantes en todo el dominio. Por ello, en esta version simplificada
% los coeficientes c_const, alpha_const y rho_const se evaluan a traves del
% helper EvaluaConstantePto, que simplemente devuelve matrices llenas con el
% valor constante correspondiente. Esto es el analogo, en el caso k
% constante, del uso de EvaluaShoenbergPto (splines de Shoenberg) en el caso
% k variable. Mantener la misma firma de salida (cxy, alphaxy, rhoxy)
% permite que este evaluaQ.m sea casi identico al de la carpeta
% "Variable k Helmholtz".
%
% Entrada:
%   u            - coeficientes en la base B-spline (multipatch y producto
%                  tensorial bicubica) de la solucion aproximada de la
%                  ecuacion de Helmholtz, salida de
%                  radiationmatricesHelmholtz_3Daxial.m
%   space        - objeto que representa el espacio multipatch al que
%                  pertenece la solucion aproximada (ver sp_multipatch)
%   geometry     - estructura que describe el mapeo multipatch, parametro
%                  de salida de mp_geo_load
%   c_const      - valor constante de la velocidad de propagacion c del
%                  ultrasonido en todo el dominio (m/s)
%   alpha_const  - valor constante del coeficiente de atenuacion alpha en
%                  todo el dominio (Np/m)
%   rho_const    - valor constante de la densidad rho del tejido en todo
%                  el dominio (kg/m^3)
%
% Salida:
%   XYQ          - matriz [N,3] cuyas columnas son, respectivamente, la
%                  coordenada X, la coordenada Y y el valor de Q en cada
%                  punto de evaluacion. Se acumula parche por parche.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('--------------------------------------------------------------\n')
fprintf(' Calculos para graficar Q ...\n')

% Numero de puntos por direccion parametrica para la malla de evaluacion
if nargin < 7 || isempty(output_dir)
    output_dir = pwd;
end
if nargin < 8 || isempty(npoints)
    npoints = 250;
end
pts = {linspace(0, 1, npoints), linspace(0, 1, npoints)};

% Evaluamos la funcion Q para graficarla parche por parche
fig351 = figure(351); clf(fig351); hold on, grid on
axis([-0.15 0.15 -0.02 0.15])
colorbar
axis tight, axis equal
ax = gca;
ax.FontSize = 14;

XYQ = [];
for iptc = 1:space.npatch
    % Ciclo para cada parche
    % Evaluacion de la solucion aproximada u sobre la malla parametrica del
    % parche iptc. Zval contiene los valores de u y Fval las coordenadas
    % fisicas (Xval,Yval) en el dominio Omega'.
    [Zval, Fval] = sp_eval(u(space.gnum{iptc}), space.sp_patch{iptc}, geometry(iptc), pts);
    [Xval, Yval] = deal(squeeze(Fval(1,:,:)), squeeze(Fval(2,:,:)));

    % Evaluacion de las magnitudes constantes c, alpha y rho sobre la malla
    % de puntos del parche. Como en este caso k es constante, estas
    % magnitudes no dependen de (x,y); el helper devuelve matrices llenas
    % con los valores constantes.
    [cxy, alphaxy, rhoxy] = EvaluaConstantePto(Xval, Yval, c_const, alpha_const, rho_const);

    % Calculo de Q = (alpha / (rho*c)) * |u|^2
    deno   = rhoxy .* cxy;
    coefQ  = alphaxy ./ deno;
    Qval   = coefQ .* (real(Zval).^2 + imag(Zval).^2);

    % Grafico del parche actual de Q
    surf(Xval, Yval, Qval, 'EdgeColor', 'none')

    % Acumulacion de (X, Y, Q) en el array de salida
    XYQnew = [Xval(:), Yval(:), Qval(:)];
    XYQ    = [XYQ; XYQnew];
end
save_figure(fig351, output_dir, 'Q_distribution');
end
