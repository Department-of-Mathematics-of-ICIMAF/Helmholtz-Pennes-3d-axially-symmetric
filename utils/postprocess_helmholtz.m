function [last_max, preim, quarter_points] = postprocess_helmholtz(u, space, geometry, params, output_dir)
%POSTPROCESS_HELMHOLTZ Postprocesa la solucion y guarda todas las figuras.
%   QUARTER_POINTS contiene, por filas, los puntos horizontal, vertical
%   superior y vertical inferior cercanos a |u(0,z_f)|/4. Sus columnas son
%   [rho, z, |u|, numero_de_parche, xi, eta].

fprintf('Postprocesamiento grafico\n')

% Solucion a lo largo del eje z.
zp = linspace(0, 1, params.npoints_axis);
pts = {0, zp};
[up1, Fzp1] = sp_eval(u(space.gnum{1}), space.sp_patch{1}, geometry(1), pts);
Zp1 = deal(squeeze(Fzp1(2,:,:)));
[up2, Fzp2] = sp_eval(u(space.gnum{2}), space.sp_patch{2}, geometry(2), pts);
Zp2 = deal(squeeze(Fzp2(2,:,:)));

fig400 = figure(400); clf(fig400); grid on; hold on;
plot(Zp1, real(up1), 'LineWidth', 2, 'MarkerSize', 6);
plot(Zp2, real(up2), 'LineWidth', 2, 'MarkerSize', 6);
title('Parte Real de la solucion de Helmholtz sobre el eje z');

fig500 = figure(500); clf(fig500); grid on; hold on;
plot(Zp1, imag(up1), 'LineWidth', 2, 'MarkerSize', 6);
plot(Zp2, imag(up2), 'LineWidth', 2, 'MarkerSize', 6);
title('Parte imaginaria de la solucion de Helmholtz sobre el eje z');

fig600 = figure(600); clf(fig600); grid on; hold on;
plot(Zp1, abs(up1), 'LineWidth', 2, 'MarkerSize', 6);
plot(Zp2, abs(up2), 'LineWidth', 2, 'MarkerSize', 6);
title('Valor absoluto de la solucion de Helmholtz sobre el eje z');

save_figure(fig400, output_dir, 'axis_real');
save_figure(fig500, output_dir, 'axis_imag');

% Punto focal mediante la heuristica del ultimo maximo local.
Z = [Zp1', Zp2'];
up = [abs(up1), abs(up2)];
ind = ultimo_maximo_local(up);
last_max = [Z(ind), up(ind)];

figure(fig600); hold on;
plot(last_max(1), last_max(2), 'ko', 'MarkerSize', 10);
previous = max(1, ind-10):ind-1;
next_indices = ind+1:min(numel(up), ind+10);
plot(Z(previous), up(previous), 'o', 'Color', 'green', ...
     'MarkerFaceColor', 'green');
plot(Z(next_indices), up(next_indices), 'o', 'Color', 'green', ...
     'MarkerFaceColor', 'green');

fprintf('Ultimo maximo local: (%e, %e)\n', last_max(1), last_max(2));
fprintf('Indice: %d\n', ind);
save_figure(fig600, output_dir, 'axis_abs');

% Perfil sobre la linea horizontal z = z_f (altura del foco).
z_f = Z(ind);
fprintf('Altura del foco z_f = %e\n', z_f);
preim = preimagen_linea_focal(geometry, z_f, params.npoints_focal);

fig700 = figure(700); clf(fig700); grid on; hold on;
fig800 = figure(800); clf(fig800); grid on; hold on;
fig900 = figure(900); clf(fig900); grid on; hold on;

for iptc = 1:numel(geometry)
    if preim(iptc).n == 0
        continue
    end
    % sp_eval trabaja con rejillas tensoriales: evaluar punto a punto.
    upf = zeros(1, preim(iptc).n);
    for k = 1:preim(iptc).n
        upk = sp_eval(u(space.gnum{iptc}), space.sp_patch{iptc}, ...
                      geometry(iptc), {preim(iptc).xi(k), preim(iptc).eta(k)});
        upf(k) = upk;
    end
    rhof = preim(iptc).rho';

    figure(fig700);
    plot(rhof, real(upf), 'LineWidth', 2, 'MarkerSize', 6); grid on;
    figure(fig800);
    plot(rhof, imag(upf), 'LineWidth', 2, 'MarkerSize', 6); grid on;
    figure(fig900);
    plot(rhof, abs(upf), 'LineWidth', 2, 'MarkerSize', 6); grid on;
end

figure(fig700);
xlabel('x');
title('Parte Real de la solucion de Helmholtz sobre la linea z = z_f');
figure(fig800);
xlabel('x');
title('Parte imaginaria de la solucion de Helmholtz sobre la linea z = z_f');
figure(fig900);
xlabel('x');
title('Valor absoluto de la solucion de Helmholtz sobre la linea z = z_f');

save_figure(fig700, output_dir, 'focal_real');
save_figure(fig800, output_dir, 'focal_imag');
save_figure(fig900, output_dir, 'focal_abs');

% Primeros puntos cercanos al cuarto de la presion focal en tres direcciones.
[horizontal, vertical_up, vertical_down] = ...
    quarter_pressure_points(z_f, preim, u, space, geometry, ...
                            params.npoints_axis);
quarter_points = [horizontal; vertical_up; vertical_down];
horizontal_semi_axis = abs(horizontal(1));
vertical_semi_axis = mean(abs([vertical_up(2), vertical_down(2)]-z_f));

% Solido de revolucion con una escala de color global.
pts_surface = linspace(0, 1, params.npoints_surface);
pts = {pts_surface, pts_surface};
[up1, Fzp1] = sp_eval(u(space.gnum{1}), space.sp_patch{1}, geometry(1), pts);
Rhop1 = deal(squeeze(Fzp1(1,:,:)));
Zp1 = deal(squeeze(Fzp1(2,:,:)));
[up2, Fzp2] = sp_eval(u(space.gnum{2}), space.sp_patch{2}, geometry(2), pts);
Rhop2 = deal(squeeze(Fzp2(1,:,:)));
Zp2 = deal(squeeze(Fzp2(2,:,:)));
[up3, Fzp3] = sp_eval(u(space.gnum{3}), space.sp_patch{3}, geometry(3), pts);
Rhop3 = deal(squeeze(Fzp3(1,:,:)));
Zp3 = deal(squeeze(Fzp3(2,:,:)));

sol1 = abs(up1);
sol2 = abs(up2);
sol3 = abs(up3);
global_clim = [min([sol1(:); sol2(:); sol3(:)]), ...
               max([sol1(:); sol2(:); sol3(:)])];

% Solucion 2D original y duplicado con la semi-elipse y los tres puntos.
rho_data = {Rhop1, Rhop2, Rhop3};
z_data = {Zp1, Zp2, Zp3};
sol_data = {sol1, sol2, sol3};
fig_2d = solution_2d_figure(rho_data, z_data, sol_data, global_clim, ...
    'Valor absoluto de la solucion de Helmholtz en 2D');
save_figure(fig_2d, output_dir, 'solution_2d');

[fig_ellipse, ax_ellipse] = solution_2d_figure(rho_data, z_data, sol_data, ...
    global_clim, 'Solucion 2D y semi-elipse a |u(0,z_f)|/4');
hold(ax_ellipse, 'on');
theta = linspace(-pi/2, pi/2, 400);
ellipse_rho = horizontal_semi_axis*cos(theta);
ellipse_z = z_f + vertical_semi_axis*sin(theta);
ellipse_handle = plot(ax_ellipse, ellipse_rho, ellipse_z, 'm-', ...
                      'LineWidth', 2);
points_handle = plot(ax_ellipse, quarter_points(:,1), quarter_points(:,2), ...
                     'ro', 'MarkerSize', 9, 'LineWidth', 2);
legend(ax_ellipse, [ellipse_handle, points_handle], ...
       {'Semi-elipse', 'Puntos cercanos a |u_f|/4'}, 'Location', 'best');
save_figure(fig_ellipse, output_dir, 'solution_2d_quarter_ellipse');

fig_surface = figure();
ax = axes('Parent', fig_surface);
view(ax, 3);
axis(ax, 'equal');
grid(ax, 'on');
xlabel(ax, 'X'); ylabel(ax, 'Y'); zlabel(ax, 'Z');
title(ax, 'abs(u_h)');
hold(ax, 'on');

crear_solido_revolucion2(ax, Rhop1, Zp1, sol1, params.ang, global_clim);
crear_solido_revolucion2(ax, Rhop2, Zp2, sol2, params.ang, global_clim);
crear_solido_revolucion2(ax, Rhop3, Zp3, sol3, params.ang, global_clim);
% camlight(ax, 'left');
axes(ax);           % Hace que 'ax' sea el eje actual
camlight('left');   % Ahora camlight act?a sobre el eje actual
lighting(ax, 'gouraud');
material(ax, 'dull');
colorbar(ax);
save_figure(fig_surface, output_dir, 'solution_revolution');
end

function [fig, ax] = solution_2d_figure(rho_data, z_data, sol_data, ...
                                        global_clim, plot_title)
%SOLUTION_2D_FIGURE Crea la misma vista 2D para ambas figuras guardadas.
    fig = figure();
    ax = axes('Parent', fig);
    hold(ax, 'on');
    for iptc = 1:numel(rho_data)
        surf(ax, rho_data{iptc}, z_data{iptc}, sol_data{iptc}, ...
             'EdgeColor', 'none', 'FaceColor', 'interp');
    end
    view(ax, 2);
    axis(ax, 'equal');
    axis(ax, 'tight');
    grid(ax, 'on');
    xlabel(ax, '\rho');
    ylabel(ax, 'z');
    title(ax, plot_title);
    if global_clim(2) > global_clim(1)
        caxis(ax, global_clim);
    end
    colorbar(ax);
end
