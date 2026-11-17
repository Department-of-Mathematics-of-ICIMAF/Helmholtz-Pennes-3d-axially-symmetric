# Helmholtz-Pennes 3D con simetría axial

Resolución IgA de las ecuaciones de Helmholtz y Pennes en un dominio 3D con simetría axial.

## Ejecución automatizada de Helmholtz

Los orquestadores están en `scripts/` y leen los juegos de parámetros de `inputs/`. Desde la raíz del repositorio, en MATLAB:

```matlab
addpath('scripts')

% k constante
resultados_constante = run_helmholtz_k_constant();

% k variable por capas
resultados_variable = run_helmholtz_k_variable();
```

También se puede indicar otro archivo JSON:

```matlab
resultados = run_helmholtz_k_constant('ruta/al/archivo.json');
```

Cada archivo de entrada contiene un arreglo JSON. Para ejecutar varios casos, basta con añadir más objetos al arreglo: se procesan **en orden, uno después de otro**. Cada elemento de `resultados` contiene:

- `last_max`: posición y magnitud del último máximo local.
- `preim`: preimágenes de la línea focal.
- `XYQ`: salida de `evaluaQ`, con columnas `[X, Y, Q]`.
- `quarter_points`: puntos horizontal, vertical superior y vertical inferior cercanos a `|u(0,z_f)|/4`.
- `params`, `metadata` y `output_dir`: parámetros y ubicación de la corrida.

Las filas de `quarter_points` siguen el orden horizontal, vertical superior y vertical inferior. Sus columnas son `[rho, z, |u|, parche, xi, eta]`. La estructura `params` incluye esa misma matriz en `params.quarter_points` y los semiejes de la elipse en `params.quarter_semi_axes`. El semieje vertical es la distancia media entre `z_f` y los dos puntos verticales. Estas salidas quedan disponibles para la posterior integración con Pennes.

Los valores originales de los `main` están incluidos en los JSON; `rl` y `zmin` se calculan a partir de la geometría.

### Presets de tejido para k constante

El campo `tissue` admite `phantom`, `agua`, `piel`, `grasa` y `musculo`. El preset completa `c`, `aten` y `rho` con los valores usados en el código. Esos tres campos pueden incluirse explícitamente en un caso para sobrescribir el preset. Para un material distinto, usar `"tissue": "personalizado"` y proporcionar los tres valores.

### Carpetas de salida

Cada ejecución crea una subcarpeta independiente en:

- `Helmholtz/Homogeneous Helmholtz/outputs/`
- `Helmholtz/Variable k Helmholtz/outputs/`

El nombre combina índice, identificador del caso y refinamiento `n1`; una repetición añade un sufijo para no sobrescribir resultados. Dentro de cada subcarpeta se guardan:

- `parametros.txt`: resumen completo de parámetros, incluidos los puntos al cuarto de presión y los semiejes calculados.
- `helmholtz_results.mat`: `last_max`, `preim`, `XYQ`, `quarter_points` y `params`.
- `Q_values.mat`: valores de `XYQ`.
- `solution_2d.fig/.png`: valor absoluto de la solución en el plano `(rho,z)`.
- `solution_2d_quarter_ellipse.fig/.png`: duplicado de la figura 2D con la semi-elipse y los tres puntos destacados.
- Todas las demás figuras de la cadena de ejecución, en formatos `.fig` y `.png`. Las figuras se cierran después de guardarse.

## Reutilización con otra solución IgA, incluida Pennes

`utils/quarter_pressure_points.m` es independiente del modelo físico. Solo requiere un vector de coeficientes, su espacio multiparche y su geometría, todos compatibles con `sp_eval`. No utiliza variables globales ni funciones internas de Helmholtz o Pennes.

Para una solución temporal de Pennes se debe seleccionar **un instante** mediante `alpha(j,:)`; no se debe pasar la matriz temporal completa. Si se pasa `preim = []`, la función calcula las preimágenes con la misma geometría recibida, por lo que no depende de una parametrización de Helmholtz:

```matlab
% Desde la raiz del repositorio, con space, geometry y alpha ya calculados.
addpath(fullfile(pwd, 'utils'));
j = 1; % Indice del instante que se desea evaluar.
[h, vertical_up, vertical_down] = quarter_pressure_points( ...
    z_f, [], alpha(j,:), space, geometry, 800, 1000);
quarter_points_pennes = [h; vertical_up; vertical_down];
```

La operación siempre se aplica al **valor absoluto del campo recibido**. En Pennes, si se pasa directamente `alpha(j,:)`, el objetivo es `abs(T(0,z_f))/4`. La selección de un campo derivado, como un incremento de temperatura, debe hacerse antes de llamar a la función y estar representado por sus propios coeficientes IgA.

La utilidad `utils/preimagen_linea_focal.m` también es independiente del modelo: invierte únicamente los mapas NURBS de `geometry` y no usa el espacio ni la solución.

## Requisitos

MATLAB, las dependencias NURBS/GeoPDEs utilizadas por el proyecto y el Spline Toolbox deben estar disponibles en el `path`, igual que en la ejecución original de los `main`.
