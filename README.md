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
- `params`, `metadata` y `output_dir`: parámetros y ubicación de la corrida.

Estas salidas quedan disponibles para la posterior integración con Pennes. Los valores originales de los `main` están incluidos en los JSON; `rl` y `zmin` se calculan a partir de la geometría.

### Presets de tejido para k constante

El campo `tissue` admite `phantom`, `agua`, `piel`, `grasa` y `musculo`. El preset completa `c`, `aten` y `rho` con los valores usados en el código. Esos tres campos pueden incluirse explícitamente en un caso para sobrescribir el preset. Para un material distinto, usar `"tissue": "personalizado"` y proporcionar los tres valores.

### Carpetas de salida

Cada ejecución crea una subcarpeta independiente en:

- `Helmholtz/Homogeneous Helmholtz/outputs/`
- `Helmholtz/Variable k Helmholtz/outputs/`

El nombre combina índice, identificador del caso y refinamiento `n1`; una repetición añade un sufijo para no sobrescribir resultados. Dentro de cada subcarpeta se guardan:

- `parametros.txt`: resumen completo de parámetros.
- `helmholtz_results.mat`: `last_max`, `preim`, `XYQ` y `params`.
- `Q_values.mat`: valores de `XYQ`.
- Todas las figuras de la cadena de ejecución, en formatos `.fig` y `.png`. Las figuras se cierran después de guardarse.

## Requisitos

MATLAB, las dependencias NURBS/GeoPDEs utilizadas por el proyecto y el Spline Toolbox deben estar disponibles en el `path`, igual que en la ejecución original de los `main`.
