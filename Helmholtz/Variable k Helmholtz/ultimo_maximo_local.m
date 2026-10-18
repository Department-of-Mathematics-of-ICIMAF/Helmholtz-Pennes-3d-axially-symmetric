function idx = ultimo_maximo_local(vector)
% Entrada: vector de números (muestra densa de una función)
% Salida: índice del último máximo local que supera a sus 4 vecinos
%         a cada lado. Si no existe, devuelve el índice del mayor
%         entre el primer y el último elemento.

    % este preprocesamiento es necesario para evitar comparar ceros
    % numericos del orden e-16 donde puede dar maximos que realmente no lo
    % son
    
    umbral = 1e-10;
    vector(abs(vector) < umbral) = 0;
    
    n = length(vector);
    
    comp = 10; %cantidad de datos en cada direccion a comparar
    
    %si el vector es muy pequenho no se puede procesar
    if n < 2*comp+1
        [~, idx_max] = max([vector(1), vector(n)]);
        if idx_max == 1
            idx = 1;
        else
            idx = n;
        end
        return;
    end
    
    %Procesar el vector desde el final hasta el incio
    for i = n-comp : -1 : comp+1
        if vector(i) > max(vector(i-comp:i-1)) && vector(i) > max(vector(i+1:i+comp))            
            idx = i;
            return;
        end
    end

    % No se encontró máximo local → función monótona, devolver maximo entre
    % el punto inicial y final.
    [~, idx_max] = max([vector(1), vector(n)]);
    if idx_max == 1
        idx = 1;
    else
        idx = n;
    end
end

