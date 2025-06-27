function [f_sup, c_sup] = encontrar_esquina(IDmap, id_actual)
    [filas_id, columnas_id] = find(IDmap == id_actual);
    if numel(filas_id) ~= 4
        f_sup = NaN; c_sup = NaN;
        return
    end
    % Search for possible combinations that form a 2x2 block
    for i = 1:4
        for j = 1:4
            if i ~= j
                f0 = min(filas_id([i j]));
                c0 = min(columnas_id([i j]));
                try
                    bloque = IDmap(f0:f0+1, c0:c0+1);
                    if all(bloque(:) == id_actual)
                        f_sup = f0;
                        c_sup = c0;
                        return
                    end
                catch
                    continue
                end
            end
        end
    end
    f_sup = NaN; c_sup = NaN;
end
