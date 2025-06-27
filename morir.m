function muerto=morir(cargatoxicactumoral,ant,p_muerte,min,bajo)
    % It receives the antigen level that each tumor cell has and applies a
    % Bernoulli based on the assigned death probabilities (if a CAR-T cell
    % have shoot before).
    if cargatoxicactumoral==0
        muerto=0;
    else
        muerto=Bernu(min+(p_muerte-min)/bajo*(ant-2));
        %muerto=Bernu(0.95);
    end
end