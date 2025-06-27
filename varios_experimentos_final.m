% ----------------------------
% COMPLETE CAR-T SIMULATION
% ----------------------------
clear
cd(fileparts(mfilename('fullpath')))
delete(gcp('nocreate'));          % cierra pools anteriores
p = parpool;                      % abre nuevo pool

% *** Adjunta SÓLO los archivos que usa el parfor ***
addAttachedFiles(p, { ...
    mfilename('fullpath'), ...          % este mismo script
    which('disparar_ordenado'), ...     % funciones llamadas dentro del parfor
    which('proliferar'), ...
    which('proliferar_hueco'), ...
    which('duplicar'), ...
    which('morir') });

minimo=[0.1];
bajo=[9000];
% LOAD THE REAL DATA ONLY ONCE
ruta_base = '/Users/pablosanz/Desktop/TFM/Datos/Nivel_Antigeno_Pasadas/Nalm6 FCS/p4';
archivo = fullfile(ruta_base, 'CD19ab_Nalm6+CD28.fcs');
[dat, txt] = fca_readfcs(archivo);

idx_pecy7 = -1;
for k = 1:txt.NumOfPar
    if isfield(txt.par(k), 'name')
        canal = txt.par(k).name;
    elseif isfield(txt.par(k), 'PnN')
        canal = txt.par(k).PnN;
    else
        canal = '';
    end
    if strcmpi(strtrim(canal), 'PE-Cy7-A')
        idx_pecy7 = k;
        break;
    end
end
if idx_pecy7 < 0
    error('Canal PE-Cy7-A no encontrado.');
end

datos = dat(:, idx_pecy7);
datos=datos+abs(min(datos))+1;
% Convert to constant for parfor
datos = datos(:);  
for reintentos=1:length(minimo)
    num_reps = 100;
    num_pasadas = 19;
    nx = 200;
    ny = 200;
    ntum = 9750;
    type = 1;
    all_life=zeros(num_reps,num_pasadas);
    all_max_life=zeros(num_reps,num_pasadas);
    all_min_life=zeros(num_reps,num_pasadas);
    % ----------------------------
    % PARALLEL EXECUTION
    % ----------------------------
    parfor rep = 1:num_reps
        fprintf("Simulación %d de %d\n", rep, num_reps);
        life_local=zeros(1,num_pasadas);
        max_life_local=zeros(1,num_pasadas);
        min_life_local=ones(1,num_pasadas)*Inf;
        A = zeros(nx, ny);
        IDmap = zeros(nx, ny);
        filas = zeros(1, ntum);
        columnas = zeros(1, ntum);
        k_local = 0;
        id = 1;
        % Create the tumor and load values from FCS files
        while k_local < ntum
            fila = randi([1, nx - 1]);
            columna = randi([1, ny - 1]);
            bloque_IDs = IDmap(fila:fila+1, columna:columna+1);
            bloque_ocupado = any(bloque_IDs(:) ~= 0);
            if ~bloque_ocupado
                antigen_level = datos(randi(length(datos)));
                while antigen_level <= 0
                    antigen_level = datos(randi(length(datos)));
                end
                valor = antigen_level + 2;
                A(fila:fila+1, columna:columna+1) = valor;
                IDmap(fila:fila+1, columna:columna+1) = id;
                filas(k_local+1) = fila;
                columnas(k_local+1) = columna;
                k_local = k_local + 1;
                id = id + 1;
            else
                IDs_bloque = unique(bloque_IDs(bloque_IDs ~= 0));
                recolocadas = true;
                for idxb = 1:length(IDs_bloque)
                    id_actual = IDs_bloque(idxb);
                    [f_ant, c_ant] = find(IDmap == id_actual, 1);
                    if isempty(f_ant)
                        recolocadas = false;
                        break;
                    end
                    direcciones = [0 1; 0 -1; 1 0; -1 0];
                    movido = false;
                    for d = 1:4
                        df = direcciones(d, 1);
                        dc = direcciones(d, 2);
                        f_nueva = f_ant + df;
                        c_nueva = c_ant + dc;
                        if f_nueva >= 1 && f_nueva <= nx-1 && c_nueva >= 1 && c_nueva <= ny-1
                            id_block = IDmap(f_nueva:f_nueva+1, c_nueva:c_nueva+1);
                            if all(id_block(:) == 0 | id_block(:) == id_actual)
                                antigen_val = A(f_ant, c_ant);
                                A(f_ant:f_ant+1, c_ant:c_ant+1) = 0;
                                IDmap(f_ant:f_ant+1, c_ant:c_ant+1) = 0;
                                A(f_nueva:f_nueva+1, c_nueva:c_nueva+1) = antigen_val;
                                IDmap(f_nueva:f_nueva+1, c_nueva:c_nueva+1) = id_actual;
                                filas(id_actual) = f_nueva;
                                columnas(id_actual) = c_nueva;
                                movido = true;
                                break;
                            end
                        end
                    end
                    if ~movido
                        recolocadas = false;
                        break;
                    end
                end
                if recolocadas
                    bloque_IDs = IDmap(fila:fila+1, columna:columna+1);
                    bloque_ocupado = any(bloque_IDs(:) ~= 0);
                    if ~bloque_ocupado
                        antigen_level = datos(randi(length(datos)));
                        while antigen_level <= 0
                            antigen_level = datos(randi(length(datos)));
                        end
                        valor = antigen_level + 2;
                        A(fila:fila+1, columna:columna+1) = valor;
                        IDmap(fila:fila+1, columna:columna+1) = id;
                        filas(k_local+1) = fila;
                        columnas(k_local+1) = columna;
                        k_local = k_local + 1;
                        id = id + 1;
                    end
                end
            end
        end
        B=zeros(nx,ny);
        m_vida=zeros(nx,ny);
        m_duplicar=zeros(nx,ny);
        tiempo_dup=3;
        p=1/36;
        p_muerte=0.95;
        proporcion=[0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.25,0.1,0.1,0.25,0.25,0.25,0.1,0.25,0.1,0.1,0.25,0.25];
        media_antigeno_local = zeros(1, num_pasadas+1);
        valores_antigeno_local=cell(1,num_pasadas+1);
        for pasada = 1:num_pasadas
            % FILL MATRIX WITH CAR-T
            for fi = 1:nx
                for ci = 1:ny
                    if A(fi,ci)==0
                        A(fi,ci)=1;
                        B(fi,ci)=1;
                        m_vida(fi,ci)=336;
                    end
                end
            end
            % CALCULATE ANTIGEN LEVELS
            niveles_antigeno = A(A >= 2) - 2;
            media_antigeno_local(pasada) = mean(niveles_antigeno);
            valores_antigeno_local{pasada}=niveles_antigeno(:);
            % RUN SIMULATION
            T = 24;
            for tiempo = 1:T
                carga_toxica_ctumoral = zeros(1,length(filas));
                % First, CAR-T cells fire
                [filas,columnas,carga_toxica_ctumoral,B,m_duplicar,A,m_vida] = disparar_ordenado(filas, columnas, carga_toxica_ctumoral, B, m_duplicar, tiempo_dup,A,m_vida);
                contador=0;
                contador_muerte=zeros(1,length(carga_toxica_ctumoral));
                % Decide which tumor cells die
                for i=1:length(carga_toxica_ctumoral)
                    matamos=morir(carga_toxica_ctumoral(i),A(filas(i),columnas(i)),p_muerte,minimo(reintentos),bajo(reintentos));
                    contador_muerte(i)=matamos;
                    if matamos==1
                        A(filas(i),columnas(i))=0; A(filas(i)+1,columnas(i)+1)=0;
                        A(filas(i),columnas(i)+1)=0; A(filas(i)+1,columnas(i))=0;
                    else
                        contador=contador+1;
                    end
                end
                fila_new=zeros(1,contador); 
                columna_new=zeros(1,contador); 
                contador=1;
                for i=1:length(contador_muerte)
                    if contador_muerte(i)==0
                        fila_new(contador)=filas(i); 
                        columna_new(contador)=columnas(i); 
                        contador=contador+1;
                    end
                end
                filas=fila_new; 
                columnas=columna_new;
                %Now the tumor proliferate
                [A,B,m_vida,m_duplicar,filas,columnas] = proliferar(filas,columnas,A,B,m_vida,m_duplicar,p,type);
                %And CAR-T cells duplicate (fir tgise scheuled to do so)
                [A,B,m_vida,m_duplicar] = duplicar(A,B,m_vida,m_duplicar,type);
                [f_dup2,c_dup2]=find(m_duplicar>1);
                % Subtract one hour from the CAR-T duplication process
                for r=1:length(f_dup2) 
                    m_duplicar(f_dup2(r),c_dup2(r))=m_duplicar(f_dup2(r),c_dup2(r))-1; 
                end
                % Kill CAR-T cells with one hour of life remaining (irrelevant in this case)
                [f_muerte,c_muerte]=find(m_vida==1);
                for r=1:length(f_muerte)
                    m_vida(f_muerte(r),c_muerte(r))=0;
                    A(f_muerte(r),c_muerte(r))=0; 
                    B(f_muerte(r),c_muerte(r))=0;
                    m_duplicar(f_muerte(r),c_muerte(r))=0; 
                end
                % Sustract one hour of life from the CAR-T cells (irrelevant in this case)
                [f_vida,c_vida]=find(m_vida>1);
                for r=1:length(f_vida) 
                    m_vida(f_vida(r),c_vida(r))=m_vida(f_vida(r),c_vida(r))-1; 
                end
            end
            %Compute the tumor cells
            numero_tum_pasada=nnz(A>=2)/4;
            life_local(pasada)=numero_tum_pasada;
            max_life_local(pasada)=max(max_life_local(pasada),numero_tum_pasada);
            min_life_local(pasada)=min(min_life_local(pasada),numero_tum_pasada);
            % PREPARE NEXT ROUND
            if pasada < num_pasadas
                B=zeros(nx,ny); 
                m_vida=zeros(nx,ny); 
                m_duplicar=zeros(nx,ny); 
                objetivo=proporcion(pasada+1); 
                objetivo = 9750*(objetivo==0.1) + 9408*(objetivo~=0.1);
                for fi=1:nx
                    for ci=1:ny
                        if A(fi,ci)==1
                            A(fi,ci)=0; B(fi,ci)=0; m_vida(fi,ci)=0;
                        end
                    end
                end
                IDmap = zeros(nx,ny); id = 1;
                while(id<=length(filas))
                    IDmap(filas(id),columnas(id))=id;
                    IDmap(filas(id)+1,columnas(id))=id;
                    IDmap(filas(id),columnas(id)+1)=id;
                    IDmap(filas(id)+1,columnas(id)+1)=id;
                    id=id+1;
                end
                numero_tum = nnz(A>=2)/4;
                %We let the tumor proliferate until the goal size
                while numero_tum<objetivo
                    [A,B,m_vida,m_duplicar,filas,columnas,IDmap,id] = proliferar_hueco(filas,columnas,A,B,m_vida,m_duplicar,p,type,IDmap,id);
                    numero_tum=nnz(A>=2)/4;
                end
                % CALCULATE ANTIGEN LEVELS
                niveles_antigeno = A(A >= 2) - 2;
                media_antigeno_local(pasada+1) = mean(niveles_antigeno);
                valores_antigeno_local{pasada+1}=niveles_antigeno(:);
            end
        end
        % PREPARE LAST PROLIFERATION FOR MEASUREMENT
        B=zeros(nx,ny); 
        m_vida=zeros(nx,ny); 
        m_duplicar=zeros(nx,ny); 
        objetivo=9750;
        for fi=1:nx
            for ci=1:ny
                if A(fi,ci)==1
                    A(fi,ci)=0; 
                    B(fi,ci)=0; 
                    m_vida(fi,ci)=0;
                end
            end
        end
        IDmap = zeros(nx,ny); id = 1;
        while(id<=length(filas))
            IDmap(filas(id),columnas(id))=id;
            IDmap(filas(id)+1,columnas(id))=id;
            IDmap(filas(id),columnas(id)+1)=id;
            IDmap(filas(id)+1,columnas(id)+1)=id;
            id=id+1;
        end
        numero_tum = nnz(A>=2)/4;
        while numero_tum<objetivo
            [A,B,m_vida,m_duplicar,filas,columnas,IDmap,id] = proliferar_hueco(filas,columnas,A,B,m_vida,m_duplicar,p,type,IDmap,id);
            numero_tum=nnz(A>=2)/4;
        end
        % COMPUTE ANTIGEN LEVELS
        niveles_antigeno = A(A >= 2) - 2;
        media_antigeno_local(num_pasadas+1) = mean(niveles_antigeno);
        valores_antigeno_local{num_pasadas+1}=niveles_antigeno(:);
        guardar_resultado(rep, media_antigeno_local,valores_antigeno_local);
        all_life(rep,:)=life_local;
        all_max_life(rep,:)=max_life_local;
        all_min_life(rep,:)=min_life_local;
    end
    media_tum=mean(all_life,1);
    max_tum=max(all_max_life,[],1);
    min_tum=min(all_min_life,[],1);
    save('summary_results_1_01_9000.mat','media_tum','max_tum','min_tum');
    delete(gcp('nocreate'));
    p = parpool('Processes');              % ← nuevo pool
    addAttachedFiles(p, { mfilename('fullpath') });   % ← ¡esta línea!
end
%%
% COLLECT ALL RESULTS
num_reps = 100;
num_pasadas = 19;
bins = 10.^((-100:500)/100);
media_total = zeros(num_pasadas+1, 1);
valores_todos=cell(num_pasadas+1,1);
todas_las_medias = zeros(num_reps, num_pasadas+1);
for rep=1:num_reps
    datos=load(sprintf('temp_result_rep%d_3_0_9000.mat',rep));
    media_total=media_total+datos.media';
    todas_las_medias(rep, :) = datos.media;
    for p=1:num_pasadas+1
        valores_todos{p}=[valores_todos{p};datos.valores_antigeno{p}];
    end
end

max_y = 0;
for i = 1:num_pasadas+1
    [~, edges] = histcounts(valores_todos{i}, bins);
    counts = histcounts(valores_todos{i}, edges);
    max_y = max(max_y, max(counts));
end

media_prom = media_total / num_reps;
media_minima = min(todas_las_medias, [], 1)'; 
media_maxima = max(todas_las_medias, [], 1)'; 
figure('Name', 'Curvas histogramas CD28', 'Position', [100, 100, 1600, 800]);
tiledlayout(2, ceil(num_pasadas / 2), 'TileSpacing', 'compact', 'Padding', 'compact');

color_rgb = [128, 176, 141]/255;

for i = 1:num_pasadas+1
    nexttile;
    if i==1
        h = histogram(valores_todos{i}, 'BinEdges', bins);
        h.FaceColor = color_rgb;
        h.EdgeColor = 'none';
        ylim([0 max_y]);
        set(gca, ...
            'XScale',     'log', ...
            'XTick',      [1 10 100 1000 10000], ...
            'XTickLabel', {'10^0','10^1','10^2','10^3','10^4'});

        title(sprintf('Pob. inicial'));
        xlabel('Nivel de antígeno');
    else
        h = histogram(valores_todos{i}, 'BinEdges', bins);
        h.FaceColor = color_rgb;
        h.EdgeColor = 'none';
        ylim([0 max_y]);
        set(gca, ...
            'XScale',     'log', ...
            'XTick',      [1 10 100 1000 10000], ...
            'XTickLabel', {'10^0','10^1','10^2','10^3','10^4'});

        title(sprintf('Pasada %d',i-1));
        xlabel('Nivel de antígeno');
    end
end
% ----------------------------
% MEAN ANTIGEN LEVEL
% ----------------------------
color2=[178, 208, 182]/255;
figure('Name','Media del antígeno', 'Position', [100, 100, 800, 400]);
fill([0:num_pasadas, fliplr(0:num_pasadas)], [media_maxima', fliplr(media_minima')],color2, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
hold on;
hold on;
plot(0:num_pasadas, media_prom, '-o', 'LineWidth', 2, 'Color', color_rgb);
hold on;
plot(0:num_pasadas, media_maxima, '--', 'LineWidth', 1.5, 'Color', color_rgb);
plot(0:num_pasadas, media_minima, '--', 'LineWidth', 1.5, 'Color', color_rgb);
xlim([0,19])
legend({'Media'}, ...
       'Location','best');
ylim([0 10200])
yticks([0 2500 5000 7500 10000])
set(gca,'Fontsize',20)
xlabel('Pasada');
ylabel('Nivel antígeno');
%title('Evolución media del antígeno');
grid on;
% Cargar datos
load('summary_results_3_0_9000.mat');  % media_tum, max_tum, min_tum

% Datos de configuración
proporcion = [0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.25,0.1,0.1,0.25,0.25,0.25,0.1,0.25,0.1,0.1,0.25,0.25];
num_pasadas = length(proporcion);
celulas_iniciales = 9750 * (proporcion == 0.1) + 9408 * (proporcion ~= 0.1);

% Proporciones de muerte
proporcion_muertes_media = (celulas_iniciales - media_tum) ./ celulas_iniciales;
proporcion_muertes_min = (celulas_iniciales - max_tum) ./ celulas_iniciales;
proporcion_muertes_max = (celulas_iniciales - min_tum) ./ celulas_iniciales;

proporcion_muertes_min = max(0, min(1, proporcion_muertes_min));
proporcion_muertes_max = max(0, min(1, proporcion_muertes_max));

% Calcular desviación
err_low = proporcion_muertes_media - proporcion_muertes_min;
err_high = proporcion_muertes_max - proporcion_muertes_media;

% Crear gráfico
figure('Name','Proporción de células muertas por pasada con rango', 'Position', [100, 100, 900, 400]);
hold on;

% Dibujar las barras una a una según la proporción de CAR-T
for p = 1:num_pasadas
    if proporcion(p) == 0.25
        color_barra = [197, 122, 120] / 255;  % color más oscuro
    else
        color_barra = [235, 170, 169] / 255;  % color original
    end
    bar(p, proporcion_muertes_media(p), 'FaceColor', color_barra);
end
% Barras de error vertical tipo desviación estándar
errorbar(1:num_pasadas, proporcion_muertes_media, err_low, err_high, 'Color',  [160, 80, 78]/255, 'LineStyle', 'none', 'LineWidth', 1.5, 'CapSize', 8);

dummy_x = 1;
dummy_y = proporcion_muertes_media(dummy_x);
dummy_err_low = dummy_y - proporcion_muertes_min(dummy_x);
dummy_err_high = proporcion_muertes_max(dummy_x) - dummy_y;

hErr = errorbar(dummy_x, dummy_y, dummy_err_low, dummy_err_high, ...
    'Color', [160, 80, 78]/255, ...
    'LineStyle', 'none', ...
    'LineWidth', 1.5, ...
    'CapSize', 8, ...
    'DisplayName', 'Rango [mín–máx]');

% Barras invisibles para la leyenda
b10 = bar(NaN, NaN, 'FaceColor', [235, 170, 169]/255, 'DisplayName', '10% CAR-T');
b25 = bar(NaN, NaN, 'FaceColor', [197, 122, 120]/255, 'DisplayName', '25% CAR-T');

% Leyenda
legend([b10, b25, hErr], {'E:T=0.1:1', 'E:T=0.25:1', 'Rango [mín–máx]'}, 'Location', 'northeast');
xlabel('Pasada');
ylabel('Proporción');
%title('Células tumorales eliminadas por pasada');
ylim([0 1]);
xlim([0.5 19.5]);
xticks(1:19)
set(gca,'Fontsize',20)
grid on;

function guardar_resultado(rep, media, valores_antigeno)
nombre = sprintf('temp_result_rep%d_1_01_9000.mat', rep);
save(nombre, 'media', 'valores_antigeno');
end
