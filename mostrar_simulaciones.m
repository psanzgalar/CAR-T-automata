% ----------------------------
% SIMULACIÓN CAR-T PARA MOSTRAR
% ----------------------------

num_pasadas = 19;
nx = 200;
ny = 200;
ntum = 9750;

% CARGA UNA SOLA VEZ LOS DATOS REALES CD28
ruta_base = '/Users/pablosanz/Desktop/Articulo_CART_dominioActivacion/Datos/Nivel_Antigeno_Pasadas/Nalm6 FCS/p4';
ruta_base='C:\Users\Pablo\Desktop\TFM_definitivo\Datos\Nivel_Antigeno_Pasadas\Nalm6 FCS\p4';
archivo_CD28 = fullfile(ruta_base, 'CD19ab_Nalm6+CD28.fcs');
[dat, txt] = fca_readfcs(archivo_CD28);

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

datos_CD28 = dat(:, idx_pecy7);
datos_CD28=datos_CD28+abs(min(datos_CD28))+1;
% Convertir a constante para parfor
datos_CD28 = datos_CD28(:);  % aseguramos vector columna
A = zeros(nx, ny);
IDmap = zeros(nx, ny);
filas = zeros(1, ntum);
columnas = zeros(1, ntum);
k_local = 0;
id = 1;

%Creamos el tumor inicial

while k_local < ntum
    fila = randi([1, nx - 1]);
    columna = randi([1, ny - 1]);
    bloque_IDs = IDmap(fila:fila+1, columna:columna+1);
    bloque_ocupado = any(bloque_IDs(:) ~= 0);
    if ~bloque_ocupado
        antigen_level = datos_CD28(randi(length(datos_CD28)));
        while antigen_level <= 0
            antigen_level = datos_CD28(randi(length(datos_CD28)));
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
                antigen_level = datos_CD28(randi(length(datos_CD28)));
                while antigen_level <= 0
                    antigen_level = datos_CD28(randi(length(datos_CD28)));
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
%%
B=zeros(nx,ny);
m_vida=zeros(nx,ny);
m_duplicar=zeros(nx,ny);
tiempo_dup=3;
p=1/36;
p_muerte=0.95;
type=1;
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
    % RUN SIMULATION
    T = 24;
    for tiempo = 1:T
        % Visualización de la matriz A
        figure(100); clf;
        imagesc(A);

        % Límite máximo del valor tumoral (corte al percentil 95 para evitar extremos)
        niveles_tumorales = A(A >= 2) - 2;
        if isempty(niveles_tumorales)
            max_visual = 1;
        else
            max_visual = prctile(niveles_tumorales, 95);  % percentil 95
        end
        n_tum_vals = round(max_visual);  % número de colores del degradado

        % Evitar que sea menor que 1
        n_tum_vals = max(n_tum_vals, 1);

        % Crear degradado de rojo suave a oscuro
        rojo_suave = [1 0.8 0.8];
        rojo_oscuro = [0.5 0 0];
        rojos = [linspace(rojo_suave(1), rojo_oscuro(1), n_tum_vals)', ...
            linspace(rojo_suave(2), rojo_oscuro(2), n_tum_vals)', ...
            linspace(rojo_suave(3), rojo_oscuro(3), n_tum_vals)'];

        % Colormap completo: blanco (0), verde (1), rojo suave→oscuro (2+)
        mapa_colores = [
            1 1 1;        % blanco
            128/255 176/255 141/255;        % azul
            rojos         % tumorales
            ];

        colormap(mapa_colores);
        caxis([0, n_tum_vals + 1]);  % forzar rango visual hasta el corte
        title(sprintf('Pasada %d - Hora %d', pasada, tiempo));
        axis equal tight;
        axis off;
        drawnow;
        if tiempo == 1 && ismember(pasada, [1, 6, 12, 19])
            nombre_archivo = sprintf('malla_p%d_t1_3_01_9000.png', pasada);
            exportgraphics(gcf, nombre_archivo, 'Resolution', 300);
        end
        carga_toxica_ctumoral = zeros(1,length(filas));
        % First, CAR-T cells fire
        [filas,columnas,carga_toxica_ctumoral,B,m_duplicar,A,m_vida] = disparar_ordenado(filas, columnas, carga_toxica_ctumoral, B, m_duplicar, tiempo_dup,A,m_vida);
        contador=0;
        contador_muerte=zeros(1,length(carga_toxica_ctumoral));
        % Decide which tumor cells die
        for i=1:length(carga_toxica_ctumoral)
            matamos=morir(carga_toxica_ctumoral(i),A(filas(i),columnas(i)),p_muerte,0.1,9000);
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
    end
end



