function [new_A,new_B,vida_new,duplicar_new,FILA,COLUMNA,IDMAP,ID]=proliferar_hueco(filas,columnas,A,B,m_vida,m_duplicar,rho,type,IDmap,id)
[nx,ny]=size(A);
% So that the proliferete order is not orderly
if(type==1)
    indices_random = randperm(length(filas));
    fila_random = filas(indices_random);
    columna_random = columnas(indices_random);
else
    indices_random = 1:length(filas); 
    fila_random=filas;
    columna_random=columnas;
end
relative_positions = [
    -2, -2; -2, -1; -2, 0; -2, 1; -2,2;
    -1, -2; -1, 2; 
    0, -2;  0,2;
    1, -2; 1,2; 
    2, -2; 2, -1; 2, 0; 2, 1; 2,2;
];
%We go through all the tumor cells to proliferate (or not, them)
for j=1:length(filas)
    prolifera=Bernu(rho);
    if prolifera==1
        count=0;
        F=zeros(1,16);
        C=zeros(1,16);
        ant_level=A(fila_random(j),columna_random(j));
        if ant_level<2
            disp("Error")
        end
        if ant_level~=2
            des1=rand;
            if des1<0.15
                desv = 0;
                ant_level2=ant_level+desv;
            else
                a=-ant_level/3.75;
                b=ant_level/3.75;
                %a=-2933;
                %b=2933;
                desv = a + (b-a)*rand;
                ant_level2=ant_level+desv;
            end
        else
            ant_level2=2;
        end
        while ant_level2<2
            a=-ant_level/3.75;
            b=ant_level/3.75;
            %a=-2933;
            %b=2933;
            desv = a + (b-a)*rand;
            ant_level2=ant_level+desv;
        end
        for k = 1:size(relative_positions, 1)
            fila = fila_random(j) + relative_positions(k, 1);
            columna = columna_random(j) + relative_positions(k, 2);
            if(fila>0 && columna>0 && fila<=nx-1 && columna<=ny-1)
                if(A(fila,columna)<1 && A(fila+1,columna)<1 && A(fila,columna+1)<1 && A(fila+1,columna+1)<1)
                    count=count+1;
                    F(count)=fila;
                    C(count)=columna;
                end
            end
        end
        if count>0
            k=round((count-1)*rand(1,1))+1;
            IDmap(F(k),C(k))=id;
            IDmap(F(k)+1,C(k))=id;
            IDmap(F(k),C(k)+1)=id;
            IDmap(F(k)+1,C(k)+1)=id;
            id=id+1;
            A(F(k),C(k))=ant_level2;
            A(F(k)+1,C(k))=ant_level2;
            A(F(k),C(k)+1)=ant_level2;
            A(F(k)+1,C(k)+1)=ant_level2;
            filas(end+1)=F(k);
            columnas(end+1)=C(k);
            fila_random(end+1)=F(k);
            columna_random(end+1)=C(k);
        else
            for k = 1:size(relative_positions, 1)
                fila = fila_random(j) + relative_positions(k, 1);
                columna = columna_random(j) + relative_positions(k, 2);
                if(fila>0 && columna>0 && fila<=nx-1 && columna<=ny-1)
                    bloque_IDs = IDmap(fila:fila+1, columna:columna+1);
                    IDs_bloque = unique(bloque_IDs(bloque_IDs ~= 0));
                    recolocadas = true;
                    for idx = 1:length(IDs_bloque)
                        id_actual = IDs_bloque(idx);
                        [f_ant, c_ant] = encontrar_esquina(IDmap, id_actual);
                        if isnan(f_ant)
                            recolocadas = false;
                            break;
                        end
                        for pos_fila=1:length(filas)
                            if(filas(pos_fila)==f_ant && columnas(pos_fila)==c_ant)
                                break;
                            end
                        end
                        for pos_fila_random=1:length(fila_random)
                            if(fila_random(pos_fila_random)==f_ant && columna_random(pos_fila_random)==c_ant)
                                break;
                            end
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
                                    filas(pos_fila)=f_nueva;
                                    columnas(pos_fila)=c_nueva;
                                    fila_random(pos_fila_random)=f_nueva;
                                    columna_random(pos_fila_random)=c_nueva;
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
                end
            end
            if recolocadas==true
                count=0;
                F=zeros(1,16);
                C=zeros(1,16);
                for k = 1:size(relative_positions, 1)
                    fila = fila_random(j) + relative_positions(k, 1);
                    columna = columna_random(j) + relative_positions(k, 2);
                    if(fila>0 && columna>0 && fila<=nx-1 && columna<=ny-1)
                        if(A(fila,columna)<1 && A(fila+1,columna)<1 && A(fila,columna+1)<1 && A(fila+1,columna+1)<1)
                            count=count+1;
                            F(count)=fila;
                            C(count)=columna;
                        end
                    end
                end
                if count>0
                    k=round((count-1)*rand(1,1))+1;
                    IDmap(F(k),C(k))=id;
                    IDmap(F(k)+1,C(k))=id;
                    IDmap(F(k),C(k)+1)=id;
                    IDmap(F(k)+1,C(k)+1)=id;
                    id=id+1;
                    A(F(k),C(k))=ant_level2;
                    A(F(k)+1,C(k))=ant_level2;
                    A(F(k),C(k)+1)=ant_level2;
                    A(F(k)+1,C(k)+1)=ant_level2;
                    filas(end+1)=F(k);
                    columnas(end+1)=C(k);
                    fila_random(end+1)=F(k);
                    columna_random(end+1)=C(k);
                end
            end
        end
    end
end
new_A=A;
new_B=B;
vida_new=m_vida;
duplicar_new=m_duplicar;
FILA=filas;
COLUMNA=columnas;
IDMAP=IDmap;
ID=id;
end