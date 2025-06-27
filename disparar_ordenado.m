function [FILAS,COLUMNAS,pium, B2, m_duplicar_new,A2,m_vida_new]=disparar_ordenado(filas,columnas,cargatoxicactumoral,B,m_duplicar,tiempo_dup,A,m_vida)
% The parameter tiempo_dup refers to the time it takes for CAR-T cells to duplicate once they have been triggered

% With this loop, we will go through all the tumor cells.
% Since each tumor occupies 1 spaces, it can be in contact with 8 possible CAR-T cells.
% We will check these 12 possible positions and add the charge they have.
% If there is no CAR-T in a position, nothing will be added. If there is a CAR-T, its available charge will be added.

% If it enters the "if" statement, it's because there is a CAR-T cell,
% and a shot has occurred. If it doesn't enter the "if,"
% it's because either there is nothing, or the CAR-T cell has been depleted

%When a CAR-T cell fires, the tumor cell can kill it with a probability
% Shots are carried out in order, 
% giving tumor cells with higher antigen levels a greater chance of being targeted by a CAR-T cell
    [nx,ny]=size(B);
    i=1;
    cargaminima=0;
    niveles_antigeno = A(sub2ind(size(A), filas, columnas));
    [~, orden] = sort(niveles_antigeno, 'descend');
    filas = filas(orden);
    columnas = columnas(orden);
    cargatoxicactumoral = cargatoxicactumoral(orden);
    relative_positions = [
    -1, -1; -1, 0; -1, 1; -1, 2;
     0, -1;  0, 2;
     1, -1;  1, 2;
     2, -1;  2, 0;  2, 1;  2, 2
];
    while i<=length(filas)
        if (A(filas(i),columnas(i))>2)
        for k = 1:size(relative_positions, 1)
            fila = filas(i) + relative_positions(k, 1);
            columna = columnas(i) + relative_positions(k, 2);
            if(fila>0 && columna>0 && fila<=nx && columna<=ny)
                if B(fila, columna) > cargaminima && cargatoxicactumoral(i)==0 
                    cargatoxicactumoral(i) = 1;
                    m_duplicar(fila, columna) = tiempo_dup;
                    %We set the charge of the CAR-T cells to 0. Since they have
                    %fired, they are now depleted of charge
                    B(fila,columna)=0;
                end
            end
        end
        end
        i=i+1;
    end
pium = cargatoxicactumoral;
B2 = B;
m_duplicar_new = m_duplicar;
A2=A;
m_vida_new=m_vida;
FILAS=filas;
COLUMNAS=columnas;
end