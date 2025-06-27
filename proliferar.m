function [new_A,new_B,vida_new,duplicar_new,FILA,COLUMNA]=proliferar(filas,columnas,A,B,m_vida,m_duplicar,rho,type)
[nx,ny]=size(A);
% So that the proliferete order is not orderly
if(type==1)
    indices_random = randperm(length(filas));
    fila_random = filas(indices_random);
    columna_random = columnas(indices_random);
else
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
    %First we decided if the tumor cell can proliferate
    prolifera=Bernu(rho);
    %If prolifera=1, it means it proliferate
    if prolifera==1
        %First we must check the available positions, the tumor cell will
        % be able to duplicate itself in a box if it is empty or has a CAR-T cell,
        % in which case it will move it
        count=0;
        F=zeros(1,16);
        C=zeros(1,16);
        ant_level=A(fila_random(j),columna_random(j));
        for k = 1:size(relative_positions, 1)
            fila = fila_random(j) + relative_positions(k, 1);
            columna = columna_random(j) + relative_positions(k, 2);
            %We need to check if the tumor cell have enought space in
            %the automata
            if(fila>0 && columna>0 && fila<=nx-1 && columna<=ny-1)
                %Now we check if there isn't a tumor cell in this space
                if(A(fila,columna)<2 && A(fila+1,columna)<2 && A(fila,columna+1)<2 && A(fila+1,columna+1)<2)
                    %We can put a tumor cell in that space
                    count=count+1;
                    F(count)=fila;
                    C(count)=columna;
                end
            end
        end
        %Once we have checked where, the tumor CELL can proliferate, we
        %generate a random to decided where we will put the tumor CELL
        if count>0
            k=round((count-1)*rand(1,1))+1; %We generate the random position to which the cell moves
            %Now, we need to check if the tumor cell will move a CAR-T
            %cell
            %We need to check if there is a CAR-T cell in the place
            %the tumor will take place
            %We need to check if we have a CAR-T cell
            posible_positions=[0,0;1,0;0,1;1,1];
            for K = 1:size(posible_positions, 1)
                fila = F(k) + posible_positions(K, 1);
                columna = C(k) + posible_positions(K, 2);
                %We need to check if we have a CAR-T cell
                if(fila>0 && columna>0 && fila<=nx && columna<=ny && A(fila,columna)==1)
                    %We have a CAR-T cell in that place, so we try to move it
                    cont_moverT=0; %Counter for ther CAR-T cells that need to be moved to make space for the tumor cells
                    FF=zeros(1,12); %We are going to check 12 places, if there isn't any place, CAR-T cell
                    CC=zeros(1,12); %will desapeared
                    cart_move=[-1,-1;-1,0;-1,1;-1,2;0,-1;0,2;1,-1;1,2;
                        2,-1;2,0;2,1;2,2];
                    for k2=1:size(cart_move,1)
                        fila2 = F(k) + cart_move(k2, 1);
                        columna2 = C(k) + cart_move(k2, 2);
                        if(fila2>0 && columna2>0 && fila2<=nx && columna2<=ny && A(fila2,columna2)==0)
                            %There is a place where my CAR-T can move
                            cont_moverT=cont_moverT+1;
                            FF(cont_moverT)=fila2;
                            CC(cont_moverT)=columna2;
                        end
                    end
                    if cont_moverT>0
                        %The CAR-T cell have space to move
                        k2=round((cont_moverT-1))*rand(-1,1)+1;
                        A(FF(k2),CC(k2))=A(fila,columna);
                        B(FF(k2),CC(k2))=B(fila,columna);
                        B(fila,columna)=0;
                        m_vida(FF(k2),CC(k2))=m_vida(fila,columna);
                        m_duplicar(FF(k2),CC(k2))=m_duplicar(fila,columna);
                        m_duplicar(fila,columna)=0;
                        m_vida(fila,columna)=0;
                    end
                end
            end
            %We add the new cell that has proliferated (with a different ant
            %level)
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
            A(F(k),C(k))=ant_level2;
            A(F(k)+1,C(k))=ant_level2;
            A(F(k),C(k)+1)=ant_level2;
            A(F(k)+1,C(k)+1)=ant_level2;
            B(F(k),C(k))=0;
            B(F(k)+1,C(k))=0;
            B(F(k),C(k)+1)=0;
            B(F(k)+1,C(k)+1)=0;
            m_duplicar(F(k),C(k))=0;
            m_duplicar(F(k)+1,C(k))=0;
            m_duplicar(F(k),C(k)+1)=0;
            m_duplicar(F(k)+1,C(k)+1)=0;
            filas(end+1)=F(k); %we store the new tumor cell
            columnas(end+1)=C(k);
        end
    end
end
new_A=A;
new_B=B;
vida_new=m_vida;
duplicar_new=m_duplicar;
FILA=filas;
COLUMNA=columnas;
end
