function [new_A,new_B,vida_new,duplicar_new]=duplicar(A,B,m_vida,m_duplicar,type)
    [nx,ny]=size(A);
    recarga=1;
    %We search for the positions of the CAR-T cells that are in the last
    %hour of the duplication phase. These are the cells that will duplicate
    %in next hour
    [filaD,columnaD]=find(m_duplicar==1);
    % So that the duplicate is not orderly
    if type==1
    indices_random = randperm(length(filaD));
    filaD = filaD(indices_random);
    columnaD = columnaD(indices_random);
    end
    %We go through all the CAR-T cells that are going to duplicate
    for j=1:length(filaD)
        %A cell can duplicate if it has space for it 
        %To save the places where a CAR-T cell can duplicate
        F=zeros(1,8);
        C=zeros(1,8);
        %To count the number of places where a CAR-T cell can duplicate
        count=0;
        relative_positions=[-1,-1;-1,0;-1,1;0,-1;0,1;1,-1;1,0;1,1];
        for k = 1:size(relative_positions, 1)
            fila = filaD(j) + relative_positions(k, 1);
            columna = columnaD(j) + relative_positions(k, 2);
            if(fila>0 && columna>0 && fila<=nx && columna<=ny && A(fila,columna)==0)
                %The cell has place to duplicate
                count=count+1;
                F(count)=fila;
                C(count)=columna;
            end
        end
        if count>0
            %CAR-T cell can duplicate to some location
            k=round((count-1)*rand(1,1))+1;
            A(F(k),C(k))=1; %we set the position of the new CAR-T to 1
            %We update the martix B. Both CAR-T cells are fully recharged
            B(filaD(j),columnaD(j))=recarga;
            B(F(k),C(k))=recarga;
            %Full life in the new CAR-T
            m_vida(F(k),C(k))=336;
            %We update m_duplicar since the original CAR-T exits the
            %duplication phase
            m_duplicar(filaD(j),columnaD(j))=0;
        end
    end
    new_A=A;
    new_B=B;
    vida_new=m_vida;
    duplicar_new=m_duplicar;
end