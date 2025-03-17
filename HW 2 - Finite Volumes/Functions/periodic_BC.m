% -------------------------------------------------------------------
% Patrick Heng
% 08 Mar. 2025
% Imposes periodic boundary conditions on the A matrix.
% -------------------------------------------------------------------

function A = periodic_BC(A,u,v,n,m)

    NN = @(i,j) n*(j-1) + i;

    for i = 1:n
        if v(i) < 0
            A(NN(i,m),:) = 0;
            A(NN(i,m),:) = A(NN(i,1),:);
        else
            A(NN(i,1),:) = 0;
            A(NN(i,1),:) = A(NN(i,m),:);
        end
    end

    for i = 1:m
        if u(i) < 0
            A(NN(n,i),:) = 0;
            A(NN(n,i),:) = A(NN(1,i),:);
        else
            A(NN(1,i),:) = 0;
            A(NN(1,i),:) = A(NN(n,i),:);
        end
    end

end