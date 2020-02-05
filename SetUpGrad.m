function [GX,GY] = SetUpGrad(K,J,C,dx,dy)

% construct a matrix operator that will give the gradient using
% second-order centered spatial differencing

N = K*J; % size of the operator matrix

X=[-1/(2*dx) 1/(2*dx)];
Y=[-1/(2*dy) 1/(2*dy)];

GX = sparse(N,N);
GY = sparse(N,N);

for k=1:N % loop through ROWS
    if (k-K) < 1 % the current element of U requires a "wraparound" row due to periodic boundary conditions
        GX(k,[k-K+J*K k+K]) = X;
        if k==1
            % need to build in support for periodic y-boundaries
        else
            GY(k,[k-1 k+1]) = Y;
        end
    elseif (k+K) > (J*K) % the current element of U requires a "wraparound" in the opposite direction
        GX(k,[k-K k+K-J*K]) = X;
        if k==N
            % need to build in support for periodic y-boundaries
        else
            GY(k,[k-1 k+1]) = Y;
        end        
    else % it's just a regular old row
        GX(k,[k-K k+K]) = X;
        GY(k,[k-1 k+1]) = Y;
    end
end


% loop through elements in the boundary matrix C and zero out operator rows
% corresponding to cells for which C==1
for i=1:K
    for j=1:J
        if C(i,j)==1
            k = K*(j-1)+i;
            GX(k,:)=0;
            GY(k,:)=0;
        end
    end
end

