function [Aplus, Aminus, Bplus, Bminus] = SetUpADI_bdry(K,J,D,dx,dy,dt,C)

% in this variation of the function, an argument C specifies cells in
% addition to the y boundaries for which the diffusion term should not
% alter elevations (& therefore it will be assigned zero curvature by the
% operators constructed here)

N=J*K;

X=D*dt/2*(1/dx^2); X=[X -2*X X];
Y=D*dt/2*(1/dy^2); Y=[Y -2*Y Y];

A = sparse(N,N);
B = sparse(N,N);

for k=1:N % loop through ROWS
    if (k-K) < 1 % the current element of U requires a "wraparound" row due
            % to periodic boundary conditions
        A(k,[k-K+J*K k k+K]) = X;
        if k==1
            % need to build in support for periodic y-boundaries
        else
            B(k,(k-1):(k+1)) = Y;
        end
    elseif (k+K) > (J*K) % the current element of U requires a "wraparound"
            % in the opposite direction
        A(k,[k-K k k+K-J*K]) = X;
        if k==N
            % need to build in support for periodic y-boundaries
        else
            B(k,(k-1):(k+1)) = Y;
        end        
    else % it's just a regular old row
        A(k,[k-K k k+K]) = X;
        B(k,(k-1):(k+1)) = Y;
    end
end


% loop through elements in the boundary matrix C and zero out operator rows
% corresponding to cells for which C==1
for i=1:K
    for j=1:J
        if C(i,j)==1
            k = K*(j-1)+i;
            A(k,:)=0;
            B(k,:)=0;
        end
    end
end


I = sparse(1:N,1:N,1);

% The following 4 sparse arrays are used in the Crank-Nicolson ADI steps
Aplus = I + A;
Aminus = I - A;
Bplus = I + B;
Bminus = I - B;
