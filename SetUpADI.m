function [Aplus, Aminus, Bplus, Bminus] = SetUpADI(K,J,D,dx,dy,dt)

N=J*K;

X=D*dt/2*(1/dx^2); X=[X -2*X X];
Y=D*dt/2*(1/dy^2); Y=[Y -2*Y Y];

A = sparse(N,N);
B = sparse(N,N);

for k=1:N % loop through ROWS
    if (rem(k-1,K)==0) || (rem(k,K)==0)
        % If the current element is on the upper or lower boundary, 
        % it is defined to have zero curvature, so A(k,:) and B(k,:) remain zero. Do nothing.
    elseif (k-K) < 1 % the current element of U requires a "wraparound" row due to periodic boundary conditions
%         A(k,k) = min2X; 
%         A(k,(k-K+J*K)) = X; 
%         A(k,k+K) = X;
        A(k,[k-K+J*K k k+K]) = X;
        B(k,(k-1):(k+1)) = Y;
    elseif (k+K) > (J*K) % the current element of U requires a "wraparound" in the opposite direction
%         A(k,k) = min2X; 
%         A(k,k-K) = X; 
%         A(k,(k+K-J*K)) = X;
        A(k,[k-K k k+K-J*K]) = X;
        B(k,(k-1):(k+1)) = Y;        
    else % it's just a regular old row
%         A(k,k) = min2X; 
%         A(k,k-K) = X; 
%         A(k,k+K) = X;
        A(k,[k-K k k+K]) = X;
        B(k,(k-1):(k+1)) = Y;
    end
end

I = sparse(1:N,1:N,1);

% The following 4 sparse arrays are used in the Crank-Nicolson ADI steps
Aplus = I + A;
Aminus = I - A;
Bplus = I + B;
Bminus = I - B;
