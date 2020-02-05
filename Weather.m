function [H, X] = Weather(H,X,C,nX,kx,Ax,sx,wx,rhos,xr,dt)
%
% Update soil thickness H and mineral abundances X according to the 
% perturbations to these quantities by mineral weathering. This is part of
% a splitting method, whereby this function computes a portion of the 
% temporal changes in H and X, and other functions compute the other 
% portions of those changes.


%% Update mineral abundances X

Xnplus1 = zeros(size(X));

% calculate the summation term
BigSigma = 0; 
for i=1:nX
    BigSigma = BigSigma + (kx(i)*Ax(i)*X(:,:,i) - sx(i)*wx(i)/rhos);
end

% calculate Xn+1 for each mineral species
    
for i=1:nX
    Xtemp2 = X(:,:,i);
    Xtemp = X(:,:,i);
    DeltaX = dt * (-kx(i)*Ax(i)*X(:,:,i) + sx(i)*wx(i)/rhos + ...
        X(:,:,i).*BigSigma);
    Xtemp = Xtemp + DeltaX;

    % Apply boundary conditions.
    Xtemp(C==1) = Xtemp2(C==1); % enforce boundary condition: no change
        % in X at boundaries

    % Finalize X.
    Xnplus1(:,:,i) = Xtemp;
end


%% Update soil thickness H

% calculate the summation term at time n+1/2
BigSigma = 0; 
for i=1:nX
    BigSigma = BigSigma + (kx(i)*Ax(i)*0.5*(X(:,:,i)+Xnplus1(:,:,i)) - ...
        sx(i)*wx(i)/rhos);
end

Hnplus1 = H.*exp(-BigSigma*dt);

% Boundary condition: no change in H due to weathering at boundaries.
Hnplus1(C==1) = H(C==1);


%% Finish

% assign soil thickness and mineral concentrations at time t + dt to output
% arguments
H = Hnplus1;
X = Xnplus1;

