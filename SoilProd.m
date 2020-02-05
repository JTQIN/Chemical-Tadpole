function [H, B, X] = SoilProd(H,B,X,C,nX,alpha,P,rhos,rhor,xr,dt)
%
% Update soil thickness H, bedrock elevation B, and mineral abundances X
% according to the perturbations to these quantities by soil production.
% This is part of a splitting method, whereby this function computes a
% portion of the temporal changes in H, B, and X, and other functions
% compute the other portions of those changes.


%% Update soil thickness H

% Use exact solution for change in soil thickness over an interval dt.
% This solution can be derived via separation of variables.
Hnplus1 = (1/alpha)*log(alpha*P*dt/rhos + exp(alpha*H));

% boundary condition: no soil is produced, and therefore no bedrock
% lowering due to soil production, at boundaries.
Hnplus1(C==1) = H(C==1);


%% Update bedrock elevation B

% Conservation of mass gives new bedrock surface elevation
B = B + (rhos/rhor)*(H - Hnplus1);


%% Update soil mineral abundances X

for i=1:nX
    
    Xtemp = X(:,:,i).*H./Hnplus1 + xr(i)*(Hnplus1-H)./Hnplus1;
    
    % Apply boundary conditions.
    Xtemp2 = X(:,:,i);
    Xtemp(C==1) = Xtemp2(C==1); % enforce boundary condition: no change
        % in X at boundaries
    X(:,:,i) = Xtemp;
end


%% Finish

% assign soil thickness at time t + dt to output argument
H = Hnplus1;

