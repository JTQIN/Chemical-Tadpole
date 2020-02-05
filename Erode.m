function [H, B, X] = Erode(H,Hnminus1,B,Bnminus1,X,Xnminus1,C,nX,K,J, ...
    Aplus,Aminus,Bplus,Bminus,D,dt,dx,dy,xr,n)
%
% Update soil thickness H, bedrock elevation B, and mineral abundances X
% according to the perturbations to these quantities by soil erosion.
% This is part of a splitting method, whereby this function computes a
% portion of the temporal changes in H, B, and X, and other functions
% compute the other portions of those changes.


%% Update soil thickness H

% DIFFUSION
S = B + H; % calculate surface elevation
S = ADI(S,K,J,Aplus,Aminus,Bplus,Bminus); % step surface elevation forward 
    % in time due to hillslope soil transport
Hnplus1 = S - B; % calculate new soil thickness


%% Update mineral abundances X

if n==1  % for the first timestep (this n is timestep, not stream power n)

    Xnplus1 = ForwardEuler(H,B,X,C,nX,D,dt,dx,dy);
    
else  % for all other timesteps
            
    Xnplus1 = RK4(H,S,X,nX,C,D,dx,dy,dt);

end
    
% Apply boundary conditions.
for i = 1:nX
    Xtemp = Xnplus1(:,:,i);
    Xtemp2 = X(:,:,i);
    Xtemp(C==1) = Xtemp2(C==1); % enforce boundary condition: no change
        % in X at boundaries.
    Xnplus1(:,:,i) = Xtemp;
end


%% Finish

% assign soil thickness and mineral concentrations at time t + dt to output
% arguments
H = Hnplus1;
X = Xnplus1;


