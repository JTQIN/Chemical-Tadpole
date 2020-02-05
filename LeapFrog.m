function X = LeapFrog(H,B,X,Xnminus1,C,nX,K,J,GradX,GradY,D,dt,dx,dy)

S = H+B; % surface elevations

% Gradients from upwind differencing
[GradXS GradYS GradXX GradYX] = UpwindGrad(S,X,nX,dx,dy);

for i=1:nX

    Xtemp2 = Xnminus1(:,:,i); % don't need to keep this if we use Xn to correct boundaries
    Xtemp3 = X(:,:,i);
    Xtemp = Xtemp2 + 2*dt*D./H .* ( GradXS.*GradXX(:,:,i) + GradYS.*GradYX(:,:,i) ); % careful about dividing by H==0 here
    
    % boundary condition: no change in X at boundaries (also to correct
    % locations where H==0 at boundaries)
    Xtemp(C==1) = Xtemp3(C==1);
    X(:,:,i)=Xtemp;    
end


%%%%%%%%%%%%%%%%%%%%%%%%%

function V = ViscTerm(X,visc)

% construct numerical viscosity term as per Numerical Recipes, section 19-1

% pad X periodically in both directions
X = [X(:,end) X X(:,1)]; % x direction
X = [X(end,:); X; X(1,:)]; % y direction

Xdiff = X(2:end-1,3:end) - 2*X(2:end-1,2:end-1) + X(2:end-1,1:end-2); % this is just the numerator of a second-order second derivative in the x direction
Ydiff = X(3:end,2:end-1) - 2*X(2:end-1,2:end-1) + X(1:end-2,2:end-1); % ...same thing in the y direction

V = visc*(Xdiff + Ydiff);