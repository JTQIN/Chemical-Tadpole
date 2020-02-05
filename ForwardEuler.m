function X = ForwardEuler(H,B,X,C,nX,D,dt,dx,dy)

S = H + B;
[GradXS GradYS GradXX GradYX] = UpwindGrad(S,X,nX,C,dx,dy);

for i=1:nX
    Xtemp2 = X(:,:,i);
    Xtemp = Xtemp2 + dt*D./H .* ( GradXS.*GradXX(:,:,i) + GradYS.*GradYX(:,:,i) ); ...
        % careful about dividing by H==0 here
    
    % boundary condition: no change in X at boundaries (also to correct
    % locations where H==0 at boundaries)
    X(:,:,i)=Xtemp;    
end
