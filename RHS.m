function X = RHS(H,S,X,nX,C,D,dx,dy)

% function to evaluate the RHS of a nonlinear PDE: dX/dt = f(X,t), where f
% is a nonlinear operator

[GradXS GradYS GradXX GradYX] = UpwindGrad(S,X,nX,C,dx,dy);

for i=1:nX
    
    X(:,:,i) = D./H.*(GradXS.*GradXX(:,:,i) + GradYS.*GradYX(:,:,i));

end