function X = AdamsBashforth(H,Hnminus1,S,Snminus1,X,Xnminus1,C,nX,D,dt,dx,dy)

[GradXS GradYS GradXX GradYX] = UpwindGrad(S,X,nX,C,dx,dy);
[GradXSnminus1 GradYSnminus1 GradXXnminus1 GradYXnminus1] = UpwindGrad(Snminus1,Xnminus1,nX,C,dx,dy);

for i=1:nX

    Xtemp2 = X(:,:,i);

    Xtemp = Xtemp2 + 1.5*dt*D./H.*(GradXS.*GradXX(:,:,i) + GradYS.*GradYX(:,:,i)) - 0.5*dt*D./Hnminus1.*(GradXSnminus1.*GradXXnminus1(:,:,i) + GradYSnminus1.*GradYXnminus1(:,:,i)); % where A is D*gradS/H dot gradX
        
    % boundary condition: no change in X at boundaries (also to correct
    % locations where H==0 at boundaries)
    Xtemp(C==1) = Xtemp2(C==1);
    X(:,:,i)=Xtemp;
    
    if min(min(X(:,:,i))) < 0
        disp(['Instability with mineral' num2str(i)])
    end
end
