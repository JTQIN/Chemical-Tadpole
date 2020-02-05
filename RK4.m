function X = RK4(H,S,X,nX,C,D,dx,dy,dt)

% The function RHS evaluates the right-hand side of an equation with the
% form dX/dt = f(x,t), where f is a nonlinear operator
    
    % 4th-order Runge-Kutta
    q1 = dt*RHS(H,S,X,nX,C,D,dx,dy);
    q2 = dt*RHS(H,S,X+0.5*q1,nX,C,D,dx,dy);
    q3 = dt*RHS(H,S,X+0.5*q2,nX,C,D,dx,dy);
    q4 = dt*RHS(H,S,X+q3,nX,C,D,dx,dy);
    X = X + (q1+2*q2+2*q3+q4)/6;

