function U = RK2(U,p,dt)
    
% 2nd-order Runge-Kutta
q1 = dt*Wave(U,p);
q2 = dt*Wave(U+q1,p)-q1;
U = U + q1 + q2/2;

