function U = ADI(U,K,J,Aplus,Aminus,Bplus,Bminus)

% Crank-Nicolson (2nd-order accurate and unconditionally stable), 
% using an alternating-direction implicit (ADI) scheme

U = Aminus\(Bplus*U(:));
U = Bminus\(Aplus*U(:));
U = reshape(U,K,J);

