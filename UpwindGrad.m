function [GradXS GradYS GradXX GradYX] = UpwindGrad(S,X,nX,C,dx,dy)

% calculates gradients in X and Y directions using first-order upwind
% differencing

[K J] = size(S);
invdx=1/dx;
invdy=1/dy;


%% Gradients in topography

% mirror boundary conditions
S = [S(:,2) S S(:,end-1)]; % x direction
S = [S(2,:); S; S(end-1,:)]; % y direction

% Compute centered difference grids in x- and y-directions.  Same size as
% initial S grid.
Xdiff = S(2:end-1,3:end) - S(2:end-1,1:end-2); % centered differences in positive x direction
Ydiff = S(3:end,2:end-1) - S(1:end-2,2:end-1); % ...same thing in the y direction

% Set very small values of Xdiff and Ydiff to zero.
diffthreshold = 1e-10;
Xdiff(abs(Xdiff) < diffthreshold) = 0;
Ydiff(abs(Ydiff) < diffthreshold) = 0;

% centered surface differencing in interior
GradXS = 0.5*invdx.*Xdiff;
GradYS = 0.5*invdy.*Ydiff;

% first-order differencing at boundaries
% On the left and right edges, only GradXS (not GradYS) is corrected.
% On the top and bottom edges, only GradYS (not GradXS) is corrected.
GradXS(:,1) = invdx*(S(2:end-1,3) - S(2:end-1,2));  % left edge
GradXS(:,end) = invdx*(S(2:end-1,end-1) - S(2:end-1,end-2));  % right edge
GradYS(1,:) = invdy*(S(3,2:end-1) - S(2,2:end-1));  % top edge
GradYS(end,:) = invdy*(S(end-1,2:end-1) - S(end-2,2:end-1));  % bottom edge

% diagonal differencing at corners - note this assumes dx==dy
GradXS([1 end],end) = (invdx/sqrt(2))*(S([2 end-1],end-1) - S([3 end-2],end-2));  % top and bottom right corners
GradXS([1 end],1) = (invdx/sqrt(2))*(S([3 end-2],3) - S([2 end-1],2));  % top and bottom left corners
GradYS(end,[1 end]) = (invdy/sqrt(2))*(S(end-1,[2 end-1]) - S(end-2,[3 end-2]));  % left and right bottom corners
GradYS(1,[1 end]) = (invdy/sqrt(2))*(S(3,[3 end-2]) - S(2,[2 end-1]));  % left and right top corners


%% Gradients in soil mineral abundances

GradXX = zeros(K,J,nX);
GradYX = zeros(K,J,nX);

for i=1:nX

    Xtemp = X(:,:,i);

    % mirror boundary conditions
    Xtemp = [Xtemp(:,2) Xtemp Xtemp(:,end-1)]; % x direction
    Xtemp = [Xtemp(2,:); Xtemp; Xtemp(end-1,:)]; % y direction
    
    % upwind differencing in interior.
    GXtemp = invdx*(Xdiff>0).*(Xtemp(2:end-1,3:end)-Xtemp(2:end-1,2:end-1)) + ...
        invdx*(Xdiff<=0).*(Xtemp(2:end-1,2:end-1)-Xtemp(2:end-1,1:end-2));
    GYtemp = invdy*(Ydiff>0).*(Xtemp(3:end,2:end-1)-Xtemp(2:end-1,2:end-1)) + ...
        invdy*(Ydiff<=0).*(Xtemp(2:end-1,2:end-1)-Xtemp(1:end-2,2:end-1));
 
    % first-order differencing at boundaries
    GXtemp(:,1) = invdx*(Xtemp(2:end-1,3) - Xtemp(2:end-1,2));  % left edge
    GXtemp(:,end) = invdx*(Xtemp(2:end-1,end-1) - Xtemp(2:end-1,end-2));  % right edge
    GYtemp(1,:) = invdy*(Xtemp(3,2:end-1) - Xtemp(2,2:end-1));  % top edge
    GYtemp(end,:) = invdy*(Xtemp(end-1,2:end-1) - Xtemp(end-2,2:end-1));  % bottom edge
     
    % diagonal differencing at corners - note this assumes dx==dy
    GXtemp([1 end],end) = (invdx/sqrt(2))*(Xtemp([2 end-1],end-1) - Xtemp([3 end-2],end-2)); % left and right bottom corners
    GXtemp([1 end],1) = (invdx/sqrt(2))*(Xtemp([3 end-2],3) - Xtemp([2 end-1],2)); % left and right top corners
    GYtemp(end,[1 end]) = (invdy/sqrt(2))*(Xtemp(end-1,[2 end-1]) - Xtemp(end-2,[3 end-2])); % top and bottom right corners
    GYtemp(1,[1 end]) = (invdy/sqrt(2))*(Xtemp(3,[3 end-2]) - Xtemp(2,[2 end-1])); % top and bottom left corners

    % finalize updated grids for GradXX and GradYX
    GradXX(:,:,i) = GXtemp;
    GradYX(:,:,i) = GYtemp;

end

