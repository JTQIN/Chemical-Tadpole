% Script for computing coevolution of topography and soil chemical
% weathering, by Ken Ferrier and Taylor Perron.
%
% Input grids and parameters:
%
%   B: Matrix of initial bedrock surface elevation
%   H: Matrix of initial soil thickness
%   X: K x J x nX array of initial mineral concentrations (where n is the
%      number of mineral species)
%   C: Matrix specifying location of 'boundary' cells where soil erosion, 
%      soil production, and weathering are zero
%   run_name: Character string naming the run. If specified, the parameters 
%             and elevations at each timestep will be saved in a binary 
%             file called <run_name>.mat
%   p: A struct array of parameters that includes:
% 
%     p.dt             Timestep (yr)
%     p.tf             Total time of the simulation (yr)
%     p.dx, p.dy       Grid spacing in the x and y directions (m)
%     p.D              Hillslope diffusivity (m^2/yr)
%     p.K              Coefficient in stream power-type incision law 
%                       (kg m^(1+2m) yr^-2)
%     p.tauc           Threshold for fluvial incision
%     p.m, p.w         Exponents on drainage area and slope, respectively
%     p.Kw             Coefficient in channel width power law: 
%                       width = Kw*A^wexp
%     p.wexp           Exponent in channel width power law
%     p.E              Rate of relative uplift or base level lowering (m/yr)
%     p.saveint        Grids are saved every saveint iterations. 
%                      No output will be saved if save_interval==0.
%     p.plotint        Plot will be redrawn every plotint iterations
%     p.zrange         Range of vertical axis of plots
%     p.vexag          Vertical exaggeration of output plots
%     p.P              Soil production rate at zero soil depth 
%                       (kg m^-2 yr^-1)
%     p.alpha          Decay constant of soil production rate with soil 
%                       depth (m^-1)
%     p.rhor           Rock bulk density (kg m^-3)
%     p.rhos           Soil bulk density (kg m^-3)
%     p.wx             Mineral molar mass (kg mol^-1)
%     p.kx             Mineral reactivity (mol m^-2 yr^-1)
%     p.sx             Clay production rate (mol m^-3 yr^-1)
%     p.Ax             Mineral specific surface area (m^2 mol^-1)
%     p.xr             nX-length vector of mineral concentrations in 
%                       bedrock (dimensionless)
%     p.xname          Cell array of mineral names, e.g. {'quartz,
%                       'plagioclase','K-feldspar'}
%     p.x2plot         Index of mineral to plot (1 <= x2plot <= nX)
%     p.Hchannelvalue  Soil thickness assigned to channels (m).  A very
%                       small value greater than zero is chosen here to 
%                       avoid divide-by-zero issues.
%     p.Estepfactor    Factor by which rock uplift rate increases in
%                       experiments with imposed changes in rock uplift
%                       rate
%     p.Pstepfactor    Factor by which soil production rate increases in
%                       experiments with imposed changes in soil production
%                       rate
%     p.Dstepfactor    Factor by which soil transport rate increases in
%                       experiments with imposed changes in soil production
%                       rate
%     p.kxstepfactor   Factor by which dissolution rates increase in
%                       experiments with imposed changes in dissolution
%                       rate
%     p.mark           criterion for approaching new steady state after
%                       step
%     p.n_atstep       Time step at which to impose perturbation
%     p.experiment_type   String indicating the type of imposed
%                          perturbation.  Must be one of 'tectonic',
%                          'climatic', 'D', 'P', 'kx', or 'steady_state'.
%
%
% Important note: The boundaries are treated in a way that is specific to
% one grid, particularly in the UpwindGrad function.  This will need to be
% modified for grids with other sizes, shapes, or channel networks.


%% Load input grids and parameter values.

clear
close all
load('input_twophase101_10yr_P0.5mmyr_E0.1mmyr.mat')  % input grids


%% Create perturbation.

% Here is an example for a step perturbation in rock uplift rate.
p.Estepfactor = 1.5;  % factor by which this changes in the imposed step.
    % E.g., a factor of 1.5 would be a 50% increase.
p.experiment_type = 'tectonic';
p.tf = 1e3;  % duration of simulation (years)
run_name = 'test';

% Save inputs as struct for the output.
inputs.B = B;
inputs.C = C;
inputs.H = H;
inputs.X = X;
inputs.p = p;
inputs.run_name = run_name;

% Save outputs for do_output = 1.  Don't save outputs for do_output = 0.
do_output = 1;


%% PARAMETERS

% Parse the parameters input vector
dx = p.dx; 
dy = p.dy;
D = p.D;
E = p.E;
Kf = p.K;
P = p.P;
alpha = p.alpha;
rhor = p.rhor;
rhos = p.rhos;
kx = p.kx;
Ax = p.Ax;
sx = p.sx;
wx = p.wx;
xr = p.xr;
xname = p.xname;
x2plot = p.x2plot;
mark = p.mark;  % criterion for approaching new steady state after step
experiment_type = p.experiment_type;

% Define other parameters.
Cnonchannels = C~=1;  % grid cells not in channels = 1, otherwise = 0
Couteredge = zeros(size(C));  % grid cells on outer edge; this avoids the
    % channels in the interior of the grid
Couteredge(1,:) = 1;
Couteredge(end,:) = 1;
Couteredge(:,1) = 1;
Couteredge(:,end) = 1;
Anonchannels = sum(sum(Cnonchannels * dx * dy));  % landscape area not in
    % channels
    
% Define x-y grids and save as parts of inputs.
[x, y] = meshgrid(0:dx:dx*(size(C,1)-1), 0:dy:dy*(size(C,1)-1));
inputs.x = x;
inputs.y = y;

% Don't save output if p.saveint <= 0.
if round(p.saveint) <= 0
    do_output = 0;
else
    save_interval = round(p.saveint);
end


%% TIME VARIABLES

t = 0; % time in yr; this is the time at the first timestep
dt = p.dt; % delta t
N = round(p.tf/dt); % number of iterations to do                           
dt_initial = p.dt;  % save this in case you change dt mid-simulation


%% INITIAL CONDITIONS

[K, J, nX] = size(X); % nX is the number of mineral species
B = B - min(B(:)); % make min bedrock elevation (including boundaries) zero


%% Conditions for imposed step change.

% Important note: If you change D, you need to re-run SetUpADI_bdry each
% time D is changed, since the Aplus, Aminus, Bplus, and Bminus terms it
% generates depend directly on D.

switch experiment_type
    
    case 'tectonic'  % consider a perturbation in channel lowering rate
        E_afterstep = E * p.Estepfactor;
        
    case 'climatic'  % consider simultaneous perturbations in D, P, kX
        D_afterstep = D * p.Dstepfactor;
        P_afterstep = P * p.Pstepfactor;
        kx_afterstep = kx * p.kxstepfactor;
        [Aplus_afterstep, Aminus_afterstep, Bplus_afterstep, ...
            Bminus_afterstep] = SetUpADI_bdry(K,J,D_afterstep,dx,dy,dt,C);
        
    case 'D'  % only change soil transport coefficient
        D_afterstep = D * p.Dstepfactor;
        [Aplus_afterstep, Aminus_afterstep, Bplus_afterstep, ...
            Bminus_afterstep] = SetUpADI_bdry(K,J,D_afterstep,dx,dy,dt,C);
        
    case 'P'  % only change soil production coefficient
        P_afterstep = P * p.Pstepfactor;
        
    case 'kx'  % only change dissolution kinetic rate constants
        kx_afterstep = kx * p.kxstepfactor;
        
    case 'steady_state'  % change nothing
        
    otherwise
        error('Choose one of the valid experiment types')

end


%% MATRIX OPERATORS

% Construct matrix operators that calculate curvature in x and y directions
% for use in the Alternating-Direction Implicit (ADI) solution for linear
% diffusion
[Aplus, Aminus, Bplus, Bminus] = SetUpADI_bdry(K,J,D,dx,dy,dt,C);


%% MEMORY ALLOCATION

Xnminus1 = zeros(K,J,nX); % value of X at time n-1 (used in leapfrog scheme)
Hnminus1 = zeros(K,J); % used in Adams-Bashforth
Bnminus1 = zeros(K,J); % used in Adams-Bashforth

Xtemp = zeros(K,J,nX); % these three are used to keep track of X,H,B at time n-1
Htemp = zeros(K,J); 
Btemp = zeros(K,J); 

% Assign values to first timestep
if do_output % if run_name was specified, initialize grids
    
	output.topo_m = zeros(K,J,ceil(N/save_interval)+1); % each 'layer' of
        % this data cube holds the matrix of elevations at time n.
	output.H_m = zeros(K,J,ceil(N/save_interval)+1);  % soil thickness
	output.X = zeros(K,J,nX,ceil(N/save_interval)+1);  % mineral abundances
    output.W_mperyr = zeros(K,J,ceil(N/save_interval)+1);  % chemical erosion rate
    output.R_mperyr = zeros(K,J,ceil(N/save_interval)+1);  % physical erosion rate
    
    % Save values at first save timestep.
    output.topo_m(:,:,1) = B + H;  % topography (m)
    output.H_m(:,:,1) = H;  % soil thickness (m)
    output.X(:,:,:,1) = X;  % soil mineral abundances ()
    output.t_yr(1) = t;  % time since onset of simulation (yr)

    % Compute W at first timestep as sum of chemical erosion fluxes from 
    % all mineral phases.
    temp = zeros(size(X));
    for i = 1:length(kx)
        temp(:,:,i) = kx(i)*Ax(i)*X(:,:,i) - sx(i)*wx(i)/rhos;
    end
    output.W_mperyr(:,:,1) = H .* sum(temp,3);  % Sum over all mineral
        % phases (m/yr)
    clear temp
    
    % Compute R at first timestep using Erode function.
    Hbeforeerosion = H;  % use to infer erosion rate grid
    [Haftererosion, Baftererosion, Xaftererosion] = Erode(H,Hnminus1,B, ...
        Bnminus1,X,Xnminus1,C,nX,K,J,Aplus,Aminus,Bplus,Bminus,D,dt,dx, ...
        dy,xr,1);
    output.R_mperyr(:,:,1) = (Hbeforeerosion - Haftererosion)/dt;  % m/yr 
    clear Hbeforeerosion Haftererosion Baftererosion Xaftererosion

end  


%% MAIN ITERATION LOOP

n_atstep = p.n_atstep;  % timestep at which step occurs
tstep = dt * n_atstep;  % time (yr) at the step iteration
inputs.timeatstep_yr = tstep;  % save time of step (yr)

for n=1:N  % timesteps (N = total number of timesteps)
    
    % This iterative loop computes changes in H, B, and X in a splitting
    % method, which computes portions of the changes in each term in
    % sequential steps rather than all at once. For example, part of the 
    % change in H from one timestep to the next is computed in one function
    % (SoilProd), and the updated value of H is fed into a second function
    % (Weather), which updates H again and feeds the newly updated H into a
    % third function (Erode).  Such splitting methods are favored here
    % because of their higher accuracy over non-splitting methods.
    
    % Impose step perturbation.  Note that if you change D here, this is
    % also where you'll need to reassign values for Aplus, Aminus, Bplus,
    % and Bminus too.
    
    if n == n_atstep
        
        switch experiment_type
            
            case 'tectonic'  % change only channel lowering rate E
                E = E_afterstep;
            
            case 'climatic'  % simultaneously change D, P, and kX
                P = P_afterstep;
                kx = kx_afterstep;
                D = D_afterstep;  % requires changing Aplus Aminus Bplus
                    % Bminus too
                Aplus = Aplus_afterstep;  % required if D changes
                Aminus = Aminus_afterstep;  % required if D changes
                Bplus = Bplus_afterstep;  % required if D changes
                Bminus = Bminus_afterstep;  % required if D changes
            
            case 'D'  % only change D
                D = D_afterstep;  % requires changing Aplus Aminus Bplus
                    % Bminus too
                Aplus = Aplus_afterstep;  % required if D changes
                Aminus = Aminus_afterstep;  % required if D changes
                Bplus = Bplus_afterstep;  % required if D changes
                Bminus = Bminus_afterstep;  % required if D changes

            case 'P'  % only change P
                P = P_afterstep;

            case 'kx'  % only change kx
                kx = kx_afterstep;
                
            case 'steady_state'  % change nothing
        end
        
    end

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Time %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    t = t + dt;
    Htemp = H;
    Btemp = B;
    Xtemp = X;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%%%%%% Soil Production %%%%%%%%%%%%%%%%%%%%%%%%
    
    [H, B, X] = SoilProd(H,B,X,C,nX,alpha,P,rhos,rhor,xr,dt);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        
    %%%%%%%%%%%%%%%%%%%%%%%%%% Weathering %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [H, X] = Weather(H,X,C,nX,kx,Ax,sx,wx,rhos,xr,dt);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    %%%%%%%%%%%%%% Boundary conditions/perturbations %%%%%%%%%%%%%%

    B(C==0) = B(C==0) + E*dt; % bedrock uplift relative to boundaries

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    %%%%%%%%%%%%%%%%%%%%%% Surface Evolution %%%%%%%%%%%%%%%%%%%%%%%

    Hbeforeerosion = H;  % use to infer erosion rate grid
    [H, B, X] = Erode(H,Hnminus1,B,Bnminus1,X,Xnminus1,C,nX,K,J,Aplus, ...
        Aminus,Bplus,Bminus,D,dt,dx,dy,xr,n);
    Haftererosion = H;  % use to infer erosion rate grid
 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    %%%%%%%%%%%%%%%%%%% Store previous timestep %%%%%%%%%%%%%%%%%%%%

    Hnminus1 = Htemp;
    Bnminus1 = Btemp;
    Xnminus1 = Xtemp;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
    
    %%%%%%%%%% Physical and chemical erosion rates %%%%%%%%%%%%%%%%%
    
    % Compute grids of physical (R) and chemical (W) erosion rates (m/yr).  
    % Note that these are thinning rates of the soil, so both would have to
    % be multiplied by rhos to obtain the mass flux out of the soil.

    % Compute W as sum of chemical erosion fluxes from all mineral phases.
    temp = zeros(size(X));
    for i = 1:length(kx)
        temp(:,:,i) = kx(i)*Ax(i)*X(:,:,i) - sx(i)*wx(i)/rhos;
    end
    W = H .* sum(temp,3);  % Sum over all mineral phases (m/yr)
    clear temp
    
    % Physical erosion rate R.
    R = (Hbeforeerosion - Haftererosion)/dt;  % m/yr 
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    %%%%%%%%%%%%%%%%%%%%%%%%% SAVE DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if do_output  % only if we're saving results
        if rem(n,save_interval)==0
            lastsave = n/save_interval+1;  % Start at second save timestep.
            output.topo_m(:,:,lastsave) = B + H;  % topography (m)
            output.H_m(:,:,lastsave) = H;  % soil thickness (m)
            output.X(:,:,:,lastsave) = X;  % soil mineral abundances ()
            output.W_mperyr(:,:,lastsave) = W;  % chemical erosion rate (m/yr)
            output.R_mperyr(:,:,lastsave) = R;  % physical erosion rate (m/yr)
            output.t_yr(lastsave) = n*dt;  % time since onset of simulation (yr)
            output.t_yr(lastsave) = t;  % time since onset of simulation (yr)
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END MAIN ITERATION LOOP % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Compute CDF

CDF = zeros(size(output.H_m));
zirconindex = find(string(p.xname)=='zircon');
for i = 1:size(output.X,4)
    CDF(:,:,i) = 1 - p.xr(zirconindex)./output.X(:,:,zirconindex,i);
end

% Set edges to NaNs for certain variables.
CDF(C==1) = NaN;
output.H_m(C==1) = NaN;
output.R_mperyr(C==1) = NaN;
output.W_mperyr(C==1) = NaN;
output.X(C==1) = NaN;


%% Compute time lags for W, CDF, H, and topography

% Allocate space for time lag grids.
tlag_CDF = zeros(K,J);  % time lag for CDF (years)
tlag_H = zeros(K,J);  % time lag for soil thickness (years)
tlag_topo = zeros(K,J);  % time lag for topography (B+H) (years)
tlag_W = zeros(K,J);  % time lag for chemical erosion rate W (years)
CDF_at_tlag = zeros(K,J);  % CDF at time of time lag ()
H_at_tlag = zeros(K,J);  % soil thickness at time of lag (m)
topo_at_tlag = zeros(K,J);
W_at_tlag = zeros(K,J);

% Compute time lag grids.
for i = 1:size(output.W_mperyr,1)
    for j = 1:size(output.W_mperyr,2)
        tempCDF = squeeze(CDF(i,j,:));  % CDF history at one site (vector)
        tempH = squeeze(output.H_m(i,j,:));  % history of W
        temptopo = squeeze(output.topo_m(i,j,:));  % history of topography
        tempW = squeeze(output.W_mperyr(i,j,:));  % history of W
        [tlag_CDF(i,j), CDF_at_tlag(i,j)] = lagafterstep(output.t_yr, ...
            tempCDF, mark, tstep);
        [tlag_H(i,j), H_at_tlag(i,j)] = lagafterstep(output.t_yr, ...
            tempH, mark, tstep);
        [tlag_topo(i,j), topo_at_tlag(i,j)] = lagafterstep(output.t_yr, ...
            temptopo, mark, tstep);
        [tlag_W(i,j), W_at_tlag(i,j)] = lagafterstep(output.t_yr, ...
            tempW, mark, tstep);
    end
end

% Save time lags
output.tlag_CDF_yr = tlag_CDF;  % years
output.tlag_H_yr = tlag_H;  % years
output.tlag_topo_yr = tlag_topo;  % years
output.tlag_W_yr = tlag_W;  % years


%% Save inputs and outputs.

if do_output % only if we're saving results
    outputname = ['output_' run_name '_' clocktime2text '.mat'];
    save(outputname, 'inputs', 'p', 'output');
end


