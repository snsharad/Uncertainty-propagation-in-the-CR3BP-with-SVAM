clear;close all;clc

addpath('CUT_fcns')
addpath('CR3BP_fcns')
addpath('SVAM_fcns')

mu = 0.012150463292236;
moon = 1 - mu ;
L1 = 0.836915727639612; 
L2 = 1.155681694999380;

%% Nominal trajectory

% load('TRANSFERS1.mat')

options = odeset('RelTol',1e-12,'AbsTol',1e-14);

%% Halo case

load L1_Halo_Fam_1_new.mat

z0 = PO_states(600,:)';
Tp = PO_time(600);

C0 = jacobiEnergy(z0, mu);
[beta0, gamma0] = vel2angles(z0);
[th0, ph0, r0] = cart2sph(z0(1), z0(2), z0(3));
y0 = [r0; th0; ph0; gamma0; beta0]; % 5-dim vector (r, th, ph, gamma, beta)

% options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
tspan = linspace(0, 2*Tp, 1000);

[~,dx] = ode45(@(t,x) CR3BP(t,x,mu),tspan, z0, options);
[~,yRef] = ode45(@(t,yRef)cr3bpEOMs_SVAM(t,yRef,C0,mu), tspan, y0, options);


%% Define uncertainty in S-VAM
tspan1 = tspan;

N_sample = 5;

delSI = [1e-3 1e-3 1e-3 2 2]'; % [km, deg x 4]
del = [delSI(1)/384400, delSI(2)*pi/180, delSI(3)*pi/180, delSI(4)*pi/180, delSI(5)*pi/180]';

dev = 2 * del .* rand(5, N_sample) - repmat(del, [1,N_sample]); % uniform disbn


y0_samples = zeros(5, N_sample);

for j = 1:N_sample

    y0_samples(:,j) = y0 + dev(:,j);
    
end

parfor j = 1:N_sample
   
   [T,yy] = ode45(@(t,y)cr3bpEOMs_SVAM(t,y,C0,mu), tspan1, y0_samples(:,j), options);
   
   for ct = 1:length(T)
   
       Y{j}(:,ct) = yy(ct,:);
       YC{j}(:,ct) = svam2cart(yy(ct,:), C0, mu);
       
   end    
   
end


%% Reachability
dim = 5; % dimension of the problem
NT = 4; %Total Order of Basis Functions (multidimensional)
N = repmat(4,[dim,1]); %Order of 1-D Basis Functions (4 in each dimension)
Type = {repmat(-1,[dim,1])}; %-1=uniform in all 5 dimensions
Order = 8; % CUT Point accuracy
ind = true; %Include the mean (this is for centralizing during optimization)
Reachability = GetReachability(NT,N,Type,Order,ind);

n = length(Reachability.W);

%% Calculate coefficients

% Scale CUT points

for j = 1:dim
    
    if j==1
        ub = max(y0_samples(j,:));
        lb = min(y0_samples(j,:));
    else
        ub = [ub, max(y0_samples(j,:))];
        lb = [lb, min(y0_samples(j,:))];
    end
end

tmpA = (lb + ub)/2;
tmpB = (ub - lb)/2;
tmp_cut = ((ub - lb)/2)' .* Reachability.Z';

X = repmat(((lb + ub)/2)', [1,n]) + tmp_cut; % X0 in Cartesian

% Evaluate function at every CUT point
for ct = 1:n
    
    xi = X(:,ct);
    [T,Ytmp] = ode45(@(t,y)cr3bpEOMs_SVAM(t,y,C0,mu), tspan1, xi, options);
    
    for i = 1:length(T)
    
        y{i}(:,ct) = Ytmp(i,:)'; % cell index - timestep; each column inside the cell corresponds to a CUT point
        
    end

end


for i = 1:length(T)
   
    Coeff{i} = findCoeff(Reachability, y{i}); % find coefficients for each timestep
        
end

%% SMM

% Scale the perturbed ICs S.T. they belong to [-1,1]. 

parfor j = 1:N_sample
   
    tmpC = y0_samples(:,j)' - tmpA;
    tmpD = tmpB .^ -1;
    MC_phi(j,:) = tmpD .* tmpC;  

end

phi = EvaluateBasis(Reachability.BASIS,Reachability.index, MC_phi);


for ct = 1:length(T)
    
    tmp = Coeff{ct} * phi;
    
    for j = 1:N_sample
        y_cut{j}(:,ct) = tmp(:,j);
        y_cutC{j}(:,ct) = svam2cart(tmp(:,j),C0, mu);
        
        y_cutTimeCell{ct}(:,j) = y_cut{j}(:,ct);
        YTimeCell{ct}(:,j) = Y{j}(:,ct); 
        
        y_cutTimeCellC{ct}(:,j) = y_cutC{j}(:,ct);
        YTimeCellC{ct}(:,j) = YC{j}(:,ct); 
    end

end

% Mean and Covariance - CUT
for ct = 1:length(T)
    [mu_CUT{ct}, cov_CUT{ct}] = MUCOV_CUTpts(y{ct}', Reachability.W);
end

%% Mahalanobis Distance calculation

for ct = 1:length(T)
    
    mcStates = YTimeCell{ct};
    cutStates = y_cutTimeCell{ct};
    
    mcStatesC = YTimeCellC{ct};
    cutStatesC = y_cutTimeCellC{ct};
    
    mcMu(:,ct) = mean(mcStates, 2);
    
    tmp2 = zeros(5,5); % initialize tmp

    for j = 1:N_sample

        tmp2 = tmp2 + (mcStates(:,j) - mu_CUT{ct}') * (mcStates(:,j) -  mu_CUT{ct}')';

    end

    mcCov{ct} = tmp2 ./ (N_sample-1);
    
    % Mahalanobis distance

    for j = 1:N_sample

        tmpMd = mcStates(:,j) - mu_CUT{ct}';
        MahalD(ct,j) = abs(tmpMd' * inv(cov_CUT{ct}) * tmpMd); % Md_squared
        tmpDiff = abs(mcStates(:,j) - cutStates(:,j));
        ErrNorm(ct,j) = norm(tmpDiff);
        ErrPos(ct,j) = sqrt((mcStatesC(1,j) - cutStatesC(1,j)).^2 + (mcStatesC(2,j) - cutStatesC(2,j)).^2 + (mcStatesC(3,j) - cutStatesC(3,j)).^2 );
        ErrVel(ct,j) = sqrt((mcStatesC(4,j) - cutStatesC(4,j)).^2 + (mcStatesC(5,j) - cutStatesC(5,j)).^2 + (mcStatesC(6,j) - cutStatesC(6,j)).^2 );

    end
        
end

