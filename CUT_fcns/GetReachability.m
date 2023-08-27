%% GetReachability.m
%
% Author: Zach Hall
% Date: 4/28/2020
% Email: zyh5059@psu.edu
%
% This function compiles a REACHABILITY structure which contains the
% necessary information to apply the higher order sensitivity matrix (HOSM)
% reachability set computation method
%
% INPUTS:
%   NT: [scalar] Maximum Polynomial Order used in basis (any orders higher than NT
%           will be removed)
%
%   N: [(n x 1) vector] 1D basis function orders. Can construct polynomial
%           basis with varying order for different state dimensions if desired
%
%   Type: [(b x 1) cell vector] where b is the number of different
%           distributions. Type indicates the distribution type associated with 
%           each dimension of the state. -1=Uniform, -2=Gaussian, -3=Spherical
%
%   Order: [(b x 1) vector] Indicates the order of quadrature points (i.e.
%           number of moments replicated) for each distribution type. 
%   
%   mean_ind: [logical] true indicates that the mean will be included in
%           the basis functions, false means that the mean is removed
%
% Example Input: { NT=3,N=[2,2,2],Type=[{[-2,-2]},{[-1]}],Order=[8,4],true }
%         This input  retrieves the reachability set information for a 3
%         dimensional system where quad points x1 and x2 are gaussian
%         distributed and replicate 8th order moments, and x3 is uniformly
%         distributed and replicates 4th order moments.
%         The polynomial basis for this set is 3rd order and composed of
%         2nd order basis functions in each variable
%
%   
%   OUTPUT:
%       REACHABILITY.Phi: [(M x N) matrix] basis functions evaluated at
%        all quadrature points
%
%       REACHABILITY.Basis: [(n x 1) cell vector] Each cell contains a
%       vector of 1D basis function handles for each dimension
%
%       REACHABILITY.index: [(M x n) matrix] Contains dimension indexes for
%       each basis function. An index of [2,4,3] corresponds to a basis
%       function of phi_2(x_1)phi_4(x_2)phi_3(x_3) where phi are 1D basis
%       functions
%
%       REACHABILITY.Z: [(N x n) matrix] All normalized quad points. Will
%       be constructed via tensor product if there are different
%       distribution types
%
%       REACHABILITY.W: [(N x 1) vector] weights for all quad points
%
%       REACHABILITY.Bi [(M x M) matrix] inverted B matrix. B=<Phi,Phi>.
%       Will output error if matrix is singular due to choice of basis
%       fucntions
%
%% Begin Code
function REACHABILITY=GetReachability(NT,N,Type,Order,mean_ind)

%Check inputs
% [~,ind]=unique(Type);
% TensorType=Type(sort(ind));
if length(Type)~=length(Order)
    'Error! The number of different distribution types must equal the number of statistical order specifications (Order)'
    return
end

% Compute basis function index matrix
ND=length(N);
index=GenerateIndex(ND,N+1); 
total_deg=sum(index,2)-ND;
index=index(total_deg<=NT,:);
% Remove mean if desired
if ~mean_ind
    index(1,:)=[];
end
REACHABILITY.index=index;

% Main Reachability Loop
ind=1;
for ct=1:length(Type)
    tempType=Type{ct};
    type=tempType(1); % current distribution type
    order=Order(ct);  % quadrature point order
    type_dim=length(tempType); %dimension of the current distribution
    
    %loop over dimensions in current distribution and get basis functions
    for ct1=1:type_dim
        if type==-1
            BASIS{ind}=GetBasisHandle('legendre',N(ct));
            ind=ind+1;
        elseif type==-2
            BASIS{ind}=GetBasisHandle('hermite',N(ct));
            ind=ind+1;
        elseif type==-3
            %Orthogonal polynomials for spherical distribution not yet included
            BASIS{ind}=GetBasisHandle('monomial',N(ct));
            ind=ind+1;
        end
    end
    
    %Get quad points and weights for current distribution
    if type==-1
        if type_dim==1
            [z{ct},w{ct}]=lgwt(order,-1,1);
        else
            Zi=cut_points_uniform(type_dim,order);
            w{ct}=Zi(:,end);
            z{ct}=Zi(:,1:type_dim);
        end
    elseif type==-2
        if type_dim==1
            [z{ct},w{ct}]=ghwa(order);
        else
            if order==2
                [z{ct},w{ct}]=unscented_points(type_dim);
            else
                Zi=cut_points_gaussian(type_dim,order);
                w{ct}=Zi(:,end);
                z{ct}=Zi(:,1:type_dim);
            end
        end
    elseif type==-3
        [z{ct},w{ct}]=cut_points_spherical(order);
    end
    
    % Tensor product of distributions if necessary
    if ct==1
        Z=z{1};
        W=w{1};
    else
        zz=z{ct};
        ww=w{ct};
        zm=Z;
        wm=W;
        Z=[];W=[];
        for ct1=1:length(ww)
            Z=[Z;[zm,ones(length(wm),1)*zz(ct1,:)]];
            W=[W;wm*ww(ct1)];
        end
    end
end
REACHABILITY.BASIS=BASIS;
REACHABILITY.W=W;
REACHABILITY.Z=Z;

%Evaluate Basis Functions
phi = EvaluateBasis(BASIS,index,Z);
REACHABILITY.phi=phi;

%Determine normal matrix B
for ct=1:length(phi(:,1))
    for ct1=1:length(phi(:,1))
    B(ct,ct1)=sum(W'.*phi(ct,:).*phi(ct1,:));
    B(ct1,ct)=B(ct,ct1);
    end
end
Bi=inv(B);
if any(isnan(Bi))
   'Bi is singular. possible that quad points specified cannot integrate basis functions'
    return
end
REACHABILITY.B = B;
REACHABILITY.Bi=Bi;
end