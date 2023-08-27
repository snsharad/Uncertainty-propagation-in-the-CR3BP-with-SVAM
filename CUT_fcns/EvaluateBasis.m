%% EvaluateBasis.m
%This function evaluates a list of polynomial basis functions at points x
%
% BASIS: (1xn) cell of 1D polynomial basis for each state dimension.
% index: (Mxn) array of basis function orders
% x: (Nxn) array of points to evaluate the basis functions at
%
% PHI (MxN) array of every basis function evaluated at x

function PHI=EvaluateBasis(BASIS,index,x)
n=length(BASIS);
M=length(index(:,1));
N=length(x);
for ct=1:n
    Basis=BASIS{ct};
    for ct1=1:length(Basis)
        basis=Basis{ct1};
        PHI1(ct1,:,ct)=basis(x(:,ct)');
    end
end
PHI=ones(M,N);
for ct=1:M
    for ct1=1:n
        PHI(ct,:)=PHI(ct,:).*PHI1(index(ct,ct1),:,ct1);
    end
end
end