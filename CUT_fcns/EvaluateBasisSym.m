%% EvaluateBasis.m
%This function evaluates a list of polynomial basis functions at points x
%
% BASIS: (1xn) cell of 1D polynomial basis for each state dimension.
% index: (Mxn) array of basis function orders
% x: (Nxn) array of points to evaluate the basis functions at
%
% PHI (MxN) array of every basis function evaluated at x

function phiSym = EvaluateBasisSym(BASIS,index,xSym)
n=length(BASIS);
M=length(index(:,1));
N=size(xSym,1);
for ct=1:n
    Basis=BASIS{ct};
    for ct1=1:length(Basis)
        basis=Basis{ct1};
%         PHI1(ct1,:,ct)=basis(x(:,ct)');
        PHI1sym(ct1,:,ct) = basis(xSym(ct));
    end
end
PHI=sym(ones(M,N));
for ct=1:M
    for ct1=1:n
        PHI(ct,:)=PHI(ct,:).*PHI1sym(index(ct,ct1),:,ct1);
    end
end
phiSym = PHI;
end