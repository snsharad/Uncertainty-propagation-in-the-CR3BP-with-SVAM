function [phi,bindex,InnerProd,xind,wind] = ComputeInnerProd(N,m,ND,Np)

%
% Author: Puneet Singla.
%Modified by Zach Hall
% Last Modified: Mar 3, 2017.
%
% This function computes the orthogonal polynomials and their inner
% products in n-dim parameter space.
%
% Input Variables:
% N is the maximum order of polynomials along each direction and m detrmines the pdf.
% m = -2 corresponds to Gaussin pdf of Zero mean and variance 1
% m = -1 corresponds to uniform pdf over [-1,1]
% m >=0 corresponds to GLOMAP pdfs
% ND is number of uncertain parameters, i.e., our basis function will lie in
% ND dim space.
% Np order of quadrature points
%
% Output Variables:
% phip is Ninteg x numbasis matrix of polynomial values 
% bindex is numbasis x ND matrix of all combination of 1D poly
% InnerProd is numbasis x numbasis x numbasis tensor of basis function inner product.
% xint is Ninteg x ND matrix of quadrature points
% wint,pw are Ninteg x 1 vectors of quadrature weights and pdf evaluations, respectively.
%

%%
% Generate indicies matrix for ND basis functions
%%
numbasis = repmat(N+1, 1, ND); % number of basis functions along each direction.

index = GenerateIndex(ND,numbasis); % permutation of 1-D basis functions

total_degree = sum(index-ones(size(index)),2);

ind = find(total_degree <= max(N));

bindex = index(ind,:); % bindex contains indicies for N-D basis functions.

%%
% Generate Inner Product Matrices
%%
NB = size(bindex,1); % total number of basis functions in ND space
if m == -2
    Xw =cut_points_gaussian(ND,Np); 
    xind=Xw(:,1:ND);wind=Xw(:,end);
else if m == -1
       Xw =cut_points_uniform(ND,Np); 
       xind=Xw(:,1:ND);wind=Xw(:,end);
    end
end

[phia,phi,IP] = GramSchmidt(max(N),m,ND,Np,bindex,xind,true); % inner products etc for 1-D basis functions

p3ind = GenerateIndex(3,[NB NB NB]); % p3ind contain indicies for product of 3 basis fcns.

for ct  = 1:size(p3ind,1)
    i = p3ind(ct,1); j = p3ind(ct,2); k = p3ind(ct,3);
    b1 = bindex(i,:); b2 = bindex(j,:); b3 = bindex(k,:);
    
    InnerProd(i,j,k) = 1;
    
    for ct2 = 1:length(b1)
        dprod = IP(b1(ct2),b2(ct2),b3(ct2));
        if dprod == 0
            InnerProd(i,j,k) = 0;
            break;
        end
        InnerProd(i,j,k) = InnerProd(i,j,k)*dprod;
    end
end


    