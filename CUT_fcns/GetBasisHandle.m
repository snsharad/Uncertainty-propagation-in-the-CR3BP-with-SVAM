%% GetBasisHandle.m
% This function returns a cell vector containing function handles. The
% basis functions are 1D and dth order. The type can be specified as;
%   
%type='monomial': returns basis functions 1,x,x^2...,x^d
%type='hermite': reutrns dth order hermite basis functions
%type='legendre': returns dth order legendre basis functions

function basis=GetBasisHandle(type,d)
basis(1)={@(x)x.^0};
basis(2)={@(x)x.^1};
if d>=2
    if strcmp(type,'monomial')
        for ct=3:d+1
            basis{ct}=str2func(strcat('@(x)',sprintf('x.^%i',ct-1)));
        end
    elseif strcmp(type,'hermite')
        syms x
        temp{1}=1;
        temp{2}=x;
        for ct=3:d+1
            temp{ct}=x*temp{ct-1}-(ct-2)*temp{ct-2};
            temp{ct}=collect(temp{ct},x);
            basis{ct}=matlabFunction(temp{ct});
        end
    elseif strcmp(type,'legendre')
        syms x
        temp{1}=1;
        temp{2}=x;
        for ct=3:d+1
            temp{ct}=1/(ct-1)*((2*ct-3)*x*temp{ct-1}-(ct-2)*temp{ct-2});
            temp{ct}=collect(temp{ct},x);
            basis{ct}=matlabFunction(temp{ct});
        end
    elseif strcmp(type,'spherical')
        syms x
        if d>=2
            basis{3}=matlabFunction(1/5*(5*x^2-1));
        end
        if d>=3
            basis{4}=matlabFunction(1/7*(7*x^3-3*x));
        end
        if d>=5
            basis{5}=matlabFunction(1/21*(21*x^4-14*x^2+1));
        end
    else
        disp('Error: Type must be monomial,spherical, hermite or legendre')
    end
end
end
