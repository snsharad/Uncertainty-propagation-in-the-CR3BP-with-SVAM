%This function gives the weights and abscissa values for the gauss hermite
%polynomials of order n
%
%Zach Hall
%5/16/2018
%
function [xi,w]=ghwa(n,Case)

if ~exist('Case','var')
    Case = 'Prob';
end

numel=ceil((n+1)/2);
COEFF=zeros(1,n+1);
for m=1:numel
    if strcmp(Case,'Phys')
        if mod(n,2)==0
            COEFF(n+1-2*(m-1))=factorial(n)*(-1)^(m-1)/(factorial(m-1)*factorial(n-2*(m-1)))*2^(n-2*(m-1));
        elseif mod(n,2)==1
            COEFF(n+1-2*(m-1))=factorial(n)*(-1)^(m-1)/(factorial(m-1)*factorial(n-2*(m-1)))*2^(n-2*(m-1));
        end
    elseif strcmp(Case,'Prob')
        if mod(n,2)==0
            COEFF(n+1-2*(m-1))=factorial(n)*(-1)^(m-1)/(factorial(m-1)*factorial(n-2*(m-1)))/2^(m-1);
        elseif mod(n,2)==1
            COEFF(n+1-2*(m-1))=factorial(n)*(-1)^(m-1)/(factorial(m-1)*factorial(n-2*(m-1)))/2^(m-1);
        end
    end
end
COEFF=fliplr(COEFF);
xi=roots(COEFF);
w=2.^(n-1)*factorial(n)*sqrt(pi)./(n^2*HermiteP(n-1,xi,Case).^2);
w=w/sum(w);
end
