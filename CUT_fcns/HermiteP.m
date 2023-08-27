function X=HermiteP(n,x,Case)

if ~exist('Case','var')
    Case = 'Prob';
end

lim2=floor(n/2);
sum=0;
for m=1:lim2+1
    if strcmp(Case,'Prob')
    sum=sum+(-1)^(m-1)/(factorial(m-1)*factorial(n-2*(m-1))).*(x.^(n-2*(m-1)))/2^(m-1);
    elseif strcmp(Case,'Phys')
    sum=sum+(-1)^(m-1)/(factorial(m-1)*factorial(n-2*(m-1))).*(2*x).^(n-2*(m-1));
    end
end
X=factorial(n)*sum;
end