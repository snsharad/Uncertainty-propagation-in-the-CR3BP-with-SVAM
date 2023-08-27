function [zk,wk]=unscented_points(n)
k=0;
zk=zeros(2*n+1,n);
wk=zeros(2*n+1,1);
I=eye(n);
wk(1)=k/(n+k);
for ct=1:n
    zk(ct+1,:)=sqrt(n+k)*I(ct,:);
    zk(ct+n+1,:)=-sqrt(n+k)*I(ct,:);
    wk(ct+1)=1/(2*(n+k));
    wk(ct+1+n)=1/(2*(n+k));
end
end