function [L] = Lagrange_Pts(mu)

syms X 

d1 = norm(X + mu);
d2 = norm(X + mu - 1);

f = X - (1 - mu) * (X + mu)/d1^3 - (mu) * (X - 1 + mu)/d2^3;


% Finding colinear points

assume((X + mu) > 0 & (X - (1 - mu)) < 0) 
L1 = vpasolve(f == 0, X);

assume((X + mu) > 0 & (X - (1 - mu)) > 0)
L2 = vpasolve(f == 0, X,1); % For some reason, was not able to find L2 without initial guess, even with the assume function prescribing the ranges

assume((X + mu) < 0 & (X - (1 - mu)) < 0)
L3 = vpasolve(f == 0, X);

% L4 & L5 - from geometry

L4x = -mu + 1/2;
L5x = L4x;
L4y = sqrt(3)/2;
L5y = -L4y;

L = [L1, 0, 0; L2, 0, 0; L3, 0, 0; L4x, L4y, 0; L5x, L5y,0];
    
end