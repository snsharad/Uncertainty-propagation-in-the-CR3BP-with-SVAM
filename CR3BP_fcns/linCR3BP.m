function xdot = linCR3BP(~, z, x, mu)

%     xdot = zeros(6,1);

    r1 = norm([(x(1) + mu), x(2), x(3)]);
    r2 = norm([(x(1) + mu - 1), x(2), x(3)]);
    
    Uxx = 1 - ((1-mu)/r1^3) - mu/r2^3 + (3*(1-mu)*(x(1)+mu)^2)/r1^5 + (3*mu*(x(1)+mu-1)^2)/r2^5;
    Uxy = (3*(1-mu)*(x(1)+mu)*x(2))/r1^5 + (3*mu*(x(1) + mu -1)*x(2))/r2^5;
    Uxz = (3*x(3)*(1-mu)*(x(1)+mu))/r1^5 + (3*mu*x(3)*(x(1)+mu-1))/r2^5;
    Uyx = Uxy;
    Uyy = 1 - (1-mu)/r1^3 - mu/r2^3 + (3*(1-mu)*x(2)^2)/r1^5 + (3*mu*x(2)^2)/r2^5;
    Uyz = (3*x(3)*(1-mu)*x(2))/r1^5 + (3*mu*x(3)*x(2))/r2^5;
    Uzx = Uxz;
    Uzy = Uyz;
    Uzz = -(1-mu)/r1^3 - mu/r2^3 + (3*(1-mu)*x(3)^2)/r1^5 + (3*mu*x(3)^2)/r2^5;
    
    C = [Uxx Uxy Uxz; Uyx Uyy Uyz; Uzx Uzy Uzz];
    a = zeros(3);
    b = eye(3);
    K = [0 2 0;-2 0 0; 0 0 0];

    A = [a b;C K];
    
    xdot = A * z;

end