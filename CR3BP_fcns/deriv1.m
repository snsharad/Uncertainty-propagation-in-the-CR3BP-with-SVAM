function dx = deriv1(~,X, mu)
  
    x = X(1:6);
    Phi = X(7:42);
    Phi = reshape(Phi,6,6);

    Oxy = partialxy(x, mu);

        Ox = Oxy(1);
        Oy = Oxy(2);
        Oz = Oxy(3);

    %     Equation (1) in KC Howell's paper
        dx(1) = x(4);
        dx(2) = x(5);
        dx(3) = x(6);
        dx(4) = Ox + 2 * x(5);
        dx(5) = Oy - 2 * x(4);
        dx(6) = Oz;

    r1 = sqrt((x(1)+mu)^2 + x(2)^2 + x(3)^2);
    r2 = sqrt((x(1)+mu-1)^2 + x(2)^2 + x(3)^2);

    C11 = 1 - ((1-mu)/r1^3) - mu/r2^3 + (3*(1-mu)*(x(1)+mu)^2)/r1^5 + (3*mu*(x(1)+mu-1)^2)/r2^5;
    C12 = (3*(1-mu)*(x(1)+mu)*x(2))/r1^5 + (3*mu*(x(1) + mu -1)*x(2))/r2^5;
    C13 = (3*x(3)*(1-mu)*(x(1)+mu))/r1^5 + (3*mu*x(3)*(x(1)+mu-1))/r2^5;
    C21 = C12;
    C22 = 1 - (1-mu)/r1^3 - mu/r2^3 + (3*(1-mu)*x(2)^2)/r1^5 + (3*mu*x(2)^2)/r2^5;
    C23 = (3*x(3)*(1-mu)*x(2))/r1^5 + (3*mu*x(3)*x(2))/r2^5;
    C31 = C13;
    C32 = C23;
    C33 = -(1-mu)/r1^3 - mu/r2^3 + (3*(1-mu)*x(3)^2)/r1^5 + (3*mu*x(3)^2)/r2^5;

    C = [C11 C12 C13; C21 C22 C23; C31 C32 C33];
    a = zeros(3);
    b = eye(3);
    K = [0 2 0;-2 0 0; 0 0 0];

    A = [a b;C K];



    Phi_dot = A*Phi;
    Phi_dot = reshape(Phi_dot,36,1);



    for i = 7:42
        dx(i) = Phi_dot(i-6);

    end
    % dx = [dx;Phi_dot];
     dx = dx(:);
end

