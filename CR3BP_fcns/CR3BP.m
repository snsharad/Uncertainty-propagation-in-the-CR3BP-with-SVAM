function dx = CR3BP(~,x,mu)

    dx = zeros(6,1);

    X = x(1);
    Y = x(2);
    Z = x(3);
    Xd = x(4);
    Yd = x(5);
    Zd = x(6);

    r1 = norm([(X+mu),Y,Z]);
    r2 = norm([(X+mu-1),Y,Z]);

    dx(1) = Xd;
    dx(2) = Yd;
    dx(3) = Zd;
    dx(4) = 2*Yd + X - (1-mu)*(X+mu)/r1^3 - (mu)*(X-1+mu)/r2^3;
    dx(5) = -2*Xd + Y - (1-mu)*(Y)/r1^3 - (mu)*(Y)/r2^3;
    dx(6) =  - (1-mu)*(Z)/r1^3 - (mu)*(Z)/r2^3;


end