function C = jacobiEnergy(x,mu)

    X = x(1);
    Y = x(2);
    Z = x(3);
    Xd = x(4);
    Yd = x(5);
    Zd = x(6);

    d1 = norm([(X + mu),Y,Z]);
    d2 = norm([(X + mu - 1),Y,Z]);

    Omega = 0.5 * (X^2 + Y^2) + (1 - mu)/d1 + mu/d2; 
    C = 2 * Omega - (Xd^2 + Yd^2 + Zd^2);


end