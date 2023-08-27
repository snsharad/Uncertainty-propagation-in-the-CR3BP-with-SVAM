function J = jacobianXYZ(pos, mu)

    syms X Y Z Xd Yd Zd 
    
    d1 = norm([(X + mu), Y, Z]);
    d2 = norm([(X + mu - 1), Y, Z]);
    
    F(1) = Xd;
    F(2) = Yd;
    F(3) = Zd;
    F(4) = 2 * Yd + X - (1 - mu) * (X + mu)/d1^3 - (mu) * (X - 1 +mu)/d2^3;
    F(5) = -2 * Xd + Y - (1 - mu) * Y/d1^3 - mu * Y/d2^3;
    F(6) =  - (1-mu)*(Z)/d1^3 - (mu)*(Z)/d2^3;
    x = [X,Y,Z,Xd,Yd,Zd];

    Jtmp = jacobian(F,x);

    J = double(subs(Jtmp,[X,Y,Z],[pos(1), pos(2), pos(3)]));


end