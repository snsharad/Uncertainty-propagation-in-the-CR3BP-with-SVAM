function Y = svam2cart(X, C, mu)

   Y = zeros(6,1); 

   [Y(1), Y(2), Y(3)] = sph2cart(X(2), X(3), X(1));
    
    pot = pseudoPot_SVAM(X, mu);
    
    Y(6) = sqrt(2 * pot - C) * sin(X(5));
    Y(5) = sqrt(2 * pot - C - Y(6)^2) * sin(X(4));
    Y(4) = sqrt(2 * pot - C - Y(6)^2) * cos(X(4));
    

end

