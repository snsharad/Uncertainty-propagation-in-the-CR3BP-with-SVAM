function v = angles2vel(z, C, mu)
    
%    z(4): gamma
%    z(5): beta
     
    v = zeros(3,1);
    
    pot = pseudoPot(z, mu);
    
    v(3) = sqrt(2 * pot - C) * sin(z(5));
    v(2) = sqrt(2 * pot - C - v(3)^2) * sin(z(4));
    v(1) = sqrt(2 * pot - C - v(3)^2) * cos(z(4));
    
end
