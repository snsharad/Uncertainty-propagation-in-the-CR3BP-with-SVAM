function [beta,gamma] = vel2angles(z)
    

    v_xy = sqrt(z(4)^2 + z(5)^2);
    
    beta = atan2(z(6), v_xy);
    gamma = atan2(z(5), z(4));
    
end