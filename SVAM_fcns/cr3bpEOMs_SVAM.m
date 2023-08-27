function zdot = cr3bpEOMs_SVAM(~,z,C,mu) 
    
    pot = pseudoPot_SVAM(z, mu);
    r = z(1);
    theta = z(2);
    phi = z(3);
    gamma = z(4);
    beta = z(5);
    
    d1 = sqrt(r^2 + mu^2 + 2 * mu * r * cos(phi) * cos(theta));
    d2 = sqrt(r^2 + (mu - 1)^2 + 2 * r * cos(phi) * cos(theta) * (mu - 1));

    zdot = zeros(5,1);
    
    % Define some constants
    A = (1 - mu) / d1^3;
    B = mu / d2^3;
    tmp = sqrt(2 * pot - C);
    tmp1 = r * cos(phi) * cos(theta) - A * (r * cos(phi) * cos(theta) + mu) - B * (r * cos(phi) * cos(theta) - 1 + mu);
    
    zdot(1) = tmp * (cos(phi) * cos(beta) * cos(gamma - theta) + sin(phi) * sin(beta));
    zdot(2) = tmp / (r * cos(phi)) * cos(beta) * sin(gamma - theta);
    zdot(3) = tmp / r * (sin(beta) * cos(phi) - sin(phi) * cos(beta) * cos(gamma - theta));
    zdot(4) = 1 / (tmp * cos(beta)) * (r * cos(phi) * sin(theta) * cos(gamma) * (1 - A - B) - sin(gamma) * tmp1) - 2;
    zdot(5) = 1 / tmp * (-r * sin(phi) * cos(beta) * (A + B) - cos(gamma) * sin(beta) * tmp1 - r * cos(phi) * sin(theta) * sin(gamma) * sin(beta) * (1 - A - B));

end