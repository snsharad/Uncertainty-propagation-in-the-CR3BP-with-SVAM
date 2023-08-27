function omega = pseudoPot_SVAM(z, mu)

    r = z(1);
    theta = z(2);
    phi = z(3);
%     gamma = z(4);
%     beta = z(5);
    
    d1 = sqrt(r^2 + mu^2 + 2 * mu * r * cos(phi) * cos(theta));
    d2 = sqrt(r^2 + (mu - 1)^2 + 2 * r * cos(phi) * cos(theta) * (mu - 1));
    
    omega = 0.5 * r^2 * cos(phi)^2 + (1 - mu)/d1 + mu/d2;
    
end