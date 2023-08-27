function [C, A] = findCoeff(Reachability, y)
    
    %Compute Coefficients
    A = y*diag(Reachability.W)*Reachability.phi';
    C = A*Reachability.Bi; %Note B matrix is computed and inverted automatically inside GetReachability Function

end