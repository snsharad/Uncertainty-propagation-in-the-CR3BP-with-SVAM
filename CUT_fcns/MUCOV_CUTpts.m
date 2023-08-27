function [mu, cov] = MUCOV_CUTpts( x1, x_W)

    n = length(x_W);

    mu = x_W' * x1;
    
    d = length(mu);

    X = x1 - mu;

    temp = x_W.*X;

    cov = temp'*X;

    

end


    


