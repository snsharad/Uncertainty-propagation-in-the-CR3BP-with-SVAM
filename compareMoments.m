%% Run this section directly to generate moment comparison plot. "statMomPlot.mat" contains data for various cases with different sample sizes.
clear

load statMomPlot.mat

% Define system parameters
LU = 384400;
TU = 375190;

for j = 1:size(mu_MC,2)
    
    muMCnorm(j) = norm(mu_MC(:,j)); 
    covNorm(j) = norm(cov_MC1(:,:,j));
    errCOV(j,1) = abs(covNorm(j) - norm(cov_CUT))';
    errMU(j) = abs(muMCnorm(j) - norm(mu_CUT));
    Emean{1}(:,j) = mu_MC(:,j) - mu_CUT';
    EmeanNorm(j,1) = norm(Emean{1}(:,j));
    
end

figure
hold on
loglog(N_sample, EmeanNorm(:,1), 'b')
loglog(N_sample, errCOV, 'r')
loglog(N_sample, skewErr', 'k')
set(gca, 'XScale', 'log', 'YScale', 'log');
grid on;
xlabel('Number of MC samples')
ylabel('Error between MC and CUT')
legend('$|E[x]|$', '$|E[x^2]|$', '$|E[x^3]|$', 'interpreter', 'latex')
hold off