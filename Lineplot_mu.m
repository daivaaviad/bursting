%
% A plot of absolute values of the parameter mu_k 
% of the fractional difference logistic map of matrices 
% with the nilpotent matrix of initial conditions.
%
% If you find this code useful, please cite:
%
% D. Petkevičiūtė-Gerlach, R. Šmidtaitė and M. Ragulskis. "Intermittent bursting in the 
% fractional difference logistic map of matrices", Int. J. Bifurcation and Chaos 32 (2022).
%

nu = 0.8;
a = 3.36;
n = 5000;
lambda0 = 0.12;

[lambda, mu] = seqmu(lambda0,a,nu,n);
plot(abs(mu),'k-', 'linewidth',1.2);
set(gca, 'TickLabelInterpreter', 'latex');
set(gca,'fontsize',18);
txi = [1 1000 2000 3000 4000 4999];
xticks(txi);
xticklabels({'0', '1000','2000','3000','4000','$k$'});
