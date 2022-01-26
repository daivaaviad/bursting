%
% A color plot of absolute values of the parameter sequence mu_k 
% of the fractional difference logistic map of matrices 
% with the nilpotent matrix of initial conditions,
% as a function of the parameter lambda_0.
%
% If you find this code useful, please cite:
%
% D. Petkevičiūtė-Gerlach, R. Šmidtaitė and M. Ragulskis. "Intermittent bursting in the 
% fractional difference logistic map of matrices", Int. J. Bifurcation and Chaos 32 (2022).
%

close all;
clear all;

i = 0; j = 0;

n =  5000;
xx = zeros(100, n);

a = 3.378;
nu = 0.8;


for lambda0 = 0.001:0.001:1 
    
   i = i+1; 
   j = 0;
   [lambda, mu] = seqmu(lambda0, a, nu, n+102);
   
   for ind1 = 1:n-100 
       
     j=j+1;
     xx(i,j) = max(log10(abs(mu(ind1:ind1+9)))); 
    
   end
end    


figure; hold on;
imagesc(xx)
cb = colorbar;
colormap(jet);
caxis([-4 6]);
set(cb, 'TickLabelInterpreter', 'latex');
set(gca, 'TickLabelInterpreter', 'latex');

xlim([1,4900])
txi = [1 1000 2000 3000 4000 4900];
xticks(txi);
xticklabels({'0', '1000','2000','3000','4000','$k$'});

ylim([1,1000]);
yticks([1 500 999]);
yticklabels({'0','0.5','$\lambda_0^{(0)}$'});

set(gca,'fontsize',18);

