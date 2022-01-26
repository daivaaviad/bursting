%
% A function for computing parameters lambda_k and mu_k 
% of the classical logistic map of matrices 
% with the nilpotent matrix of initial conditions.
%
% Here x is lambda, y is mu (see formula (12) of the article
% D. Petkevičiūtė-Gerlach, R. Šmidtaitė and M. Ragulskis. "Intermittent 
% bursting in the fractional difference logistic map of matrices", 
% Int. J. Bifurcation and Chaos 32 (2022).)
%
% Parameter x0 can be set from 0 to 1 and 0 < a < 4.
%
% If you find this code useful, please cite:
%
% D. Petkevičiūtė-Gerlach, R. Šmidtaitė and M. Ragulskis. "Intermittent 
% bursting in the fractional difference logistic map of matrices", 
% Int. J. Bifurcation and Chaos 32 (2022).
%
% Also see the article for more detailed explanations.
%


function [x, y] = seqmu_classical(x0,a,n)
  
  y0 = 1;
  x = zeros(n,1);
  y = zeros(n,1);
  
  x(1) = a * x0 *(1-x0);
  y(1) = a * y0 *(1-2*x0); 
  
  for k = 2:n
  
    x(k) = a * x(k-1) *(1-x(k-1));
    y(k) = a * y(k-1) *(1-2*x(k-1)); 
      
  end

end


