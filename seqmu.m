%
% A function for computing parameters lambda_k and mu_k 
% of the fractional difference logistic map of matrices 
% with the nilpotent matrix of initial conditions.
%
% Here x is lambda, y is mu (see formula (24) of the article
% D. Petkevičiūtė-Gerlach, R. Šmidtaitė and M. Ragulskis. "Intermittent bursting in the 
% fractional difference logistic map of matrices", Int. J. Bifurcation and Chaos 32 (2022).
%
% Parameters x0 and nu can be set from 0 to 1, 0 < a < 4 
% and n is an integer greater than 1.
%
% If you find this code useful, please cite:
%
% D. Petkevičiūtė-Gerlach, R. Šmidtaitė and M. Ragulskis. "Intermittent bursting in the 
% fractional difference logistic map of matrices", Int. J. Bifurcation and Chaos 32 (2022).
%
% Also see the article for more detailed explanations.
%

function [x, y] = seqmu(x0,a,nu,n)
  
  g = coeff2(nu,n);  
  g0 = 1;
  y0 = 1;
  x = zeros(n,1);
  y = zeros(n,1);
  
  x(1) = a * x0 *(1-x0);
  y(1) = a * y0 *(1-2*x0); 
  
  for k = 2:n
  
      gx = g0 * (a * x(k-1)*(1-x(k-1)) - x(k-1));
      gy = g0 * (a * y(k-1)*(1-2*x(k-1)) - y(k-1));
     
      for j = 2:k-1
         gx = gx + g(j-1) * (a*x(k-j)*(1-x(k-j)) - x(k-j)); 
         gy = gy + g(j-1) * (a*y(k-j)*(1-2*x(k-j)) - y(k-j));
      end    
      gx = gx + g(k-1) * (a*x0*(1-x0) - x0);
      gy = gy + g(k-1) * (a*y0*(1-2*x0) - y0); 

      x(k)= x0 + gx;
      y(k)= y0 + gy;
      
  end

end


