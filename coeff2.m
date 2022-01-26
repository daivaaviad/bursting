%
% Coefficients G for the (scalar) fractional difference 
% logistic map.
%
% If you find this code useful, please cite:
%
% D. Petkevičiūtė-Gerlach, R. Šmidtaitė and M. Ragulskis. "Intermittent bursting in the 
% fractional difference logistic map of matrices", Int. J. Bifurcation and Chaos 32 (2022).
%


function g = coeff2(nu,n)
  
  g = zeros(n,1);
  g(1) = nu;
  
  for j=2:n
    g(j) = (1-(1-nu)/j)*g(j-1);  
  end
end
