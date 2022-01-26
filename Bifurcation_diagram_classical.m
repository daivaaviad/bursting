%
% Bifurcation diagram for the (scalar) classical logistic map.
%
% If you find this code useful, please cite:
%
% D. Petkevičiūtė-Gerlach, R. Šmidtaitė and M. Ragulskis. "Intermittent bursting in the 
% fractional difference logistic map of matrices", Int. J. Bifurcation and Chaos 32 (2022).
%

clear all; 
close all;

N = 2000; % number of iterations
nn = 600; % number of iterations omitted for transient processes
lambda0  = 0.25;

x_bounds = [0, 1];
x_numpoints = 800;
x_step = (x_bounds(2)-x_bounds(1))/x_numpoints; 
xx = x_bounds(1):x_step:x_bounds(2);

a_bounds = exp([0.01, 4]);
a_numpoints = 2000; 
a_step = (a_bounds(2)-a_bounds(1))/a_numpoints;
aa = a_bounds(1)+a_step:a_step:a_bounds(2);

plotvalues = zeros(x_numpoints+1,a_numpoints+1);

i = 1;

 for a = aa
     
      A = log(a);
      i = i+1; 
      
      [lambda, ~] =  seqmu_classical(lambda0,A,N);
      
      for j = nn+1:N
         
          [M,index] = min(abs(xx - lambda(j))) ;
          plotvalues(index,i) = plotvalues(index,i)+1;
          
      end   
 end

figure('Units','normalized','Position',[0.05 0.05 0.95 0.4],'Color',[1 1 1]);   
imshow(ones(x_numpoints+1,a_numpoints+1) - plotvalues,[0 1]);

set(gca,'YDir','normal');  axis on;
set(gca, 'TickLabelInterpreter', 'latex');
set(gca,'FontSize',18); 

tyi  = [0 0.2 0.4 0.6 0.8 1]*x_numpoints;  tyi(1,1) = 1;
yticklabels({'0','0.2','0.4','0.6','0.8','$\lambda_0^{\left(k\right)}$'});
set(gca,'YTick',tyi);

txi = ((exp([0.02 1 2 3 3.6 3.999]) - exp(0))/(exp(4) - exp(0)))*(a_numpoints); 
xticklabels({'0', '1','2','3','3.6','$a$'});
set(gca,'XTick',txi);
                      
