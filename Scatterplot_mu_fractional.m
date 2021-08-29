%
% Scatter plot of the parameter mu for the fractional difference logistic map
% of matrices with a nilpotent matrix of initial conditions,
% as a function of the parameter "a".
%
% If you find this code useful, please cite:
%
% D. Petkevičiūtė-Gerlach, R. Šmidtaitė and M. Ragulskis. "Intermittent bursting in the 
% fractional difference logistic map of matrices", 2021.
%

clear all;
close all;

N = 1000; 
lambda0  = 0.1;
nu = 0.8;

x_bounds = [-100, 100]; 
x_numpoints = 600; 
x_step = (x_bounds(2)-x_bounds(1))/x_numpoints;
xx = x_bounds(1):x_step:x_bounds(2); 

a_bounds = exp([0.01, 3.8]);
a_numpoints = 2000; 
a_step = (a_bounds(2)-a_bounds(1))/a_numpoints;
aa = a_bounds(1)+a_step:a_step:a_bounds(2);

plotvalues = zeros(x_numpoints+1,a_numpoints+1); 

i = 1;
 for a = aa 
     
      A = log(a);
      i = i+1;
      
      [~, mu] =  seqmu(lambda0,A,nu,N);
 
      for j=1:N
           
          [M,index] = min(abs(xx-mu(j))) ;
          plotvalues(index,i) = plotvalues(index,i)+1;
          
      end   
 end

figure('Units','normalized','Position',[0.05 0.05 0.95 0.4],'Color',[1 1 1]);
imshow(ones(x_numpoints+1,a_numpoints+1)-plotvalues,[0 1]);

set(gca,'YDir','normal');  axis on;
set(gca, 'TickLabelInterpreter', 'latex');
set(gca,'FontSize',18); 

tyi  = ([-100 -50 0 50 100]+[100 100 100 100 100])/200*x_numpoints;  
tyi(1,1) = 1;
yticklabels({'-100','-50','0','50','$\mu^{\left(k\right)}$'});
set(gca,'YTick',tyi);

txi = ((exp([0 1 2 3 3.38 3.7999])-exp(0))/(exp(3.7999)-exp(0)))*(a_numpoints); % txi(1,1) = 1; % - sita reiksme vaizduoja 0 ant x asies
xticklabels({'0', '1','2','3','3.38','$a$'});
set(gca,'XTick',txi);

