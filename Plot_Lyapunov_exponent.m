%
% Plots of the Lyapunov exponent as a function of the parameter "a"
% for the (scalar) classical and the (scalar) fractional difference 
% logistic maps. 
%
% If you find this code useful, please cite:
%
% D. Petkevičiūtė-Gerlach, R. Šmidtaitė and M. Ragulskis. "Intermittent bursting in the 
% fractional difference logistic map of matrices", Int. J. Bifurcation and Chaos 32 (2022).
%
% Also see the article for more detailed explanations.
%


clear all; 
close all;

classical = 0;
fractional = 1;

%----------------------------------------------------------
%------ Classical logistic map ----------------------------
if classical

N = 500; % number or iterations

f = @(x,a) a*x*(1- x); 
f_grad = @(x,a) a*(1-2*x); 

a_numpoints = 1000; 
a_bounds = exp([.01, 4]); 
a_step = (a_bounds(2)-a_bounds(1))/a_numpoints;
aa = a_bounds(1):a_step:a_bounds(2);

jj = 0;

  for a = aa
      
      A = log(a);
      jj = jj+1;
      lambda0  = 0.1;
      L(jj) = 0;
      for ii = 1:N
          
          T(ii) = ii-1;
          L(jj) = L(jj) + log(abs(f_grad(lambda0,A))); 
          lambda0 = f(lambda0,A);  
          
      end     
      L(jj) = L(jj)/N;
      
  end

figure('Units','normalized','Position',[0.05 0.05 0.8 0.4],'Color',[1 1 1]); %[left bottom width height]                       

belowCutoff = L;
belowCutoff(L >= 0) = NaN; 
aboveCutoff = L;

h = plot(a_bounds(1):a_step:a_bounds(2),aboveCutoff,'r',...
       a_bounds(1):a_step:a_bounds(2),belowCutoff,'k');
set(h,{'LineWidth'},{1.2;1.2}) 

hold on
plot(aa, zeros(a_numpoints+1,1),'k', 'LineStyle','--');

set(gca, 'TickLabelInterpreter', 'latex');
set(gca,'FontSize',18);  

yLimits = [-6, 2]; 
ylim(yLimits);
xlim(a_bounds);
    
yticklabels({'-6','-4','-2','0','$L$'});

txi = (exp([0.02 1 2 3 3.6 3.999])); 
xticklabels({'0', '1','2','3','3.6','$a$'});
set(gca,'XTick',txi);
                   
clear L yLimits;  

end

%-----------------------------------------------------------
%------ Fractional difference logistic map -----------------
if fractional

N = 1000;  % number or iterations

a_numpoints = 1000;  
a_bounds = exp([.01, 3.8]); 
a_step = (a_bounds(2) - a_bounds(1))/(a_numpoints-1);
aa = a_bounds(1):a_step:a_bounds(2);

jj = 0;

  for a = aa
      
      A=log(a);
      jj=jj+1;
      lambda0 = 0.1;
      L(jj) = seqxa(lambda0,A,0.8,N);
   
  end


figure('Units','normalized','Position',[0.05 0.05 0.8 0.4],'Color',[1 1 1]); %[left bottom width height]                       
 
belowCutoff = L;
belowCutoff(L >= 0.001) = NaN;  
aboveCutoff = L;

h=plot(aa,aboveCutoff,'r',...
       aa,belowCutoff,'k');
set(h,{'LineWidth'},{1.2;1.2}) 

hold on
plot(aa, zeros(a_numpoints,1),'k', 'LineStyle','--');

yLimits = [-0.2, 0.6]; 
ylim(yLimits);
xlim(a_bounds);

set(gca,'FontSize',18);   
set(gca, 'TickLabelInterpreter', 'latex');

yticklabels({'-0.2','0','0.2','0.4','$L$'});

txi = (exp([0.02 1 2 3 3.38 3.8])); 
xticklabels({'0', '1','2','3','3.38','$a$'});
set(gca,'XTick',txi);

end
%----------------------------------------------------------


function l = seqxa(x0,a,nu,n)
  
  g = coeff2(nu,n);  
  g0 = 1;
  b0 = 1;
  x = zeros(n,1);
  b = zeros(n,1);
  
  x(1) = a * x0*(1-x0);
  b(1) = a * b0*(1-2*x0);
  
  for k = 2:n
  
      gx = g0 * (a * x(k-1)*(1-x(k-1)) - x(k-1));
      gy = g0 * (a * b(k-1)*(1-2*x(k-1)) - b(k-1));
      
      for j = 2:k-1
         gx = gx + g(j-1) * (a*x(k-j)*(1-x(k-j)) - x(k-j));
         gy = gy + g(j-1) * (a*b(k-j)*(1-2*x(k-j)) - b(k-j));
      end    
      gx = gx + g(k-1) * (a*x0*(1-x0) - x0); 
      gy = gy + g(k-1) * (a*b0*(1-2*x0) - b0); 

      x(k)= x0 + gx;
      b(k)= b0 + gy;
          
  end
  
  c = b(end);
  l = 1/(n+1)*log(abs(c));
 

end

