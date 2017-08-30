% Homework 1
clear all; close all;

% meshing parameters
dx = 1;                     % grid size

x = [0 : dx : 100]';        % coordinates of grid points: 0-100
nx = length(x);             % number of grid points

% model parameters
D = 1.0;             % thermal diffusivity

% timing parameters
FACTOR = 0.3;

dt = FACTOR*dx^2/D;         % time step
t = [0 : dt : 25]';         % simulation time: 0-25
nt = length(t);

% initial condition
%T=zeros(nx,1);  T(51)=1;  

% analytical solution
sigma = 1.0;
K = 1.0;
T_analytical = zeros(nx,1);
Tmax = 1.0;

% initial analytical solution
T = zeros(nx,1);
time = 0.0;
for i=1:nx
   T(i) = Tmax/sqrt(1+4*time*K/sigma^2) * exp( - ( x(i) - 50.0)^2/( sigma^2 + 4*time*K) );
end    

%plotting
iplot = 0;
imod = round(nt/10);
ylim('manual');

% time marching
for it = 1:nt               

    % update the temperature field
    T_old = T;              
    
    %> TODO: implement your finite-difference scheme
    % e.g. explicit scheme: 
    %      forward difference in time, centre differences in space 
    for i = 2:nx-1
        T(i) = T_old(i) + ...
    end
    %<TODO
    
    % boundary conditions
    T(1) = 0; T(nx) = 0;      
         
     % analytical solution
     time = (it-1) * dt;
     for i = 1:nx
       T_analytical(i) = Tmax/sqrt(1+4*time*K/sigma^2) * exp( - ( x(i) - 50.0)^2/( sigma^2 + 4*time*K) );
     end    
     
    % figures
     if (it==1 || mod(it-1,imod)==0 )
        iplot = iplot + 1;
        
        if( iplot <= 5)
          subplot(6,2, 2*iplot-1);
          
          axis([0 100 0 1]);
          ylim([0 1]);
          %plot(x,T_old,'b');          
          plot(x,T_old,'b',x,T_analytical,'r');
         
          text(20, 0.5, ['t=' num2str(dt*(it-1))]);
          if (it==1)        
            text(60,0.8,'Finite-difference Scheme');
            text(60,0.4,['dt = ', num2str(dt/(dx^2/D)), '* dx^2/D']);
          end
          
          subplot(6,2,2*iplot);
          axis([0 100 -0.1 0.1]);
          ylim([-0.1 0.1]);
          err = T_old - T_analytical;
          err = err / max(T_analytical);
          plot(x,err,'r');
                    
          %figure(gcf);
          %pause;
        end
    end
end

% end snapshot figure
subplot(6,2, 11);
ylim([0 1]);
plot(x,T_old,'b',x,T_analytical,'r'); 
text(20, 0.5, ['t=' num2str(dt*(nt-1))]);          

subplot(6,2,12);
ylim([-0.001 0.001]);
plot(x,T_old-T_analytical,'r');

figure_number = input('input figure number: \n','s');

%%% pdf format
filename = ['./figures/figure_',figure_number,'.pdf'];
saveas(gcf,filename,'pdf');
%%% eps format
%print(gcf, '-depsc', ['./figures/figure_',figure_number,'.eps']);
disp(['plotted file: ',filename]);
%!rm *.m~

