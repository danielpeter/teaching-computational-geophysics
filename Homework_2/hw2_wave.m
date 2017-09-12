% Homework 2
clear all; close all;

% meshing parameters
dx = 0.1;                   % grid size
x = [0 : dx : 100]';        % coordinates of grid points: 0-100
nx = length(x);             % number of grid points

% model parameters
rho = ones(nx,1);
kappa = ones(nx,1);

% material contrast
%x_discon = 60;
%rho(fix(x_discon/dx):nx) = 1.0;    % rho contrast
%kappa(fix(x_discon/dx):nx) = 4.0;  % kappa contrast

c2 = kappa./rho;                   % wavespeed square: c^2

% timing parameters
FACTOR = 0.25;
dt = FACTOR * min(dx^2./c2);        % time step
t = [0 : dt : 200]';            % simulation time: 0-200
nt = length(t);

disp(['time step dt = ',num2str(dt)]);
disp(['number of time steps nt = ',num2str(nt)]);

% initial condition
sigma = 0.1;
% second-order system
u = exp(-sigma*(x-50).^2);        % displacement
v = zeros(size(u));               % velocity 
u_old1 = u; 

% first-order system
T = kappa.*(-2).*(x-50).*u*sigma; % stress
V = zeros(size(T));               % velocity 
T_old = T; V_old = V;

% time marching
for it = 1:nt               
  % movie
  if mod(it-1,200) == 0
    v = (u-u_old1)/dt;
    subplot(2,2,1);
    plot(x,u,'b');
    ylim([-1 1]); 
    title('displacement ( 2nd-order equation  )');    
    text(20, 0.5, [num2str(it) 't=' num2str(dt*(it-1))]);
    
    subplot(2,2,2);
    hold off;
    plot(x,v,'b');
    ylim([-1 1]); title('velocity');
    hold on;
    plot(x,V,'r');
    legend('2nd','1st');
    
    subplot(2,2,3);
    plot(x,T,'b');
    ylim([-1 1]); 
    title('stress       ( 1st-order equations )');

    subplot(2,2,4);
    plot(x,V-v,'b');
    ylim([-1 1]); 
    title('velocity difference');
    pause(0.01);
  end
  
  %%%%%%%%%%%%%%%%%%%%%%% 2nd order equation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %> TODO: implement your finite-difference scheme
  ..
  
  % Dirichlet boundary condition
  ..
  
  % Neumann boundary condition
  ..
  %<TODO
  
  %%%%%%%%%%%%%%%%%%%%%%% 1st order equations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %> TODO: implement your finite-difference scheme
  ..
  
  % Dirichlet boundary condition
  ..
  % Neumann boundary condition
  ..
  %<TODO
  
  % figures
  %if (it==1 || mod(it-1,4000)==0 )
  %    it
  %    if( ((it-1)/4000*2+2) <13)
  %        figure(1);
  %        v=(u-u_old1)/dt;
  %        subplot(6,2, (it-1)/4000*2+1);hold off;
  %        plot(x,v,'b');ylim([-1 1]); title('velocity');
  %        hold on;plot(x,V,'r');line([60 60],[-1 1],'Color',[1 0 0]);
  %        text(20, 0.5, ['blue = 2nd and red = 1st']);
  %        subplot(6,2,(it-1)/4000*2+2);plot(x,V-v,'b');ylim([-0.5 0.5]); title('velocity difference');
  %        line([60 60],[-0.5 0.5],'Color',[1 0 0]);
  %        text(20, 0.25, ['t=' num2str(dt*(it-1))]);
  %    elseif( ((it-1)/4000*2+2) <25)
  %        figure(2);
  %        v=(u-u_old1)/dt;
  %        subplot(6,2, (it-1)/4000*2+1-12);hold off;
  %        plot(x,v,'b');ylim([-1 1]); title('velocity');
  %        hold on;plot(x,V,'r');line([60 60],[-1 1],'Color',[1 0 0]);
  %        text(20, 0.5, ['blue = 2nd and red = 1st']);
  %        subplot(6,2,(it-1)/4000*2+2-12);plot(x,V-v,'b');ylim([-0.5 0.5]); title('velocity difference');
  %        line([60 60],[-0.5 0.5],'Color',[1 0 0]);
  %        text(20, 0.25, ['t=' num2str(dt*(it-1))]);
  %    elseif( ((it-1)/4000*2+2) <37)
  %        figure(3);
  %        v=(u-u_old1)/dt;
  %        subplot(6,2, (it-1)/4000*2+1-24);hold off;
  %        plot(x,v,'b');ylim([-1 1]); title('velocity');
  %        hold on;plot(x,V,'r');line([60 60],[-1 1],'Color',[1 0 0]);
  %        text(20, 0.5, ['blue = 2nd and red = 1st']);
  %        subplot(6,2,(it-1)/4000*2+2-24);plot(x,V-v,'b');ylim([-0.5 0.5]); title('velocity difference');
  %        line([60 60],[-0.5 0.5],'Color',[1 0 0]);
  %        text(20, 0.25, ['t=' num2str(dt*(it-1))]);
  %    end
  %end

end

%figure(1);saveas( gcf, ['./figures/figure_1.pdf'],'pdf');
%figure(2);saveas( gcf, ['./figures/figure_2.pdf'],'pdf');
%figure(3);saveas( gcf, ['./figures/figure_1.pdf'],'pdf');
% disp(['plotted files in figures/ folder']);
% !rm *.m~



