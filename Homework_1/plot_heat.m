% Homework 1: heat equation
% plots output files from hw1_heat.f90
clear all; close all;

% meshing parameters
dx = 1;                     % grid size
x = [0 : dx : 100]';        % coordinates of grid points: 0-100
nx = length(x);             % number of grid points

% model parameters
D = 1.0;                    % thermal diffusivity

% timing parameters
%FACTOR = 0.6;
%dt = FACTOR * dx^2/D;       % time step
%t = [0 : dt : 25]';         % simulation time: 0-25
%nt = length(t);

% last output time step and nout from hw1_heat.f90 
nmax = 11;
nout = 1;

it = 0;

for i = 1:nout:nmax

  filename = ['figures/',sprintf('HW1_%6.6i',it),'.dat']; 
  disp(['reading file: ',filename]);
  A = load (filename); 
  
  x = A(:,1);   % position
  y1 = A(:,2);  % numerical result
  y2 = A(:,3);  % analytical solution
  y3 = A(:,4);  % error

  if ( (it*2+1)<12)
    disp(['  subplot: ',num2str(it)]);
    % temperature plot
    subplot(6,2, it*2+1);
    plot(x,y1,'b');
    ylim([0 1]);
    %text(20, 0.5, ['t=' num2str(dt*(i-1))]);
    text(20, 0.5, ['it=' num2str(it)]);
    if (it==0)
      text(60,0.8,'Finite-difference Scheme');
    end
    % error plot
    subplot(6,2, it*2+2);
    plot(x,y3,'b');
    if (it==0)
      text(60,0.8,'Error');
    end
  end
      
  %if ( (it*2+1)<13)
  %  subplot(6,2, it*2+2);
  %  plot(x,y2,'b');
  %  ylim([0 1]);
  %  text(20, 0.5, ['t=' num2str(dt*(i-1))]);
  %  if (it==0)
  %    text(60,0.8,'Implicit Scheme');
  %  end
  %end

  it = it + 1;
end

figure_number=input('input figure number: \n','s');

%%% pdf format
filename = ['./figures/figure_',figure_number,'.pdf'];
saveas(gcf,filename,'pdf');
%%% eps format
%print(gcf, '-depsc', ['./figures/figure_',figure_number,'.eps']);
disp(['plotted file: ',filename]);

