% Homework 2: wave equation
% plots output files from hw2_wave.f90
clear all; close all;

% meshing parameters
Nel = 100;                  % number of elements
dx = 10;                    % grid size
x = [0 : dx : dx*(Nel-1)]'; % coordinates of grid points: 0-100
nx = length(x);             % number of grid points

dt = 0.5e-6;                % time step

% program starts here
it = 0;
for i = 0:3000:300000
  it = it+1
  %seismograms
  filename = ['figures/',sprintf('S_%6.6i',it),'.dat'];
  disp(['reading file: ',filename]);
  A = load (filename); 
  
  disp(['  subplot: ',num2str(it)]);
  subplot(7,1, it);
  plot(x,A,'r');
  ylim([0 1]);
  text(20, 0.5, ['t=' num2str(dt*(it-1))]);
end

figure_number=input('input figure number: \n','s');

%%% pdf format
filename = ['./figures/figure_',figure_number,'.pdf'];
saveas(gcf,filename,'pdf');
%%% eps format
%print(gcf, '-depsc', ['./figures/figure_',figure_number,'.eps']);
disp(['plotted file: ',filename]);
