% Homework 4
clear all; close all;

% meshing parameters
dx = 0.1;                   % grid size
x = [0 : dx : 100-dx]';     % coordinates of grid points: 0-100
nx = length(x);             % number of grid points

% wavenumber increment
dk = 2.0 * pi / nx / dx;

% model parameters
rho = ones(nx,1);
kappa = ones(nx,1);

x_discon = 60;
rho(fix(x_discon/dx):nx) = 1.0;     % rho contrast
kappa(fix(x_discon/dx):nx) = 1.0;   % kappa contrast

% wavespeed square: c^2
c2 = kappa ./ rho;            

% timing parameters
% time step size
% TODO: use different factors to find good stability
dt = 0.025 * min(dx ./ sqrt(c2) / (2.0 * pi));     

% simulation time: 0-50
t = [0 : dt : 50]';        
nt = length(t);

disp(['time step dt = ',num2str(dt)]);
disp(['number of time steps nt = ',num2str(nt)]);

% initial condition
sigma = 1e-1;
u = exp(-sigma*(x-50).^2);
T = kappa .* (-2) .* (x-50) .* u * sigma;
V = zeros(size(T));

V_old1 = V;
V_old2 = V;

T_old1 = T; 
T_old2 = T; 

% time marching
for it = 1:nt               
    % movie
    if mod(it-1,5000) == 0
        disp(['it = ',num2str(it),' out of ',num2str(nt)]);
        subplot(1,1,1);
        plot(x,V,'r');ylim([-1 1]); title('velocity');
        text(20, 0.8, ['t = ' num2str(dt*(it-1))]);
        pause(0.01);
    end

    % velocity
    V_old1 = V;       
    V_old2 = V_old1;  
    
    % stress
    T_old1 = T;
    T_old2 = T_old1;     
    
    %>TODO: implement your pseudo-spectral method
    %       based on a Fourier transform (fft)
    %       (note: if X is real, then Y = fft(X) is conjugate symmetric
    %        and the number of unique points is ceil((n+1)/2)
    ..
    
    % Dirichlet boundary condition
    ..
    
    % Neumann boundary condition
    ..
    
    %<TODO
    
    % figures
    %if (it==1 || mod(it-1,4000)==0 )
    %   it
    %   if( ((it-1)/4000*2+2) <13)
    %       v = (u-u_old1)/dt;
    %       subplot(6,2, (it-1)/4000*2+1);
    %       hold off;
    %       plot(x,v,'b');ylim([-1 1]); title('velocity');
    %       hold on;
    %       plot(x,V,'r');
    %       line([60 60],[-1 1],'Color',[1 0 0]);
    %       text(20, 0.5, ['blue = 2nd and red = 1st']);
    %       subplot(6,2,(it-1)/4000*2+2);
    %       plot(x,V-v,'b');ylim([-0.5 0.5]); title('velocity difference');
    %       line([60 60],[-0.5 0.5],'Color',[1 0 0]);
    %       text(20, 0.25, ['t=' num2str(dt*(it-1))]);
    %   end
    %end

    % stability check
    if max(V) > 1.e3
        error('simulation became unstable and blew up, please reduce time step size');
    end
end

%% pdf figure
figure_number = input('input figure number: \n','s');
filename = ['./figures/figure_',figure_number,'.pdf'];
saveas(gcf,filename,'pdf');
disp(['plotted figure: ',filename]);

% cleanup
% !rm *.m~

%% derivative of example functions to evaluate pseudo-spectral accuracy
if 1 == 1
    disp('derivative comparison example:');
    clf; clear all;
    n = 10000; 
    L = 2.0;
    x = [0 : L*pi*2/n : L*pi*2*(n-1)/n];
    dx = x(2) - x(1);
    dk = 2.0*pi / dx / n;

    disp(['  grid points n = ',num2str(n)]);
    disp(['  dx = ',num2str(dx),' -> dk = ',num2str(dk)]);
    
    % sin-function
    y = sin(x);
    % step-function
    % y(1:n/3) = 0; y(n/3+1:2*n/3) = 1; y(2*n/3+1:n) = 0;
    % triangle-function
    % y(1:n/2) = 1:n/2; y(n/2+1:n) = n/2-1:-1:0;

    % Fourier transform
    fft_y = fft(y); 

    % plot([0:n-1]*dk,abs(fft_y),'b.-');pause;

    % multiplication with i * k
    fft_y(1:n/2+1) = fft_y(1:n/2+1) * i .* [0:n/2] * dk;
    fft_y(n:-1:n/2+2) = conj(fft_y(2:n/2));

    % inverse Fourier transform
    dy_FFT = real(ifft(fft_y));

    % comparison with finite-difference solution (center difference)
    dy_FD = zeros(size(y));
    dy_FD(2:n-1) = (y(3:n)-y(1:n-2))/ (2*dx);
    dy_FD(1) = (y(2)-y(1))/dx;
    dy_FD(n) = (y(n)-y(n-1))/dx;

    % maximum errors
    ymax_FFT = max(dy_FFT-cos(x));
    ymax_FD = max(dy_FD-cos(x));
    disp(['  maximum error: FFT = ',num2str(ymax_FFT),'  FD = ',num2str(ymax_FD)]);

    % initial function plot
    subplot(3,1,1); plot(x,y);title('f(x)');

    % compares pseudo-spectral and finite-difference solution
    %subplot(3,1,2); plot(x,dy_FFT,'r',x,dy_FD,'g');
    %legend('FFT','FD');title('$\frac{df}{dx}(x)$','interpreter','latex');

    % sine-function example
    subplot(3,1,2); plot(x,dy_FFT,'r',x,dy_FD,'g',x,cos(x),'b');
    legend('FFT','FD','ANA');title('$\frac{df}{dx}(x)$','interpreter','latex');
    
    subplot(3,1,3); plot(x,dy_FFT-cos(x),'r',x,dy_FD-cos(x),'g');
    legend('FFT','FD');title('error');
    ymax = max(ymax_FFT,ymax_FD);    
    ylim([-abs(ymax),abs(ymax)]); 
end
