 % Homework 6: unsteady-state heat equation
 clear all; close all;

 %% meshing parameters
 % number of elements
 Nel = 10;

 % start/end point location
 X1 = 0; X2 = pi/2;

 % regular grid
 He = (X2-X1)/Nel;
 % point locations
 x1(1:Nel) = X1 + ((1:Nel)-1) * He;
 x2(1:Nel) = x1 + He;
 % element size
 he(1:Nel) = x2(1:Nel) - x1(1:Nel);

 % irregular grid
 %Nel2 = Nel/2;
 %He = 2*(X2-X1)/(3*Nel);
 % point locations
 %x1(1:Nel2) = X1 + ((1:Nel2)-1) * 2 * He;
 %x1(Nel2+1:Nel) = x1(Nel2) +2 * He+ ((1:Nel2)-1) * He;
 % element size
 %he(1:Nel2) = 2 * He; he(Nel2+1:Nel) = He;
 %x2(1:Nel) = x1(1:Nel) + he(1:Nel);

 % number of points
 Np = Nel + 1;

 %% boundary condition
 T1 = 1;  % initial temperature at point x = L
 q0 = 0;  % heat flux at point x = 0

 % discretized material parameters
 rhoc(1:Nel) = 1.0;
 ka(1:Nel/2) = 1.0; ka(Nel/2+1:Nel) = 1.0; % constant conductivity

 % constant force factor
 f_const = 0.0;
 % given force
 f = f_const * ones(Np,1);


 %% time marching parameters
 dt = 0.1 * min(he(1:Nel).^2.*rhoc(1:Nel)./ka(1:Nel));
 Ntime = 1000;
 alpha = 0.5;

 disp(['time step size:             dt    = ',num2str(dt)]);
 disp(['number of time steps:       Ntime = ',num2str(Ntime)]);
 disp(['predictor-corrector scheme: alpha = ',num2str(alpha)]);

 %% initialization
 M = zeros(Nel,Nel); % mass (capacity) matrix
 K = zeros(Nel,Nel); % stiffness matrix
 F = zeros(Nel,1);   % force (right-hand-side) vector

 d_dot = zeros(Nel,1);

 % single harmonic initial
 xgrid = x1';
 d = cos(xgrid(1:Nel))+1;

 % time loop
 for itime = 1 : Ntime
   %%>TODO: predictor
   d = ..
   d_dot = ..

   %% solver
   for e = 1: Nel
      %% number of local shape functions
      Nen = 2;

      %% local to global view
      % setup local to global and equation numbering
      ID(1:Nel) = 1:Nel;
      ID(Np) = 0;
      IEN(1:Nen) = [e,e+1];

      %> TODO: Location matrix: setup local to global numbering
      % (entries for this element)
      LM(1:Nen) = ..

      if itime == 1
        disp(['Location matrix: ',num2str(LM)]);
      end
      
      %>TODO: setup local mass, stiffness matrix and rhs vector
      me = ..
      ke = .. % for const. rhoc,ka over element
      fe = ..

      % boundaries
      if e == 1
          fe(1) = ..
      elseif e == Nel
          fe(1) = ..
      end

      %% assembly
      % add local contribution to the global matrix and rhs
      ind = find(LM); % note: find(X) returns only the non-zero elements from vector X
      %> TODO: construct global mass matrix
      M(LM(ind),LM(ind)) = ..
      % global stiffness matrix
      K(LM(ind),LM(ind)) = ..
      % global force vector
      F(LM(ind)) = ..
   end

   %> TODO: solve linear system: M d_dot + K d = F  for unknowns in d
   R = ..
   delta_d_dot = ..

   %%>TODO: corrector
   d = ..
   d_dot = ..

   %% exact solution
   % exact solution setup (problem A)
   N_ex = 200;
   x_ex = linspace(X1,X2,N_ex);
   T_ex = cos(x_ex) * exp(-dt * itime) + T1;

   %% plotting
   dtime = floor((Ntime+1)/4);
   if mod(itime-1,dtime) == 0      
      ca = subplot(2,2,(itime-1)/dtime+1);
      set(ca,'fontsize',12,'linewidth',2);
      plot([xgrid(1:Nel);X2],[d',T1],'r*-')
      hold on;
      plot(x_ex,T_ex);
      xlabel('x'); ylabel('T'); axis([0,pi/2,1,2]);
      title(['time = ',num2str(itime*dt)]); legend('FEM solution', 'Exact solution');
   end
 end % itime
