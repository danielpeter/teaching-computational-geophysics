% Homework 5: steady-state heat equation
clear all; close all;

%% meshing parameters
% number of elements
Nel = 10;   % or: 2

% start/end point location
x1 = 0.0; x2 = 1.0; 

% element size
he = (x2-x1)/Nel; 

% number of points
Np = Nel + 1; 
% point locations
x = linspace(x1,x2,Np);

%% boundary condition
T1 = 1;                 % initial temperature at point x = 1
q0 = 1;                 % heat flux at point x = 0

% constant force factor
f_const = 0;      % or: 1
% given force
f = f_const * ones(Np,1);

%% initialization
K = zeros(Nel,Nel);  % stiffness matrix
F = zeros(Nel,1);    % force vector

for e = 1 : Nel
    %% number of local shape functions
    Nen = 2;
    
    %% local to global view
    % sets up global to equation numbering
    ID(1:Nel) = 1:Nel; 
    ID(Np) = 0;
    IEN(1:Nen) = [e,e+1];
    
    %> TODO: Location matrix: setup local to global numbering
    % (entries for this element)
    LM(1:Nen) = ..    
    
    disp(['Location matrix: ',num2str(LM)]);
    
    %> TODO: setup local stiffness matrix and rhs vector
    ke = ..
    fe = ..
    
    %% assembly
    % (add local contribution to the global matrix)    
    ind = find(LM);  % note: find(X) returns only the non-zero elements from vector X

    %> TODO: construct global stiffness matrix
    K(LM(ind),LM(ind)) = ..
    
    % global force vector
    F(LM(ind)) = ..
end

%> TODO: solve linear system: K d = F  for unknowns in d
d = ..

%% exact solution
N = 200;
x_ex = linspace(x1,x2,N);
T_ex = T1 + (1-x_ex)*q0 + (1-x_ex.^2)*f_const/2;

%% plot result
plot([x(1:Nel),x2],[d',T1],'r*-',x_ex,T_ex); 
xlabel('x'); ylabel('T'); 
title(['Nel = ',num2str(Nel),' , f = ',num2str(f_const)]);
legend('FEM solution', 'Exact solution');
