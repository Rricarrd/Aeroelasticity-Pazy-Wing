
%% PAZY WING PROBLEM
clc
clear
%% Calculate Wing Stiffness
% K Matrix Parameters
if exist('EI','var') ~= 1
    [GJ,EI,xs] = CalcStiffness();
end
p.EI = EI;
p.GJ = GJ;
p.xs = xs;

% Material properties
p.rhoAl = 2795; % [kg/m^3] Density of Aluminium
p.rhoNyl = 930; % [kg/m^3] Density of Nylon 12
p.z = 0.18;           % [m] Maximum thickness
p.xi= 0.0075;         % [m] Initial x position
p.xf = 0.095;         % [m] Final x position
p.c = 0.1;            % [m] Chord
p.b = 0.06;           % [m] Aluminium bar width
p.t = 0.00225;        % [m] Aluminium bar thickness 
p.bi = 0.02;          % [m] Aluminium bar front edge position
p.h = 0.004;          % [m] Depth
p.a_slot = 0.005;     % [m] Plastic slot width
p.x1 = 0.035;         % [m] Plastic slot start position
p.x2 = 0.065;         % [m] Plastic slot end position

% Fluid properties
p.rho = 1.225;


%% Numerical discretization
% Parameters
p.n = 14;              % Panel node number
p.Rs = 0.0384;         % [m] Rib spacing
p.Lw = 0.550;          % [m] Wing lenght
p.Ls = p.Lw - 13*p.Rs; % Last wing segment space
p.nj = 3;              % DOFs per node (theta, gamma i eta)
p.ng = 1;              % Gauss points per element
p.en = 2;              % Element nodes
p.Rt = 0.004;          % [m] Rib thickness

p.sparse = 1;

% Computing beam mass parameters
[p.mu,p.r0,p.rS,p.d] = computeBeamMassParameters(p);

% Discretization and conectivities
[y,Tn,Tr,p.Nnode,p.Nelem,p.Ndof] = DiscretizeWing(p);


%% Element mass matrices

M = AssemblyM(y,Tn,Tr,p);

K = AssemblyK(y,Tn,Tr,p);


%% Free vibrations
% Fixing first nodes
K = K(4:end,4:end);
M = M(4:end,4:end);



k = 10;
[Q,W] = eigs(K,M,k,'sm'); % Solve for eigenvalues

f = sqrt(diag(W))/(2*pi);


%% Aerodynamic elemental matrices
% Elemental matrices
[A0k, A1kc, A1knc] = SymbolicAerodynamicMatrices(p); % Matrices divided by U, hk and Ck

%% Aerodynamic global matrices
% Assembly of matrices (Without U and Ck)
A0 = AssemblyA(A0k, y, Tn, p);
A1c = AssemblyA(A1kc, y, Tn, p);
A1nc = AssemblyA(A1knc, y, Tn, p);

% Excluding BCs
A0 = A0(4:end,4:end);
A1c = A1c(4:end,4:end);
A1nc = A1nc(4:end,4:end);



%%%%%%%%%%%%%%%%%%%%%% FALTA REDUCCIÓ D'ORDRE %%%%%%%%%%%%%%%%%%%%%%



%% Divergence
[Ud] = Divergence(K,A0,50);


%% Flutter p method
[U_,p_,Vp_] = FlutterPMethod(p,K,M,A0,A1c,A1nc);

figure
for j = 1:6
    subplot(2,1,1)
    hold on
    plot(U_,real(p_(j,:).*p.c./(2*U_)));
    xlabel('U')
    ylabel('p_R c / 2U\infty')
    grid on
    subplot(2,1,2)
    hold on
    plot(U_,abs(imag(p_(j,:))/(2*pi)));
    xlabel('U')
    ylabel('p_I / 2\pi [Hz]')
    grid on
end


