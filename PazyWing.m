
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

p.rhoAl = 2795; % [kg/m^3] Density of Aluminium
p.rhoNyl = 930; % [kg/m^3] Density of Nylon 12
p.z = 0.018;          % [m] Maximum thickness
p.xi= 0.0075;         % [m] Initial x position
p.xf = 0.095;         % [m] Final x position
p.c = 0.1;            % [m] Chord
p.b = 0.06;           % [m] Aluminium bar width
p.t = 0.00225;        % [m] Aluminium bar thickness 
p.bi = 0.02;           % [m] Aluminium bar front edge position
p.h = 0.004;          % [m] Depth
p.a_slot = 0.005;     % [m] Plastic slot width



%% Numerical discretization
% Parameters
p.n = 14;              % Panel node number
p.Rs = 0.0384;         % [m] Rib spacing
p.Lw = 0.550;          % [m] Wing lenght
p.Ls = p.Lw - 13*p.Rs; % Last wing segment space
p.nj = 3;              % DOFs per node (theta, eta i gamma)
p.ng = 1;              % Gauss points per element
p.en = 2;              % Element nodes
p.Rt = 0.004;          % [m] Rib thickness

p.sparse = 0;

% Computing beam mass parameters
[p.mu,p.r0,p.rS,p.d] = computeBeamMassParameters(p);

% Discretization and conectivities
[y,Tn,Tr,p.Nnode,p.Nelem,p.Ndof] = DiscretizeWing(p);


%% Element mass matrices

M = AssemblyM(y,Tn,Tr,p);

K = AssemblyK(y,Tn,Tr,p);


%% Free vibrations
k = 10;
[Q,W] = eig(K,M); % Solve for eigenvalues

f = sqrt(diag(W))/(2*pi);
