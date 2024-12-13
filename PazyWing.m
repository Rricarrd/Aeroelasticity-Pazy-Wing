
%% PAZY WING PROBLEM
clc
%% Calculate Wing Stiffness
% K Matrix Parameters
if exist('EI','var') ~= 1
    [GJ,EI,xs] = CalcStiffness()
end
p.EI = EI
p.GJ = GJ
p.xs = xs


% M Matrix parameters
p.mu = 1;
p.rS = 1;
p.r0 = 0;
p.d = 1;

%% Numerical discretization
% Parameters
p.n = 14;         % Panel node number
p.Rs = 38.4;       % Rib spacing
p.Lw = 550;        % Wing lenght
p.Ls = p.Rn*p.Rs - p.Lw; % Last wing segment space
p.nj = 3;          % DOFs per node (theta, eta i gamma)
p.ng = 1;          % Gauss points per element
p.en = 2;          % Element nodes
 

% Discretization and conectivities
[y,Tn,Tr,p.Nnode,p.Nelem,p.Ndof] = DiscretizeWing(p);


%% Element mass matrices

M = AssemblyM(y,Tn,Tr,p);

K = AssemblyK(y,Tn,Tr,p);




