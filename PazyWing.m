
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
KLL = K(4:end,4:end);
MLL = M(4:end,4:end);

k = 10;
[Q,W] = eigs(KLL,MLL,k,'sm'); % Solve for eigenvalues

f = sqrt(diag(W))/(2*pi);

%% Aerodynamic matrices

% Definition of symbolic variables
syms xi hk y t C c rho Uinf xs a b
syms eta1(t) eta2(t) gamma1(t) gamma2(t) theta1(t) theta2(t) theta_xi(t) eta_xi(t)
syms deltaeta1(t) deltaeta2(t) deltagamma1(t) deltagamma2(t) deltatheta1(t) deltatheta2(t) deltatheta_xi(t) deltaeta_xi(t)
syms l ms % Aerodynamic functions

% Definition of nodal variable vector
q1 = {eta1; gamma1; theta1; eta2; gamma2; theta2};
q2 = {deltaeta1; deltagamma1; deltatheta1; deltaeta2; deltagamma2; deltatheta2};

% Definition of y (according to P1-14), b and a
y = p.Lw*(xi+1)/4;
b = xs - 3*c/4;
a = xs - c/4;

% Definition of shape functions (according to P1-1)
psi1 = (1 + xi)/2;
psi2 = (1 - xi)/2;

phi1 = (2 + 3*xi - xi^3)/4;
phi2 = (2 - 3*xi + xi^3)/4;

phi_bar1 = diff(y,xi)*(-1 - xi + xi^2 + xi^3)/4;
phi_bar2 = diff(y,xi)*(1 - xi - xi^2 + xi^3)/4;

% Substitute expressions for theta(xi) and eta(xi)
theta_xi = psi1*theta1 + psi2*theta2;
eta_xi = phi1*eta1 + phi_bar1*gamma1 + phi2*eta2 + phi_bar2*gamma2;

deltatheta_xi = psi1*deltatheta1 + psi2*deltatheta2;
deltaeta_xi = phi1*deltaeta1 + phi_bar1*deltagamma1 + phi2*deltaeta2 + phi_bar2*deltagamma2;

% Definition of lift and moment expressions of our model
lc = pi*rho*Uinf^2*c*C*(theta_xi - b*diff(theta_xi,t)/(2*c) - diff(eta_xi,t)/(2*c));
msc = pi*rho*Uinf^2*c*a*C*(theta_xi - b*diff(theta_xi,t)/(2*c) - diff(eta_xi,t)/(2*c));
lnc = 0.5*pi*rho*Uinf^2*c*(diff(theta_xi,t) - (2*xs/c - 1)*diff(diff(theta_xi,t),t) - 2*(diff(diff(eta_xi,t),t))/c);
msnc = -0.5*pi*rho*Uinf^2*c^2/2*(3/2*diff(theta_xi,t) - (2*xs/c - 9/8)*diff((diff(theta_xi,t)),t) - 2*(diff(diff(eta_xi,t),t))/c);

% Definition of the virtual work term deltaW
deltaWk = int((lc+lnc)*deltaeta_xi*hk/2,xi,-1,1) + ...
          int((msc+msnc+xs*lnc)*deltatheta_xi*hk/2,xi,-1,1);
deltaWkc = int(lc*deltaeta_xi*hk/2,xi,-1,1) + ...
          int(msc*deltatheta_xi*hk/2,xi,-1,1);
deltaWknc = int(lnc*deltaeta_xi*hk/2,xi,-1,1) + ...
          int(msnc*deltatheta_xi*hk/2,xi,-1,1);

% Computation of element matrix coefficients symbolically
A0k = sym(zeros(6, 6)); 
A1kc = sym(zeros(6, 6));
A1knc = sym(zeros(6, 6));

for i=1:6
    for j=1:6
        A0k(i,j) = simplify(diff(diff(deltaWk, q2(i)), q1(j)));
        A0kc(i,j) = simplify(diff(diff(deltaWkc, q2(i)), diff(q1(j),t)));
        A0knc(i,j) = simplify(diff(diff(deltaWknc, q2(i)), diff(q1(j),t)));
    end
end

