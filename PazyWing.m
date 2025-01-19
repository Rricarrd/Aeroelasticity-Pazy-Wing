
%% PAZY WING PROBLEM
clc
clear
%% Calculate Wing Stiffness
% K Matrix Parameters
[GJ,EI,xs] = CalcStiffness();

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
p.a_Al = 0.00225;     % [m] Aluminium slot width
p.x1 = 0.035;         % [m] Plastic slot start position
p.x2 = 0.065;         % [m] Plastic slot end position
p.xAl1 = 0.02;        % [m] Plastic slot start position
p.xAl2 = 0.08;        % [m] Plastic slot end position

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

M_ = AssemblyM(y,Tn,Tr,p);

K_ = AssemblyK(y,Tn,Tr,p);


%% Free vibrations
% Fixing first nodes
K_ = K_(4:end,4:end);
M_ = M_(4:end,4:end);

Nm = 10;
[Q,W] = eigs(K_,M_,Nm,'smallestabs'); % Solve for eigenvalues

f = sqrt(diag(W))/(2*pi);


%% Aerodynamic elemental matrices
% Elemental matrices
[A0k, A1kc, A1knc] = SymbolicAerodynamicMatrices(p); % Matrices divided by U, hk and Ck

%% Aerodynamic global matrices
% Assembly of matrices (Without U and Ck)
A0_ = AssemblyA(A0k, y, Tn, p);
A1c_ = AssemblyA(A1kc, y, Tn, p);
A1nc_ = AssemblyA(A1knc, y, Tn, p);

% Excluding BCs
A0_ = A0_(4:end,4:end);
A1c_ = A1c_(4:end,4:end);
A1nc_ = A1nc_(4:end,4:end);


%% Order reduction
M = Q'*M_*Q;
K = Q'*K_*Q;
A0 = Q'*A0_*Q;
A1c = Q'*A1c_*Q;
A1nc = Q'*A1nc_*Q;


%% Divergence
[Ud] = Divergence(p,K,A0,50);


%% Flutter p method
Umax = 100; %[m/s]
[U_,p_,Vp_] = FlutterPMethod(Umax,Nm,p,K,M,A0,A1c,A1nc);

figure
for j = 1:Nm
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


%% Flutter k method
[Uk_,gk_,wk_] = FlutterKMethod(Umax,Nm,p,K,M,A0,A1c,A1nc);
figure
subplot(2,1,1)
plot(Uk_',gk_');
xlabel("U_\infty");
ylabel("g");
xlim([0 20])
grid on
subplot(2,1,2)
plot(Uk_',wk_'/(2*pi));
xlabel("U_\infty ");
ylabel("\omega / 2\pi [Hz]");
xlim([0 20])
grid on

%% Flutter pk method
[U_, gam_,w_] = FlutterPKMethod(Umax,Nm,p,K,M,A0,A1c,A1nc);
figure
subplot(2,1,1)
plot(U_, gam_.*p.c./(2*U_))
xlabel("U_\infty")
ylabel("\gammac/2U\infty")
subplot(2,1,2)
plot(U_, w_/(2*pi))
xlabel("U_\infty")
ylabel("\omega/ 2\pi [Hz]")
grid on



%% Plot modes
for i = 1:Nm
    modes_legend{i} = sprintf("Mode %i, $f = %.2f Hz$",i,f(i));
end


%Eta
figure
for i = 1:Nm
    plot(y(2:end),Q(1:3:end,i))
    hold on
end
grid minor;
title(sprintf("First %i modal displacements",Nm))
xlabel("y [m]", 'Interpreter', 'latex');
ylabel("Modal displacements $\Phi(\eta)$", 'Interpreter', 'latex');
legend(modes_legend{1:Nm},'Interpreter',"latex");


% Gamma
figure
for i = 1:Nm
    plot(y(2:end),Q(2:3:end,i))
    hold on
end
grid minor;
title(sprintf("First %i modal displacements",Nm))
xlabel("y [m]", 'Interpreter', 'latex');
ylabel("Modal displacements $\Phi(\gamma)$", 'Interpreter', 'latex');
legend(modes_legend{1:Nm},'Interpreter',"latex");


% Theta
figure
for i = 1:Nm
    plot(y(2:end),Q(3:3:end,i))
    hold on
end
grid minor;
title(sprintf("First %i modal displacements",Nm))
xlabel("y [m]", 'Interpreter', 'latex');
ylabel("Modal displacements $\Phi(\theta)$", 'Interpreter', 'latex');
legend(modes_legend{1:Nm},'Interpreter',"latex");
