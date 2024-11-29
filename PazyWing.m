
%% PAZY WING PROBLEM
clc
%% prova

%% Calculate Wing Stiffness
[GJ,EI,xs] = CalcStiffness()

%% Numerical discretization
% Rib positions
yj = 0:38.4:550
yj(15) = 550

% Nodal coordinates
n = 14
y = 0
for i = 1:length(yj)-1
    yi = linspace(yj(i),yj(i+1),n);
    y = [y,yi(2:end)];
end


% Element connectivities
Tn = [1:length(y); 2:length(y)+1]' % node conectivities
Tr = (n:n:length(y))' 
Tr(end) = length(y)+1; % rib conectivities




%% Element mass matrices
Li = 1
mu = 1
rs = 1
r0 = 0
d = 1
% Element mass matrix [M]_i
M1 = (mu*Li/420) * [156 22*Li 0 54 -13*Li 0;
                   22*Li 4*Li^2 0 13*Li -3*Li^2 0;
                   0 0 0 0 0 0;
                   54 13*Li 0 156 -22*Li 0;
                   -13*Li -3*Li^2 0 -22*Li 4*Li^2 0;
                   0 0 0 0 0 0];

M2 = (mu*r0^2/(30*Li)) * [36 3*Li 0 -36 -3*Li 0;
                          3*Li Li^2 0 -3*Li -Li^2 0;
                          0 0 0 0 0 0;
                          -36 -3*Li 0 36 3*Li 0;
                          -3*Li -Li^2 0 3*Li Li^2 0;
                          0 0 0 0 0 0];

M3 = (mu*rS^2*Li/6) * [0 0 0 0 0 0;
                       0 2 0 0 0 1;
                       0 0 0 0 0 0;
                       0 0 0 0 0 0;
                       0 0 0 0 2 0;
                       0 1 0 0 0 2];

M4 = (mu*d*Li/60) * [0 21 0 0 9 0;
                    21 3*Li 0 0 0 2*Li;
                    0 0 0 0 0 0;
                    0 0 0 0 21 3*Li;
                    9 0 0 21 -2*Li -3*Li;
                    0 2*Li 0 -3*Li 0 0];

M_i = M1 + M2 + M3 + M4;

% Element stiffness matrix [K]_i
K1 = (EI/Li^3) * [12 6*Li 0 -12 6*Li 0;
                 6*Li 4*Li^2 0 -6*Li 2*Li^2 0;
                 0 0 0 0 0 0;
                 -12 -6*Li 0 12 -6*Li 0;
                 6*Li 2*Li^2 0 -6*Li 4*Li^2 0;
                 0 0 0 0 0 0];

K2 = (GJ/Li) * [0 0 0 0 0 0;
               0 0 0 0 0 0;
               0 0 1 0 0 -1;
               0 0 0 0 0 0;
               0 0 0 0 0 0;
               0 0 -1 0 0 1];

K_i = K1 + K2;

