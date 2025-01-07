function [A0k_num, A1kc_num, A1knc_num] = SymbolicAerodynamicMatrices(p)
%% Aerodynamic matrices
% Definition of symbolic variables
syms xi hk y t Ck c rho Uinf xs 
syms eta1(t) eta2(t) gamma1(t) gamma2(t) theta1(t) theta2(t) 
syms deltaeta1(t) deltaeta2(t) deltagamma1(t) deltagamma2(t) deltatheta1(t) deltatheta2(t) 

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
lc = pi*rho*Uinf^2*c*Ck*(theta_xi - b*diff(theta_xi,t)/(2*c) - diff(eta_xi,t)/(2*c));
msc = pi*rho*Uinf^2*c*a*Ck*(theta_xi - b*diff(theta_xi,t)/(2*c) - diff(eta_xi,t)/(2*c));
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
        A1kc(i,j) = simplify(diff(diff(deltaWkc, q2(i)), diff(q1(j),t)));
        A1knc(i,j) = simplify(diff(diff(deltaWknc, q2(i)), diff(q1(j),t)));
    end
end

% Substituting constant parameters (And dividing by hk and Uinf (And Ck where applicable))
A0k_num = double(subs(A0k,[xs,c,rho],[p.xs,p.c,p.rho])/(Uinf^2*hk*Ck));
A1kc_num = double(subs(A1kc,[xs,c,rho],[p.xs,p.c,p.rho])/(Uinf^2*hk*Ck));
A1knc_num = double(subs(A1knc,[xs,c,rho],[p.xs,p.c,p.rho])/(Uinf^2*hk));




end
