
syms eta theta rho xs

% Expressió analítica del NACA 0018
https://en.wikipedia.org/wiki/NACA_airfoil

airfoil = importdata('NACA0018.txt'); % Geometry of the airfoil

h = 0.04;
x = flip(airfoil(:,1))/10; % [m]
y = flip(airfoil(:,2))/10; % [m]
A = double(trapz(x(5:17),y(5:17)))*2 - 0.005*0.03 - 0.06*0.00225;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
syms x y t % variables
syms theta(t) eta(t) delta_theta(t) delta_eta(t) % functions

% Parameters
xs = sym("x_s")


% Kinetics
w(x,y,t) = eta + theta*(xs-x)

% Kinetic energy
T = int(int((diff(w,t))^2,x,0,c),y,0,h)

 