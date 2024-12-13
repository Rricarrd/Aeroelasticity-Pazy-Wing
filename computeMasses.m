function [mu mrib] = computeMasses(p)

airfoil = importdata('NACA0018.txt'); % Geometry of the airfoil

x = flip(airfoil(:,1))/10; % [m]
y = flip(airfoil(:,2))/10; % [m]

A1 = double(trapz(x(1:5),y(1:5)))*2; % [m^2] Surface of the first plastic part of the rib
mu1 = A1*p.rhoNyl;

A2 = double(trapz(x(5:17),y(5:17)))*2 - 0.005*0.03 - 0.06*0.00225; % [m^2] Surface of the second plastic part of the rib
mrib = A2*p.rhoNyl*p.Rt; % [kg] Mass of the rib

A3 = double(trapz(x(17:18),y(17:18)))*2; % [m^2] Surface of the third plastic part of the rib
mu3 = A3*p.rhoNyl;

A4 = 0.06*0.00225; % [m^2] Surface of the Aluminium bar
mu4 = A4*p.rhoAl;

mu = mu1 + mu3 + mu4; % [kg/m] Linear beam density