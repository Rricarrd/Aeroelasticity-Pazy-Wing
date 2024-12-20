function [mu,r0,rS,d] = computeBeamMassParameters(p)

syms x y t % variables

% Numeric parameters
xs = p.xs;         % [m] Shear center
z = p.z;           % Maximum thickness
xi = p.xi;         % [m] Initial x position
xf = p.xf;         % [m] Final x position
c = p.c;            % [m] Chord
rhoNyl = p.rhoNyl;       % [m] Nylon density
rhoAl = p.rhoAl;       % [m] Aluminium density
b = p.b;           % [m] Aluminium bar width
t = p.t;        % [m] Aluminium bar thickness 
bi = p.bi;      % [m] Aluminium bar starting position

% NACA 0018 expression
yt = 5 * z * c * ( ...
    0.2969 * sqrt(x/c) - ...
    0.1260 * (x/c) - ...
    0.3516 * (x/c).^2 + ...
    0.2843 * (x/c).^3 - ...
    0.1015 * (x/c).^4); %[m]

%%%%%%% Linear density
mu = double(2*int(int(rhoNyl,y,0,yt),x,0,xi) + ...
    b*t*rhoAl + ...
    int(2*int(rhoNyl,y,0,yt),x,xf,c));


%%%%%%% Parameter r0
r0 = double(sqrt(2*((int(int(rhoNyl*y^2,y,0,yt),x,0,xi)) + ...
           (int(int(rhoAl*y^2,y,0,t/2),x,bi,b+bi)) + ...
           int(int(rhoNyl*y^2,y,0,yt),x,xf,c))) ...
    /mu);

%%%%%%% Parameter rS
rS = double(sqrt(2*((int(int(rhoNyl*((xs-x)^2 + y^2),y,0,yt),x,0,xi)) + ...
           (int(int(rhoAl*((xs-x)^2 + y^2),y,0,t/2),x,bi,b+bi)) + ...
           int(int(rhoNyl*((xs-x)^2 + y^2),y,0,yt),x,xf,c))) ...
    /mu);

%%%%%%% Parameter d
d = double(2*((int(int(rhoNyl*(xs-x),y,0,yt),x,0,xi)) + ...
           (int(int(rhoAl*(xs-x),y,0,t/2),x,bi,b+bi)) + ...
           int(int(rhoNyl*(xs-x),y,0,yt),x,xf,c)) ...
    /mu);


% mu = 0.445;
% r0 = 0.0012;
% rS = 0.0234;
% d = -0.0061;



end

