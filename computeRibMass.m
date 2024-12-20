function [Mribnum] = computeRibMass(xsnum,hnum,znum,xinum,xfnum,cnum,a_slotnum,rhonum)

syms x y t % variables
syms theta(t) eta(t) gamma_(t) % functions

% Variables
xs = sym("x_s");         % [m] Shear center
h = sym("h");            % [m] Depth
z = sym("z");            % Maximum thickness
xi = sym("xi");          % [m] Initial x position
xf = sym("xf");          % [m] Final x position
c = sym("c");            % [m] Chord
a_slot = sym("a_slot");  % [m] Plastic slot width
rho = sym("rho");        % [m] Nylon density



% NACA 0018 expression
yt = 5 * z * c * ( ...
    0.2969 * sqrt(x/c) - ...
    0.1260 * (x/c) - ...
    0.3516 * (x/c).^2 + ...
    0.2843 * (x/c).^3 - ...
    0.1015 * (x/c).^4); %[m]


% Slot expression POTSER FALTA EL FORAT DEL NYLON PER L'ALUMINI
ys = piecewise(x >= 0 & x <= 0.035, 0, ...
               x > 0.035 & x <= 0.065, a_slot/2, ...
               x > 0.065 & x <= 0.1, 0);

% Kinetics
w(x,y,t) = eta + theta*(xs-x);

% Kinetic energy

T = h*int(2*int(0.5*rho*diff(w,t)^2,y,0,yt),x,xi,xf);
q = {theta(t) gamma_(t) eta(t)};

% Find Me rib 

for i=1:3
    for j=1:3
        Mrib(i,j) = simplify(diff(diff(diff(T,diff(q{i},t)),t),diff(diff(q{j},t),t)));
    end
end

Mribnum = double(subs(Mrib,[xs,h,z,xi,xf,c,a_slot,rho],[xsnum,hnum,znum,xinum,xfnum,cnum,a_slotnum,rhonum]));
end

