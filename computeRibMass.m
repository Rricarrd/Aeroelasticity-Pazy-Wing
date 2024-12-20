function [Mribnum] = computeRibMass(xsnum,hnum,znum,xinum,xfnum,cnum,a_slotnum,rhonum,x1num,x2num)

syms x y t % variables
syms theta(t) eta(t) gamma_(t) % functions

% Variables
xs = sym("x_s");         % [m] Shear center
h = sym("h");            % [m] Depth
z = sym("z");            % Maximum thickness
xi = sym("xi");          % [m] Initial x position
x1= sym("x1");           % [m] Start of slot
x2= sym("x1");           % [m] End of slot
xf = sym("xf");          % [m] Final x position
c = sym("c");            % [m] Chord
a_slot = sym("a_slot");  % [m] Plastic slot width
rho = sym("rho");        % [m] Nylon density



% NACA 0018 expression
yt(x) = 5 * z * c * ( ...
    0.2969 * sqrt(x/c) - ...
    0.1260 * (x/c) - ...
    0.3516 * (x/c).^2 + ...
    0.2843 * (x/c).^3 - ...
    0.1015 * (x/c).^4); %[m]


% Slot expression POTSER FALTA EL FORAT DEL NYLON PER L'ALUMINI
ys = a_slot

% Kinetics
w(x,y,t) = eta + theta*(xs-x);

% Kinetic energy

T = h*int(2*int(0.5*rho*diff(w,t)^2,y,0,yt-ys),x,xi,xf)...
     - h*int(2*int(0.5*rho*diff(w,t)^2,y,0,a_slot),x,x1,x2);
q = {theta(t) gamma_(t) eta(t)};

% Find Me rib 

for i=1:3
    for j=1:3
        Mrib(i,j) = simplify(diff(diff(diff(T,diff(q{i},t)),t),diff(diff(q{j},t),t)));
    end
end

Mribnum = double(subs(Mrib,[xs,h,z,xi,xf,c,a_slot,rho,x1,x2],[xsnum,hnum,znum,xinum,xfnum,cnum,a_slotnum,rhonum,x1num,x2num]));
end

