function [Mribnum] = computeRibMass(p)

% Parameters
xsnum = p.xs
hnum = p.h;
znum = p.z;
xinum = p.xi;
xfnum = p.xf;
cnum = p.c;
a_slotnum = p.a_slot;
a_Alnum = p.a_Al;
rhonum = p.rhoNyl;
x1num = p.x1;
x2num = p.x2;
xAl1num = p.xAl1;
xAl2num = p.xAl2;

% Variables
syms x y t % variables
syms theta(t) eta(t) gamma_(t) % functions

xs = sym("x_s");         % [m] Shear center
h = sym("h");            % [m] Depth
z = sym("z");            % Maximum thickness
xi = sym("xi");          % [m] Initial x position
x1= sym("x1");           % [m] Start of slot
x2= sym("x2");           % [m] End of slot
xAl1= sym("xAl1");           % [m] Start of slot for aluminium
xAl2= sym("xAl2");           % [m] End of slot for aluminium
xf = sym("xf");          % [m] Final x position
c = sym("c");            % [m] Chord
a_slot = sym("a_slot");  % [m] Plastic slot width
a_Al = sym("a_Al");      % [m] Aluminium slot width
rho = sym("rho");        % [m] Nylon density



% NACA 0018 expression
yt(x) = 5 * z * c * ( ...
    0.2969 * sqrt(x/c) - ...
    0.1260 * (x/c) - ...
    0.3516 * (x/c).^2 + ...
    0.2843 * (x/c).^3 - ...
    0.1015 * (x/c).^4); %[m]


% Kinetics
w(x,y,t) = eta + theta*(xs-x);

% Kinetic energy

    % T1 = h*int(2*int(0.5*rho*diff(w,t)^2,y,0,yt-ys),x,xi,xf)...
    %      - h*int(2*int(0.5*rho*diff(w,t)^2,y,0,a_slot),x,x1,x2); % OLD

% The following two expressions for T are equivalent:

    T = h*int(2*int(0.5*rho*diff(w,t)^2,y,0,yt),x,xi,xf)...
         - h*int(2*int(0.5*rho*diff(w,t)^2,y,0,a_slot/2),x,x1,x2)...
         - h*int(2*int(0.5*rho*diff(w,t)^2,y,0,a_Al/2),x,xAl1,xAl2);

    % T2 = h*int(2*int(0.5*rho*diff(w,t)^2,y,0,yt),x,xi,xAl1)...
    %      + h*int(2*int(0.5*rho*diff(w,t)^2,y,a_Al/2,yt),x,xAl1,x1)...
    %      + h*int(2*int(0.5*rho*diff(w,t)^2,y,a_slot/2,yt),x,x1,x2)...
    %      + h*int(2*int(0.5*rho*diff(w,t)^2,y,a_Al/2,yt),x,x2,xAl2)...
    %      + h*int(2*int(0.5*rho*diff(w,t)^2,y,0,yt),x,xAl2,xf);

    

q = {theta(t) gamma_(t) eta(t)};

% Find Me rib 
for i=1:3
    for j=1:3
        Mrib(i,j) = simplify(diff(diff(diff(T,diff(q{i},t)),t),diff(diff(q{j},t),t)));
    end
end

Mribnum = double(subs(Mrib,[xs,h,z,xi,xf,c,a_slot,a_Al,rho,x1,x2,xAl1,xAl2],[xsnum,hnum,znum,xinum,xfnum,cnum,a_slotnum,a_Alnum,rhonum,x1num,x2num,xAl1num,xAl2num]));

end

