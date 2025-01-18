function [Mrib] = computeRibMass(p)

% Parameters
xs = p.xs;
h = p.h;
z = p.z;
xi = p.xi;
xf = p.xf;
c = p.c;
a_slot = p.a_slot;
a_Al = p.a_Al;
rho = p.rhoNyl;
x1 = p.x1;
x2 = p.x2;
xAl1 = p.xAl1;
xAl2 = p.xAl2;

% Variables
syms x y t % variables
syms theta(t) eta(t) gamma_(t) % functions




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

    T = 0.5*rho*h*int(int(diff(w,t)^2,y,-yt,yt),x,xi,xf);
    
    % Testing
    % ...
    %      - 0.5*rho*h*int(int(diff(w,t)^2,y,-a_slot/2,a_slot/2),x,x1,x2)...
    %      - 0.5*rho*h*int(int(diff(w,t)^2,y,-a_Al/2,a_Al/2),x,xAl1,xAl2);

    % T2 = h*int(2*int(0.5*rho*diff(w,t)^2,y,0,yt),x,xi,xAl1)...
    %      + h*int(2*int(0.5*rho*diff(w,t)^2,y,a_Al/2,yt),x,xAl1,x1)...
    %      + h*int(2*int(0.5*rho*diff(w,t)^2,y,a_slot/2,yt),x,x1,x2)...
    %      + h*int(2*int(0.5*rho*diff(w,t)^2,y,a_Al/2,yt),x,x2,xAl2)...
    %      + h*int(2*int(0.5*rho*diff(w,t)^2,y,0,yt),x,xAl2,xf);

    

q = {eta(t) gamma_(t) theta(t)};

% Find Me rib 
for i=1:3
    for j=1:3
        
        dt_dT_dqdot = diff(diff(T,diff(q{i},t)),t);
        dqdotdot = diff(diff(q{j},t),t);
            
        Mrib(i,j) = simplify(diff(dt_dT_dqdot,dqdotdot));
    end
end

Mrib = double(Mrib);

end

