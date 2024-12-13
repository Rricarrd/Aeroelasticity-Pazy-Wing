%% Virtual test Pazy wing
clc
clear
close all

% Choose masses (in kg) to add at each slot in the wing tip rod
f = 1; % kg
fload = [
    0      % x01 = -0.10 m
    0      % x02 = -0.08 m
    0      % x03 = -0.06 m
    0      % x04 = -0.04 m
    0      % x05 = -0.02 m
    0      % x06 =  0.00 m
    0      % x07 =  0.02 m
    0      % x08 =  0.04 m
    0      % x09 =  0.06 m
    0      % x10 =  0.08 m
    0      % x11 =  0.10 m
    0      % x12 =  0.12 m
    0      % x13 =  0.14 m
    0      % x14 =  0.16 m
    0      % x15 =  0.18 m
    f      % x16 =  0.20 m
];

% Test is performed and vertical displacements are recorded in displ(:,1).
% displ(:,2) are the x-coordinates of each sensor (chordwise direction).
% displ(:,3) are the y-coordinates of each sensor (spanwise direction).
displ = PazyWingLoad(fload,false);

% Note: to prevent plot from showing set second input in the 'PazyWingLoad'
% funciton to 'false'.


%% Data input
% Vertical displacements
w1 = displ(1:2:end,1);
w2 = displ(2:2:end,1);

% Other coordinates
x2 = displ(2:2:end,2);
x1 = displ(2:2:end,1);
y = displ(1:2:end,3);
 
theta = (w1-w2)/0.6;

% Plot of the vertical displacements
figure(2)
scatter(y,w1)
hold on
scatter(y,w2)
hold off
legend("Leading edge", "Trailing edge")
xlabel("y-position [m]")
ylabel("Vertical displacement [m]")


%% Finding shear center position
x_pos_ = -0.10:0.02:0.2;
twists_ = zeros(1,length(x_pos_));

for i = 1:length(x_pos_)
  fload = zeros(1,length(x_pos_));
  fload(i) = f;
  displ = PazyWingLoad(fload,false);
  w1_27 = displ(27,1);
  w2_28 = displ(28,1);
  twists_(i) = (w1_27-w2_28)/60;
end

% Plot of twist vs applied load x position
figure(3)
scatter(x_pos_, twists_)
xlabel("Applied load X position [m]")
ylabel("Twist angle (\theta)")

p = polyfit(twists_, x_pos_,1);
xs = polyval(p,0);

%% Bending at shear
ws = w2 + theta.*(x2-xs);

%% Finding coefficients
% Fitting twist to y
a = polyfit(y, theta,1);
a1 = a(1);

% Plot of the twist
figure(4)
scatter(y,theta)
hold on
plot(y,polyval(a,y))
xlabel("y-position [m]")
ylabel("Angle of twist (\theta)")

% Fitting bending to y
b = polyfit(y, ws,3);
b3 = b(1);

% Plot of the bending at the shear center
figure(5)
scatter(y,ws)
hold on
plot(y,polyval(b,y))
xlabel("y-position [m]")
ylabel("Vertical displacement at xs (ws)")

%%  Physical parameters
F = f*9.81;
GJ = F*(xs-x_pos_)/a1;
EI = -F/(6*b3);

