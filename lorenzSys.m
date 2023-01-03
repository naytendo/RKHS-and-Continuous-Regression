function [dphidt] =lorenzSys(t,phi)
% Defining system parameters
sigma = 10;
r = 28;
b = 8/3;

% Define states
x = phi(1);
y = phi(2);
z = phi(3);

dxdt = sigma*(y - x);
dydt = x*(r- z)-y;
dzdt = x*y -  b*z;

dphidt = [dxdt;dydt;dzdt];

end