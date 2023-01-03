
n = 2;
alpha0 = [7
   20.5];
gamma = 0.85;
delta = 0.4;
tau = 1.25;

[boolean,kProb,checks] = isDiophantine(delta,alpha0,gamma,tau);

alpha = alpha0/norm(alpha0);

if boolean
    T = (1+n^2*factorial(n))^(tau+1)/gamma/delta;
    t = 1:0.01:T;
    rho1 = 4;
    rho2 = 1;
    theta1 = alpha(1)*t;
    theta2 = alpha(2)*t;
    theta1T = zeros(length(theta1),1);
    theta2T = zeros(length(theta2),1);
    for ii = 1:length(theta1)
        theta1T(ii) = mod(theta1(ii),1);
        theta2T(ii) = mod(theta2(ii),1);
    end
    figure()
    plot(theta1T,theta2T,'.','markersize',0.5)
    grid on

    w1 = alpha(1);
    w2 = alpha(2);
    x = (rho1+rho2*cos(2*pi*w2*t)).*cos(2*pi*w1*t);
    y = (rho1+rho2*cos(2*pi*w2*t)).*sin(2*pi*w1*t);
    z = -rho2*sin(2*pi*w2*t);
    figure()
    plot3(x,y,z);
end