function dx = acc(t,x)
global u v0;
% noise1 = 8*(rand() - 0.5);
% noise2 = 1.8*(rand() - 0.5);
noise1 = 0;
noise2 = 0;  %without noise
m = 1650; f0 = 0.1; f1 = 5; f2 = 0.25;
dx = zeros(4,1);
dx(1) = x(2) + noise1;
dx(2) = (1/m)*(u(1) - f0 - f1*x(2) - f2*x(2).^2) + noise2;
dx(3) = v0 - x(2) - noise1;
dx(4) = u(3);