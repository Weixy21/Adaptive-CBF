function dx = acc_racbf(t,x)
global u v0;
% noise1 = 8*(rand() - 0.5); % safety is not guaranteed with noise, unless consider the noise bound in the racbf
% noise2 = 1.8*(rand() - 0.5);
noise1 = 0;
noise2 = 0;
m = 1650; f0 = 0.1; f1 = 5; f2 = 0.25;
dx = zeros(5,1);
dx(1) = x(2) + noise1;
dx(2) = (1/m)*(u(1) - f0 - f1*x(2) - f2*x(2).^2) + noise2;
dx(3) = v0 - x(2) - noise1;
dx(4) = x(5);
dx(5) = u(3);