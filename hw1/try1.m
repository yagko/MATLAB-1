clc;clear all;close all;


y = @(t) 9600*(1-exp(-t/15))-480*t; % height
x = @(t) 2400*(1-exp(-t/15)); % range

% initial guess
t0 = 10;

% convergence criteria constants
delta = 10^(-6);
eps = 10^(-6);

% first differentials
y1 = @(t) 640*exp(-t/15) - 480;

maxI = 1000;

for i = 1:maxI
    t1 = t0 - (y(t0)/y1(t0));

    err = abs(t1-t0); 
    relerr = err/(abs(t1) + eps);
    err2 = abs(y(t1)); %required error values
    if (err < delta) && (relerr < delta) && (err2 < delta) 
        break
    end
    t0 = t1;
end

impact_time = t1
range = x(t1)

t = 0:0.01:impact_time;

plot(x(t),y(t))
xlabel('x (m)')
ylabel('y (m)')
title('Projectile Motion')
axis([0 1200, 0 400])
grid on