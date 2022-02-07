clear all
close all
clc

format short

y = @(t) 9600*(1-exp(-t/15))-480*t;
dy = @(t) 9600*(exp(-t/15)/15)-480; %derivative of y function
x = @(t) 2400*(1-exp(-t/15));
p0 = 10; %initial guess
delta = 10.^(-6);
epsilon = 10.^(-6);
count = 50; %since one may not know the approximate solution counter should be high enough

for i=1:count
    p1 = p0 - y(p0)/dy(p0); %newton-raphson equation
    err = abs(p1-p0); 
    relerr = err/(abs(p1) + epsilon);
    err2 = abs(y(p1)); %required error values
    p0 = p1; %assigning the new initial guess

    if (err<delta) && (relerr<delta) && (err2<delta)
        break
    end %error checking
end

impact_time = p0
range = x(p0)
y_value = y(p0)

t=0:0.01:impact_time; %assigning values to variable 't' as a vector such that it can plot in terms of 't'

plot(x(t),y(t))

% Please plot the projectile motion below this line.
xlabel('x (m)')
ylabel('y (m)')
title('Projectile Motion')
axis([0 1200, 0 400])
grid on
