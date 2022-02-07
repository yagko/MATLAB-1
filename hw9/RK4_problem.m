clc; clear all; close all;

% y'(t) = k(a-y(t))(b-y(t))
% y(0) = 0
% y(t) = 350(1-exp(-0.2t))/(7-5*exp(-0.2t))

% constants for the function
k = 0.01;
a = 70; % millimoles/liter
b = 50; % millimoles/liter

% time interval
ti = 0;
tf = 20;

% step sizes 
h1 = 2;
h2 = 1;

% defining the function
f = @(t,y) k*(a-y)*(b-y);

% exact solution
y_exact = @(t) 350*(1-exp(-0.2*t))./(7-5*exp(-0.2*t));


% initial value condition
y0 = 0;

% dividing the interval with step size h1 and h2
t1 = ti:h1:tf;
t2 = ti:h2:tf;

% predefining the result vectors
y1_RK4 = zeros(1,length(t1));
y2_RK4 = zeros(1,length(t2));

% initial values
y1_RK4(1) = y0;
y2_RK4(1) = y0;

% RK4 method for step size (h) = 2
for i = 1:length(t1)-1
    f1 = f(t1(i) , y1_RK4(i));
    f2 = f(t1(i) + h1/2 , y1_RK4(i) + h1/2 * f1);
    f3 = f(t1(i) + h1/2 , y1_RK4(i) + h1/2 * f2);
    f4 = f(t1(i) + h1 , y1_RK4(i) + h1 * f3);

    y1_RK4(i+1) = y1_RK4(i) + h1/6 * (f1 + 2*f2 + 2*f3 + f4);
end

% RK4 method for step size (h) = 1
for i = 1:length(t2)-1
    f1 = f(t2(i) , y2_RK4(i));
    f2 = f(t2(i) + h2/2 , y2_RK4(i) + h2/2 * f1);
    f3 = f(t2(i) + h2/2 , y2_RK4(i) + h2/2 * f2);
    f4 = f(t2(i) + h2 , y2_RK4(i) + h2 * f3);

    y2_RK4(i+1) = y2_RK4(i) + h2/6 * (f1 + 2*f2 + 2*f3 + f4);
end

% defining a time vector for the exact solution
t_exact = ti:0.001:tf;

plot(t1,y1_RK4,'g-*')
grid on
hold on
plot(t2,y2_RK4,'k-*')
plot(t_exact,y_exact(t_exact),'r')
legend( 'y_{RK4}(h=2)', 'y_{RK4}(h=1)', 'y_{exact}' )

Y1 = y1_RK4;
Y2 = y2_RK4;