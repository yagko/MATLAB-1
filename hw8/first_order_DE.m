clear all
close all
clc
% y' = ay + by^2 + cy(integral(y(t))(t = 0 to t))

% define min and max domain values
t0 = 0;
t1 = 20;

% define step size
h = 0.2;

% we will divide the domain interval into M equal subintervals
M = (t1-t0)/h;

% define coefficients
a = 1.3;
b = -0.25;
c = -0.0001;

% initial values
y0_1 = 1;  % for part a
y0_2 = 10; % for part b

% integral's value at t=0
K0 = 0;

% predefining the result vectors
y1 = zeros(1,M+1);
y2 = zeros(1,M+1);

% predefining the vectors which will be used for approximating ...
% the integrals 
K1 = zeros(1,M+1);
K2 = zeros(1,M+1);


% first terms of the vectors for part a
y1(1) = y0_1;
K1(1) = K0;

% function for part a
y_1 = @(K1,y1) a*y1 + b*y1.^2 + c*y1*K1;

% Euler's method with trapezoidal rule 
for i = 1:M
    y1(i+1) = y1(i) + h*y_1(K1(i),y1(i));
    K1(i+1) = K1(i) + h*(y1(i)+y1(i+1))/2;
end

% first terms of the vectors for part b
y2(1) = y0_2;
K2(1) = K0;

% function for part b
y_2 = @(K2,y2) a*y2 + b*y2.^2 + c*y2*K2;

% Euler's method with trapezoidal rule 
for i=1:M
    y2(i+1) = y2(i) + h*y_2(K2(i),y2(i));
    K2(i+1) = K2(i) + h*(y2(i)+y2(i+1))/2;
end

% defining the x axis
tx = t0:h:t1;

figure (1)
plot(tx,y1,'green','Linewidth',2)
hold on
plot(tx,y2,'black','Linewidth',2)
title("First-Order Integro-Ordinary Differential Equation")
legend('y, with y(0)=1','y, with y(0)=10')