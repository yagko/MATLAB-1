% ME 303 - Lab Session 8 - December 7th, 2021
% Solution of Ordinary Differential Equations
% Heun's Method

clear all
close all
clc

t0 = 0;
tf = 3;
y0 = 1;
h = 0.5;
M = (tf - t0)/h;

t = t0:h:tf;

f = @(t,y) (t-y)/2;
y = @(t) 3*exp(-t/2) + t - 2;

y_e = zeros(1,M+1);
y_h = zeros(1,M+1);

y_e(1) = y0;
y_h(1) = y0;

for j = 1:M
    % Euler's Part
    f1_e = f(t(j),y_e(j));
    y_e(j+1) = y_e(j) + h*f1_e;
    % Heun's Part
    f1_h = f(t(j),y_h(j));
    f2_h = f(t(j+1),y_h(j) + h*f1_h);
    y_h(j+1) = y_h(j) + h/2*(f1_h + f2_h);
end

t_exact = t0:0.001:tf;
figure
plot(t,y_e,'r--')
hold on
plot(t,y_h,'g*-')
hold on
plot(t_exact,y(t_exact),'k','Linewidth',3)
legend('Euler','Heun','Exact')

