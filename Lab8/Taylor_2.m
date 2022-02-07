% ME 303 - Lab Session 8 - December 7th, 2021
% Solution of Ordinary Differential Equations
% Taylor's Method 2

clear all
close all
clc

N = 4;
syms t y
f = (t-y)/2;
d(1) = f;

for k = 1:N-1
    d(k+1) = diff(d(k),t) + f*diff(d(k),y);
end
d = d';

% d = [y'(t_k)
%      y''(t_k)
%      y'''(t_k)
%      y''''(t_k)]

d = matlabFunction(d);
f = matlabFunction(f);


t0 = 0;
tf = 3;
y0 = 1;
h = 1;
M = (tf - t0)/h;
t = t0:h:tf;

y_e = zeros(1,M+1);
y_h = zeros(1,M+1);
y_t = zeros(1,M+1);

y_e(1) = y0;
y_h(1) = y0;
y_t(1) = y0;

y = @(t) 3*exp(-t/2) + t - 2;

N = 4;
C = zeros(1,N);

for i = 1:N
    C(i) = h^(i)/factorial(i);
end
% C = [h, h^2/2!, h^3/3!, h^4/4!]

for j = 1:M
    % Euler's part
    f1_e = f(t(j),y_e(j));
    y_e(j+1) = y_e(j) + h*f1_e;
    % Heun'S part
    f1_h = f(t(j),y_h(j));
    f2_h = f(t(j+1),y_h(j) + h*f1_h);
    y_h(j+1) = y_h(j) + h/2*(f1_h + f2_h);
    % Taylor's part
    y_t(j+1) = y_t(j) + C*d(t(j),y_t(j));   
end

t_exact = t0:0.001:tf;
figure
plot(t,y_e,'r--')
hold on
plot(t,y_h,'g*-')
hold on
plot(t,y_t,'bd-')
plot(t_exact,y(t_exact),'k','Linewidth',3)
legend('Euler','Heun','Taylor','Exact')

