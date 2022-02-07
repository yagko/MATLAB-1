% ME 303 - Numerical Differentiation
% 30.11.2021
% -------------------------------------------------

clear all
close all
clc

format long

% f = @(x) 60*x.^(45) - 32*x.^(33) + 233*x.^(5) - 47*x.^(2) - 77;
% x = 1/sqrt(3);
f = @(x) exp(x);
x = 0.5;

max1 = 15;

h = 1;
d = (f(x + h) - f(x - h))/(2*h);
% d = (-f(x + 2*h) + 8*f(x + h) - 8*f(x - h) + f(x - 2*h))/(12*h);
H = h;
D = d;
E = 0;

for i = 2:3
    h = h/10;
    H = [H,h];
    d = (f(x + h) - f(x - h))/(2*h);
%     d = (-f(x + 2*h) + 8*f(x + h) - 8*f(x - h) + f(x - 2*h))/(12*h);
    D = [D,d];
    e = abs(D(i) - D(i - 1));        
    E = [E,e];
end

while (E(i-1) > E(i)) && (i < max1)
    i = i + 1;
    h = h/10;
    H = [H,h];
    d = (f(x + h) - f(x - h))/(2*h);
%     d = (-f(x + 2*h) + 8*f(x + h) - 8*f(x - h) + f(x - 2*h))/(12*h);
    D = [D,d];
    e = abs(D(i) - D(i - 1));     
    E = [E,e];
end

L = [H' D' E']
n_best = length(D)-1
H_best = H(n_best)
D_best = D(n_best)
D_exact = f(x)

