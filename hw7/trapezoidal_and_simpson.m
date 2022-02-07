clear all
clc
% given function
f= @(x) 2*pi*exp(-x).*sqrt(exp(-2*x)+1);
% given constants
a=0; 
b=1;
k=8; % constant for iteration (should be even)
h=0.125; % step size
s=f(a)+f(b); 
l=0;

% iteration for trapezoidal rule A = h*(f(a)+f(b))/2+h*l
% I calculated l through it
for i = 1:(k-1)
    x = a + h*i;
    l = l + f(x);
end

% iterations for simpsons rule A = 1/3*h*s
% I calculated s through them
for i = 1:2:k-1;
    s = s + 4*f(a+i*h);
end
for i=2:2:k-2;
    s = s + 2*f(a+i*h);
end

A_matlab = integral(f,a,b) % built-in integration function which takes the funciton and the limit values
A_trapez = h*(f(a)+f(b))/2+h*l 
A_simpson = h/3 * s % simpson's rule is closer to the exact value