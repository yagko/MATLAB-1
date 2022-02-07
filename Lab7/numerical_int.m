clear all
clc


f = @(x) 2 + sin(2*sqrt(x));

a = 1;
b = 6;
M = 10;
h = (b-a)/M;

s = 0;

for k = 1:(M-1)
    x = a + h*k;
    s = s + f(x);
end

T = h*(f(a)+f(b))/2+h*s
I = integral(f,a,b)