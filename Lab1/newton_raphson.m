clear all; close all; clc;

% p_k = p_(k-1) - f(p(k_1))/f'(p_(k-1))

A = 100;
f = @(x) x^3 - A;
df = @(x) 3*x^2;


delta = 10^(-6);
eps = 10^(-8);
p0 = 30;
maxi = 50;
for k=1:maxi
    p1 = p0(k) - f(p0(k))/df(p0(k));
    err = abs(p1-p0(k));
    relerr = err/(abs(p1) + eps);
    p0 = [p0 p1];
    if (err < delta) && (relerr < delta)
        break
    end
end

p0'