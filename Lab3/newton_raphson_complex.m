clear all; close all; clc;

% p_k = p_(k-1) - f(p(k_1))/f'(p_(k-1))
p = [1,0,0,-1];
r = roots(p)


f = @(x) x^3 - 1;
df = @(x) 3*x^2;

delta = 10^(-6);
p0 = complex(1,-6);
maxi = 50;
for k=1:maxi
    p1 = p0(k) - f(p0(k))/df(p0(k));
    err = norm(p1-p0(k));
    relerr = err/(norm(p1) + delta);
    p0 = [p0;p1];
    if (err < delta) && (relerr < delta)
        break
    end
end

err
relerr
p0