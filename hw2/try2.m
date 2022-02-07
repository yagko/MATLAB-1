clc;clear all;close all;

l = 24;
w = 16;
root0 = 2;
root1 = 3;
delta = 10^(-6);

f= @(x) x.^3-20*x.^2+96*x-125;

root = [root0,root1];
Y = [f(root0), f(root1)];

maxI = 50;

for k = 3:maxI
    p_approx = root(k-1) - f(root(k-1))*(root(k-1)-root(k-2))/(f(root(k-1)) - f(root(k-2)));
    root = [root,p_approx];
    y = f(p_approx);
    Y = [Y,y];
    err = abs(root(k) - root(k-1));
    if (err < delta)
        break
    end
end

relerr = (w-(2*root(end)))/2;
if root(end) > relerr
    disp('Try different initial guesses')
end

root_end = root(end)
vol = (l-2*root_end)*(w-2*root_end)*root_end