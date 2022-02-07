clear all
close all
clc

format long

width = 16;
length = 24;


f= @(x) x.^3-20*x.^2+96*x-125;

p0 = 1;
p1 = 2;

P = [p0,p1];
Y = [f(p0),f(p1)];
delta = 10^(-6);
maxI = 100;

for k= 3:maxI
    
    p_approx = P(k-1) - f(P(k-1))*(P(k-1)-P(k-2))/(f(P(k-1))-f(P(k-2)));
    
    P = [P,p_approx];
    
    err = abs(P(k) - P(k-1));
    if (err < delta)
        break
    end
end

if p_approx>= width-2*p_approx
    disp('Try different initial guesses')
else
    root = p_approx
    vol = (length-2*root)*(16-2*root)*root
end
