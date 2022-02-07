clear all; close all; clc;

% f = -exp(x)-2-x

a = -3; 
b=1;
delta = 10^(-3);

f = @(x) -exp(x)-2-x;

x = a:0.2:b;
y = f(x);
plot(x,y)

ya = f(a);
yb = f(b);

if ya*yb > 0
    disp('ya*yb > 0:')
    return
end

maxI = 100;

for i =1:maxI

    c = (a+b)/2;
    yc =f(c);

    if yc == 0

        a = c;
        b = c;

    elseif yb*yc > 0

        b = c;
        yb = yc;

    else

        a = c;
        ya = yc;

    end

    err = abs(b-a);

    if err < delta
        break
    end
end

root = c
err