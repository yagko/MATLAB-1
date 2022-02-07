clear all; clc;

% there are two roots (1 and 2)
% by changing the initial guesses the other root can be found

f = @(x) x.^2 - 3*x + 2;

p0 = 2.5;
p1 = 5;

P = [p0,p1];
Y = [f(p0),f(p1)];
delta = 10^(-6);
maxI = 100;

for k = 3:maxI
    p_approx = P(k-1) - f(P(k-1))*(P(k-1)-P(k-2))/(f(P(k-1)) - f(P(k-2)));
    P = [P,p_approx];
    y = f(p_approx);
    Y = [Y,y];
    err = abs(P(k) - P(k-1));

    if (err < delta) && (abs(y) < delta)
        break
    end
end

P = P';
Y = Y';
index = (1:1:k)';

Data = [index,P,Y];
Varnames = {'k','p','y'};
T = table(Data(:,1),Data(:,2),Data(:,3),'VariableNames',Varnames)