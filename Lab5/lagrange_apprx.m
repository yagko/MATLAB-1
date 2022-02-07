clc;clear all;close all;

% CHAPTER 4 - LINEAR INTERPOLATION

% define the function 
f = @(x) 2*sin(pi/6*x);

% data points
X = [0,1,3,5];

% calculate
% N+1 data points -> length(X)
N = length(X) - 1;

% evaluate our function at X
Y = f(X);

% assign an initial value
p_n = 0;
syms x

% calculate L_k(x)
for k=1:N+1
    L_k = 1;
    for j=1:N+1
        if k ~= j
            L_k = L_k*(x-X(j)) / (X(k)-X(j));
        end
    end
    % calculate p_n(x) = f(x_k)L_k(x)
    p_n = p_n + Y(k)*L_k;
end

p_n = matlabFunction(p_n);


% compare the polynomial app. and the exact function
% this interval should contain all the initial points scarcly
x = 0: 0.01 : 5;
plot(x,f(x))
grid on
hold on
plot(x,p_n(x))
legend('f(x)','p_n(x)')