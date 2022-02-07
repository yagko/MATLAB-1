clc;clear all;close all;

% CHAPTER 4 - LINEAR INTERPOLATION

% define the function 
f = @(x) 2*sin(pi/6*x);

% define the 4th derivative wrt x
f4 = @(x) 2*(pi/6)^4*sin(pi/6*x);

% define spacing, h
h = 2;
% data points
X = 0 : h : 6;

% calculate
% N+1 data points -> length(X)
N = length(X) - 1;

% evaluate our function at X
Y = f(X);

% calculate L_k(x)
for k=1:N+1
    product = 1;
    for j=1:N+1
        if k ~= j
            product = conv(product,poly(X(j))/(X(k) - X(j)));
% conv function gives the coefficients of the polynomial
% created by the multiplication of the terms inside
% poly function is x-x(j) polynomial
        end
    end
    C(k,:) = product;
end

% calculate A 
A = Y*C;

p_n = 0;
syms x
for i = 0:N
    p_n = p_n + A(i+1)*x.^(N-i);
end
p_n = matlabFunction(p_n);

% compare the polynomial app. and the exact function
% this interval should contain all the initial points scarcly
x = 0: 0.01 : 6;
subplot(2,1,1)
plot(x,f(x))
grid on
hold on
plot(x,p_n(x))
legend('f(x)','p_n(x)')

% define error
E = abs(f(x) - p_n(x));
subplot(2,1,2)
plot(x,E)

% error calculations
% maximum error should be less than or equal to the error bound
max_error = max(E)
M4 = max(f4(x));
error_bound = h^4*M4/24