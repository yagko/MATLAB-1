clc;clear all;close all

data = [1 0.6
        2 1.9
        3 4.3
        4 7.6
        5 12.6];

N = length(data);

% define RMSE function
E_2 = @(x) sum(1/N*(x(1).*exp(x(2).*data(:,1)) - data(:,2)).^2);

% specify the initial guess for fminsearch
x0 = [0.363381 , 0.747534];

[coef,E_2] = fminsearch(E_2 , x0);

% calculate the RMSE
E_rms = sqrt(E_2);

% calculate the coefficients
C = coef(1);
A = coef(2);

expfit = @(x) C.*exp(A.*x);

fprintf('C = %10f\nA = %10f\nRMSE = %10f\n',C,A,E_rms);

x = 0:0.01:6;

plot(x,expfit(x))
grid on
hold on
scatter(data(:,1),data(:,2),'rs','filled')
legend('Non-linear Least Square Fit','Data')