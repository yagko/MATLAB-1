clc;clear all;close all

data = [1 0.6
        2 1.9
        3 4.3
        4 7.6
        5 12.6];

% conduct variable change
X = data(:,1);
Y = log(data(:,2));

% calculate the coefficients A and B 
coef = lsline_linear(X,Y);
A = coef(1);
B = coef(2);

C = exp(B);

expfit = @(x) C.*exp(A.*x);

N = length(X);
E_rms = sqrt (sum(1/N*(expfit(X)-Y).^2));

fprintf('C = %10f\nA = %10f\nRMSE = %10f\n',C,A,E_rms);

x = 0:0.01:6;

plot(x,expfit(x))
grid on
hold on
scatter(data(:,1),data(:,2),'rs','filled')
legend('Least Square Line','Data')