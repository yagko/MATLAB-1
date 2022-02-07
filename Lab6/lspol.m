function [E_rms , C] = lspol(X,Y,M)

% calculate N
N = length(X);

% allocate space for the matrix F
F = zeros(M+1 , N);

% calculate F
for i = 1:M+1
    F(i,:) = X.^(M-i+1);
end

% calculate the coefficient matrix
A = F*F';

% calculate the outpu vector
B = F*Y;

C = A\B;

f = @(x) 0;
for i = 1:M+1
    f = @(x) f(x) + C(i).*x.^(M-i+1);
end

% calculate RMSE
E_rms = sqrt(sum(1/N*(f(X)-Y).^2));

% interval for the plots
x = -2 : 0.01 : 10;

plot(x,f(x))
grid on
hold on
scatter(X,Y,'rs', 'filled')
legend('Least Squares Polynomial', 'Data')

end