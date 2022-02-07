function C = lsline_linear(X,Y)
% write a function to find the constants
% A and B for the least squares line 

% calculate N
N = length(X);

% calculate the coefficient matrix
% [sum(x^2) sum(x)
% sum(x)   N]

A = [sum(X.^2) , sum(X)
     sum(X)    ,  N];

% calculate the output vector
B = [sum(X.*Y) ; sum(Y)];

C = A\B;

end