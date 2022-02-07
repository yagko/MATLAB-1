clear all; clc;

% system of linear equations

% 4x - y + z = 7
% 4x -8y + z = -21
% -2x + y + 5z = 15
% AX = B

% A matrix should be strictly diagonally dominant

A = [4 -1 1
     4 -8 1
     -2 1 5];

B = [7;-21;15];
N = length(B);

maxI = 100;
delta = 1e-6;

P0 = [1;2;2];
P = P0;
X = zeros(N,maxI);
for k=1:maxI % until it converges
    for j = 1:N % number of variables
        X(j,k) = (B(j) - A(j,[1:j-1,j+1:N])*P([1:j-1,j+1:N]))/A(j,j);
    end
    err = norm(X(:,k)-P);
    relerr = err/(norm(X(:,k) + delta));
    P = X(:,k);
    if (err < delta) &&  (relerr < delta)
        break
    end
end

k = (0:1:k)';
x_k = [P0(1) X(1,1:k(end))]';
y_k = [P0(2) X(2,1:k(end))]';
z_k = [P0(3) X(3,1:k(end))]';

Data = [k,x_k,y_k,z_k]