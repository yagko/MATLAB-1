clc;clear all;close all;

A = [   3  7  16  2   2
        2  4   8 24   8
        1 15   5  6   2
        2  1 -10  7 -21
      -13  1  -3 -4  -2];

K = [0 0 0 0 1
     0 0 1 0 0
     1 0 0 0 0
     0 1 0 0 0
     0 0 0 1 0];

B = [14;23;4;18;11];

Aug = [A B];

S = K*Aug;


N = length(B);

maxI = 100;
delta = 10^(-6);

X = zeros(N,1);
P = [1;2;3;4;5];

for i = 1:maxI
    for j = 1:N
        X_2 = X(j);
        X(j) = S(j,end) - S(j,1:j-1)*X(1:j-1) - S(j,j+1:N)*X(j+1:N) / S(j,j);
        err(j) = abs(X(j) - X_2);
        if (err < delta) %tolerance condition checking
            break
        end
    end
end

X(:,end)
X_star = [-1.322644867279853
          -0.355330111009888
           1.295315796809742
           1.111456782021438
          -1.246360677443388]