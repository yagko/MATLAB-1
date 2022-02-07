% Please do not use the commands below in your code.
clear all
close all
clc
% Matrix form of equations
A = [3 7 16 2 2
     2 4 8 24 8 
     1 15 5 6 2 
     2 1 -10 7 -21 
    -13 1 -3 -4 -2];

B = [14;23;4;18;11];
Aug = [A B]; % Augmented matrix

% Permutation matrix
% In order to correcctly use the gauss-seidel method augmented matrix
% should be rearrenged by starting from the left hand side each column's 
% highest number (as absolute value) should be at the diagonal
Per = [0 0 0 0 1; 0 0 1 0 0; 1 0 0 0 0; 0 1 0 0 0; 0 0 0 1 0 ];

%Diagonally Dominant Matrix
S = Per*Aug;

tolerance_value = 10^-6; % error value
maxI = 20; % maximum number of iterations

% Predefining the matrices
X_star = zeros(5,1); %since the matrix A is a 5-5 matrix
err = zeros(5,1);

%Gauss-Seidel Method
for i = 1:maxI
    for j = 1:5
        X_2 = X_star(j);
        k = S(j,end) - S(j,1:j-1)*X_star(1:j-1) - S(j,j+1:5)*X_star(j+1:5);
        X_star(j) = k / S(j,j);
        err(j) = abs(X_star(j) - X_2);
        if (err < tolerance_value) %tolerance condition checking
            break
        end
    end
end
format long

X_star