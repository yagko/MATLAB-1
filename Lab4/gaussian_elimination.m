clc;clear all;close all;

% solution of system of linear equations

% coefficient matrix
A = [1 2 1 4
     2 0 4 3
     4 2 2 1
    -3 1 3 2];

% output vector
B = [13 ; 28 ; 20 ; 6];

% elimination process
Aug_prime = elimination(A,B)

% backtracing process
X = backtracing(Aug_prime)

% ----------- built-in function ---------------
X_matlab = A\B % gives soln to AX=B
% ---------------------------------------------