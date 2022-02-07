clc;clear all;close all;

data = [ -1 10
          0  9
          1  7
          2  5
          3  4
          4  3
          5  0
          6 -1];

% define the degree of the least squares polynomial
M = 7; %cubic function if M = 3

[E_rms,  C] = lspol(data(:,1),data(:,2), M);

for i = 1:M+1
    fprintf('C%d = %10f\n' , i , C(i));
end

fprintf('RMSE = %10f\n',E_rms);

% built in function
% cftool(data(:,1),data(:,2))
