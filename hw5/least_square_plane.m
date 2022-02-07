% clear all
% close all
% clc
K   =  [1.0 1.3  8.9
        1.2 2.2  9.8
        2.1 1.7 10.7
        2.3 2.4 11.4
        2.7 3.3 13.1
        3.7 2.5 13.6
        2.9 0.9 11.2
        3.1 1.8 12.4
        0.8 2.9 10.9
        1.4 1.6 10.4];

X = K(:,1); % x values
Y = K(:,2); % y values
Z = K(:,3); % z values

N = length(X);

A_sqd = [sum(X.^2) sum(X.*Y) sum(X)
         sum(X.*Y) sum(Y.^2) sum(Y)
         sum(X)    sum(Y)      N   ]; % calculating the coefficient A 

B_sqd = [sum(X.*Z)
         sum(Z.*Y)
         sum(Z)]; % calculating the coefficient B 

C = A_sqd\B_sqd; % getting the resultant matrix

% attend every column of matrix C to A,B,C vectors
A = C(1) 
B = C(2)
C = C(3)

% plane plotting
cfit=@(x,y) A*x+B*y+C;
[x,y]=meshgrid(0:0.1:5,0:0.1:5);
surf(x,y,cfit(x,y))
grid on
hold on
plot3(X,Y,Z,'R*')