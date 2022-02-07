clc;clear all;close all;

data = [ -1 10
          0  9
          1  7
          2  5
          3  4
          4  3
          5  0
          6 -1];

C = lsline(data(:,1),data(:,2));

A = C(1);
B = C(2);

cfit = @(x) A.*x + B;

x= -2:0.01:10;
plot(x,cfit(x))
grid on
hold on
scatter(data(:,1),data(:,2),'rs' , 'filled')
legend('Least Square Line', 'Data')