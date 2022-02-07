clear all; close all; clc;

x = 0:1:20;

[Y] = my_function(x);

y1 = Y(1,:);
y2 = Y(2,:);

plot(x,y1)
hold on
plot(x,y2)
legend('y1', 'y2')
xlabel('x')
ylabel('y')