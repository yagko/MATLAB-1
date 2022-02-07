% ME 303 - Lab Session 8 - December 7th, 2021
% Solution of Ordinary Differential Equations
% Heun's Method 2

clear all
close all
clc

t0 = 0;
tf = 3;
y0 = 1;
h1 = 1;
h2 = 0.25;
M1 = (tf - t0)/h1;
M2 = (tf - t0)/h2;

t1 = t0:h1:tf;
t2 = t0:h2:tf;

f = @(t,y) (t-y)/2;
y = @(t) 3*exp(-t/2) + t - 2;

y_h1 = zeros(1,M1+1);
y_h2 = zeros(1,M2+1);

y_h1(1) = y0;
y_h2(1) = y0;


for j = 1:M1
    f1_h = f(t1(j),y_h1(j));
    f2_h = f(t1(j+1),y_h1(j) + h1*f1_h);
    y_h1(j+1) = y_h1(j) + h1/2*(f1_h + f2_h);
end

for j = 1:M2
    f1_h = f(t2(j),y_h2(j));
    f2_h = f(t2(j+1),y_h2(j) + h2*f1_h);
    y_h2(j+1) = y_h2(j) + h2/2*(f1_h + f2_h);
end

t = t0:0.001:tf;
figure
plot(t1,y_h1,'r--')
hold on
plot(t2,y_h2,'g*-')
hold on
plot(t,y(t),'k','Linewidth',3)
legend('h1=1s','h2=0.5s','Exact')

FGE_1 = abs(y(tf) - y_h1(end))
FGE_2 = abs(y(tf) - y_h2(end))

ratio = FGE_1/FGE_2



