% ME 303 - Lab Session 8 - December 7th, 2021
% Solution of Ordinary Differential Equations
% Euler's Method

clear all
close all
clc

t0 = 0;
tf = 3;
y0 = 1;
h1 = 1;
h2 = h1/4;
M1 = (tf - t0)/h1;
M2 = (tf - t0)/h2;

f = @(t,y) (t-y)/2;
y = @(t) 3*exp(-t/2) + t - 2;

t1 = t0:h1:tf;
t2 = t0:h2:tf;

y1 = zeros(1,M1+1);
y2 = zeros(1,M2+1);

y1(1) = y0;
y2(1) = y0;


for j = 1:M2
    if j <=M1
        y1(j+1) = y1(j) + h1*f(t1(j),y1(j));
        y2(j+1) = y2(j) + h2*f(t2(j),y2(j));
    else
        y2(j+1) = y2(j) + h2*f(t2(j),y2(j));
    end
end

t = t0:0.001:tf;
figure
plot(t1,y1,'r--')
hold on
plot(t2,y2,'g*-')
hold on
plot(t,y(t),'k','Linewidth',3)
legend('h1=1s','h2=0.5s','Exact')

FGE_1 = abs(y(tf) - y1(end))
FGE_2 = abs(y(tf) - y2(end))
