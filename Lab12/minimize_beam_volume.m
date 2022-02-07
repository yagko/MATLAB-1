clear all
close all
clc

global L E P sigma_max delta_max

L = 1;
E = 200e9;

P = 100;
sigma_max = 300e6;
delta_max = 0.01;

t0 = L/15;
w0 = L/200;

fun = @beam_volume;
x0 = [t0,w0];
% A = [1,-20];
% b = 0;
A = [];
b = [];
Aeq = [];
beq = [];
lb = [L/100,L/100];
ub = [L/10,L/10];
nonlcon = @constraints;

x = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon);

t = 1000*x(1)
w = 1000*x(2)

% I = w*t^3/12;
% EI = E*I;
% sigma_b = P*L*t/(2*I)
% delta = P*L^3/(3*EI)
% 
% 
% c1 = sigma_b - sigma_max
% c2 = delta - delta_max

