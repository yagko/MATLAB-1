clear all; clc;

% solution for nonlinear equations

p_s = 0.710;
p_w = 1.000;

R = 0.15;
V_s = 4/3*pi*R^3;
m_s = p_s*V_s;

V_out = @(x) pi*((R-x)^2)*(2*R+x)/3;
V_in = @(x) V_s - V_out(x);
f = @(x) m_s - p_w*V_in(x);

a = 0.01;
b = 0.14;

delta = 10^(-6);
maxI = 100;
c = b - f(b)*(b-a)/(f(b)-f(a));
dx = min(abs(c-a),abs(c-b));
yc = f(c);

while (dx > delta) || (abs(yc) > delta)
    c = b - f(b)*(b-a)/(f(b)-f(a));
    dx = min(abs(c-a),abs(c-b));
    yc = f(c);

    if f(b)*f(c) > 0

        b = c;

    elseif f(a)*f(c) > 0

        a = c;
    end
end

dx
abs(yc)
x = c
V_sub = V_in(c)
ratio = V_sub/V_s