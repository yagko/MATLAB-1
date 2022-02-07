function [c,ceq] = constraints(x)

global L E P sigma_max delta_max

t = x(1);
w = x(2);

I_zz = w*t^3/12;
M = P*L;
sigma_b = M*(t/2)/I_zz;
delta = (P*L^3)/(3*E*I_zz);

c1 = sigma_b - sigma_max;
c2 = delta - delta_max;

c = [c1,c2];
ceq = [];

end