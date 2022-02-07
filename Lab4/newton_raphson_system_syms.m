clc;clear all;close all;

syms x y z
f1 = x^2 + y^2 + z^2 - x - 5;
f2 = x^2 + y^2 + z^2 - y - 4;
f3 = x^2 + y^2 + z^2 + z - 6;

F = [ f1;f2;f3];

f1_x = diff(f1,x);
f1_y = diff(f1,y);
f1_z = diff(f1,z);
f2_x = diff(f2,x);
f2_y = diff(f2,y);
f2_z = diff(f2,z);
f3_x = diff(f3,x);
f3_y = diff(f3,y);
f3_z = diff(f3,z);

J = [ f1_x f1_y f1_z
      f2_x f2_y f2_z
      f3_x f3_y f3_z];

F = matlabFunction(F);
J = matlabFunction(J);

% specify our initial guess
P0 = [1;1;10];
maxI = 1000;
eps = 1E-6;

for i = 1:maxI
    P1 = P0 - J( P0(1),P0(2),P0(3)) \ F(P0(1),P0(2),P0(3));
    err = norm(P1-P0);
    relerr = err / (norm(P1) + eps);
    Y = norm(F(P1(1),P1(2),P1(3)));
    P0 = P1;
    if (err < eps) && (relerr < eps) && ( Y< eps)
        break
    end
end

P1

f1 = matlabFunction(f1);
f2 = matlabFunction(f2);
f3 = matlabFunction(f3);

f_1 = f1(P1(1),P1(2),P1(3))
f_2 = f2(P1(1),P1(2),P1(3))
f_3 = f3(P1(1),P1(2),P1(3))