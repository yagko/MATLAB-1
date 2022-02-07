clc;clear all;close all;

f1 = @(x,y,z) x.^2 + y.^2 + z.^2 - x - 5;
f2 = @(x,y,z) x.^2 + y.^2 + z.^2 - y - 4;
f3 = @(x,y,z) x.^2 + y.^2 + z.^2 + z - 6;

F = @(x,y,z) [f1(x,y,z);f2(x,y,z);f3(x,y,z)];

% partial derivatives of fi wrt x,y,z
f1_x = @(x,y,z) 2*x-1;
f1_y = @(x,y,z) 2*y;
f1_z = @(x,y,z) 2*z;
f2_x = @(x,y,z) 2*x;
f2_y = @(x,y,z) 2*y-1;
f2_z = @(x,y,z) 2*z;
f3_x = @(x,y,z) 2*x;
f3_y = @(x,y,z) 2*y;
f3_z = @(x,y,z) 2*z+1;

% define jacobian matrix
J = @(x,y,z) [f1_x(x,y,z) f1_y(x,y,z) f1_z(x,y,z)
              f2_x(x,y,z) f2_y(x,y,z) f2_z(x,y,z)
              f3_x(x,y,z) f3_y(x,y,z) f3_z(x,y,z)];

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

f_1 = f1(P1(1),P1(2),P1(3))
f_2 = f2(P1(1),P1(2),P1(3))
f_3 = f3(P1(1),P1(2),P1(3))