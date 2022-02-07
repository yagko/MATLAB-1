% ME 303 - Numerical Differentiation (Richardson's Extrapolation)
% 30.11.2021
% -------------------------------------------------

% D = [D_11 0    0    0
       D_21 D_22 0    0
       D_31 D_32 D_33 0
       D_41 D_42 D_43 D_44]

D_22 = (4^(k)*D_21 - D_11)/(4^(k) - 1)


clear all
clc

% f=@(x) 60*x.^(45)-32*x.^(33)+233*x.^(5)-47*x.^(2)-77;
% x=1/sqrt(3);

f = @(x) exp(x);
x = 0.5;

delta = 10^-6;
err = 1;
relerr = 1;
h = 1;
j = 1;
D(1,1) = (f(x + h) - f(x - h))/(2*h);

maxI = 50;

while (relerr > delta) && (err > delta) && (j < maxI)
    
    h = h/2;
    D(j+1,1) = (f(x + h) - f(x - h))/(2*h);
    
    for k = 1:j
        D(j+1,k+1) = D(j+1,k) + (D(j+1,k) - D(j,k))/((4^k) - 1);
    end
    
    err = abs(D(j+1,j+1) - D(j,j));
    relerr = err/(abs(D(j+1,j+1)) + eps);
    
    j = j + 1;
    
end

[n,n]=size(D);
n
D
D_approx=D(n,n)