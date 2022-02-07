% clear all
% close all
% clc
% The units of x and y are in meters.
x=[ 0.9958    0.9978    0.9040    0.7209    0.4571    0.1240   -0.2638 ...
   -0.6891   -1.1324   -1.5728   -1.9888   -2.3590   -2.6634   -2.8838 ...
   -3.0053   -3.0162   -2.9094   -2.6823   -2.3373   -1.8820   -1.3286 ...
   -0.6943   -0.0000    0.7299    1.4686    2.1878    2.8588    3.4536 ...
    3.9461    4.3133    4.5359    4.5998    4.4962    4.2227    3.7830 ...
    3.1874    2.4527    1.6013    0.6608   -0.3368   -1.3564   -2.3608 ...
   -3.3124   -4.1743   -4.9116   -5.4933   -5.8931   -6.0907   -6.0726 ...
   -5.8331   -5.3746   -4.7073   -3.8498   -2.8281   -1.6748   -0.4287    0.8674    2.1675];

y=[ 0.3160    0.4685    0.6221    0.7660    0.8901    0.9855    1.0442    1.0605 ...
    1.0299    0.9506    0.8226    0.6484    0.4327    0.1822   -0.0945   -0.3872 ...
   -0.6845   -0.9744   -1.2445   -1.4828   -1.6779   -1.8199   -1.9003   -1.9132 ...
   -1.8547   -1.7237   -1.5221   -1.2546   -0.9284   -0.5537   -0.1427    0.2905 ...
    0.7305    1.1607    1.5648    1.9265    2.2308    2.4642    2.6155    2.6764 ...
    2.6413    2.5085    2.2796    1.9600    1.5585    1.0875    0.5621    0.0000 ...
   -0.5792   -1.1548   -1.7054   -2.2102   -2.6494   -3.0050   -3.2616   -3.4070   -3.4331   -3.3354];

mass = 0.025; 
m = length(x); 
% defining and assigning zeros to my vectors beforehand
vx = zeros(1,m);       
vy = zeros(1,m);
ax = zeros(1,m);
ay = zeros(1,m);
fx = zeros(1,m);
fy = zeros(1,m);

% forward difference formula for the first term of ax
ax(1) = (2*x(1)-5*x(2)+4*x(3)-x(4))/(0.05^2);
ay(1) = (2*y(1)-5*y(2)+4*y(3)-y(4))/(0.05^2);

% central formula for ax
for i = 2:1:m-1
    ax(i) = (x(i+1)-2*x(i)+x(i-1))/(0.05^2);
    ay(i) = (y(i+1)-2*y(i)+y(i-1))/(0.05^2);
end

% backward difference formula for the last term of ax
ax(m) = (2*x(m)-5*x(m-1)+4*x(m-2)-x(m-3))/(0.05^2);
ay(m) = (2*y(m)-5*y(m-1)+4*y(m-2)-y(m-3))/(0.05^2);

% force = (accelaration)x(mass)
fx = ax.*mass;
fy = ay.*mass;

% forward difference formula for the first term of vx
vx(1) = (-3*x(1)+4*x(2)-x(3))/0.1; 
vy(1) = (-3*y(1)+4*y(2)-y(3))/0.1;

% central formula for vx
for i = 2:1:m-1                
    vx(i) = (x(i+1)-x(i-1))/0.1;
    vy(i) = (y(i+1)-y(i-1))/0.1;
end

% backward difference formula for the last term of vx
vx(m) = (3*x(m)-4*x(m-1)+x(m-2))/0.1;
vy(m) = (3*y(m)-4*y(m-1)+y(m-2))/0.1;


figure(1);
scatter(x,y,20,"black","filled"); % (respectively) first data, second data, size of the circle, color of the circle, filled or not option
hold on
quiver(x,y,vx,vy,0); % built-in quiver func for vectors
title("Particle Velocity Field")
xlabel("x");
ylabel("y");

figure(2);
scatter(x,y,20,"magenta","filled"); % (respectively) first data, second data, size of the circle, color of the circle, filled or not option
hold on
quiver(x,y,fx,fy,0); % built-in quiver func for vectors
title("Particle Force Field")
xlabel("x");
ylabel("y");