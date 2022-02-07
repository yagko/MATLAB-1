clc;clear all;close all;

% DEFINE CONSTANTS
% -------------------------------------------------------------------------
% Masses (in kg)
m1 = 1;
m2 = 1;
% Lengths (in cm) 
L1 = 0.1;
L2 = 0.1;
g = 9.81;
u2 = m2/(m2+m1);

% Define initial conditions for theta
theta1_init = pi/2;
theta2_init = pi/2;
diff_theta1_init = 0;
diff_theta2_init = 0;

%Initial conditions for z
z1_init = theta1_init;
z2_init = diff_theta1_init;
z3_init = theta2_init;
z4_init = diff_theta2_init;

% Define time values and step size for rk4 function
ti=0;
tf=2;
h=1/200;


% RK4 METHOD
% -------------------------------------------------------------------------

% z1 = theta1
% z2 = first derivative of theta1
% z3 = theta2
% z4 = first derivative of theta2

% Define matrix A (all of its elements are functions)
A = { @(t, z1, z2, z3, z4) L1              , @(t, z1, z2, z3, z4) u2*L2*cos(z1-z3);
      @(t, z1, z2, z3, z4) L1*(cos(z1-z3)) , @(t, z1, z2, z3, z4)  L2              };

% Define matrix B (all of its elements are functions)
B = { @(t, z1, z2, z3, z4) -g*sin(z1) - u2*L2*sin(z1-z3)*z4^2;
      @(t, z1, z2, z3, z4) -g*sin(z3) + L1*sin(z1-z3)*z2^2    };

% First order differentials
d_z1 = @(t, z1, z2, z3, z4) z2;
d_z2 = @(t, z1, z2, z3, z4) ( A{2,2}(t, z1, z2, z3, z4) .* B{1}(t, z1, z2, z3, z4) -   ...
                              A{1,2}(t, z1, z2, z3, z4) .* B{2}(t, z1, z2, z3, z4) ) / ...
                            ( A{1,1}(t, z1, z2, z3, z4) .* A{2,2}(t, z1, z2, z3, z4) - ...
                              A{1,2}(t, z1, z2, z3, z4) .* A{2,1}(t, z1, z2, z3, z4) );
d_z3 = @(t, z1, z2, z3, z4) z4;
d_z4 = @(t, z1, z2, z3, z4) ( -A{2,1}(t, z1, z2, z3, z4) * B{1}(t, z1, z2, z3, z4) +   ...
                               A{1,1}(t, z1, z2, z3, z4) * B{2}(t, z1, z2, z3, z4) ) / ...
                            (  A{1,1}(t, z1, z2, z3, z4) * A{2,2}(t, z1, z2, z3, z4) - ...
                               A{1,2}(t, z1, z2, z3, z4) * A{2,1}(t, z1, z2, z3, z4) );

% Calling "RK4_fun" and appending the corresponding values
y = RK4_fun({d_z1;d_z2;d_z4;d_z4} , [z1_init,z2_init,z3_init,z4_init], [ti tf], h);

theta1_rk4 = y(:,1);
theta2_rk4 = y(:,3);

% Position of the first bob (x and y coordinates)
x1_rk4 = sin(theta1_rk4)*L1;
y1_rk4 = -cos(theta1_rk4)*L1;

% Position of the second bob (x and y coordinates)
x2_rk4 = L1*sin(theta1_rk4) + sin(theta2_rk4)*L2;
y2_rk4 = -L1*cos(theta1_rk4) - cos(theta2_rk4)*L2;


% ODE45 METHOD
% -------------------------------------------------------------------------

[t1_ode45 , theta_ode45] = ode45(@( t , x) odefun(t,x,A,B) , [ti tf] , [z1_init,z2_init,z3_init,z4_init] );

theta1_ode = theta_ode45(:,1);
theta2_ode = theta_ode45(:,3);

% Position of the first bob (x and y coordinates)
x1_ode = sin(theta1_ode)*L1;
y1_ode = -cos(theta1_ode)*L1;

% Position of the second bob (x and y coordinates)
x2_ode = L1*sin(theta1_ode) + sin(theta2_ode)*L2;
y2_ode = -L1*cos(theta1_ode) - cos(theta2_ode)*L2;


% PLOT
% -------------------------------------------------------------------------
figure("Name", "RK4 Method")
scatter(x1_rk4,y1_rk4,"filled")
grid on
hold on
scatter(x2_rk4, y2_rk4, "filled")
title("Double Pendulum via RK4 Method")
xlabel("x(meters)")
ylabel("y(meters)")
xlim([-0.25 0.25])
ylim([-0.25 0.25])
axis equal
legend("First Bob","Second Bob")

figure("Name", "ode45 Method")
scatter(x1_ode,y1_ode,"filled")
grid on
hold on
scatter(x2_ode,y2_ode,"filled")
title("Double Pendulum via ode45 Method")
xlabel("x(meters)")
ylabel("y(meters)")
xlim([-0.25 0.25])
ylim([-0.25 0.25])
axis equal
legend("First Bob","Second Bob")

% Animation
figure(3)
  for i = 1:1:length(x1_ode)
    
    clf % clear the figure
    
    plot(0, 0,'m.',x1_ode(i),y1_ode(i),'k.',x2_ode(i),y2_ode(i),'k.','markersize',40);
    xlim([-0.25 0.25])
    ylim([-0.25 0.25])
    axis equal
    grid on
    line([0 x1_ode(i)], [0 y1_ode(i)],'Linewidth',2);
    line([x1_ode(i) x2_ode(i)], [y1_ode(i) y2_ode(i)],'linewidth',2,'color','r');
    drawnow

  end






function dx_dt = odefun( t , x ,A ,B)



dx_dt = zeros( 4 , 1 );
dx_dt( 1 ) = x( 2 );
dx_dt( 2 ) = (A{2,2}(t, x(1), x(2), x(3), x(4)) .* B{1}(t, x(1), x(2), x(3), x(4)) - ...
                            A{1,2}(t, x(1), x(2), x(3), x(4)) .* B{2}(t, x(1), x(2), x(3), x(4)) )/...
                            ...
                          ( A{1,1}(t, x(1), x(2), x(3), x(4)) .* A{2,2}(t, x(1), x(2), x(3), x(4)) ...
                          - A{1,2}(t, x(1), x(2), x(3), x(4)) .* A{2,1}(t, x(1), x(2), x(3), x(4)) );

dx_dt( 3 ) = x( 4 );
dx_dt( 4 ) = ( -A{2,1}(t, x(1), x(2), x(3), x(4)) * B{1}(t, x(1), x(2), x(3), x(4)) + ...
                            A{1,1}(t, x(1), x(2), x(3), x(4)) * B{2}(t, x(1), x(2), x(3), x(4)) )/...
                            ...
                          ( A{1,1}(t, x(1), x(2), x(3), x(4))*A{2,2}(t, x(1), x(2), x(3), x(4)) ...
                          - A{1,2}(t, x(1), x(2), x(3), x(4))*A{2,1}(t, x(1), x(2), x(3), x(4)));

end