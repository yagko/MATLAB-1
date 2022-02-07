clear all; close all; clc

% DEFINE CONSTANTS
% -------------------------------------------------------------------------
% initial and final time
ti = 0;
tf = 17.06521656015796;
%tf = 40;
% initial conditions
y1_init = 0.994;
dy1_init = 0;
y2_init = 0;
dy2_init = -2.0015851063790825;
% normalized mass values 
u1 = 0.012277471; % moon
u2 = 1-u1; % earth


% EULER'S METHOD
% -------------------------------------------------------------------------
tic 
M1 = 24000; % subintervals
h_e = (tf-ti)/M1; % step size
t_e = ti:h_e:tf ; % time array

% Define differential equations (all given)
D1_e = @(t_e, y1, y2) ((y1+u1)^2+y2^2)^(3/2);
D2_e = @(t_e, y1, y2) ((y1-u2)^2+y2^2)^(3/2);
diff1_y1_e = @(t_e, y1, dy1, y2, dy2) dy1;
diff1_y2_e = @(t_e, y1, dy1, y2, dy2) dy2;
diff2_y1_e = @(t_e, y1, dy1, y2, dy2) y1 + 2*dy2 - u2*((y1+u1)/D1_e(t_e,y1,y2)) - u1*((y1-u2)/D2_e(t_e,y1,y2));
diff2_y2_e = @(t_e, y1, dy1, y2, dy2) y2 - 2*dy1 - u2*y2/D1_e(t_e,y1,y2) - u1*y2/D2_e(t_e,y1,y2);

% Preallocate space for solution matrices
y1_euler = zeros(1 , M1+1);
dy1_euler = zeros(1 , M1+1);
y2_euler = zeros(1 , M1+1);
dy2_euler = zeros(1 , M1+1);

% Append the initial values to solution matrices
y1_euler(1) = y1_init;
dy1_euler(1) = dy1_init;
y2_euler(1) = y2_init;
dy2_euler(1) = dy2_init;

% Apply the Euler's method 
for i=1:M1
    
    y1_euler(i+1) = y1_euler(i) + h_e * diff1_y1_e(t_e(i), y1_euler(i), dy1_euler(i), y2_euler(i), dy2_euler(i));
    dy1_euler(i+1) = dy1_euler(i) + h_e * diff2_y1_e(t_e(i), y1_euler(i), dy1_euler(i), y2_euler(i), dy2_euler(i));
    y2_euler(i+1) = y2_euler(i) + h_e * diff1_y2_e(t_e(i), y1_euler(i), dy1_euler(i), y2_euler(i), dy2_euler(i));
    dy2_euler(i+1) = dy2_euler(i) + h_e * diff2_y2_e(t_e(i), y1_euler(i), dy1_euler(i), y2_euler(i), dy2_euler(i)); 
    
end

% Plot for Euler
figure('Name','Euler')
plot(y1_euler,y2_euler,'m');
grid on
legend('Orbit')
title("Orbit by Euler's Method")
xlabel('y1')
ylabel('y2')

time_euler = toc

% RUNGE-KUTTE 4 METHOD
% -------------------------------------------------------------------------
tic
M2 = 6000; % subintervals
h_r=(tf-ti)/M2; % step size
t_r=ti:h_r:tf ; % time array

% Define differential equations (all given)
D1_r = @(t_r, y1, y2) ((y1+u1)^2+y2^2)^(3/2);
D2_r = @(t_r, y1, y2) ((y1-u2)^2+y2^2)^(3/2);
diff1_y1_r = @(t_r, y1, dy1, y2, dy2) dy1;
diff1_y2_r = @(t_r, y1, dy1, y2, dy2) dy2;
diff2_y1_r = @(t_r, y1, dy1, y2, dy2) y1 + 2*dy2 - u2*((y1+u1)/D1_r(t_r,y1,y2)) - u1*((y1-u2)/D2_r(t_r,y1,y2));
diff2_y2_r = @(t_r, y1, dy1, y2, dy2) y2 - 2*dy1 - u2*y2/D1_r(t_r,y1,y2) - u1*y2/D2_r(t_r,y1,y2);

% Preallocate space for solution matrices
y1_rk4 = zeros(1 , M2+1);
dy1_rk4 = zeros(1 , M2+1);
y2_rk4 = zeros(1 , M2+1);
dy2_rk4 = zeros(1 , M2+1);

% Append the initial values to solution matrices
y1_rk4(1) = 0.994;
dy1_rk4(1) = 0;
y2_rk4(1) = 0;
dy2_rk4(1) = -2.0015851063790825;

% Apply RK4 Method
for i=1:M2 
     % first order derivatives
  p1_1 = diff1_y1_r(t_r(i), y1_rk4(i), dy1_rk4(i), y2_rk4(i), dy2_rk4(i)); %defining p1 in terms of partition vector, y1, dy1, y2, and dy2
  q1_1 = diff2_y1_r(t_r(i), y1_rk4(i), dy1_rk4(i), y2_rk4(i), dy2_rk4(i)); %defining q1 in terms of partition vector, y1, dy1, y2, and dy2
  p2_1 = diff1_y2_r(t_r(i), y1_rk4(i), dy1_rk4(i), y2_rk4(i), dy2_rk4(i)); %defining p2 in terms of partition vector, y1, dy1, y2, and dy2
  q2_1 = diff2_y2_r(t_r(i), y1_rk4(i), dy1_rk4(i), y2_rk4(i), dy2_rk4(i)); %defining q2 in terms of partition vector, y1, dy1, y2, and dy2
 
     % second order derivatives
  p1_2 = diff1_y1_r(t_r(i) + 0.5*h_r , y1_rk4(i) + 0.5*h_r*p1_1 , dy1_rk4(i) + 0.5*h_r*q1_1 , y2_rk4(i) + 0.5*h_r*p2_1 , dy2_rk4(i) + 0.5*h_r*q2_1);
  q1_2 = diff2_y1_r(t_r(i) + 0.5*h_r , y1_rk4(i) + 0.5*h_r*p1_1 , dy1_rk4(i) + 0.5*h_r*q1_1 , y2_rk4(i) + 0.5*h_r*p2_1 , dy2_rk4(i) + 0.5*h_r*q2_1);
  p2_2 = diff1_y2_r(t_r(i) + 0.5*h_r , y1_rk4(i) + 0.5*h_r*p1_1 , dy1_rk4(i) + 0.5*h_r*q1_1 , y2_rk4(i) + 0.5*h_r*p2_1 , dy2_rk4(i) + 0.5*h_r*q2_1);
  q2_2 = diff2_y2_r(t_r(i) + 0.5*h_r , y1_rk4(i) + 0.5*h_r*p1_1 , dy1_rk4(i) + 0.5*h_r*q1_1 , y2_rk4(i) + 0.5*h_r*p2_1 , dy2_rk4(i) + 0.5*h_r*q2_1);
 
     % third order derivatives
  p1_3 = diff1_y1_r(t_r(i) + 0.5*h_r , y1_rk4(i) + 0.5*h_r*p1_2 , dy1_rk4(i) + 0.5*h_r*q1_2 , y2_rk4(i) + 0.5*h_r*p2_2 , dy2_rk4(i) + 0.5*h_r*q2_2);
  q1_3 = diff2_y1_r(t_r(i) + 0.5*h_r , y1_rk4(i) + 0.5*h_r*p1_2 , dy1_rk4(i) + 0.5*h_r*q1_2 , y2_rk4(i) + 0.5*h_r*p2_2 , dy2_rk4(i) + 0.5*h_r*q2_2);
  p2_3 = diff1_y2_r(t_r(i) + 0.5*h_r , y1_rk4(i) + 0.5*h_r*p1_2 , dy1_rk4(i) + 0.5*h_r*q1_2 , y2_rk4(i) + 0.5*h_r*p2_2 , dy2_rk4(i) + 0.5*h_r*q2_2);
  q2_3 = diff2_y2_r(t_r(i) + 0.5*h_r , y1_rk4(i) + 0.5*h_r*p1_2 , dy1_rk4(i) + 0.5*h_r*q1_2 , y2_rk4(i) + 0.5*h_r*p2_2 , dy2_rk4(i) + 0.5*h_r*q2_2);
 
     % fourth order derivatives
  p1_4 = diff1_y1_r(t_r(i) + h_r , y1_rk4(i) + h_r*p1_3 , dy1_rk4(i) + h_r*q1_3 , y2_rk4(i) + h_r*p2_3 , dy2_rk4(i) + h_r*q2_3);
  q1_4 = diff2_y1_r(t_r(i) + h_r , y1_rk4(i) + h_r*p1_3 , dy1_rk4(i) + h_r*q1_3 , y2_rk4(i) + h_r*p2_3 , dy2_rk4(i) + h_r*q2_3);
  p2_4 = diff1_y2_r(t_r(i) + h_r , y1_rk4(i) + h_r*p1_3 , dy1_rk4(i) + h_r*q1_3 , y2_rk4(i) + h_r*p2_3 , dy2_rk4(i) + h_r*q2_3);
  q2_4 = diff2_y2_r(t_r(i) + h_r , y1_rk4(i) + h_r*p1_3 , dy1_rk4(i) + h_r*q1_3 , y2_rk4(i) + h_r*p2_3 , dy2_rk4(i) + h_r*q2_3);
    
     % Append the approximated values at (i+1)
  y1_rk4(i+1)= y1_rk4(i) + h_r/6*(p1_1+2*p1_2+2*p1_3+p1_4);
  dy1_rk4(i+1)= dy1_rk4(i) + h_r/6*(q1_1+2*q1_2+2*q1_3+q1_4);
  y2_rk4(i+1)= y2_rk4(i) + h_r/6*(p2_1+2*p2_2+2*p2_3+p2_4);
  dy2_rk4(i+1)= dy2_rk4(i) + h_r/6*(q2_1+2*q2_2+2*q2_3+q2_4);

end

% Plot for RK4
figure("Name","RK4")
plot(y1_rk4,y2_rk4,'g');
grid on
legend('Orbit')
title('Orbit by Runge-Kutta4 Method')
xlabel('y1')
ylabel('y2')

time_rk4 = toc

% ODE45
% -------------------------------------------------------------------------
tic
% ode45(differential equations, time span, initial conditions)
% "odefun" is defined at the bottom
[t_ode45 , x_ode45] = ode45(@(t,x) odefun(t,x), [ti tf] ,...
                      [y1_init dy1_init y2_init dy2_init]);

% Plot for ode45
figure("Name", "ode45")
plot(x_ode45(:,1),x_ode45(:,3),'b') % first and third column stores the location values
grid on
legend('Orbit')
title('Orbit by ode45')
xlabel('y1')
ylabel('y2')

time_ode45 = toc

% TIME COMPARISON FOR ALL METHODS
% -------------------------------------------------------------------------
if (time_euler < time_ode45) 
    if (time_euler < time_rk4)
        disp('Euler solution is the fastest solution.')
    end
end

if (time_ode45 < time_rk4) 
    if( time_ode45 < time_euler)
        disp('ode45 solution is the fastest solution.')
    end
end
        
if (time_rk4 < time_rk4) 
    if(time_rk4 < time_ode45)
        disp('RK4 solution is the fastest solution.')
    end
end

% COMPARISON OF ALL CURVES
% -------------------------------------------------------------------------

figure ('Name','Comparing All Methods for Restricted 3 Body Problem')
plot(y1_euler, y2_euler, 'm') % euler
grid on
hold on
plot(y1_rk4, y2_rk4, 'g') % rk4
plot (x_ode45(:,1),x_ode45(:,3), 'b') % ode45
title("Comparison of All Methods")
legend('Euler Solution','RK4 Solution','ode45 Solution')


n=zeros(length(x_ode45));


%animation
figure("Name","Animation")
for k=1:3:length(t_ode45)
    
    clf % clear the figure
    
    t_k = t_ode45(k); % time
    x_k = x_ode45(k,1); % y1
    y_k = x_ode45(k,3); % y2
    z_k = 0; % since the motion is planar motion        

    % location of earth
    plot3(-u1,0,0,'b.','LineWidth',10,'MarkerSize',40)
    hold on
    %location of moon
    plot3(u2,0,0,'m.','LineWidth',5,'MarkerSize',25)
    hold on
    % location of satellite
    plot3(x_k,y_k,z_k,'-^r','LineWidth',2,'MarkerSize',10)
    % plot entire orbit
    hold on
    plot3(x_ode45(:,1),x_ode45(:,3),n,'k-','LineWidth',2)
    grid on
    xlabel('y1')
    ylabel('y2')
    zlabel('z')
    title(['spacecraft at t= ',num2str(t_k)])
    view([0 90]) % angle of view
    drawnow
    legend('earth','moon','spacecraft')
    
end



% Define "odefun"
function dx_dt=odefun(t,x)

    u1 = 0.012277471;
    u2 = 1-u1;
    dx_dt = zeros(4,1);
    dx_dt(1) = x(2);
    dx_dt(2) = x(1) + 2*x(4) - u2*(x(1)+u1)/((x(1)+u1)^2+x(3)^2)^(3/2)-u1*(x(1)-u2)/((x(1)-u2)^2+x(3)^2)^(3/2);
    dx_dt(3) = x(4);
    dx_dt(4) = x(3) - 2*x(2) - u2*x(3)/((x(1)+u1)^2+x(3)^2)^(3/2)-u1*x(3)/((x(1)-u2)^2+x(3)^2)^(3/2);

end


