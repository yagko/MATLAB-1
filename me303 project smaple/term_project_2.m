%TERM PROJECT QUESTION 2
clc
clear all
close all

% Euler's Method

% Initial conditions

moon_mass=0.012277471; %mass of the moon (given)
earth_mass=1-moon_mass; %mass of the earth (given)
initial_time=0 %starting time (given)
final_time=17.06521656015796; %ending time (given)
h_euler= (final_time-initial_time)/24000;  % step size (given)
N_euler = initial_time:h_euler:final_time;  % Time division 
M_euler=(final_time-initial_time)/h_euler % Number of intervals

D1_euler=@(N_euler, y1, y2) ((y1+moon_mass)^2+y2^2)^(3/2) % D1 formula(given)
D2_euler=@(N_euler, y1, y2) ((y1-earth_mass)^2+y2^2)^(3/2) % D2 formula(given)
p2_euler=@(N_euler, y1, dy1, y2, dy2) dy2; %the differential of y2
q2_euler=@(N_euler, y1, dy1, y2, dy2) -2*dy1 + (1-earth_mass/D1_euler(N_euler, y1, y2) -moon_mass/D2_euler(N_euler, y1,y2))*y2; % y''2(t), the second derivative formula of y2
p1_euler=@(N_euler, y1, dy1, y2, dy2) dy1; %the differential of y1
q1_euler=@(N_euler, y1, dy1, y2, dy2) (1-earth_mass/D1_euler(N_euler, y1, y2)-moon_mass/D2_euler(N_euler, y1, y2))*y1 +2*dy2-moon_mass*earth_mass*(1/D1_euler(N_euler, y1, y2) - 1/D2_euler(N_euler, y1, y2)); % y''1(t), the second derivative formula of y1


%Defining the size of variables. There will be M_euler + 1 points for S intervals.
y1_euler=zeros(1,M_euler+1);
dy1_euler=zeros(1,M_euler+1);
y2_euler=zeros(1,M_euler+1);
dy2_euler=zeros(1,M_euler+1);

%Initial conditions, given
y1_euler(1)=0.994; %given and here is written y1(1) instead of y1(0) although it was given as y1(0)=0.994 in the question. The reason is that here y1 is a matrix and the number inside the parantheses is the index of the matrix, so the initial element of the matrix here has the index of 1, not 0. In the question, it was meant y1 value at the time t=0 by the expression y1(0)=0.994
dy1_euler=0; %given, which is p1 stated above 
y2_euler=0; %given %the reason why y2(1) is written instead of y2(0) is stated above
dy2_euler= -2.0015851063790825;  %given, which is p2 stated above 

tic; %in order to measure the operating time for this Euler's loop, tic & toc is used. The measurement starts with tic and ends with toc
% The for loop to solve the differential equation by using the Euler's method:
for i=1:M_euler
    y1_euler(i+1)=y1_euler(i)+h_euler*p1_euler(N_euler(i), y1_euler(i), dy1_euler(i), y2_euler(i), dy2_euler(i));
    dy1_euler(i+1)=dy1_euler(i)+h_euler*q1_euler(N_euler(i), y1_euler(i), dy1_euler(i), y2_euler(i), dy2_euler(i));
    y2_euler(i+1)=y2_euler(i)+h_euler*p2_euler(N_euler(i), y1_euler(i), dy1_euler(i), y2_euler(i), dy2_euler(i));
    dy2_euler(i+1)=dy2_euler(i)+h_euler*q2_euler(N_euler(i), y1_euler(i), dy1_euler(i), y2_euler(i), dy2_euler(i));
end
toc; %The measurement starts with tic and ends with toc, so the toc value gives the time required to run this loop for Euler's method
euler_time=toc; %the operating time data is stored in this variable

figure(1)
plot(y1_euler, y2_euler)
grid on
title("Orbit by Euler Method")
xlabel("y1(t)")
ylabel("y2(t)")

%Runge Kutta Method

D1_rk4= @(N_rk4, y1, y2) ((y1+moon_mass)^2 + y2^2)^(3/2); %given
D2_rk4= @(N_rk4, y1, y2) ((y1-earth_mass)^2 + y2^2)^(3/2); %given
p1_rk4= @(N_rk4, y1, dy1, y2, dy2) dy1;
q1_rk4= @(N_rk4, y1, dy1, y2, dy2) (1-earth_mass/D1_rk4(N_rk4, y1, y2)-moon_mass/D2_rk4(N_rk4, y1, y2))*y1 +2*dy2-moon_mass*earth_mass*(1/D1_rk4(N_rk4, y1, y2) - 1/D2_rk4(N_rk4, y1, y2)); % y''1(t)
p2_rk4= @(N_rk4, y1, dy1, y2, dy2) dy2;
q2_rk4= @(N_rk4, y1, dy1, y2, dy2) -2*dy1 + (1-earth_mass/D1_rk4(N_rk4, y1, y2) -moon_mass/D2_rk4(N_rk4, y1,y2))*y2; % y''2(t)

h_rk4= (final_time-initial_time)/6000;% step size
N_rk4 = initial_time:h_rk4:final_time;% Time range;
M_rk4=(final_time-initial_time)/h_rk4; % number of intervals
%Defining the size of variables. There will be R+1 points for R intervals.
y1_rk4=zeros(1,M_rk4+1);
dy1_rk4=zeros(1,M_rk4+1);
y2_rk4=zeros(1,M_rk4+1);
dy2_rk4=zeros(1,M_rk4+1);

%Initial conditions, given
y1_rk4(1)=0.994;
dy1_rk4=0;
y2_rk4=0;
dy2_rk4= -2.0015851063790825;

tStart = tic;  %to start measuring the operating time of the RK4 loop 
for i=1:M_rk4  % calculation loop
    %the starting step of RK4 method, so the 1st order formulations/equations of RK4 method
  p1_1= p1_rk4(N_rk4(i), y1_rk4(i), dy1_rk4(i), y2_rk4(i), dy2_rk4(i)); %defining p1 in terms of partition vector, y1, dy1, y2, and dy2
  q1_1= q1_rk4(N_rk4(i), y1_rk4(i), dy1_rk4(i), y2_rk4(i), dy2_rk4(i)); %defining q1 in terms of partition vector, y1, dy1, y2, and dy2
  p2_1= p2_rk4(N_rk4(i), y1_rk4(i), dy1_rk4(i), y2_rk4(i), dy2_rk4(i)); %defining p2 in terms of partition vector, y1, dy1, y2, and dy2
  q2_1= q2_rk4(N_rk4(i), y1_rk4(i), dy1_rk4(i), y2_rk4(i), dy2_rk4(i)); %defining q2 in terms of partition vector, y1, dy1, y2, and dy2
 
     %2nd order formulations/equations of RK4 method
  p1_2=p1_rk4(N_rk4(i)+0.5*h_rk4,y1_rk4(i)+0.5*h_rk4*p1_1,dy1_rk4(i)+0.5*h_rk4*q1_1,y2_rk4(i)+0.5*h_rk4*p2_1,dy2_rk4(i)+0.5*h_rk4*q2_1);
  q1_2=q1_rk4(N_rk4(i)+0.5*h_rk4,y1_rk4(i)+0.5*h_rk4*p1_1,dy1_rk4(i)+0.5*h_rk4*q1_1,y2_rk4(i)+0.5*h_rk4*p2_1,dy2_rk4(i)+0.5*h_rk4*q2_1);
  p2_2=p2_rk4(N_rk4(i)+0.5*h_rk4,y1_rk4(i)+0.5*h_rk4*p1_1,dy1_rk4(i)+0.5*h_rk4*q1_1,y2_rk4(i)+0.5*h_rk4*p2_1,dy2_rk4(i)+0.5*h_rk4*q2_1);
  q2_2=q2_rk4(N_rk4(i)+0.5*h_rk4,y1_rk4(i)+0.5*h_rk4*p1_1,dy1_rk4(i)+0.5*h_rk4*q1_1,y2_rk4(i)+0.5*h_rk4*p2_1,dy2_rk4(i)+0.5*h_rk4*q2_1);
 
      %3rd order formulations/equations of RK4 method
  p1_3=p1_rk4(N_rk4(i)+0.5*h_rk4,y1_rk4(i)+0.5*h_rk4*p1_2,dy1_rk4(i)+0.5*h_rk4*q1_2,y2_rk4(i)+0.5*h_rk4*p2_2,dy2_rk4(i)+0.5*h_rk4*q2_2);
  q1_3=q1_rk4(N_rk4(i)+0.5*h_rk4,y1_rk4(i)+0.5*h_rk4*p1_2,dy1_rk4(i)+0.5*h_rk4*q1_2,y2_rk4(i)+0.5*h_rk4*p2_2,dy2_rk4(i)+0.5*h_rk4*q2_2);
  p2_3=p2_rk4(N_rk4(i)+0.5*h_rk4,y1_rk4(i)+0.5*h_rk4*p1_2,dy1_rk4(i)+0.5*h_rk4*q1_2,y2_rk4(i)+0.5*h_rk4*p2_2,dy2_rk4(i)+0.5*h_rk4*q2_2);
  q2_3=q2_rk4(N_rk4(i)+0.5*h_rk4,y1_rk4(i)+0.5*h_rk4*p1_2,dy1_rk4(i)+0.5*h_rk4*q1_2,y2_rk4(i)+0.5*h_rk4*p2_2,dy2_rk4(i)+0.5*h_rk4*q2_2);

      %4th order formulations/equations of RK4 method
  p1_4=p1_rk4(N_rk4(i)+h_rk4,y1_rk4(i)+h_rk4*p1_3,dy1_rk4(i)+h_rk4*q1_3,y2_rk4(i)+h_rk4*p2_3,dy2_rk4(i)+h_rk4*q2_3);
  q1_4=q1_rk4(N_rk4(i)+h_rk4,y1_rk4(i)+h_rk4*p1_3,dy1_rk4(i)+h_rk4*q1_3,y2_rk4(i)+h_rk4*p2_3,dy2_rk4(i)+h_rk4*q2_3);
  p2_4=p2_rk4(N_rk4(i)+h_rk4,y1_rk4(i)+h_rk4*p1_3,dy1_rk4(i)+h_rk4*q1_3,y2_rk4(i)+h_rk4*p2_3,dy2_rk4(i)+h_rk4*q2_3);
  q2_4=q2_rk4(N_rk4(i)+h_rk4,y1_rk4(i)+h_rk4*p1_3,dy1_rk4(i)+h_rk4*q1_3,y2_rk4(i)+h_rk4*p2_3,dy2_rk4(i)+h_rk4*q2_3);

  y1_rk4(i+1)= y1_rk4(i) + h_rk4*(1/6)*(p1_1+2*p1_2+2*p1_3+p1_4); % main equation of RK4 for y1
  dy1_rk4(i+1)= dy1_rk4(i) + h_rk4*(1/6)*(q1_1+2*q1_2+2*q1_3+q1_4);% main equation of RK4 for dy1, namely p1
  y2_rk4(i+1)= y2_rk4(i) + h_rk4*(1/6)*(p2_1+2*p2_2+2*p2_3+p2_4); % main equation of RK4 for y2
  dy2_rk4(i+1)= dy2_rk4(i) + h_rk4*(1/6)*(q2_1+2*q2_2+2*q2_3+q2_4); % main equation of RK4 for dy2, namely p2

end
tEnd = toc(tStart); %here the tic toc pair is used with another syntax for the seek of variety
rk4_time=tEnd;  %the toc value, so the time required for RK4 method, is stored in the variable rk4_time

figure(2)
plot(y1_rk4, y2_rk4)
grid on
title("Orbit by Runge Kutta 4th Order Method")
xlabel("y1(t)")
ylabel("y2(t)")


%ODE45 SOLUTION 
tic; 
ode45method=@(t,x)[x(3);x(4);
    x(1) + 2*x(4) - earth_mass*(x(1)+moon_mass)/((((x(1)+moon_mass)^2)+(x(2)^2))^1.5) - moon_mass*(x(1)-earth_mass)/((((x(1)-earth_mass)^2)+(x(2)^2))^1.5);
    x(2) - 2*x(3) - earth_mass*x(2)/((((x(1)+moon_mass)^2)+(x(2)^2))^1.5) - moon_mass*x(2)/((((x(1)-earth_mass)^2)+(x(2)^2))^1.5)
     ]
 x0=[0.994 0 0 -2.0015851063790825]; %the initial condition matrix
[t,x]=ode45(ode45method,[initial_time final_time],x0); %the timespan and the x operated with the ode45 function
toc;
ode45_time=toc; %measured the time required
display(ode45_time); %written into command window


%plotting ode45 results:
figure('Name','ode45 method')
plot(x(:,1),x(:,2), 'b')
grid on

%comparing the 3 methods in terms of operation speed to select the fastest
%method by comparing the operating time datas taken from each method
if(euler_time<ode45_time) 
    if(euler_time<rk4_time)
    display('Euler solution is faster')
    end
end
 if (ode45_time<rk4_time) 
        if(ode45_time<euler_time)
         display('ode45 solution is faster')
        end
    end
        
if (rk4_time<euler_time) 
            if(rk4_time<ode45_time)
            display('rk4 solution is faster')
            end
end

% comparison of all curves
figure ('Name','Comparing Methods for Restricted 3 Body Problem')
plot(y1_euler, y2_euler, 'm') %euler method solution plotting
hold on
plot(y1_rk4, y2_rk4, 'g') %k1=y1, k2=dy1, k3=y2, k4=dy2 as stated in hand calculation attached to the report %in this line the rk4 solution is plotted
hold on
plot (x(:,1),x(:,2), 'b')% ode45 method solution plotting
hold on
grid on
title("Comparison of All Methods")
legend('Euler Solution','RK4 Solution','ode45 Solution')
