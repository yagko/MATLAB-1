clear all
close all
clc
%Süleyman Enes ŞEKER
%arenstorf orbit
% 3 method(euler, rungekutta 4 and ode45) are used to solve problem
% firstly ı change 2, 2nd degree differential equation into 4, 1st degree
% differential equation calculations given on the pdf in the project zip.

%euler
% initially ı defined my constant parameters
tic
ti=0;
tf=17.06521656015796;
h=(tf-ti)/24000;
m=ti:h:tf ;
u1=0.012277471;
u2 = 1-u1;
% ı create matrices to send solution values of each differential equation
y1=zeros(1,length(m));
x1=zeros(1,length(m));
y2=zeros(1,length(m));
x2=zeros(1,length(m));
D1=zeros(1,length(m));
D2=zeros(1,length(m));
% ı added initial values to 1st element of the matrices 
y1(1)=0.994;
x1(1)=0;
y2(1)=0;
x2(1)=-2.0015851063790825;

% ı define a for loop in order to solve ODEs by using euler formula as we
% did in the lab and homeworks.General formula y(i+1)=y(i)+h*y'(i)
for i=1:length(m)
    
    D1(i)=((y1(i)+u1)^2+y2(i)^2)^(3/2);
    D2(i)=((y1(i)-u2)^2+y2(i)^2)^(3/2);
    y1(i+1)=y1(i)+h*x1(i);
    x1(i+1)=x1(i)+h*(y1(i)+2*x2(i)-u2*(y1(i)+u1)/D1(i)-u1*(y1(i)-u2)/D2(i));
    y2(i+1)=y2(i)+h*x2(i);
    x2(i+1)=x2(i)+ h*(y2(i)-2*x1(i)-u2*y2(i)/D1(i)-u1*y2(i)/D2(i));
    
    
end
% then ı plot y1,y2 together and ı calculate the time by tic-toc function
% to compare methods
plot(y1,y2);
grid on
legend('orbit')
title('Euler Method')
xlabel('y1')
ylabel('y2')
figure
timeeuler=toc


% runge kutta4
%as ı did in the euler method ı define constant values firstly
tic
ti=0;
tf=17.06521656015796;
h=(tf-ti)/6000;
m=ti:h:tf ;
u1=0.012277471;
u2 = 1-u1;
% then ı created matrices to assign solution values on
y1=zeros(1,length(m));
x1=zeros(1,length(m));
y2=zeros(1,length(m));
x2=zeros(1,length(m));
D1=zeros(1,length(m));
D2=zeros(1,length(m));
% ı added initial values to 1st element of the matrices 
y1(1)=0.994;
x1(1)=0;
y2(1)=0;
x2(1)=-2.0015851063790825;
% ı write my differential equations in function form. By doing so applying 
% runge kutta 4 method in the for loop is easier.
fy1=@(x1) x1;
fx1=@(y1,y2,x2) y1+2*x2-u2*(y1+u1)/((y1+u1)^2+y2^2)^(3/2)-u1*(y1-u2)/((y1-u2)^2+y2^2)^(3/2);
fy2=@(x2) x2;
fx2=@(x1,y1,y2) y2-2*x1-u2*y2/((y1+u1)^2+y2^2)^(3/2)-u1*y2/((y1-u2)^2+y2^2)^(3/2);
% ı define a for loop in order to solve ODEs by using runge kutta4 formula 
%as we did in the lab and homeworks.
%General formula:
%f1=y'(t(i),y(i)),f2=(t(i)+h/2,y(i)+h/2*f1),f3=(t(i)+h/2,y(i)+h/2*f2),
%f4=(t(i),y(i)+h*f3); y(i+1)=y(i)+(h/6)*(f1+2*f2+2*f3+f4)
for i=1:length(m)
    
    
    f1y1=fy1(x1(i));
    f1x1=fx1(y1(i),y2(i),x2(i));
    f1y2=fy2(x2(i));
    f1x2=fx2(x1(i),y1(i),y2(i));
    
    f2y1=fy1(x1(i)+h*f1x1/2);
    f2x1=fx1(y1(i)+h*f1y1/2,y2(i)+h*f1y2/2,x2(i)+h*f1x2/2);
    f2y2=fy2(x2(i)+h*f1x2/2);
    f2x2=fx2(x1(i)+h*f1x1/2,y1(i)+h*f1y1/2,y2(i)+h*f1y2/2);
    
    f3y1=fy1(x1(i)+h*f2x1/2);
    f3x1=fx1(y1(i)+h*f2y1/2,y2(i)+h*f2y2/2,x2(i)+h*f2x2/2);
    f3y2=fy2(x2(i)+h*f2x2/2);
    f3x2=fx2(x1(i)+h*f2x1/2,y1(i)+h*f2y1/2,y2(i)+h*f2y2/2);
    
    f4y1=fy1(x1(i)+h*f3x1);
    f4x1=fx1(y1(i)+h*f3y1,y2(i)+h*f3y2,x2(i)+h*f3x2);
    f4y2=fy2(x2(i)+h*f3x2);
    f4x2=fx2(x1(i)+h*f3x1,y1(i)+h*f3y1,y2(i)+h*f3y2);
    
    y1(i+1)=y1(i)+(h*(f1y1+2*f2y1+2*f3y1+f4y1)/6);
    x1(i+1)=x1(i)+(h*(f1x1+2*f2x1+2*f3x1+f4x1)/6);
    y2(i+1)=y2(i)+(h*(f1y2+2*f2y2+2*f3y2+f4y2)/6);
    x2(i+1)=x2(i)+(h*(f1x2+2*f2x2+2*f3x2+f4x2)/6);
    
end
% then ı plot y1,y2 together and ı calculate the time by tic-toc function
% to compare methods
plot(y1,y2);
grid on
legend('orbit')
title('Runge-Kutta4 Method')
xlabel('y1')
ylabel('y2')
figure

timerk4=toc

% ode45
tic
% in ode45 method first ı define may time and initial values as matrices
% since we have 4 initial value it is 1 by 4 matrices
tSpan = linspace(0,17.06521656015796);
x0 = [0.994 0 0 -2.0015851063790825];
% secondly ı call the ode 45 function to solve my differential equations
% which ı define as odefun(t,x) at the bottom.
[t1 ,x1] = ode45(@(t,x) odefun(t,x), tSpan, x0);
%ince 2nd and 4th thh column of the x1 matrix hold the solutions of x1 and x1
% only 1st and 3rd column of the x1 matrix is plotted in order to obtain y1
% and y2.
plot(x1(:,1),x1(:,3))
grid on
legend('orbit')
title('ode45 solution')
xlabel('y1')
ylabel('y2')
figure
n=zeros(length(x1));
timeode45=toc

%animation
for k=1:length(t1)
   clf
    
    
    t_k=t1(k); % time
    x_k=x1(k,1); %y1
    y_k=x1(k,3); %y2
    z_k=0;        
    % ın order to  show in 3d ı took z as 0 since it is a 2d problme
    % location of earth
    plot3(-0.012277471,0,0,'b.','LineWidth',10,'MarkerSize',40)
    hold on
    %location of moon
    plot3(1-0.012277471,0,0,'m.','LineWidth',5,'MarkerSize',25)
    hold on
    % location of satellite
    plot3(x_k,y_k,z_k,'ro','LineWidth',2,'MarkerSize',10)
    % plot entire orbit
    hold on
    plot3(x1(:,1),x1(:,3),n,'k-','LineWidth',2)
    grid on
    xlabel('y1')
    ylabel('y2')
    zlabel('z')
    title(['spacecraft at t= ',num2str(t_k)])
    view([-30 35]) %view angle
    drawnow
    legend('earth','moon','spacecraft')
    
end


% here ı define odefun(t,x) to send this column matrix to ode45 function
% in the ode45, functions ,that will be solved, must be defined as column
% matrix. since we have 4, 1st order differential equation i create 4 by 1
% matrix.Every row contain one differential equation
function dx_dt=odefun(t,x);

u1=0.012277471;
u2 = 1-u1;
dx_dt=zeros(4,1);
dx_dt(1)=x(2);
dx_dt(2)=x(1)+2*x(4)-u2*(x(1)+u1)/((x(1)+u1)^2+x(3)^2)^(3/2)-u1*(x(1)-u2)/((x(1)-u2)^2+x(3)^2)^(3/2);
dx_dt(3)=x(4);
dx_dt(4)=x(3)-2*x(2)-u2*x(3)/((x(1)+u1)^2+x(3)^2)^(3/2)-u1*x(3)/((x(1)-u2)^2+x(3)^2)^(3/2);

end


