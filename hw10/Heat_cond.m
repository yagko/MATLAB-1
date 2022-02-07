clc; clear all; close all;

% given constants
h = 1;
k = 1;
p = 1;
c = 1;
q = 10;

% specify the position vector
xi = 0;
xf = 1;
dx = 0.1;
x = xi:dx:xf;

% specify the time step and time vector
dt = 0.001;
maxTime = 5;
t=0:dt:maxTime;

% calculate the number of nodes
N = length(x);

% define temperature array, T
T = zeros(N,1);

% specify the initial condition
T_air = 20;
T0 = 0;

% define the temperature distribution at the previous time step
T_prev = T;

% create the coefficient matrices, A and B
A = zeros(N,N);
B = zeros(N,1);

% calculate the Fourier Number
Fo = (k*dt)/dx^2;

% calculate the elements of A
for i = 1 : N
    if i == 1
        % Left Boundary
         A(i,i)=1-2*Fo;
        A(i,i+1)=2*Fo;

        B(i)=Fo*q*2*dx/k;

    elseif i == N
        % Right Boundary
        A(i,i)=(1 - Fo*(h*2*dx/k + 2));
        A(i,i-1)=2*Fo;

        B(i)=T_air*h*2*dx*Fo/k;

    else
        % Interior Nodes
        A(i,i-1) = Fo;
        A(i,i+1) = Fo;
        A(i,i) = 1 - 2*Fo;

    end
end

% plot the temperature distributions
figure
hold on
xlabel('x (m)')
ylabel('Temperature (C)')


tol = 1E-6;

for time_step = 0: dt : maxTime
    T = A*T_prev + B;
    if mod(time_step , 0.1) == 0
        plot(x,T)
    end
    if max(abs(T-T_prev)) < tol
        break
    end
    T_prev = T;
end