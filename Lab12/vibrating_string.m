clear all
close all
clc

% String parameters
L = 1;              % String length [m]
M = 0.1;           % String mass [kg]
pA = M/L;          % Linear mass density [kg/m]

T = 500;           % Tension [N]
c = sqrt(T/pA);     % Wave speed [m/s]

% Initial Conds
A = 0.05;               % Maximum displacement amplitude [m]
f = @(x) A*sin(pi*x/L); % Initial displacement of the string [0,L]
g = @(x) 0;             % Initial velocity of the string [0,L]

n = 11;
dx = L/(n-1);
x = 0:dx:L;

tf = 0.1;
dt = 0.0005;
t = 0:dt:tf;
m = tf/dt + 1;

r = c*dt/dx;

u = zeros(n,m);

for i=2:n-1
    u_tt = (c^2)*(f(x(i+1))-2*f(x(i))+f(x(i-1)))/dx^2;
    u_t = g(x(i));
    u(i,1) = f(x(i));
    u(i,2) = u(i,1) + u_t*dt + u_tt*dt^2/2;
end

for j=3:m
    for i=2:n-1
        u(i,j) = (2-2*r^2)*u(i,j-1) + r^2*(u(i+1,j-1)+...
            u(i-1,j-1))-u(i,j-2);
    end
end

figure
axis([0 L -2*A 2*A])
xlabel('x(m)')
ylabel('u(m)')
grid on
ax = gca;

% Set the font size, tick direction, tick length, and y-axis limits 
% for the current axes. Use gca to refer to the current axes.

ax.NextPlot = 'replaceChildren';
replacechildren

% Deletes the figure's children (axes objects and their descendants) and 
% makes this figure the current figure.

for k=1:m
    plot(x,u(:,k),'LineWidth',1.5)
    F(k) = getframe;
end


% u_tt = c*u_xx
% u_xx = (f(x+h)-2f(x)+f(x-h))/dx^2