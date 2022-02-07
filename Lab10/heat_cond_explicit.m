% ------------------------------------------
% Me 303 - Lab Session 10
% 1D Transient Heat Conduction Problem
% Date : 21.12.2021
% ------------------------------------------

clc;clear all;close all;

% Given properties
% k/(dens*cp) = 1
alpha = 1;
L = 1;

% Specify the grid for the position, x
xmin = 0;
xmax = L;
dx = ( xmax - xmin ) / 40;
x = xmin : dx : xmax;

% Specify the time step, dt
dt = 0.001;

% Calculate the number of nodes
N = length( x );

% Define the temperature array, T
T = zeros( N , 1 );

% Specify the initial condition
Ti = 20;
T_left = 40;

% Apply the initial condition
T( 1 ) = T_left;
T( 2 : N , 1 ) = Ti;

% Define the temperature distribution
% at the previous time step
T_prev = T;

% Create the coefficient matrix, A
A = zeros( N , N );

% Calculate the Fourier number
Fo = alpha * dt / dx^2;

% Calculate the elements of A
for i = 1 : N
    
    if i == 1
        % Left Boundary
        A( i , i ) = 1;
    elseif i == N
        % Right Boundary
        A( i , i - 1 ) = 2 * Fo;
        A( i , i ) = 1 - 2 * Fo;
    else
        % Interior Nodes
        A( i , i - 1 ) = Fo;
        A( i , i + 1 ) = Fo;
        A( i , i ) = 1 - 2 * Fo;
    end
    
end

% Plot the temperature distributions
figure
hold on
xlabel('x (m)')
ylabel('Temperature (C)')

maxTime = 10;
timeInc = 0.1;
maxTstep = maxTime * ( 1 / dt );
timeIncStep = timeInc * ( 1 / dt );
tol = 1E-6;

for time_step = 1 : maxTstep
    T = A * T_prev;
    if mod( time_step , timeIncStep ) == 0
        plot( x , T )
    end
    if max( abs( T - T_prev ) ) < tol
        
        tss = time_step * dt;
        break
        
    end
    T_prev = T;
end

fprintf('System reaches steady state at %8f s.\n' , tss )
fprintf('Steady state temperature at the left boundary is %8f degrees Celcius.\n' , T( 1 ) );
fprintf('Steady state temperature at the right boundary is %8f degrees Celcius.\n' , T( N ) );