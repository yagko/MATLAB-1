% ------------------------------------------
% ME 303 - LS 9
% Solution of an IVP by RK4
% Date : 14.12.2021
% ------------------------------------------

clc;clear all;close all;

% IVP
% y'(t) = exp(-2t)-2y
% with y(0) = 0.1
% in the interval [0,4]

% Define f(t,y)
f = @( t , y ) exp( -2 .* t ) - 2 .* y;

% Exact solution
y_exact = @( t ) 0.1 * exp( -2 .* t ) + ...
                 t .* exp( -2 .* t );
             
% Specify the interval
tmin = 0;
tmax = 4;

% Introduce the initial condition
y0 = 0.1;

% RK4 Method...

% Specify the step size
h1 = 0.2;
h2 = 0.1;

% Divide the interval into subintervals
% separated by the step size h
t1 = tmin : h1 : tmax;
t2 = tmin : h2 : tmax;

% Define the time vector to plot the exact
% solution
t_exact = tmin : 0.001 : tmax;

% Define the solution vector
y1_RK4 = zeros( 1 , length( t1 ) );
y2_RK4 = zeros( 1 , length( t2 ) );

% Assign the initial value
y1_RK4( 1 ) = y0;
y2_RK4( 1 ) = y0;

for k = 1 : length( t1 ) - 1
    
    f1 = f( t1( k ) , y1_RK4( k ) );
    f2 = f( t1( k ) + h1/2 , y1_RK4( k ) + ...
         h1/2 * f1 );
    f3 = f( t1( k ) + h1/2 , y1_RK4( k ) + ...
         h1/2 * f2 );
    f4 = f( t1( k ) + h1 , y1_RK4( k ) + ...
         h1 * f3 );
     
    % Calculate the approximated function
    % value at the next time step
    y1_RK4( k + 1 ) = y1_RK4( k ) + h1/6 *...
            ( f1 + 2 * f2 + 2 * f3 + f4 );
    
end

for k = 1 : length( t2 ) - 1
    
    f1 = f( t2( k ) , y2_RK4( k ) );
    f2 = f( t2( k ) + h2/2 , y2_RK4( k ) + ...
         h2/2 * f1 );
    f3 = f( t2( k ) + h2/2 , y2_RK4( k ) + ...
         h2/2 * f2 );
    f4 = f( t2( k ) + h2 , y2_RK4( k ) + ...
         h2 * f3 );
     
    % Calculate the approximated function
    % value at the next time step
    y2_RK4( k + 1 ) = y2_RK4( k ) + h2/6 *...
            ( f1 + 2 * f2 + 2 * f3 + f4 );
    
end

% Solution of the IVP via the built-in
% function ode45
[ t_ode45 , y_ode45 ] = ode45( f , ...
                     [ tmin tmax ] , y0 );

% Plot the solutions to the IVP
plot( t1 , y1_RK4 , '-.b' )
grid on
hold on
plot( t2 , y2_RK4 , '--r' )
plot( t_ode45 , y_ode45 , '-oc' )
plot( t_exact , y_exact( t_exact ) , '-k' )
legend( 'y_{RK4}(h=0.2)' , ...
        'y_{RK4}(h=0.1)' , 'y_{ode45}' ,...
        'y_{exact}' )