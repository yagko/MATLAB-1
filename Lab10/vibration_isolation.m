% ------------------------------------------
% Me 303 - Lab Session 10
% Vibration Isolation Problem
% Date : 21.12.2021
% ------------------------------------------

clc;clear all;close all;

% Specify the system parameters
m = 75;
c = 20;
k1 = 200;
k2 = 1000;
A = 25;

% Calculate the equivalent stiffnesses for
% designs 1 and 2
k_eq1 = k1 * k2 / ( k1 + k2 );
k_eq2 = k1 + k2;

% Calculate frequency of oscillations
% of the base,
w_n1 = sqrt( k_eq1 / m );
r = 2;
w = r * w_n1;

% Define the base excitation function
y = @( t ) A .* sin( w .* t );

% Define the slope functions
g_11 = @( t , z_11 , z_21 ) z_21;
g_21 = @( t , z_11 , z_21 ) -c/m * ...
    z_21 - k_eq1/m * z_11 ...
    + c/m * w*A*cos( w*t ) + ...
    k_eq1/m * A*sin( w*t );

g_12 = @( t , z_12 , z_22 ) z_22;
g_22 = @( t , z_12 , z_22 ) -c/m * ...
    z_22 - k_eq2/m * z_12 ...
    + c/m * w*A*cos( w*t ) + ...
    k_eq2/m * A*sin( w*t );

% Specify the initial points
z11_0 = 0.5;
z21_0 = -20;
z12_0 = 0.5;
z22_0 = -20;

% Specify the time interval
tmin = 0;
tmax = 50;

% Specify the step size
h = 0.2;

% Create the time array
t = tmin : h : tmax;

% Apply RK4 method...
y1 = RK4_fun( { g_11 ; g_21 } , ...
    [ z11_0 ; z21_0 ] , ...
    [ tmin ; tmax ] , h );
y2 = RK4_fun( { g_12 ; g_22 } , ...
    [ z12_0 ; z22_0 ] , ...
    [ tmin ; tmax ] , h );

x1 = y1( : , 1 );
x2 = y2( : , 1 );

% Plot the results
figure ( 1 )
plot( t , x1 , '-.b' )
grid on
hold on
plot( t , x2 , '--g' )
plot( t , y( t ) , '-r' )
legend( '1^{st} Design' ,...
    '2^{nd} Design' , 'Base' )
xlabel( 't (s)' )
ylabel( 'x(t) (mm)' )

% Calculate the displacement 
% transmissibilities
TR1 = max( abs( x1 ) ) / A
TR2 = max( abs( x2 ) ) / A

tspan = [ tmin tmax ];
z0_matlab1 = [ z11_0 , z21_0 ];
z0_matlab2 = [ z12_0 , z22_0 ];
[ t1_ode45 , x1_ode45 ] = ode45( ...
    @( t , x ) odefun( t , x , m , c ,...
    k_eq1 , A , w ) , tspan , ...
    z0_matlab1 );
[ t2_ode45 , x2_ode45 ] = ode45( ...
    @( t , x ) odefun( t , x , m , c ,...
    k_eq2 , A , w ) , tspan , ...
    z0_matlab2 );

x1_matlab = x1_ode45( : , 1 );

% Plot the results
figure ( 2 )
plot( t , x1 , '--b' , 'LineWidth' , 2 )
grid on
hold on
plot( t1_ode45 , x1_matlab , '-r' )
legend( 'RK4_fun' , 'ode45' )
xlabel( 't (s)' )
ylabel( 'x(t) (mm)' )

function dx_dt = odefun( t , x , m , c ,...
    k , A , w )

dx_dt = zeros( 2 , 1 );
dx_dt( 1 ) = x( 2 );
dx_dt( 2 ) = -c/m * ...
    x( 2 ) - k/m * x( 1 ) ...
    + c/m * w*A*cos( w*t ) + ...
    k/m * A*sin( w*t );

end