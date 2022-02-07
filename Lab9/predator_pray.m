% ------------------------------------------
% ME 303 - LS 9
% Predator-prey problem by RK4
% Date : 14.12.2021
% ------------------------------------------

clc;clear all;close all;

% Specify the coefficients characterizing
% the predator-prey model
A = 2;
B = 0.02;
C = 0.0002;
D = 0.8;

% Define the differential equations
f_rabbit = @( t , rabbit , fox ) A * rabbit...
    - B * rabbit * fox;
f_fox = @( t , rabbit , fox ) C * rabbit ...
    * fox - D * fox;

% Specify the time interval
tmin = 0;
tmax = 5;

% Introduce the initial conditions
rabbit_init = 3000;
fox_init = 120;

% Specify the step size
M = 50;
h = ( tmax - tmin ) / M;

% Create the time array
t = tmin : h : tmax;

% Define the solution vectors
rabbit = zeros( 1 , length( t ) );
fox = zeros( 1 , length( t ) );

% Call our function
y = RK4_fun( { f_rabbit ; f_fox } , [ rabbit_init ; fox_init ] , [ tmin ; tmax ] , h );
rabbit = y( : , 1 );
fox = y( : , 2 );
     
% Assign the initial values
%rabbit( 1 ) = rabbit_init;
%fox( 1 ) = fox_init;

%for k = 1 : M
%     
%     f1_rabbit = f_rabbit( t( k ) , ...
%                  rabbit( k ) , fox( k ) );
%     f1_fox = f_fox( t( k ) , ...
%                  rabbit( k ) , fox( k ) );
%     f2_rabbit = f_rabbit( t( k ) + h/2 ,...
%         rabbit( k ) + h/2 * f1_rabbit ,...
%         fox( k ) + h/2 * f1_fox );
%     f2_fox = f_fox( t( k ) + h/2 ,...
%         rabbit( k ) + h/2 * f1_rabbit ,...
%         fox( k ) + h/2 * f1_fox );
%     f3_rabbit = f_rabbit( t( k ) + h/2 ,...
%         rabbit( k ) + h/2 * f2_rabbit ,...
%         fox( k ) + h/2 * f2_fox );
%     f3_fox = f_fox( t( k ) + h/2 ,...
%         rabbit( k ) + h/2 * f2_rabbit ,...
%         fox( k ) + h/2 * f2_fox );
%     f4_rabbit = f_rabbit( t( k ) + h ,...
%         rabbit( k ) + h * f3_rabbit , ...
%         fox( k ) + h * f3_fox );
%     f4_fox = f_fox( t( k ) + h ,...
%         rabbit( k ) + h * f3_rabbit , ...
%         fox( k ) + h * f3_fox );
%      
%     % Calculate the approximated function
%     % value at the next time step
%     rabbit( k + 1 ) = rabbit( k ) + h/6 *...
%             ( f1_rabbit + 2 * f2_rabbit + ...
%             2 * f3_rabbit + f4_rabbit );
%     fox( k + 1 ) = fox( k ) + h/6 *...
%             ( f1_fox + 2 * f2_fox + ...
%             2 * f3_fox + f4_fox );
%     
% end

[ t_ode45 , x_ode45 ] = ode45( @( t , x )...
      odefun( t , x ) ,...
      [ tmin tmax ] , [ rabbit_init , ...
      fox_init ] );
  
rabbit_ode45 = x_ode45( : , 1 );
fox_ode45 = x_ode45( : , 2 );

% Plot the results
plot( t , rabbit , 'LineWidth' , 2 )
grid on
hold on
plot( t , fox , 'LineWidth' , 2 )
plot( t_ode45 , rabbit_ode45 , '-oc' )
plot( t_ode45 , fox_ode45 , '-^r' )
legend( 'Rabbit Population' , ...
        'Fox Population' ,...
        'Rabbit Population (ode45)',...
        'Fox Population (ode45)' )
    
% Print the population after 5 weeks
fprintf( 'Rabbit population after 5 weeks : %d\n' , round( rabbit( end ) ) );
fprintf( 'Fox population after 5 weeks : %d\n' , round( fox( end ) ) );

function dx_dt = odefun( t , x )

% Specify the coefficients characterizing
% the predator-prey model
A = 2;
B = 0.02;
C = 0.0002;
D = 0.8;

dx_dt = zeros( 2 , 1 );
dx_dt( 1 ) = A * x( 1 )...
             - B * x( 1 ) * x( 2 );
dx_dt( 2 ) = C * x( 1 ) ...
             * x( 2 ) - D * x( 2 );

end