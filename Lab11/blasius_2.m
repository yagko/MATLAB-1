% ------------------------------------------
% Me 303 - Lab Session 11
% Boundary Layer on a Flat Plate
% (Blasius Solution)
% Alternative Solution (via fzero)
% Date : 28.12.2021
% ------------------------------------------

clc;clear all;close all;

% Define the slope functions
g1 = @( eta , z1 , z2 , z3 ) z2;
g2 = @( eta , z1 , z2 , z3 ) z3;
g3 = @( eta , z1 , z2 , z3 ) - 1/2 * z1 * z3;

% Specify the initial conditions
z1_0 = 0;
z2_0 = 0;

% Specify the initial guess for the z3_0
z3_0 = 0.1;

% Specify the interval
eta_min = 0;
eta_max = 10;

% Specify the step size, h
h = 0.1;

% Calculate the eta array
eta = ( eta_min : h : eta_max )';

% Define the function that gives the correct
% initial guess for z3 when equated to zero
find_fdp_init = @( x ) RK4_blasius( ...
              { g1 ; g2 ; g3 } , ...
              [ z1_0 ; z2_0 ; x ] , ...
              [ eta_min ; eta_max ] , h ) - 1;
          
% Find the initial value of z3 that makes the
% above function 0
fdp_init = fzero( find_fdp_init , z3_0 );

% Calculate f, f', and f'' using the correct
% initial guess
[ fpInf , y ] = RK4_blasius( ...
              { g1 ; g2 ; g3 } , ...
              [ z1_0 ; z2_0 ; fdp_init ] , ...
              [ eta_min ; eta_max ] , h );

% Derive corresponding f, f', and f'' values
f = y( : , 1 );
f_p = y( : , 2 );
f_dp = y( : , 3 );

% Print the results in the form of a table
T = table( eta , f , f_p , f_dp )

% Calculate the index that corresponds to
% u/U=0.99 (f'=0.99)
[ fp_val , i_min ] = min( abs( f_p - 0.99 ) );

% Calculate the constant in the BL thickness
% equation
coef = eta( i_min );
fprintf( 'Constant in the BL thickness = %.1f\n',...
         coef );
     
% Plot the results
plot( eta , f_p )
grid on
xlabel( '{\eta}' )
ylabel( 'f''({\eta})' )
title( 'Dimensionless Velocity in the BL' )

function [ fpInf , y ] = RK4_blasius( myfun , ...
                          init_val , ...
                          interval , h )
                      
% Calculate the number of functions
N = length( myfun );

% Calculate the time array
t = interval( 1 ) : h : interval( 2 );

% Define the solution matrix y
% Columns -> diffferent functions 
% dependent on t
% Rows -> different time steps
y = zeros( length( t ) , N );

% Assign the initial values
y( 1 , : ) = init_val';

% Define slope approximation matrix
f = zeros( 4 , N );

% Conduct RK4 method
for k = 1 : length( t ) - 1
    
    for i = 1 : 4
        
        for j = 1 : N
            
            fun = myfun{ j };
            
            if i == 1
                pk = num2cell( y( k , : ) );
                f( i , j ) = fun( t( k ) , ...
                          pk{ : } );
            elseif i == 4
                temp = y( k , : ) + ...
                    h * f( i - 1 , : );
                pk = num2cell( temp );
                f( i , j ) = fun( t( k ) + h , ...
                          pk{ : } );
            else
                temp = y( k , : ) + ...
                    h/2 * f( i - 1 , : );
                pk = num2cell( temp );
                f( i , j ) = fun( t( k ) + h/2 , ...
                          pk{ : } );
            end
            
        end
        
        % Calculate the solutions
        y( k + 1 , : ) = y( k , : ) + ...
            h/6 .* ( f( 1 , : ) + ...
            2* f( 2 , : ) + ...
            2 * f( 3 , : ) + ...
            f( 4 , : ) );
        
    end
    
end

fpInf = y( end , 2 );

end