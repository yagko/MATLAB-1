% ------------------------------------------
% Me 303 - Lab Session 11
% Boundary Layer on a Flat Plate
% (Blasius Solution)
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

% Specify the initial guesses for the z3_0
z3_0_guess = [ 0 2 ];

% Specify the interval
eta_min = 0;
eta_max = 10;

% Specify the step size, h
h = 0.1;

% Calculate the eta array
eta = ( eta_min : h : eta_max )';

% Specify the tolerance for the error in f'(Inf)
tol = 1E-6;

% Define f'(Inf) array
fpInf = zeros( 1 , 2 );

% Solve for f, f', and f'' with the 1st initial 
% guess
y = RK4_fun( { g1 ; g2 ; g3 } , [ z1_0 ; z2_0 ;...
              z3_0_guess( 1 ) ] , ...
              [ eta_min ; eta_max ] , h );
          
% Derive f'(Inf) from the resultant matrix
fpInf( 1 ) = y( end , 2 );

% Solve for f, f', and f'' with the 2nd initial 
% guess
y = RK4_fun( { g1 ; g2 ; g3 } , [ z1_0 ; z2_0 ;...
              z3_0_guess( 2 ) ] , ...
              [ eta_min ; eta_max ] , h );
          
% Derive f'(Inf) from the resultant matrix
fpInf( 2 ) = y( end , 2 );

% Specify the maximum number of iterations
maxI = 1000;

for iter = 1 : maxI
    
    % Approximate the initial guess for f'' that
    % yields f'(Inf)=1 using Lagrange Approximation
    z3_0_interp = lagrange_app( fpInf , ...
                                      z3_0_guess );
                                  
    % Solve for f, f', and f'' with the updated
    % initial guess
    y = RK4_fun( { g1 ; g2 ; g3 } , ...
                [ z1_0 ; z2_0 ; z3_0_interp ] ,...
                [ eta_min ; eta_max ] , h );
            
    % Derive f'(Inf) and append to fpInf array
    fpInf = [ fpInf , y( end , 2 ) ];
    
    % Break the loop if the absolute error in
    % f'(Inf) is less than the predefined 
    % tolerance value
    if abs( fpInf( end ) - 1 ) < tol
        break
    end
    
    % Update the z3 initial guess array by 
    % appending the updated initial guess
    z3_0_guess = [ z3_0_guess , z3_0_interp ];
    
end

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

function z3_0_interp = lagrange_app( X , Y )

% Calculate
% N+1 data points -> length(X)
N = length( X ) - 1;

% Assign initial value for p_n
p_n = @( x ) 0;

% Calculate L_k(x)...
for k = 1 : N + 1
    
    L_k = @( x ) 1;
    
    for j = 1 : N + 1
        
        if k ~= j
            
            L_k = @( x ) L_k( x ) .* ( x - X( j ) ) ...
                  ./ ( X( k ) - X( j ) );
            
        end
        
    end
    
    p_n = @( x ) p_n( x ) + Y( k ) .* L_k( x );
    
end

z3_0_interp = p_n( 1 );

end