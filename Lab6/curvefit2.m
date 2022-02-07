clc;clear all;close all;

%-------------------------------------------
% ME 303 - Non-Linear Least Squares &
% Data Linearization Code
% 23.11.2021
%-------------------------------------------

data = [ 1  0.6
         2  1.9
         3  4.3
         4  7.6
         5 12.6 ];

% Calculate N
N = length( data );

%% Data Linearization

% Conduct variable change
X = data( : , 1 );
Y = log( data( : , 2 ) );

% Calculate the coefficients A and B
coef1 = lsline_linear( X , Y );
A1 = coef1( 1 );
B1 = coef1( 2 );

C1 = exp( B1 );

expfit1 = @( x ) C1 .* exp( A1 .* x );

E_rms1 = sqrt( sum( 1 / N * ...
         ( expfit1( X ) - Y ) .^ 2 ) );
     
%% Non-linear Least Squares Curve

% Define RMSE function
E_2 = @( x ) sum( 1 / N * ...
        ( x( 1 ) .* exp( x( 2 ) .* ...
        data( : , 1 ) ) - ...
        data( : , 2 ) ) .^ 2 );

% Specify the initial guess for fminsearch
x0 = [ C1 , A1 ];

[ coef2 , E_2 ] = fminsearch( E_2 , x0 );

% Calculate the RMSE
E_rms2 = sqrt( E_2 );

% Calculate the coefficients
C2 = coef2( 1 );
A2 = coef2( 2 );

expfit2 = @( x ) C2 .* exp( A2 .* x );

%% Non-linear Least Squares Power Fit

% A*x^B

% Define RMSE function
E_2pfit = @( x ) sum( 1 / N * ...
        ( x( 1 ) .* data( : , 1 ) .^ ...
        x( 2 ) - data( : , 2 ) ) .^ 2 );

% Specify the initial guess for fminsearch
x0 = [ C1 , A1 ];

[ coef3 , E_2pfit ] = fminsearch( E_2pfit , x0 );

% Calculate the RMSE
E_rms3 = sqrt( E_2pfit );

% Calculate the coefficients
C3 = coef3( 1 );
A3 = coef3( 2 );

powfit = @( x ) C3 .* x .^ A3;

fprintf( 'Results for the Data Linearization Method\n' );
fprintf( '-----------------------------------------------\n' );
fprintf( 'C = %10f\nA = %10f\nRMSE = %10f\n' , C1 , A1 , E_rms1 );
fprintf( '-----------------------------------------------\n' );
fprintf( 'Results for the Non-linear Least Squares Method\n' );
fprintf( '-----------------------------------------------\n' );
fprintf( 'C = %10f\nA = %10f\nRMSE = %10f\n' , C2 , A2 , E_rms2 );
fprintf( '-----------------------------------------------\n' );
fprintf( 'Results for the Power Fit\n' );
fprintf( '-----------------------------------------------\n' );
fprintf( 'A = %10f\nB = %10f\nRMSE = %10f\n' , C3 , A3 , E_rms3 );

x = 0 : 0.01 : 6;

plot( x , expfit1( x ) )
grid on
hold on
plot( x , expfit2( x ) )
plot( x , powfit( x ) )
scatter( data( : , 1 ) , data( : , 2 ) , ...
         'rs' , 'filled' )
legend( 'y=Ce^{Ax} (Data Lin.)' , ...
        'y=Ce^{Ax} (Nonlinear LS)' , ... 
        'y=Ax^B' , ...
        'Data' )