% ME 303 - Non-Linear Least Squares Code
% 23.11.2021
%-------------------------------------------

data = [ 1  28
         2  34
         3  36
         4  38
         5  39 ];

% Calculate N
N = length( data ); 
%% Data Linearization

% Conduct variable change
X = 1 ./ data( : , 1 );
Y = 1 ./ data( : , 2 );


% Calculate the coefficients A and B
coef1 = lsline( X , Y );
A1 =  coef1( 1 );
B1 =  coef1( 2 );




f1 = @( x ) x ./ (A1 .* x + B1);

     
% Define RMSE function
E_2 = @( x ) sum( 1 / N * ...
        ((data(: , 1) ./( x( 2 ) + ( x( 1 ) .* data( : , 1 ) ))) - data( : , 2 ) ) .^ 2 );

% Specify the initial guess for fminsearch
x0 = [ A1 , B1 ];

[ coef , E_2 ] = fminsearch( E_2 , x0 );

% Calculate the RMSE
E_rms = sqrt( E_2 );

% Calculate the coefficients
A = coef( 1 );
B = coef( 2 );

f2 = @( x ) x ./ (A .* x + B);

fprintf('C = %10f\nA = %10f\nRMSE = %10f\n',...
         A , B , E_rms );
     
x = 0 : 0.01 : 6;

plot( x , f2( x ) )
grid on
hold on
scatter( data( : , 1 ) , data( : , 2 ),...
         'rs','filled' )
legend( 'Non-linear Least Square Fit' , 'Data' )



% lsline function
function [ C ] = lsline( X , Y )

% Calculate N
N = length( X );

% Calculate the coefficient matrix
% [ sum(x^2) sum(x)
%   sum(x)     N   ]
A = [ sum( X.^2 ) , sum( X )
      sum( X )    ,    N    ];
  
% Calculate the output vector
B = [ sum( X .* Y ) ; sum( Y ) ];

C = A \ B;

end