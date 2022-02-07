function [ y ] = RK4_fun( myfun , init_val , ...
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

end