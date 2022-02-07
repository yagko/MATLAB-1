function [ y ] = RK4_fun(myfun , init_val , time_interval , step_size )
                      
% Calculate the number of functions
N = length(myfun);

% Calculate the time array
t = time_interval(1):step_size:time_interval(2);


y = zeros(length(t),N);

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
                f( i , j ) = fun( t( k ) , pk{ : } );
            elseif i == 4
                temp = y( k , : ) + step_size * f( i - 1 , : );
                pk = num2cell( temp );
                f( i , j ) = fun( t( k ) + step_size , pk{ : } );
            else
                temp = y( k , : ) + step_size/2 * f( i - 1 , : );
                pk = num2cell( temp );
                f( i , j ) = fun( t( k ) + step_size/2 , pk{ : } );
            end
            
        end
        
        % Calculate the solutions
        y(k+1, : ) = y(k, : ) + step_size/6.*(f(1, : ) + 2*f(2, : ) + 2* f(3, : ) + f(4, : ));
        
    end   
end
end