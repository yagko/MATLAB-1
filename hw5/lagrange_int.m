function [poly_deg , E_max] = lagrange_int(lower,upper)

% function definition
f = @(t) cos(t).*exp(-0.2*t);
degree = 1;

while 1 % declearing a while loop to iterate until error conditions are satisfied
    s = (upper-lower)/degree; % spacing
    X = lower: s : upper; % data points
    N = length(X) - 1; % length of X 
    Y = f(X); % data points
    poly = 0; % initializing the polynomial
    syms t % defining the t in the symbolic form
    if (lower >= 0) && (upper <= 40)   % bounds can not be exceeded
         for i = 1 : N+1
            L_i = 1;
            for j = 1 : N + 1
                if i ~= j
                    L_i = L_i * (t - X(j)) / (X(i)-X(j));
                end
            end
            poly = poly + Y(i) * L_i;
            poly = matlabFunction(poly); % redefining the polynomial in the function form
            %Error defining
            x = lower : 0.001 : upper;
            %Error calculation
            err = abs(f(x) - poly(x)); 
         end 
    end
    if (err < 10^-3) % breaking the while loop condititon
        break;
    end
        degree = degree + 1 ; % increment the degree until the loop is breaked
end
poly_deg = polynomialDegree(poly(t));
E_max = max(err);