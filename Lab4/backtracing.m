function X = backtracing(Aug)

% calculate N
[N, M] = size(Aug);

% define unknown vector, X
X = zeros(N,1);

X(N) = Aug(N,N+1) / Aug(N,N);

for i = N-1 : -1 : 1
    X(i) = (Aug(i,N+1) - Aug(i,i+1:N)*X(i+1:N)) / Aug(i,i);
end

% ------------------------------------------------
% FOR SYSTEM OF 3 EQNS ONLY
% X(3) = Aug(3,4) / Aug(3,3);
% X(2) = (Aug(2,4) - Aug(2,3)*X(3)) / Aug(2,2);
% X(1) = (Aug(1,4) - Aug(1,3)*X(3) - Aug(1,2)*X(2)) / Aug(1,1);
% ------------------------------------------------

end