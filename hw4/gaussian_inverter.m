function inv_A = gaussian_inverter(A)

N = length(A);
inv_A = zeros(N);
for k = 1 : N
    inv_A(k,k) = 1;
end

for i = 1 : N - 1
    for j = i + 1 : N
        C = A(j,i) / A(i,i);
        A(j,:) = A(j,:) - C*A(i,:);
        inv_A(j,:) = inv_A(j,:) - C*inv_A(i,:);
    end
end

for i = 1 : N-1
    for j = i+1 : N
        C = A(j,j) / A(i,j);
        A(i,:) = A(j,:)/C - A(i,:);
        inv_A(i,:) = inv_A(j,:)/C - inv_A(i,:);
    end
end

for i = 1 : N
    C = A(i,i);
    A(i,:) = A(i,:)/C;
    inv_A(i,:) = inv_A(i,:)/C;
end