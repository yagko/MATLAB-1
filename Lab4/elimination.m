function Aug = elimination(A,B)

% Calculate number of eqns
N = length(B);

Aug = [A B];

% perform row operations
% Step 1 = R2- R1 -> R2
% Step 2 = R3 + R1/2 -> R3
% Step 3 = R3 +R2/14 -> R3

for i_p = 1 : N-1
    for i_c = i_p+1 : N
        C = Aug(i_c,i_p) / Aug(i_p,i_p);
        Aug(i_c,:) = Aug(i_c,:) - C*Aug(i_p,:);
    end
end

% ------------------------------------------------
% FOR SYSTEM OF 3 EQNS ONLY
% % Step 1
% C = Aug(2,1)/Aug(1,1);
% Aug (2,:) = Aug(2,:) - C*Aug(1,:);
% 
% % Step 2
% C = Aug(3,1)/Aug(1,1);
% Aug (3,:) = Aug(3,:) - C*Aug(1,:);
% 
% % Step 1
% C = Aug(3,2)/Aug(2,2);
% Aug (3,:) = Aug(3,:) - C*Aug(2,:);
% -------------------------------------------------

end