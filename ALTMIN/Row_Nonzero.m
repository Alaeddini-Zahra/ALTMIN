function X = Row_Nonzero(A) 
% This function returns X that records the column number of nonzero entries
% at each corresponding row.

[m,n] = size(A);
X = zeros(m,n);
for i = 1:m
    t = 1;
    for j = 1:n
        if A(i,j) ~= 0
            X(i,t) = j;
            t = t + 1;
        end
    end
end

