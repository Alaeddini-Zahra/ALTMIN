function X = Column_Nonzero(A) 
% This function returns X that records the row number of nonzero entries
% at each corresponding column.

[m,n] = size(A);
X = zeros(m,n);
for j = 1:n
    t = 1;
    for i = 1:m
        if A(i,j) ~= 0
            X(t,j) = i;
            t = t + 1;
        end
    end
end

