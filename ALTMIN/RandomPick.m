function A = RandomPick(B,p)
% B is all-zero matrix. This function selects each entry with probability
% p.
[m,n] = size(B);
A = B;
for i = 1:m
    for j = 1:n
        if rand(1) < p
            A(i,j) = 1;
        end
    end
end
end