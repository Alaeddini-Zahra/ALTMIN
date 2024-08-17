function V = cut_zero(U); % Cutting zeros at the end of u
i = 1;
while i <= length(U) && U(i) ~= 0
    i=i+1;
end
V= U(1:i-1)
