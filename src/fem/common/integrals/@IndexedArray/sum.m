function a_sum = sum(a)
%SUM  Sum of all entries, as a plain scalar.
    a_sum = sum(double(a), "all");
end