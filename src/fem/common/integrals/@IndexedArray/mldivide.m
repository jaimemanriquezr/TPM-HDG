function a_div_b = mldivide(a, b)
%MLDIVIDE  Page-wise a\b of the underlying arrays, as a plain array.
    a_div_b = pagemldivide(double(a), double(b));
end