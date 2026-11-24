function B = hybrid_cols(A, d1)
% component-major hybrid columns [x(3d1) | y(3d1)] -> per-face columns
% [x(d1); y(d1)] with the face index in dim 7
ds = size(A, 1); n_el = size(A, 3);
Ax = A(:, 1:3*d1, :); Ay = A(:, 3*d1 + (1:3*d1), :);
B = zeros(ds, 2*d1, 1, 1, n_el, 1, 3);
for i = 1:3
    fi = (i-1)*d1 + (1:d1);
    B(:, 1:d1, 1, 1, :, 1, i)      = reshape(Ax(:, fi, :), ds, d1, 1, 1, n_el);
    B(:, d1+(1:d1), 1, 1, :, 1, i) = reshape(Ay(:, fi, :), ds, d1, 1, 1, n_el);
end
end
