function B = hybrid_rows(Ax, Ay, d1)
% component-major hybrid rows (x block / y block, each 3d1 tall) ->
% per-face rows [x(d1); y(d1)] with the face index in dim 6
nc = size(Ax, 2); n_el = size(Ax, 3);
B = zeros(2*d1, nc, 1, 1, n_el, 3);
for i = 1:3
    fi = (i-1)*d1 + (1:d1);
    B(1:d1, :, 1, 1, :, i)      = reshape(Ax(fi, :, :), d1, nc, 1, 1, n_el);
    B(d1+(1:d1), :, 1, 1, :, i) = reshape(Ay(fi, :, :), d1, nc, 1, 1, n_el);
end
end
