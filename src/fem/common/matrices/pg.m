function B = pg(A)
% [r x c x n_el] -> page layout [r, c, 1, 1, n_el]
B = reshape(A, size(A, 1), size(A, 2), 1, 1, size(A, 3));
end
