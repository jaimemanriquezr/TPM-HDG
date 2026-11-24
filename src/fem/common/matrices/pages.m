function A = pages(A)
% builders return squeezed arrays; the factories work in the elemental
% page convention [d, d, 1, 1, n_el]
    A = permute(A, [1 2 4 5 3]);
end
