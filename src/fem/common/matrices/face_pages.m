function A = face_pages(A)
% builders return squeezed arrays; the factories work in the elemental
% page convention [d, d, 1, 1, n_el, 3]
    A = permute(A, [1 2 5 6 4 3]);
end
