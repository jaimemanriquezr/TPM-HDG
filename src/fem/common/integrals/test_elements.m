function I = test_elements(mesh, f, k, options)
%TEST_ELEMENTS  Integrals of a function against a basis on each element.
%
%   I = test_elements(mesh, f, k, basis=@dubiner_basis, subset=)
%
%   f is a handle @(x, y) evaluated at the quadrature points of order k
%   (quadpoints_elements_2d). I(i,1,e) is the integral over the e-th
%   element of subset (default all elements) of f times the i-th function
%   of basis(k, ...); I is Dimension-by-1-by-numel(subset).
arguments
    mesh SimplicialGrid
    f function_handle
    k double = 0
    options.basis function_handle = @dubiner_basis
    options.subset = 1:size(mesh.elements, 2);
end
% QUADRATURE
[xy, rs, w] = quadpoints_elements_2d(mesh, k, options.subset, output="cell");

% COLLOCATION
p = options.basis(k, rs{:});

% INTEGRATION
I = reshape(mesh.volume(options.subset), 1, []) .* ((w.*p)' * f(xy{:}));

% RESHAPING
d2 = size(p, 2);
N_ELEMENTS = numel(options.subset);
I = reshape(I, [d2, 1, N_ELEMENTS]);
end
