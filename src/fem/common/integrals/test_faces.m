function I = test_faces(mesh, f, k, options)
%TEST_FACES  Integrals of a function against a basis on each face.
%
%   I = test_faces(mesh, f, k, subset=, basis=@legendre_basis)
%
%   f is a handle @(x, y) evaluated at the quadrature points of order k
%   (quadpoints_faces_2d). I(i,1,e) is the integral over the e-th face of
%   subset (default all faces; a row or column list) of f times the i-th
%   function of basis(k, ...); I is Dimension-by-1-by-numel(subset).
arguments
    mesh SimplicialGrid
    f function_handle
    k double = 0
    options.subset = 1:size(mesh.faces, 2);
    options.basis function_handle = @legendre_basis
end
options.subset = reshape(options.subset, 1, []);   % accept row or column face lists
% QUADRATURE
[xy, t, w] = quadpoints_faces_2d(mesh, k, options.subset, output="cell");

% COLLOCATION
l = options.basis(k, t);

% INTEGRATION
I = mesh.area(options.subset) .* ((w.*l)' * f(xy{:}));

% RESHAPING
Mh_dim = size(l, 2);
Nfaces = size(options.subset, 2);
I = reshape(I, [Mh_dim, 1, Nfaces]);
end
