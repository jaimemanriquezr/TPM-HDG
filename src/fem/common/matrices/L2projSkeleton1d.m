function proj_u = L2projSkeleton1d(mesh, u, k, options)
%L2PROJSKELETON1D  Face-wise L2 projection of a function onto a polynomial space.
%
%   proj_u = L2projSkeleton1d(mesh, u, k, subset=, basis=@legendre_basis)
%
%   mesh    SimplicialGrid.
%   u       function handle u(x, y), or a cell array of such handles, each
%           projected componentwise.
%   k       polynomial degree (default 0), also the quadrature order.
%   subset  faces to project (default all faces).
%   basis   basis handle of fe-spaces/bases/handles (default
%           @legendre_basis).
%
%   proj_u  dofs, n_basis-by-n_faces; for a cell u the components are
%           stacked by cell2mat.
arguments
    mesh SimplicialGrid
    u
    k double = 0
    options.subset = 1:size(mesh.faces, 2)
    options.basis = @legendre_basis
end
% QUADRATURE
[xy, t, w] = quadpoints_faces_2d(mesh, k, options.subset, output="cell");

% COLLOCATION
l = options.basis(k, t);
lw = w .* l;

% L2 PROJECTION
if isa(u, "function_handle")
    proj_u = (lw'*l) \ (lw'*u(xy{:}));
elseif isa(u, "cell")
    proj_u = cell(numel(u), 1);
    for i = 1:numel(u)
        proj_u{i} = L2projSkeleton1d(mesh, u{i}, k, subset=options.subset, basis=options.basis);
    end
    proj_u = cell2mat(proj_u);
end
end
