function [proj_u, error] = L2proj2D(mesh, u, k, options)
%L2PROJ2D  Element-wise L2 projection of a function onto a polynomial space.
%
%   [proj_u, error] = L2proj2D(mesh, u, k, basis = @dubiner_basis)
%
%   mesh   SimplicialGrid.
%   u      function handle u(x, y), or a cell array of such handles, each
%          projected componentwise.
%   k      polynomial degree (default 0).
%   basis  basis handle of fe-spaces/bases/handles (default
%          @dubiner_basis).
%
%   proj_u  dofs, n_basis-by-n_elements; for a cell u the components are
%           stacked by cell2mat.
%   error   L2 projection error; for a cell u, the Euclidean norm of the
%           component errors.
%
%See also: mass_matrix, dubiner_basis
arguments
    mesh SimplicialGrid
    u
    k double = 0
    options.basis = @dubiner_basis
end
[xy, rs, w] = quadpoints_elements_2d(mesh, k, output="cell");
p = options.basis(k, rs{:});
wp = w .* p;

N_ELEMENTS = size(mesh.elements, 2);

if isa(u, "function_handle")
    U = reshape(u(xy{:}), size(p, 1), N_ELEMENTS);

    % L2 PROJECTION: the reference mass matrix wp'*p is element-independent
    % (affine elements, the volume cancels), so one solve covers all elements
    proj_u = (wp'*p) \ (wp'*U);
    error = error2d(mesh, u, proj_u, k);
elseif isa(u, "cell")
    proj_u = cell(numel(u), 1);
    error = zeros(numel(u), 1);
    for i = 1:numel(u)
        [proj_u{i}, error(i)] = L2proj2D(mesh, u{i}, k, basis=options.basis);
    end
    proj_u = cell2mat(proj_u);
    error = norm(error, 2);
end
end
