function E = error2d(mesh, u, uh, k, options)
%ERROR2D  L2 error of a discrete field against a function.
%
%   E = error2d(mesh, u, uh, k, basis=@dubiner_basis)
%
%   u is a function handle u(x, y) and uh the dofs of the field in basis
%   of degree k (default 0), n_basis-by-n_elements. k is also the
%   quadrature order. E is the scalar L2 error over mesh.
arguments
    mesh SimplicialGrid
    u function_handle
    uh double
    k double = 0
    options.basis = @dubiner_basis
end
[xy, rs, w] = quadpoints_elements_2d(mesh, k, output="cell");
p = options.basis(k, rs{:});

% INTEGRATION
l2_error = @(u, uh) sqrt(sum(mesh.volume .* abs(w' * ((p*uh - u(xy{:})).^2) ) ));
if isequal(class(u),'cell')
    E = cell(size(u));
    for i = 1:numel(u)
        E{i} = l2_error(u{i}, uh{i});
    end
else
    E = l2_error(u, uh);
end
end
