function dp = evaluate_gradient(basis, x, mesh)
%EVALUATE_GRADIENT  Gradient of a scalar basis at reference points.
%
%   dp = evaluate_gradient(basis, x, mesh)
%
%   x is a 2-by-N array of reference points. dp is
%   N-by-Dimension-by-2-by-1-by-n_el. Without mesh the reference
%   gradient is returned; with mesh it is mapped to the physical
%   elements by the inverse affine map (grad = B^-T * grad_local).
    arguments
        basis ScalarFEBasis
        x double
        mesh AbstractGrid = SimplicialGrid.empty;
    end
    [dPx, dPy] = basis.evaluate_derivatives(1, x, mesh);
    n_points = size(x, 2);
    n_basis = size(dPx, 2);
    if isempty(mesh)
        dp = reshape(cat(3, dPx, dPy), n_points, n_basis, 2, 1, []);
    else
        iB = get_affine_map(mesh, inverse=true, output="array");
        % grad = B^-T * grad_local
        dPx_ = reshape(iB(1, 1, :).*dPx + iB(2, 1, :).*dPy, ...
                        n_points, n_basis, 1, 1, []);
        dPy_ = reshape(iB(1, 2, :).*dPx + iB(2, 2, :).*dPy, ...
                        n_points, n_basis, 1, 1, []);
        dp = reshape(cat(3, dPx_, dPy_), n_points, n_basis, 2, 1, []);
    end
end