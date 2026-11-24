function E = compute_l2_error(obj, dof, mesh, options)
%COMPUTE_L2_ERROR  h-weighted skeleton L2 error of a Legendre field.
%
%   E = compute_l2_error(obj, dof, mesh, exact_solution=)
%
%   dof holds the coefficients of the field, one column per face of mesh.
%   exact_solution is a handle @(x, y) evaluated at the face quadrature
%   points (default zero). E = (sum_F h_F ||dof - exact||_F^2)^(1/2)
%   with h_F = |F| the face length.
    arguments
        obj ScalarFEBasis
        dof double
        mesh AbstractGrid
        options.exact_solution function_handle = @(x,y) 0*(x.*y);
    end
    exact_solution = options.exact_solution;

    % QUADRATURE
    % k + 2 points, as FEField.compute_l2_error: the k + 1 Gauss points
    % are the zeros of the leading term of the projection error
    [xy, t, w] = quadpoints_faces_2d(mesh, obj.Order + 1, output="cell");

    % COLLOCATION
    l = obj.evaluate(t);

    % INTEGRATION
    % weighted by "h_e"
    ex = exact_solution(xy{:});
    E = sqrt(sum((mesh.area).^2 .* (w'*((l*dof - ex).^2))));
end