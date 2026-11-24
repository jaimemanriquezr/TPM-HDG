function E = compute_l2_error(obj, dof, mesh, options)
%COMPUTE_L2_ERROR  Element-based L2 error of a Dubiner field.
%
%   E = compute_l2_error(obj, dof, mesh, exact_solution=)
%
%   dof holds the coefficients of the field, one column per element of
%   mesh. exact_solution is a handle @(x, y) evaluated at the element
%   quadrature points of degree Order + 1 (default zero). E is the scalar
%   L2 error over mesh.
    arguments
        obj ScalarFEBasis
        dof double
        mesh SimplicialGrid
        options.exact_solution function_handle = @(x,y) 0*(x.*y);
    end
    exact_solution = options.exact_solution;

    % QUADRATURE
    % degree k + 1, as FEField.compute_l2_error
    [xy, rs, w] = quadpoints_elements_2d(mesh, obj.Order + 1, output="array");

    % COLLOCATION
    P = obj.evaluate(rs) * dof;

    % INTEGRATION
    P_ex = reshape(exact_solution(xy(1,:,:),xy(2,:,:)), size(P));
    E = sqrt(sum(mesh.volume .* (w'*(P - P_ex).^2)));
end