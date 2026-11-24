function varargout = evaluate_derivatives(basis, n, points, mesh, options)
%EVALUATE_DERIVATIVES  The n-th order partial derivatives of a scalar basis.
%
%   [d1, ..., d_{n+1}] = evaluate_derivatives(basis, n, points, mesh)
%
%   Declared on AbstractFEBasis, but defined for scalar bases only (basis
%   must be a ScalarFEBasis). n is the derivative order; the outputs are
%   the last n+1 outputs of basis.evaluate (n = 1: {dPx, dPy}). points,
%   mesh and the coordinate_system option are forwarded to basis.evaluate.
arguments
    basis ScalarFEBasis
    n double
    points (:, :) double
    mesh AbstractGrid = SimplicialGrid.empty
    options.coordinate_system string = "local"
end
op_args = namedargs2cell(options);

all_derivatives = cell(nchoosek(n + 2, 2), 1);
[all_derivatives{:}] = basis.evaluate(points, mesh, op_args{:});
% the n-th derivative block is the last n+1 outputs of evaluate
% (e.g. n = 1: {dPx, dPy} out of {P, dPx, dPy})
varargout = all_derivatives(end-n:end);
end
