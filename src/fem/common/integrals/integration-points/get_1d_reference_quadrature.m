function [s, w] = get_1d_reference_quadrature(k, element)
%GET_1D_REFERENCE_QUADRATURE  Gauss-Legendre rule on the reference interval.
%
%   [s, w] = get_1d_reference_quadrature(k, element)
%
%   "Simplex" (default): k+1 points on [0, 1]. "SymmetricInterval": k
%   points on [-1, 1]. "Quad": k+1 points on [0, 1] (as "Simplex"), one
%   factor of the "Quad" rule of get_2d_reference_quadrature. s is a row
%   of points, w a column of weights.
arguments
    k (1,1) double                      % polynomial degree
    element (1,1) string = "Simplex"
end

switch element
    case "Simplex"
        rule = k + 1;
        [s, w] = lgwt(rule, 0, 1);
    case "Quad"
        % one factor of the "Quad" rule of get_2d_reference_quadrature
        rule = k + 1;
        [s, w] = lgwt(rule, 0, 1);
    case "SymmetricInterval"
        rule = k;
        [s, w] = lgwt(rule, -1, 1);
end

s = reshape(s, 1, []);
w = w(:);
end

