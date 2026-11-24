function [x, w] = get_2d_reference_quadrature(k, element)
%GET_2D_REFERENCE_QUADRATURE  Quadrature rule on the reference element.
%
%   [x, w] = get_2d_reference_quadrature(k, element)
%
%   "Simplex" (default): the Grundmann-Moeller rule gm_rule_set(k, 2).
%   "Quad": the tensor Gauss-Legendre rule with k+1 points per direction
%   on [0, 1]^2, exact up to degree 2k+1 in each variable (the same
%   meaning of k as for "Simplex"). x is 2-by-n, w a column of weights.
arguments
    k (1,1) double                      % polynomial degree
    element (1,1) string = "Simplex"
end

switch element
    case "Simplex"
        rule = k; % exact up to 2*rule + 1
        [w, x] = gm_rule_set(rule, 2);
    case "Quad"
        rule = k + 1; % exact up to 2*k + 1 in each variable
        [x, w] = lgwt(rule, 0, 1);
        w = w*w';
        [x_1, x_2] = meshgrid(x);
        x = [x_1(:).'; x_2(:).'];
end
w = w(:);
end
