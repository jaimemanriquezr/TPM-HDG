function [P, dPx, dPy, dPxx, dPyy, dPxy] = dubiner_basis(k, r, s)
%DUBINER_BASIS  Dubiner basis of degree up to k on the reference triangle.
%
%   [P, dPx, dPy, dPxx, dPyy, dPxy] = dubiner_basis(k, r, s)
%
%   Evaluates the nchoosek(k+2, 2) functions at the points (r, s); if s is
%   omitted r is 2-by-n. The outputs requested (1, 3 or 6) select the
%   values, the first derivatives and the second derivatives.
    arguments
        k double
        r double
        s double = r(2, :);
    end
    if nargin < 3
        r = r(1, :);
    end

    d2 = nchoosek(k+2, 2);
    if nargout > 3
        [P, dPx, dPy, dPxx, dPyy, dPxy] = dubiner_eval(1:d2, r, s);
    elseif nargout > 1
        [P, dPx, dPy] = dubiner_eval(1:d2, r, s);
    else
        P = dubiner_eval(1:d2, r, s);
    end
end