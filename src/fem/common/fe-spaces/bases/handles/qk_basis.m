function [P, dPx, dPy, dPxx, dPyy, dPxy] = qk_basis(k, r, s)
%QK_BASIS  Tensor Legendre basis of Q^k on the reference square.
%
%   [P, dPx, dPy, dPxx, dPyy, dPxy] = qk_basis(k, r, s)
%
%   Evaluates the (k+1)^2 functions of QkBasis(k) at the points (r, s); if
%   s is omitted r is 2-by-n. The handle twin of dubiner_basis for the
%   function-handle (low-level) path on a CartesianGrid.
%
%See also: dubiner_basis, QkBasis
    arguments
        k double
        r double
        s double = r(2, :);
    end
    if nargin < 3
        r = r(1, :);
    end
    [P, dPx, dPy, dPxx, dPyy, dPxy] = QkBasis.evaluate_basis(k, r, s);
end
