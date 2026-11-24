function [P, dPx, dPy, dPxx, dPyy, dPxy] = dubinertilde_basis(k, r, s)
% Dubiner polynomials of degree exactly "k"
d2 = nchoosek(k+2, 2);
if nargout > 3
    [P, dPx, dPy, dPxx, dPyy, dPxy] = dubiner_eval((d2 - k):d2, r, s);
elseif nargout > 1
    [P, dPx, dPy] = dubiner_eval((d2 - k):d2, r, s);
else
    P = dubiner_eval((d2 - k):d2, r, s);
end
end