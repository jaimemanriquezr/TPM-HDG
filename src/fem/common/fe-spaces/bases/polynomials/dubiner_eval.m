function [P, dPx, dPy, dPxx, dPyy, dPxy] = dubiner_eval(j_list, r, s)
%-------------------------------------------------------------------------
%Description: Evaluates the Dubiner basis functions of indices j_list, and
%           their first and second derivatives, at the natural
%           coordinates (r, s)
%Input:     - j_list: basis function indices (one column of output each)
%           - r: first natural coordinate
%           - s: second natural coordinate
%Output:    - P: basis function values
%           - dPx, dPy: first partial derivatives w.r.t. r and s
%           - dPxx, dPyy, dPxy: second partial derivatives
%-------------------------------------------------------------------------
arguments
    j_list (1, :) double
    r (:, :, :) double
    s (:, :, :) double
end
r = reshape(r(:), [], size(r, 3));
s = reshape(s(:), [], size(s, 3));
assert(length(r) == length(s), "TPMHDG:inputSizeMismatch", "Coordinate vectors do not match in length.")

% (r, s): coordinates in the triangle [(0, 0), (1, 0), (0, 1)]
% The collapsed coordinate x = 2r/(1 - s) - 1 enters only through
% P_n(x) (1 - y)^n = 2^n Q_n(z, t), z = 2r + s - 1, t = 1 - s, with
% Q_n(z, t) = t^n P_n(z/t) a polynomial: no division, valid for all (r, s).
z = 2*r + s - 1;
t = 1 - s;
y = 2*s - 1;
dy_ds = 2;

[P, dPx, dPy, dPxx, dPxy, dPyy] = deal(zeros(size(r, 1), length(j_list), size(r, 2)));
for i = 1:length(j_list)
    j = j_list(i);
    cj = ceil(-3/2 + .5*sqrt(1 + 8*j));
    n = (cj+1)*(cj+2)/2 - j;
    m = cj - n;

    % NORMALIZATION
    np = 2^n * 1/sqrt((2*n+1)*(2*m+2*n+2)); % \phi_00 = sqrt(2)
    c = 2^n/np;

    [Q, Qz, Qt, Qzz, Qzt, Qtt] = scaled_legendre(n, z, t);
    [p2, dp2_dy, ddp2_dydy] = jacobi_eval(m, 2*n+1, 0, y);
    dp2 = dp2_dy*dy_ds;
    Qs = Qz - Qt;                               % d/ds of Q(2r + s - 1, 1 - s)

    P(:,i,:) = c*Q.*p2;
    if nargout > 1
        dPx(:,i,:) = c*2*Qz.*p2;
        dPy(:,i,:) = c*(Qs.*p2 + Q.*dp2);
    end
    if nargout > 3
        ddp2 = ddp2_dydy*dy_ds^2;
        dPxx(:,i,:) = c*4*Qzz.*p2;
        dPxy(:,i,:) = c*2*((Qzz - Qzt).*p2 + Qz.*dp2);
        dPyy(:,i,:) = c*((Qzz - 2*Qzt + Qtt).*p2 + 2*Qs.*dp2 + Q.*ddp2);
    end
end
end

function [Q, Qz, Qt, Qzz, Qzt, Qtt] = scaled_legendre(n, z, t)
% Q_n(z, t) = t^n P_n(z/t) and its first and second derivatives, by the
% Legendre recurrence (m + 1) Q_m+1 = (2m + 1) z Q_m - m t^2 Q_m-1
o = zeros(size(z));
[Q, Qz, Qt, Qzz, Qzt, Qtt] = deal(o + 1, o, o, o, o, o);
[R, Rz, Rt, Rzz, Rzt, Rtt] = deal(o, o, o, o, o, o);          % Q_m-1
for m = 0:n - 1
    Q1   = ((2*m + 1)*z.*Q - m*t.^2.*R)/(m + 1);
    Qz1  = ((2*m + 1)*(Q + z.*Qz) - m*t.^2.*Rz)/(m + 1);
    Qt1  = ((2*m + 1)*z.*Qt - m*(2*t.*R + t.^2.*Rt))/(m + 1);
    Qzz1 = ((2*m + 1)*(2*Qz + z.*Qzz) - m*t.^2.*Rzz)/(m + 1);
    Qzt1 = ((2*m + 1)*(Qt + z.*Qzt) - m*(2*t.*Rz + t.^2.*Rzt))/(m + 1);
    Qtt1 = ((2*m + 1)*z.*Qtt - m*(2*R + 4*t.*Rt + t.^2.*Rtt))/(m + 1);
    [R, Rz, Rt, Rzz, Rzt, Rtt] = deal(Q, Qz, Qt, Qzz, Qzt, Qtt);
    [Q, Qz, Qt, Qzz, Qzt, Qtt] = deal(Q1, Qz1, Qt1, Qzz1, Qzt1, Qtt1);
end
end
