function l2_error = compute_l2_error(obj, mesh, u, p)
%COMPUTE_L2_ERROR  L2 distance between a discrete field and a Datum.
%
%   l2_error = compute_l2_error(obj, mesh, u, p)
%
%   Evaluates obj and the Datum u at the quadrature points of obj.basis,
%   integrates the squared difference over each element (for a
%   one-dimensional basis: each face F, weighted by h_F = |F|) and
%   returns the p-norm (default 2) of the vector of elementwise errors.
%   u defaults to zero, which gives the norm of obj itself.
arguments
    obj FEField
    mesh AbstractGrid
    u Datum = Datum(@(x, y) 0 * (x.*y).^0);
    p double = 2;
end
[x_local, w] = get_quadrature(obj.basis);
u_h = obj.evaluate(x_local, mesh);
u_ex = u.evaluate(x_local, mesh);

in_dim = obj.basis.GeometricInputDimension;
switch in_dim
    case 2
        el_measure = reshape(mesh.volume, 1, 1, []);
    case 1
        el_measure = reshape(mesh.area, 1, 1, []).^2;
end
out_dim = obj.basis.GeometricOutputDimension;
if out_dim == 1
    out_dim = [1; 1];
end

e_h = u_h - u_ex;
e_h_sq = pagemtimes(w .* e_h , 'transpose', e_h, 'none');
e_h_sq = reshape(e_h_sq, out_dim(1), out_dim(2), []);
E = sqrt(sum(el_measure .* e_h_sq, 3));

l2_error = norm(E(:), p);
end
