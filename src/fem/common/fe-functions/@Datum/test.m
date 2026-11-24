function I = test(obj, mesh, basis, options)
%TEST  Integrals of a Datum against the functions of a basis.
%
%   I = test(obj, mesh, basis, element_indices=)
%
%   I(i,1,e) is the integral over element (or face) e of mesh of the Datum
%   times the i-th function of basis, by the basis quadrature. I is
%   Dimension-by-1-by-n_entities, where n_entities counts the elements of
%   mesh for a two-dimensional basis and its faces for a one-dimensional
%   one; element_indices selects a subset of them (default all). basis
%   defaults to DubinerBasis(0).
arguments
    obj Datum
    mesh SimplicialGrid
    basis AbstractFEBasis = DubinerBasis(0);
    options.element_indices = [];
end
[x_local, w] = get_quadrature(basis);
p = basis.evaluate(x_local, mesh);
f_physical = obj.evaluate(x_local, mesh);

el_dim = basis.GeometricInputDimension;
switch el_dim
    case 2
        el_measure = reshape(mesh.volume, 1, 1, 1, 1, []);
    case 1
        el_measure = reshape(mesh.area, 1, 1, 1, 1, []);
end
I = el_measure .* pagemtimes(w .* p, 'transpose', f_physical, 'none');
I = sum(I, 3:4);

if isempty(options.element_indices)
    el_ind = 1:size(mesh, 3 - el_dim);
else
    el_ind = options.element_indices;
end
I = reshape(I(:, :, :, :, el_ind), size(p, 2), 1, []);
end
