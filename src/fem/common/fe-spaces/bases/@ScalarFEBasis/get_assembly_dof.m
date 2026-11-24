function [rows, cols, rows_rhs] = get_assembly_dof(basis, mesh)
% Per-element assembly scatter blocks: for each element, the global dofs
% of its incident support entities, expanded to (rows x cols) index grids
% over subentity-pair slots. The slot count and incidence come from the
% mesh, the per-entity block is Dimension contiguous dofs.
arguments
    basis ScalarFEBasis
    mesh AbstractGrid
end
d = basis.Dimension;
codimension = mesh.geometric_dimension - basis.GeometricInputDimension;
switch codimension
    case 0
        E = 1:size(mesh.elements, 2);      % the element itself
    case 1
        E = mesh.boundaries;               % its faces
    otherwise
        error("ScalarFEBasis:enumeration", ...
            "no local enumeration for codimension %d support", codimension);
end
ns = size(E, 1);                           % subentity slots per element

D = (E(:) - 1)*d + (1:d);                          % (ns*n_el) x d
D = permute(reshape(D.', d, ns, []), [1 3 2]);     % d x n_el x ns

[J_d, I_d] = meshgrid(1:d);
[J_s, I_s] = meshgrid(1:ns);

rows = D(I_d(:), :, :);
rows = rows(:, :, I_s(:));
rows = reshape(rows, d, d, [], ns, ns);

cols = D(J_d(:), :, :);
cols = cols(:, :, J_s(:));
cols = reshape(cols, d, d, [], ns, ns);

rows_rhs = rows(:, 1, :, :, 1);
end
