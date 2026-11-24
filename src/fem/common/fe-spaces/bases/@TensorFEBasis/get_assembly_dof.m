function [rows, cols, rows_rhs] = get_assembly_dof(obj, mesh)
%GET_ASSEMBLY_DOF  Assembly index blocks of a tensor basis.
%
%   [rows, cols, rows_rhs] = get_assembly_dof(obj, mesh)
%
%   The scalar basis' blocks (ScalarFEBasis/get_assembly_dof) repeated for
%   every linear map, the c-th copy shifted by c-1 times the scalar
%   dof count on the mesh.
[r, stride, c] = scalar_enumeration(obj, mesh, "local");
n_comp = numel(obj.LinearMaps);
[r_stack, c_stack] = deal(cell(n_comp, 1));
for c0 = 1:n_comp
    r_stack{c0} = r + (c0-1)*stride;
    c_stack{c0} = c + (c0-1)*stride;
end
rows = repmat(cat(1, r_stack{:}), [1, n_comp]);
cols = repmat(cat(2, c_stack{:}), [n_comp, 1]);
rows_rhs = rows(:, 1, :, :, 1);
end

