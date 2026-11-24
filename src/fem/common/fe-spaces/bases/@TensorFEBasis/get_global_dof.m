function dof_list = get_global_dof(obj, mesh)
%GET_GLOBAL_DOF  Global dof numbering of a tensor basis on a mesh.
%
%   dof_list = get_global_dof(obj, mesh)
%
%   The scalar basis' numbering repeated for every linear map, component c
%   shifted by (c-1) times the scalar dof count on the mesh. Column e
%   lists the dof of entity e, component by component.
arguments
    obj TensorFEBasis
    mesh AbstractGrid
end
[single_dof, stride] = scalar_enumeration(obj, mesh, "global");
n_comp = numel(obj.LinearMaps);
dof_shift = permute((0:n_comp-1)*stride, [1 3 2]);
dof_list = single_dof + dof_shift;             % [sdim, n_sub, n_comp]
dof_list = reshape(permute(dof_list, [1 3 2]), [], size(single_dof, 2));
end
