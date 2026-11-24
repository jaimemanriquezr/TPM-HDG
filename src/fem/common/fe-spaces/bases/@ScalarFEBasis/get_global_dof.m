function dof_list = get_global_dof(basis, mesh)
%GET_GLOBAL_DOF  Global dof numbering of a scalar basis on a mesh.
%
%   dof_list = get_global_dof(basis, mesh)
%
%   Returns the Dimension-by-n_entities matrix whose column e holds the
%   Dimension consecutive global dof of entity e. Entities are the
%   elements, faces or vertices of mesh, according to the codimension of
%   the basis support (mesh.geometric_dimension - GeometricInputDimension
%   = 0, 1 or 2).
arguments
    basis ScalarFEBasis
    mesh AbstractGrid
end

codimension = mesh.geometric_dimension - basis.GeometricInputDimension;
switch codimension
    case 0
        n_entities = size(mesh.elements, 2);
    case 1
        n_entities = size(mesh.faces, 2);
    case 2
        n_entities = size(mesh.coordinates, 2);
    otherwise
        error("Negative codimension.")
end

num_dof = basis.Dimension;
dof_list = ((1:n_entities) - 1)*num_dof + (1:num_dof).';
end


