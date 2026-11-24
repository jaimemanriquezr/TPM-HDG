function [x_hat, dir_dof, g_hat] = get_boundary_conditions(obj, mesh, gD, gN)
% Dirichlet values and Neumann (reference traction) loads of û, as in
% LinearElasticityProblem.get_boundary_conditions (Data.boundary_conditions
% specs, else gD on the "Dirichlet" tag and gN on the "Neumann" tag).
arguments
    obj     (1,1)   HyperelasticGradientProblem
    mesh            SimplicialGrid
    gD      (2,1)   cell = obj.Data.datum_dirichlet
    gN      (2,1)   cell = obj.Data.datum_neumann
end
lin = LinearElasticityProblem(obj.Order);
if isfield(obj.Data, "boundary_conditions")
    lin.Data.boundary_conditions = obj.Data.boundary_conditions;
end
[x_hat, dir_dof, g_hat] = get_boundary_conditions(lin, mesh, gD, gN);
end
