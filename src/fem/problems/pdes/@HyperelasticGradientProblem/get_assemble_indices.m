function [lhs_rows, lhs_cols, rhs_rows] = get_assemble_indices(obj, mesh, k)
% Assembly indices of the hybrid system: the displacement trace û of
% LinearElasticityProblem (same layout [û_x; û_y], 2*(k+1) dofs per face).
arguments
    obj     (1,1)   HyperelasticGradientProblem
    mesh            SimplicialGrid
    k       (1,1)   double = obj.Order
end
[lhs_rows, lhs_cols, rhs_rows] = get_assemble_indices(LinearElasticityProblem(k), mesh, k);
end
