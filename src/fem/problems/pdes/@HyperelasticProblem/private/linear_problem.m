function lin = linear_problem(obj)
% The LinearElasticityProblem with the same data (the nonlinear keys removed).
lin = LinearElasticityProblem(obj.Order);
nonlinear = ["material", "newton_tolerance", "newton_max_iterations", "iterate", ...
    "stabilization", "characteristic_length", "stabilization_update", "solid_viscosity", "viscous_model"];
lin.Data = rmfield(obj.Data, intersect(nonlinear, fieldnames(obj.Data)));
end
