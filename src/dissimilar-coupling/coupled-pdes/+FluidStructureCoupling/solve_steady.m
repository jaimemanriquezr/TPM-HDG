function [fluid_solution, solid_solution, global_system] = solve_steady(...
    fluid_problem, solid_problem, fluid_mesh, solid_mesh, ...
    interface, source_neumann, source_dirichlet, ...
    options)
% Steady TPM-coupled FSI solve. Thin alias for solve(..., "steady").
arguments
    fluid_problem   IncompressibleElasticityProblem
    solid_problem   LinearElasticityProblem
    fluid_mesh      SimplicialGrid
    solid_mesh      SimplicialGrid
    interface       AbstractCurve
    source_neumann  cell = {@(x,y) 0*x.*y, @(x,y) 0*x.*y}
    source_dirichlet cell = {@(x,y) 0*x.*y, @(x,y) 0*x.*y}
    options.Order double {mustBeScalarOrEmpty} = fluid_problem.Order
    options.n_quadpoints      double {mustBeScalarOrEmpty} = []
    options.n_quadpoints_path double {mustBeScalarOrEmpty} = []
    options.quiet             (1,1) logical = false
    options.patch_weight      (1,1) string {mustBePatchWeight} = "rectangular"
end
opts = namedargs2cell(options);
if nargout > 2
    [fluid_solution, solid_solution, ~, global_system] = FluidStructureCoupling.solve( ...
        fluid_problem, solid_problem, fluid_mesh, solid_mesh, ...
        interface, source_neumann, source_dirichlet, "steady", opts{:});
else
    [fluid_solution, solid_solution] = FluidStructureCoupling.solve( ...
        fluid_problem, solid_problem, fluid_mesh, solid_mesh, ...
        interface, source_neumann, source_dirichlet, "steady", opts{:});
end
end
