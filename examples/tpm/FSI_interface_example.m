% Steady fluid-structure interaction across a flat interface, on the
% InterfaceProblem path: each subdomain is an ordinary FEProblem and the
% coupling owns the transmission rows. Fitted meshes with gap h^2; prints
% the velocity L2 error on each side.
%
% This is the no-slip coupling:
%       sigma_F n = sigma_S n + T1      traction balance
%       v         = u_S/dt   + T2       full velocity continuity
% See FSI_navier_slip_example.m for the same problem with the velocity
% continuity relaxed to impermeability plus a friction law.
clear;

k = 2;  h = 1/8;
viscosity = 1.0;  young_modulus = 1.0;  poisson_ratio = 0.25;  dt = 1.0;

interface = FlatCurveHorizontal(0.5);
[fluid_mesh, solid_mesh] = fitted_split(interface, [h, h], h^2);

% manufactured data: exact fields, volume sources, and the interface
% sources T1 (traction jump) and T2 (velocity jump)
[~, ~, interface_sources, pd] = FluidStructureCoupling.load_problems( ...
    k, interface, 0.0, dt, "x^2*y^2", "x^2 + y^2", ["x^2*y"; "x*y^2"], ...
    0.0, viscosity, 0.0, young_modulus, poisson_ratio, 0.0, ...
    pressure_mean = "pinned");

% --- fluid: incompressible elasticity (PEERS stress / velocity) ---------
Pk = DubinerBasis(k);
fluid = FEProblem(Order = k, ...
    Spaces = [FESpace([PEERSBasis(k), Pk.to_vector(), Pk, DubinerBasis(0)], Condense = true), ...
              FESpace([LegendreBasis(k).to_vector(2), DubinerBasis(0)])]);
fluid.Mesh = fluid_mesh;
[lhs, flux, rhs] = IncompressibleElasticityMatrices(k, fluid_mesh, viscosity, ...
    pd.f_F, "pressure_mean", "pinned");
fluid.Operator = AffineOperator(lhs + rhs) + flux;
fluid.BoundaryConditions = StrongEssentialBoundaryConditions(Datum(pd.u_F), ...
    "Dirichlet", projection_order = k);

% --- solid: linear elasticity ------------------------------------------
solid = FEProblem(Order = k, ...
    Spaces = [FESpace([PEERSBasis(k), Pk.to_vector(), Pk], Condense = true), ...
              FESpace(LegendreBasis(k).to_vector(2))]);
solid.Mesh = solid_mesh;
[slhs, sflux, srhs] = LinearElasticityMatrices(k, solid_mesh, young_modulus, poisson_ratio);
solid.Operator = AffineOperator(slhs) + sflux;
if ~isempty(srhs), solid.Operator = solid.Operator + srhs; end
w = TestFunction(Pk);
solid.Operator = solid.Operator + IndexedMatrix(IndexedArray(cat(1, ...
    double(Datum(pd.f_S{1}) * w * Measure(solid_mesh, 'dx')), ...
    double(Datum(pd.f_S{2}) * w * Measure(solid_mesh, 'dx'))), [-1, 2]));
solid.BoundaryConditions = StrongEssentialBoundaryConditions(Datum(pd.u_S), ...
    "Dirichlet", projection_order = k);

% --- couple and solve ---------------------------------------------------
coupled = InterfaceProblem.fsi(fluid, solid, interface, ...
    viscosity = viscosity, young_modulus = young_modulus, ...
    poisson_ratio = poisson_ratio, domain_measure = pd.measure_F, ...
    scheme = "steady", time_step = dt, pressure_mean = "pinned", ...
    source_neumann = interface_sources{1}, ...
    source_dirichlet = interface_sources{2});

[fluid_solution, solid_solution] = solve(coupled, k);

error_fluid = compute_l2_error(fluid_solution.Unknowns.field(2), fluid_mesh, Datum(pd.u_F))
error_solid = compute_l2_error(solid_solution.Unknowns.field(2), solid_mesh, Datum(pd.u_S))
