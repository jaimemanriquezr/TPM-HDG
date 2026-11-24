% Steady fluid-structure interaction across a flat interface with the
% Navier-slip transmission conditions, on the InterfaceProblem path.
%
% Instead of full velocity continuity, the kinematic coupling keeps only
% impermeability and the freed tangential rows carry a friction law
% (Cesik, Kampschulte & Schwarzacher, arXiv:2603.12030):
%
%       (v - u_S).n = T2.n                        impermeability
%       sigma_F n   = sigma_S n + T1              traction balance, unchanged
%                                                 (the interface force is internal)
%       t.(sigma_F n) + a t.(cdt*v - u_S) = t.g   fluid Navier law
%
% a = slip_coefficient. interface_law = "full" gives a = 0 (frictionless);
% a -> infinity recovers the no-slip coupling. The manufactured data here
% is the no-slip set, so the example prints the distance to the no-slip
% solution for free slip and for increasing a (expected O(1/a)). A
% manufactured solution of the slip problem also needs source_slip = g.
%
% cdt in the friction term is the scheme weight of the kinematic row (dt
% for "steady"), so the units of a depend on the time scheme.
%
% FSI_interface_example.m is this file with the slip options removed.
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
common = {"viscosity", viscosity, "young_modulus", young_modulus, ...
    "poisson_ratio", poisson_ratio, "domain_measure", pd.measure_F, ...
    "scheme", "steady", "time_step", dt, "pressure_mean", "pinned", ...
    "source_neumann", interface_sources{1}, ...
    "source_dirichlet", interface_sources{2}};

solve_with = @(law, a) solve(InterfaceProblem.fsi(fluid, solid, interface, ...
    common{:}, "interface_law", law, "slip_coefficient", a), k);

% the no-slip solution, as the reference the friction limit must reach
reference = solve_with("noslip", 0);
u_noslip = reference.Unknowns.field(2).dof;

% frictionless: the tangential velocities are free to differ
free_slip = solve_with("full", 0);
gap_full_slip = max(abs(free_slip.Unknowns.field(2).dof - u_noslip), [], "all")

% and the friction limit: O(1/a) towards no-slip
fprintf("\n%10s  %s\n", "a", "max|v - v_noslip|");
for a = [1e0 1e2 1e4 1e6]
    slip = solve_with("navier", a);
    fprintf("%10.0e  %.3e\n", a, ...
        max(abs(slip.Unknowns.field(2).dof - u_noslip), [], "all"));
end
