% Stationary FSI: time_step = Inf in InterfaceProblem.fsi (InterfaceProblem path).
%
% scheme = "steady" with a finite time_step is one backward-Euler step from
% rest: the kinematic row is < dt*u_F - u_S, mu_2 > = < dt*T2, mu_2 >, with
% the solid DISPLACEMENT in it. time_step = Inf is the stationary problem,
% the row divided by dt in the limit, < u_F, mu_2 > = < T2, mu_2 >: the
% solid velocity vanishes, so no solid term and no history. On a polygonal
% interface the two differ at O(1), since a tangential solid displacement
% has a component along the facet normals.
%
% Flat interface, k = 2, "constrained" closure, and stationary interface
% data T2 = u_F on Sigma. Checked:
%   - no-slip with time_step = Inf reproduces the manufactured fluid at
%     discretisation level (the stationary data is exact for it);
%   - time_step = Inf equals the dt -> Inf limit (no-slip and full slip);
%   - time_step = Inf is rejected outside scheme = "steady" or without the
%     "constrained" closure (the kinematic rows no longer fix the pressure
%     level).
clc; clear;
k = 2; mu_f = 1.0; E_s = 1.0; nu_s = 0.25; h = 1/8;
interface = FlatCurveHorizontal(0.5);
closure = "constrained";
Pk = DubinerBasis(k);
[~, ~, interface_sources, pd] = FluidStructureCoupling.load_problems( ...
    k, interface, 0.0, 1.0, "x^2*y^2", "x^2 + y^2", ["x^2*y"; "x*y^2"], ...
    0.0, mu_f, 0.0, E_s, nu_s, 0.0, pressure_mean = closure);
[fluid_mesh, solid_mesh] = fitted_split(interface, [h, h], h^2);

fluid_pb = FEProblem(Spaces = [FESpace([PEERSBasis(k), Pk.to_vector(), Pk, ...
    DubinerBasis(0)], Condense=true), ...
    FESpace([LegendreBasis(k).to_vector(2), DubinerBasis(0), ConstantBasis()])], Order = k);
fluid_pb.Mesh = fluid_mesh;
[lhs, flux, rhs] = IncompressibleElasticityMatrices(k, fluid_mesh, mu_f, pd.f_F, ...
    "pressure_mean", closure);
fluid_pb.Operator = AffineOperator(lhs + rhs) + flux;
fluid_pb.BoundaryConditions = StrongEssentialBoundaryConditions(Datum(pd.u_F), ...
    "Dirichlet", projection_order = k);

solid_pb = FEProblem(Spaces = [FESpace([PEERSBasis(k), Pk.to_vector(), Pk], ...
    Condense=true), FESpace(LegendreBasis(k).to_vector(2))], Order = k);
solid_pb.Mesh = solid_mesh;
[slhs, sflux, srhs] = LinearElasticityMatrices(k, solid_mesh, E_s, nu_s);
solid_pb.Operator = AffineOperator(slhs) + sflux;
if ~isempty(srhs); solid_pb.Operator = solid_pb.Operator + srhs; end
w = TestFunction(Pk);
F1 = double(Datum(pd.f_S{1}) * w * Measure(solid_mesh, 'dx'));
F2 = double(Datum(pd.f_S{2}) * w * Measure(solid_mesh, 'dx'));
solid_pb.Operator = solid_pb.Operator + IndexedMatrix(IndexedArray(cat(1, F1, F2), [-1, 2]));
solid_pb.BoundaryConditions = StrongEssentialBoundaryConditions(Datum(pd.u_S), ...
    "Dirichlet", projection_order = k);

% stationary kinematic data: the fluid velocity itself on Sigma
T2 = {pd.u_F{1}, pd.u_F{2}};
base = {'viscosity', mu_f, 'young_modulus', E_s, 'poisson_ratio', nu_s, ...
    'domain_measure', pd.measure_F, 'scheme', "steady", ...
    'source_neumann', interface_sources{1}, 'source_dirichlet', T2, ...
    'pressure_mean', closure};
solve_with = @(dt, law, extra) deal_solution(InterfaceProblem.fsi(fluid_pb, solid_pb, ...
    interface, base{:}, 'time_step', dt, 'interface_law', law, extra{:}), k, fluid_mesh, pd);

n_fail = 0;
[~, ~, eF] = solve_with(Inf, "noslip", {});
fprintf("noslip, time_step = Inf: e(u_F) %.3e\n", eF);
n_fail = n_fail + check(eF < 1e-4, sprintf("stationary no-slip e(u_F) = %.3e", eF));

for law = ["noslip", "full"]
    [uF_inf, uS_inf] = solve_with(Inf, law, {});
    [uF_big, uS_big] = solve_with(1e8, law, {});
    scale = max(abs([uF_inf(:); uS_inf(:)]));
    d = max(max(abs(uF_inf - uF_big), [], 'all'), max(abs(uS_inf - uS_big), [], 'all')) / scale;
    fprintf("%-6s: |Inf - 1e8| %.2e (relative)\n", law, d);
    n_fail = n_fail + check(d < 1e-6, sprintf("%s: time_step = Inf is the dt -> Inf limit (%.2e)", law, d));
end

% time_step = Inf is the stationary problem only
n_fail = n_fail + check(throws(@() solve_with(Inf, "noslip", {'scheme', "newmark"}), ...
    "InterfaceProblem:infiniteTimeStep"), "time_step = Inf with scheme = ""newmark"" is rejected");
n_fail = n_fail + check(throws(@() solve_with(Inf, "noslip", {'pressure_mean', "pinned"}), ...
    "InterfaceProblem:infiniteTimeStep"), "time_step = Inf with pressure_mean = ""pinned"" is rejected");

fprintf("%d check(s) failed.\n", n_fail);
if n_fail > 0
    error("test_interface_fsi_stationary:failed", "%d checks failed", n_fail);
end
