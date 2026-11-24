% Navier- / full-slip FSI interface law (InterfaceProblem path).
%
% Flat interface, steady closure, k = 2. The manufactured data is the
% no-slip set, so it is not an exact solution of the slip problem; what is
% checked is that the slip rows assemble and solve, and that a -> infinity
% recovers the no-slip solution:
%
%     t.(sigma_F n) + a t.(v - u_S) = 0   -->   t.(v - u_S) = 0
%
% The velocity difference must be taken at Sigma, on the same transferred
% reconstructions the kinematic row uses; taken on the traces on Gamma_h,
% the limit saturates instead of converging.
clc; clear;
k = 2; mu_f = 1.0; E_s = 1.0; nu_s = 0.25; dt = 1.0; h = 1/8;
interface = FlatCurveHorizontal(0.5);
closure = "pinned";
Pk = DubinerBasis(k);
[~, ~, interface_sources, pd] = FluidStructureCoupling.load_problems( ...
    k, interface, 0.0, dt, "x^2*y^2", "x^2 + y^2", ["x^2*y"; "x*y^2"], ...
    0.0, mu_f, 0.0, E_s, nu_s, 0.0, pressure_mean = closure);
[fluid_mesh, solid_mesh] = fitted_split(interface, [h, h], h^2);

fluid_pb = FEProblem(Spaces = [FESpace([PEERSBasis(k), Pk.to_vector(), Pk, ...
    DubinerBasis(0)], Condense=true), ...
    FESpace([LegendreBasis(k).to_vector(2), DubinerBasis(0)])], Order = k);
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

base = {'viscosity', mu_f, 'young_modulus', E_s, 'poisson_ratio', nu_s, ...
    'domain_measure', pd.measure_F, 'scheme', "steady", 'time_step', dt, ...
    'source_neumann', interface_sources{1}, 'source_dirichlet', interface_sources{2}, ...
    'pressure_mean', closure};
solve_law = @(law, a) deal_solution(InterfaceProblem.fsi(fluid_pb, solid_pb, ...
    interface, base{:}, 'interface_law', law, 'slip_coefficient', a), k, ...
    fluid_mesh, pd);

n_fail = 0;
[uF0, uS0, eF0] = solve_law("noslip", 0);
fprintf("noslip          : e(u_F) %.3e\n", eF0);
% The no-slip coupling is exact for this manufactured set, so the control
% sits at discretisation level; a flipped interface normal gives an O(1e-1)
% error here and nowhere else.
n_fail = n_fail + check(eF0 < 1e-4, sprintf("no-slip control e(u_F) = %.3e", eF0));

% full slip must assemble and solve, and differ from no-slip
[uF, uS] = solve_law("full", 0);
d_full = max(max(abs(uF - uF0), [], 'all'), max(abs(uS - uS0), [], 'all'));
fprintf("full slip       : max|u - u_noslip| %.3e\n", d_full);
n_fail = n_fail + check(isfinite(d_full) && d_full > 1e-4, ...
    sprintf("full slip differs from no-slip (%.3e)", d_full));

% a -> infinity recovers no-slip, at the O(1/a) rate
prev = d_full;                 % full slip is the a = 0 end of the sweep
for a = [1e0 1e2 1e4 1e6]
    [uF, uS] = solve_law("navier", a);
    d = max(max(abs(uF - uF0), [], 'all'), max(abs(uS - uS0), [], 'all'));
    fprintf("navier a %8.0e: max|u - u_noslip| %.3e\n", a, d);
    n_fail = n_fail + check(d < prev, sprintf("a = %g decreases the gap", a));
    prev = d;
end
n_fail = n_fail + check(prev < 1e-6, ...
    sprintf("a = 1e6 reaches no-slip (%.3e); saturation means the friction " + ...
            "residual is taken on Gamma_h, not at Sigma", prev));
% the option validation must reject combinations that would be silently ignored
n_fail = n_fail + check(throws(@() solve_law("full", 5), ...
    "InterfaceProblem:slipCoefficientUnused"), ...
    "slip_coefficient with interface_law = ""full"" is rejected");
n_fail = n_fail + check(throws(@() solve_law("navier", 0), ...
    "InterfaceProblem:slipCoefficientZero"), ...
    "interface_law = ""navier"" with a = 0 is rejected");

fprintf("%d check(s) failed.\n", n_fail);
if n_fail > 0
    error("test_interface_fsi_slip:failed", "%d checks failed", n_fail);
end
