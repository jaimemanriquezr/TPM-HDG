% FluidStructureCoupling.solve must agree dof-for-dof (tol 1e-9) with
% FluidStructureCouplingOracle.solve, the frozen reference implementation,
% from the same low-level problem objects: steady, onestep and newmark on the
% flat fitted problem (newmark with nonzero skeleton and stress histories,
% whose T2 transport the API forms itself and the oracle receives as its
% history_extension); a steady case with componentwise Neumann on Outflow,
% Symmetry and Top plus a Winkler foundation; a convected (Oseen) fluid,
% whose traction balance stays Cauchy (no convective interface flux). The
% oracle must reject a slip law.
clc; clear;
n_fail = 0;  tol = 1e-9;
k = 2;  d1 = k + 1;  d2 = nchoosek(k + 2, 2);  dS = PEERSBasis(k).Dimension;
mu_f = 1.0; E_s = 1.0; nu_s = 0.25; dt = 1e-2;
interface = FlatCurveHorizontal(0.5);
h = 1/8;
[fluid_mesh, solid_mesh] = fitted_split(interface, [h, h], h^2);
wave = @(n, m, c) reshape(sin(c + (1:n*m)), n, m);
nel_F = size(fluid_mesh.elements, 2);  nel_S = size(solid_mesh.elements, 2);
nF_S  = size(solid_mesh.faces, 2);

for scheme = ["steady", "onestep", "newmark"]
    [fp, sp, srcs, pd] = FluidStructureCoupling.load_problems(k, interface, 0.0, dt, ...
        "x^2*y^2", "x^2 + y^2", ["x^2*y"; "x*y^2"], ...
        0.0, mu_f, 0.0, E_s, nu_s, 0.0, pressure_mean = "pinned");
    state = struct.empty;
    switch scheme
        case "onestep"
            fp.Data.density = 1.0;  sp.Data.density = 1.0;  sp.Data.time = 0.0;
            state = struct("time", 0.0, "step", 0, ...
                "fluid", struct("velocity", wave(2*d2, nel_F, 0.1)), ...
                "solid", struct("displacement", wave(2*d2, nel_S, 0.3), ...
                                "velocity", wave(2*d2, nel_S, 0.4), ...
                                "displacement_skeleton", wave(2*d1, nF_S, 0.7)), ...
                "diagnostics", struct.empty);
        case "newmark"
            fp.Data.density = 1.0;  fp.Data.time_scheme = "cranknicolson";
            sp.Data.density = 1.0;  sp.Data.time_scheme = "newmark";
            state = struct("time", 0.0, "step", 0, ...
                "fluid", struct("velocity", wave(2*d2, nel_F, 0.1), ...
                                "acceleration", wave(2*d2, nel_F, 0.2)), ...
                "solid", struct("displacement", wave(2*d2, nel_S, 0.3), ...
                                "velocity", wave(2*d2, nel_S, 0.4), ...
                                "acceleration", wave(2*d2, nel_S, 0.5), ...
                                "displacement_skeleton", wave(2*d1, nF_S, 0.7), ...
                                "velocity_skeleton", wave(2*d1, nF_S, 0.8), ...
                                "acceleration_skeleton", wave(2*d1, nF_S, 0.9), ...
                                "stress", wave(dS, nel_S, 1.2), ...
                                "stress_velocity", wave(dS, nel_S, 1.3), ...
                                "stress_acceleration", wave(dS, nel_S, 1.4), ...
                                "rotation", wave(d2, nel_S, 1.5), ...
                                "rotation_velocity", wave(d2, nel_S, 1.6), ...
                                "rotation_acceleration", wave(d2, nel_S, 1.7)), ...
                "diagnostics", struct.empty);
    end
    [fo, so, opts] = oracle_data(fp, sp, state, scheme, fluid_mesh, solid_mesh, interface, k, dt);
    [fO, sO] = FluidStructureCouplingOracle.solve(fo, so, fluid_mesh, solid_mesh, interface, ...
        srcs{:}, scheme, "quiet", true, opts{:});
    [fI, sI] = FluidStructureCoupling.solve(fp, sp, fluid_mesh, solid_mesh, interface, ...
        srcs{:}, scheme, "quiet", true, state = state);
    d = max([max(abs(fO.Unknowns.field(2).dof - fI.Unknowns.field(2).dof), [], "all"), ...
             max(abs(sO.Unknowns.field(2).dof - sI.Unknowns.field(2).dof), [], "all"), ...
             max(abs(fO.Unknowns.field(1).dof - fI.Unknowns.field(1).dof), [], "all")]);
    fprintf("%-8s max|oracle - api| = %.3e\n", scheme, d);
    if d > tol, fprintf(2, "FAILED: %s oracle and api disagree\n", scheme); n_fail = n_fail + 1; end
end

%% componentwise Neumann + Winkler, steady
b0 = 1e2;
[fp, sp, srcs, pd] = FluidStructureCoupling.load_problems(k, interface, 0.0, dt, ...
    "x^2*y^2", "x^2 + y^2", ["x^2*y"; "x*y^2"], ...
    0.0, mu_f, 0.0, E_s, nu_s, 0.0, pressure_mean = "pinned", swap_domains = true);
% fluid below the interface: Symmetry on y = 0 (n = (0,-1)), Outflow on x = 1
fm = apply_tag(fluid_mesh, "Symmetry", @(X) isapprox(X(2, :), 0));
fm = apply_tag(fm, "Outflow",  @(X) isapprox(X(1, :), 1));
sm = apply_tag(solid_mesh, "Top", @(X) isapprox(X(2, :), 1));
sF = pd.sigma_F;  sS = pd.sigma_S;
fp.Data.boundary_conditions = [ ...
    boundary_condition("Dirichlet", x_data = pd.u_F{1}, y_data = pd.u_F{2});
    boundary_condition("Outflow",  x = "neumann", x_data = sF{1,1}, y_data = pd.u_F{2});
    boundary_condition("Symmetry", x = "neumann", x_data = @(x,y) -sF{1,2}(x,y), y_data = pd.u_F{2})];
sp.Data.hooke_modulus = b0;
sp.Data.source_term = {@(x,y) pd.f_S{1}(x,y) + b0*pd.u_S{1}(x,y); ...
                       @(x,y) pd.f_S{2}(x,y) + b0*pd.u_S{2}(x,y)};
sp.Data.boundary_conditions = [ ...
    boundary_condition("Dirichlet", x_data = pd.u_S{1}, y_data = pd.u_S{2});
    boundary_condition("Top", y = "neumann", x_data = pd.u_S{1}, y_data = sS{2,2})];
[fO, sO] = FluidStructureCouplingOracle.solve(fp, sp, fm, sm, interface, srcs{:}, "steady", "quiet", true);
[fI, sI] = FluidStructureCoupling.solve(fp, sp, fm, sm, interface, srcs{:}, "steady", ...
    "quiet", true);
d = max([max(abs(fO.Unknowns.field(2).dof - fI.Unknowns.field(2).dof), [], "all"), ...
         max(abs(sO.Unknowns.field(2).dof - sI.Unknowns.field(2).dof), [], "all")]);
eF = compute_l2_error(fI.Unknowns.field(2), fm, Datum(pd.u_F));
eS = compute_l2_error(sI.Unknowns.field(2), sm, Datum(pd.u_S));
fprintf("%-8s max|oracle - api| = %.3e   e(u_F) %.3e  e(u_S) %.3e\n", "hemo", d, eF, eS);
if d > tol, fprintf(2, "FAILED: hemodynamics-shaped oracle and api disagree\n"); n_fail = n_fail + 1; end
if max(eF, eS) > 1e-2, fprintf(2, "FAILED: mixed-BC + Winkler case not converged\n"); n_fail = n_fail + 1; end

%% Navier-Stokes fluid: convection in the fluid, Cauchy traction on Sigma
[fp, sp, srcs, pd] = FluidStructureCoupling.load_problems(k, interface, 0.0, dt, ...
    "x^2*y^2", "x^2 + y^2", ["x^2*y"; "x*y^2"], ...
    0.0, mu_f, 0.0, E_s, nu_s, 0.0, pressure_mean = "pinned");
beta = L2proj2D(fluid_mesh, {@(x,y) 1 + 0*x; @(x,y) 0.5*x}, k);      % 2*d2 x n_el
fp.Data.convective_velocity = {beta(1:d2, :); beta(d2+1:end, :)};
[fO, sO] = FluidStructureCouplingOracle.solve(fp, sp, fluid_mesh, solid_mesh, interface, srcs{:}, "steady", "quiet", true);
[fI, sI] = FluidStructureCoupling.solve(fp, sp, fluid_mesh, solid_mesh, interface, srcs{:}, "steady", ...
    "quiet", true);
d = max([max(abs(fO.Unknowns.field(2).dof - fI.Unknowns.field(2).dof), [], "all"), ...
         max(abs(sO.Unknowns.field(2).dof - sI.Unknowns.field(2).dof), [], "all")]);
fprintf("%-8s max|oracle - api| = %.3e\n", "oseen", d);
if d > tol, fprintf(2, "FAILED: convected oracle and api disagree\n"); n_fail = n_fail + 1; end

%% the oracle has no slip law: asking it for one must be a hard error
try
    FluidStructureCouplingOracle.solve(fp, sp, fm, sm, interface, srcs{:}, "steady", ...
        "quiet", true, "interface_law", "full");
    fprintf(2, "FAILED: the oracle accepted a slip law\n"); n_fail = n_fail + 1;
catch err
    % the oracle's argument validation rejects the unknown option; any
    % other error means the call failed for an unrelated reason
    if err.identifier ~= "MATLAB:TooManyInputs" || ~contains(err.message, "interface_law")
        fprintf(2, "FAILED: the slip-law call failed for another reason: %s\n", err.message);
        n_fail = n_fail + 1;
    end
end
fprintf("%d check(s) failed.\n", n_fail);
if n_fail > 0, error("test_interface_oracle_parity:failed", "%d checks failed", n_fail); end

%% ------------------------------------------------------------------------
function [fp, sp, opts] = oracle_data(fp, sp, state, scheme, fluid_mesh, solid_mesh, interface, k, dt)
% The oracle's inputs from the state: the volume histories through
% Data.state (read by the problems' get_matrices), the skeleton histories as
% Data.initial_*_skeleton, and the T2 transport the API forms from
% state.solid (same partition and arithmetic) as its history_extension.
opts = {};
if isempty(state), return; end
[fp.Data.state, sp.Data.state] = deal(state);
for f = ["displacement_skeleton", "velocity_skeleton", "acceleration_skeleton"]
    if isfield(state.solid, f), sp.Data.("initial_" + f) = state.solid.(f); end
end
if scheme == "newmark"
    [b, g] = deal(sp.Data.newmark_beta, sp.Data.newmark_gamma);
    c_v = dt*(1 - b/g);  c_a = dt^2*((1/2 - b) - b*(1 - g)/g);
    s = state.solid;
    [~, part_S] = get_partitions(interface, fluid_mesh, solid_mesh, quadrature_patch = k + 1, quadrature_path = k + 2);
    opts = {"history_extension", trace_extension_coefficients(sp, solid_mesh, part_S, ...
        s.stress + c_v*s.stress_velocity + c_a*s.stress_acceleration, ...
        s.rotation + c_v*s.rotation_velocity + c_a*s.rotation_acceleration, Order = k)};
end
end
