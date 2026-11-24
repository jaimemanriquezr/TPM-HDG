function results = FSI_turek_hron_benchmark(options)
% The Turek-Hron flag benchmark (Turek and Hron, "Proposal for numerical
% benchmarking of fluid-structure interaction between an elastic object
% and laminar incompressible flow", in Bungartz and Schaefer (eds.),
% Fluid-Structure Interaction, LNCSE 53, Springer 2006, pp. 371-385).
%
% Geometry: channel (0, 2.5) x (0, 0.41); rigid cylinder of radius 0.05
% centred at (0.2, 0.2); elastic flag [., 0.6] x [0.19, 0.21] clamped on
% the cylinder; control point A = (0.6, 0.2) at the flag tip. Fluid
% rho_f = 1000, nu_f = 1e-3; parabolic inflow of mean U_bar,
%     u_x(0, y) = 1.5 U_bar y (H - y)/(H/2)^2,
% ramped by (1 - cos(pi t/2))/2 over t < 2; traction-free outlet; no-slip
% walls and cylinder.
%
%   CSM1, CSM2  the flag alone under gravity g = (0, -2), static
%   CSM3        the flag alone under gravity, released from rest
%   CFD1-CFD3   the fluid alone, the flag rigid (Re = 20, 100, 200)
%   FSI1        Re = 20, steady state
%   FSI2, FSI3  Re = 100, 200, periodic self-induced oscillation
%
% T = 0 (the default for CSM1, CSM2, CFD1, CFD2, FSI1) solves for the
% steady state directly: Picard iteration on the coupled "steady" scheme,
% Stokes first, then the convective velocity and the mass-term history
% both set to the previous iterate, so that at convergence the backward-
% Euler mass term vanishes and the fixed point is the steady
% Navier-Stokes/elasticity solution. The kinematic row of the steady
% scheme, steady_time_step*u_F = u_S on Sigma, is the no-slip condition up
% to |u_S|/steady_time_step. T > 0 marches in time instead.
%
% Discretization: stress-velocity PEERS HDG fluid (Navier-Stokes, the
% convective velocity lagged one step), PEERS HDG linear elasticity
% solid, coupled across the flag's three free sides by the transfer path
% method (the interface is the open polygon of the flag, both meshes
% fitted to it). Time: Newmark(1/4, 1/2) solid; Crank-Nicolson fluid
% after backward-Euler startup steps for the periodic tests, backward
% Euler throughout for the steady ones (CN does not damp the march to a
% steady state: at large dt it rings odd/even about it).
%
% MODEL DIFFERENCES FROM THE BENCHMARK. The benchmark solid is St.
% Venant-Kirchhoff and the fluid domain follows the flag (ALE). Here the
% solid is linear elastic and both domains stay fixed: FSI1 and CSM2
% (small deflections) are close to the benchmark, while CSM1, CSM3, FSI2
% and FSI3 (tip deflections of 10-25% of the flag length) measure the
% linear, fixed-domain model against the nonlinear reference.
%
% Output: per-step tip displacement u(A), drag and lift (force of the
% fluid on cylinder + flag), saved to <out_dir>/<test>_k<k>_h<h>.mat (the
% default out_dir, ./results/turek-hron/, is relative to the current folder);
% steady tests report the final values, unsteady ones
% mean +- amplitude [frequency] over the last stats_window seconds, both
% beside the reference values.
%
%Usage:
%   FSI_turek_hron_benchmark(test = "FSI1")                 % steady
%   FSI_turek_hron_benchmark(test = "FSI3", k = 2, h = 0.005, dt = 5e-4, T = 10)
arguments
    options.test (1,1) string {mustBeMember(options.test, ...
        ["CSM1", "CSM2", "CSM3", "CFD1", "CFD2", "CFD3", "FSI1", "FSI2", "FSI3"])} = "FSI1"
    options.k     (1,1) double = 2
    options.h     (1,1) double = 0.01      % element size at the body
    options.h_far (1,1) double = 0.05      % element size away from the body
    options.dt    double {mustBeScalarOrEmpty} = []   % empty: per-test default
    options.T     double {mustBeScalarOrEmpty} = []   % empty: per-test default
    options.stats_window double {mustBeScalarOrEmpty} = []  % empty: min(2, T/3)
    options.rannacher_steps (1,1) double = 2  % backward-Euler startup steps (periodic tests)
    options.steady_time_step (1,1) double = 1e4
    options.picard_tolerance (1,1) double = 1e-8   % relative increment of u_h
    options.picard_max_iterations (1,1) double = 50
    options.out_dir (1,1) string = "./results/turek-hron/"   % relative to the current folder
end
cfg = test_parameters(options.test);
k = options.k;
dt = options.dt;  if isempty(dt), dt = cfg.dt; end
T  = options.T;   if isempty(T),  T  = cfg.T;  end
window = options.stats_window;  if isempty(window), window = min(2, T/3); end
is_steady = T == 0;
if is_steady && cfg.kind ~= "csm", dt = options.steady_time_step; end

roles = dictionary(["csm", "cfd", "fsi"], ["free", "rigid", "coupled"]);
[fluid_mesh, solid_mesh, flag, geo] = turek_hron_meshes(options.h, options.h_far, ...
    flag = roles(cfg.kind));
fprintf("== %s: k %d, h %g (far %g), dt %g, T %g; %d fluid / %d solid elements\n", ...
    options.test, k, options.h, options.h_far, dt, T, ...
    size(fluid_mesh.elements, 2), size(solid_mesh.elements, 2));

%% Problems -----------------------------------------------------------------
rho_f = 1000;  mu_f = 1;                         % nu_f = 1e-3
E_s = 2*cfg.mu_s*(1 + cfg.nu_s);
gravity = {@(x, y) 0*x.*y; @(x, y) -cfg.gravity*cfg.rho_s + 0*x.*y};

if cfg.kind ~= "cfd"
    solid_problem = LinearElasticityProblem(k, ...
        density = cfg.rho_s*~is_steady, young_modulus = E_s, poisson_ratio = cfg.nu_s, ...
        time_step = dt, time_scheme = "newmark", source_term = gravity);
    solid_problem.Data.boundary_conditions = [ ...
        boundary_condition("Dirichlet");                                 % clamped on the cylinder
        boundary_condition("Neumann", x = "neumann", y = "neumann")];   % traction-free (CSM)
end
if cfg.kind ~= "csm"
    fluid_problem = IncompressibleElasticityProblem(k, ...
        density = rho_f, viscosity = mu_f, time_step = dt, ...
        domain_measure = geo.measure_F, pressure_mean = "pinned", ...
        convective_stabilization = false);
    % one tau on the whole mesh, shared by the local problems and the
    % transmission rows: 1 + rho max|u|/(2 mu), the flow reaching about
    % twice the peak inflow speed past the body
    fluid_mesh = add_stabilization(fluid_mesh, 1 + rho_f*2*1.5*cfg.U_bar/(2*mu_f));
    ramp = @(t) (t < 2)*(1 - cos(pi*t/2))/2 + (t >= 2);
    H = geo.box(4);
    fluid_bcs = @(t) [ ...
        boundary_condition("Dirichlet");                                 % channel walls
        boundary_condition("Body");                                      % cylinder (+ rigid flag)
        boundary_condition("Inflow", x_data = @(x, y) ...
            ramp(t)*1.5*cfg.U_bar*y.*(H - y)/(H/2)^2 + 0*x);
        boundary_condition("Outflow", x = "neumann", y = "neumann")];
end
if cfg.kind == "fsi"
    zero_source = {@(x, y) 0*x.*y, @(x, y) 0*x.*y};   % physical transmission: no interface sources
    [partition_fluid, ~] = get_partitions(flag, fluid_mesh, solid_mesh, ...
        quadrature_patch = k + 1, quadrature_path = k + 2);
end

%% Time marching, or Picard iteration to the steady state ---------------------
% steady: each "step" is one Picard iterate at the full inflow (t = 2, the
% end of the ramp), t(n) its index
if is_steady && cfg.kind == "csm"
    n_steps = 1;
elseif is_steady
    n_steps = options.picard_max_iterations;
else
    n_steps = round(T/dt);
end
[t, ux, uy, drag, lift] = deal(nan(1, n_steps));
% the state at t = 0 (rest; the solid acceleration balances the load),
% advanced by every step; steady iterations keep only the fluid velocity,
% the pseudo-time history of the next Picard iterate
switch cfg.kind
    case "csm"
        state = struct.empty;
        if ~is_steady, state = initial_state(solid_problem, solid_mesh, k); end
    case "cfd"
        fluid_problem.Data.boundary_conditions = fluid_bcs(0);
        state = initial_state(fluid_problem, fluid_mesh, k);
    case "fsi"
        state = struct("fluid", struct("velocity", 0));
        if ~is_steady
            fluid_problem.Data.boundary_conditions = fluid_bcs(0);
            state = FluidStructureCoupling.initial_state(fluid_problem, solid_problem, ...
                fluid_mesh, solid_mesh, flag, zero_source, zero_source);
        end
end
t_start = tic;
for n = 1:n_steps
    t(n) = n*dt;
    t_bc = t(n);
    if is_steady, t(n) = n;  t_bc = 2; end
    if cfg.kind ~= "csm"
        u_previous = state.fluid.velocity;
        fluid_problem.Data.convective_velocity = {};
        if any(u_previous, "all")
            d2 = size(u_previous, 1)/2;
            fluid_problem.Data.convective_velocity = ...
                {u_previous(1:d2, :); u_previous(d2+1:end, :)};
        end
        if n <= options.rannacher_steps || ~cfg.periodic || is_steady
            fluid_problem.Data.time_scheme = "euler";
        else
            fluid_problem.Data.time_scheme = "cranknicolson";
        end
        fluid_problem.Data.boundary_conditions = fluid_bcs(t_bc);
    end

    switch cfg.kind
        case "csm"
            solid_problem.Data.state = state;
            solid_solution = solve(solid_problem, solid_mesh, k);
            if ~is_steady
                state = get_kinematic_state(solid_problem, solid_solution, state);
            end
        case "cfd"
            fluid_problem.Data.state = state;
            fluid_solution = solve(fluid_problem, fluid_mesh, k);
            state = get_kinematic_state(fluid_problem, fluid_solution, state);
        case "fsi"
            if is_steady
                [fluid_solution, solid_solution] = FluidStructureCoupling.solve( ...
                    fluid_problem, solid_problem, fluid_mesh, solid_mesh, flag, ...
                    zero_source, zero_source, "steady", quiet = true, state = state);
                state.fluid.velocity = fluid_solution.Unknowns.field(2).dof;
            else
                [fluid_solution, solid_solution, state] = FluidStructureCoupling.solve( ...
                    fluid_problem, solid_problem, fluid_mesh, solid_mesh, flag, ...
                    zero_source, zero_source, "newmark", quiet = true, state = state);
            end
    end

    if cfg.kind ~= "cfd"
        u_A = point_displacement(solid_solution, solid_mesh, geo.A);
        ux(n) = u_A(1);  uy(n) = u_A(2);
    end
    if cfg.kind ~= "csm"
        % resultant of sigma n_F over the fluid's no-slip body faces and,
        % with an elastic flag, over the interface; the force on the body
        % is its negative
        R = boundary_resultant(fluid_solution, fluid_mesh, "Body", k);
        if cfg.kind == "fsi"
            [x_i, t_i] = transfer_traction_to_interface(fluid_solution, fluid_mesh, ...
                partition_fluid, Order = k, normals = -partition_fluid.quadrature.normals);
            R = R + interface_loop_resultants(t_i, x_i, partition_fluid, [0; 0]);
        end
        drag(n) = -R(1);  lift(n) = -R(2);
    end
    if is_steady && cfg.kind ~= "csm"
        u = state.fluid.velocity;
        increment = norm(u - u_previous, "fro")/norm(u, "fro");
        fprintf("  Picard %2d: increment %.3e  ux %+.6e  uy %+.6e  drag %.6f  lift %+.6f  (%.0f s)\n", ...
            n, increment, ux(n), uy(n), drag(n), lift(n), toc(t_start));
        if increment < options.picard_tolerance, break, end
    elseif mod(n, 50) == 0 || n == n_steps
        fprintf("  step %d/%d  t = %.4f  ux %+.4e  uy %+.4e  drag %.4f  lift %+.4f  (%.0f s)\n", ...
            n, n_steps, t(n), ux(n), uy(n), drag(n), lift(n), toc(t_start));
    end
end
if is_steady && cfg.kind ~= "csm"
    keep = 1:n;
    [t, ux, uy, drag, lift] = deal(t(keep), ux(keep), uy(keep), drag(keep), lift(keep));
    if increment >= options.picard_tolerance
        warning("FSI_turek_hron_benchmark:picard", ...
            "Picard stopped at increment %.3e after %d iterations", increment, n);
    end
end

%% Report -------------------------------------------------------------------------
names = ["ux", "uy", "drag", "lift"];
series = {ux, uy, drag, lift};
results = struct("test", options.test, "k", k, "h", options.h, "h_far", options.h_far, ...
    "dt", dt, "T", T, "t", t, "ux", ux, "uy", uy, "drag", drag, "lift", lift, ...
    "n_fluid_elements", size(fluid_mesh.elements, 2), ...
    "n_solid_elements", size(solid_mesh.elements, 2), "seconds", toc(t_start));
fprintf("\n%-5s %26s %26s\n", "", "computed", "reference");
for q = 1:numel(names)
    if all(isnan(series{q})), continue, end
    ref = cfg.reference.(names(q));
    if cfg.periodic
        [m, a, fr] = periodic_statistics(t, series{q}, window);
        results.(names(q) + "_stats") = [m, a, fr];
        fprintf("%-5s %26s %26s\n", names(q), ...
            sprintf("%.4g +- %.4g [%.4g]", m, a, fr), ...
            sprintf("%.4g +- %.4g [%.4g]", ref(1), ref(2), ref(3)));
    else
        results.(names(q) + "_final") = series{q}(end);
        fprintf("%-5s %26.5g %26.5g\n", names(q), series{q}(end), ref);
    end
end
if ~isfolder(options.out_dir), mkdir(options.out_dir); end
file = fullfile(options.out_dir, sprintf("%s_k%d_h%g.mat", options.test, k, options.h));
save(file, "-struct", "results");
fprintf("wrote %s\n", file);
end

%=========================================================================%
function cfg = test_parameters(test)
% Material and flow parameters, default time stepping, and the reference
% values (Turek and Hron 2006, finest level): steady tests a scalar,
% periodic ones [mean, amplitude, frequency]; displacements in m.
csm = struct("kind", "csm", "U_bar", 0, "rho_s", 1000, "nu_s", 0.4, "gravity", 2);
cfd = struct("kind", "cfd", "rho_s", 0, "nu_s", 0, "mu_s", 0, "gravity", 0);
fsi = struct("kind", "fsi", "nu_s", 0.4, "gravity", 0);
switch test
    case "CSM1"
        cfg = setfield(csm, "mu_s", 0.5e6);
        cfg.dt = 0;  cfg.T = 0;  cfg.periodic = false;
        cfg.reference = struct("ux", -7.187e-3, "uy", -6.610e-2);
    case "CSM2"
        cfg = setfield(csm, "mu_s", 2.0e6);
        cfg.dt = 0;  cfg.T = 0;  cfg.periodic = false;
        cfg.reference = struct("ux", -4.690e-4, "uy", -1.697e-2);
    case "CSM3"
        cfg = setfield(csm, "mu_s", 0.5e6);
        cfg.dt = 5e-3;  cfg.T = 10;  cfg.periodic = true;
        cfg.reference = struct("ux", [-14.305e-3, 14.305e-3, 1.0995], ...
                               "uy", [-63.607e-3, 65.160e-3, 1.0995]);
    case "CFD1"
        cfg = setfield(cfd, "U_bar", 0.2);
        cfg.dt = 0.05;  cfg.T = 0;   cfg.periodic = false;
        cfg.reference = struct("drag", 14.29, "lift", 1.119);
    case "CFD2"
        cfg = setfield(cfd, "U_bar", 1);
        cfg.dt = 0.01;  cfg.T = 0;   cfg.periodic = false;
        cfg.reference = struct("drag", 136.7, "lift", 10.53);
    case "CFD3"
        cfg = setfield(cfd, "U_bar", 2);
        cfg.dt = 5e-3;  cfg.T = 10;  cfg.periodic = true;
        cfg.reference = struct("drag", [439.45, 5.6183, 4.3956], ...
                               "lift", [-11.893, 437.81, 4.3956]);
    case "FSI1"
        cfg = fsi;  cfg.U_bar = 0.2;  cfg.rho_s = 1000;  cfg.mu_s = 0.5e6;
        cfg.dt = 0.05;  cfg.T = 0;   cfg.periodic = false;
        cfg.reference = struct("ux", 2.2700e-5, "uy", 8.2209e-4, ...
                               "drag", 14.295, "lift", 0.7638);
    case "FSI2"
        cfg = fsi;  cfg.U_bar = 1;  cfg.rho_s = 10000;  cfg.mu_s = 0.5e6;
        cfg.dt = 1e-3;  cfg.T = 15;  cfg.periodic = true;
        cfg.reference = struct("ux", [-14.58e-3, 12.44e-3, 3.8], ...
                               "uy", [1.23e-3, 80.6e-3, 2.0], ...
                               "drag", [208.83, 73.75, 3.8], ...
                               "lift", [0.88, 234.2, 2.0]);
    case "FSI3"
        cfg = fsi;  cfg.U_bar = 2;  cfg.rho_s = 1000;  cfg.mu_s = 2.0e6;
        cfg.dt = 5e-4;  cfg.T = 10;  cfg.periodic = true;
        cfg.reference = struct("ux", [-2.69e-3, 2.53e-3, 10.9], ...
                               "uy", [1.48e-3, 34.38e-3, 5.3], ...
                               "drag", [457.3, 22.66, 10.9], ...
                               "lift", [2.22, 149.78, 5.3]);
end
end

%=========================================================================%
function [fluid_mesh, solid_mesh, flag, geo] = turek_hron_meshes(h, h_far, options)
% Both meshes fitted to the flag. The fluid is the channel minus the
% cylinder and the flag, graded from h at the body to h_far; the solid is
% the flag, whose left side is the arc of the cylinder it is clamped to.
% Tags: fluid "Dirichlet" (walls), "Inflow", "Outflow", "Body" (cylinder,
% and the flag when it is "rigid"), "Interface" (the flag when
% "coupled"); solid "Dirichlet" (arc), and on the free sides "Interface"
% ("coupled") or "Neumann" ("free", the structure-only tests).
arguments
    h (1,1) double
    h_far (1,1) double
    options.flag (1,1) string {mustBeMember(options.flag, ["coupled", "rigid", "free"])} = "coupled"
end
c = [0.2; 0.2];  r = 0.05;
box = [0 2.5 0 0.41];
A = [0.6; 0.2];
x_c = c(1) + sqrt(r^2 - 0.01^2);              % where the flag meets the cylinder
theta_0 = asin(0.01/r);
% Sigma: the flag's free sides, counter-clockwise, so that the normal
% (tangent rotated by -90 degrees) points out of the flag, into the fluid
flag = PolygonalCurve([[x_c; 0.19], [0.6; 0.19], A, [0.6; 0.21], [x_c; 0.21]], closed = false);

DIRICHLET = 1;  NEUMANN = 2;  INTERFACE = 3;  INFLOW = 4;  BODY = 5;  OUTFLOW = 6;

% fluid: the box counter-clockwise, then the hole (cylinder arc, flag)
n_arc = max(16, ceil(2*pi*r/(h/2)));         % arc points dense enough that mesh2d keeps them
th = linspace(theta_0, 2*pi - theta_0, n_arc);
arc = c + r*[cos(th); sin(th)];
tip = flag.vertices(:, 2:4);
V = [[box(1) box(2) box(2) box(1); box(3) box(3) box(4) box(4)], arc, tip];
i_arc = 4 + (1:n_arc);  i_tip = 4 + n_arc + (1:3);
flag_marker = INTERFACE;
if options.flag == "rigid", flag_marker = BODY; end
C = [[1; 2; DIRICHLET], [2; 3; OUTFLOW], [3; 4; DIRICHLET], [4; 1; INFLOW], ...
     [i_arc(1:end-1); i_arc(2:end); BODY*ones(1, n_arc - 1)], ...
     [[i_arc(end), i_tip]; [i_tip, i_arc(1)]; flag_marker*ones(1, 4)]];
% distance to the box [0.15, 0.6] x [0.15, 0.25] around the body
dist = @(P) hypot(max(0, max(0.15 - P(:, 1), P(:, 1) - 0.6)), max(0, abs(P(:, 2) - 0.2) - 0.05));
fluid_mesh = SimplicialGrid(V, C, mesh_generator = "mesh2d", element_diameter = h, ...
    size_function = @(P) min(h_far, h + 0.25*dist(P)));
fluid_mesh.boundary_tags("Body") = BODY;
fluid_mesh.boundary_tags("Outflow") = OUTFLOW;

% solid: the flag counter-clockwise, closed by the arc from (x_c, 0.21)
% back to (x_c, 0.19)
free_marker = INTERFACE;
if options.flag == "free", free_marker = NEUMANN; end
n_clamp = max(3, ceil(2*r*theta_0/(h/2)));
th = linspace(theta_0, -theta_0, n_clamp + 1);
Vs = [flag.vertices, c + r*[cos(th(2:end-1)); sin(th(2:end-1))]];
n = size(Vs, 2);
Cs = [[1:4; 2:5; free_marker*ones(1, 4)], [5:n; [6:n, 1]; DIRICHLET*ones(1, n - 4)]];
solid_mesh = SimplicialGrid(Vs, Cs, mesh_generator = "mesh2d", element_diameter = h/2);

geo = struct("box", box, "A", A, "measure_F", sum(fluid_mesh.volume));
end

%=========================================================================%
function u = point_displacement(solution, mesh, X)
% The displacement at a mesh vertex X, averaged over the elements
% sharing it (the HDG field is discontinuous).
v = find(vecnorm(mesh.coordinates - X, 2, 1) < 1e-12, 1);
assert(~isempty(v), "FSI_turek_hron_benchmark:noVertex", "no mesh vertex at A");
[i, K] = find(mesh.elements(1:3, :) == v);
U = solution.Unknowns.field(2).evaluate(reference_element(2), coordinate_system = "vertex");
u = zeros(2, 1);
for j = 1:numel(K)
    u = u + reshape(U(i(j), 1, :, 1, K(j)), 2, 1)/numel(K);
end
end

%=========================================================================%
function R = boundary_resultant(solution, mesh, tag, k)
% int sigma_h n ds over the boundary faces of one tag, n the outward
% normal of the fluid, by Gauss-Legendre on each face.
q = k + 2;
[s, w] = lgwt(q, 0, 1);
ref = reference_element(2);
% the stress at the q points of each reference edge (i, j), all elements
edges = [1 2; 2 3; 1 3];
pts = zeros(2, 3*q);
for e = 1:3
    a = ref(:, edges(e, 1));  b = ref(:, edges(e, 2));
    pts(:, (e - 1)*q + (1:q)) = a + (b - a)*s.';
end
S = solution.Unknowns.field(1).evaluate(pts, coordinate_system = "local");   % points x 1 x 2 x 2 x elements
R = zeros(2, 1);
for f = get_face_by_tag(mesh, tag)
    ab = mesh.faces(1:2, f);
    [~, K] = find(any(mesh.elements(1:3, :) == ab(1)) & any(mesh.elements(1:3, :) == ab(2)));
    K = K(1);
    loc = [find(mesh.elements(1:3, K) == ab(1)), find(mesh.elements(1:3, K) == ab(2))];
    e = find(all(sort(loc) == edges, 2));
    V = mesh.coordinates(:, mesh.elements(1:3, K));
    a = V(:, edges(e, 1));  b = V(:, edges(e, 2));
    c = V(:, setdiff(1:3, edges(e, :)));
    n = [b(2) - a(2); a(1) - b(1)];             % |n| = face length
    if dot(n, a - c) < 0, n = -n; end           % away from the element
    for p = 1:q
        sigma = reshape(S((e - 1)*q + p, 1, :, :, K), 2, 2);
        R = R + w(p)*sigma*n;
    end
end
end

%=========================================================================%
function [m, a, f] = periodic_statistics(t, x, window)
% mean = (max + min)/2, amplitude = (max - min)/2 over the last window,
% frequency from the up-crossings of the mean
keep = t >= t(end) - window & ~isnan(x);
t = t(keep);  x = x(keep);
m = (max(x) + min(x))/2;  a = (max(x) - min(x))/2;
up = find(x(1:end-1) < m & x(2:end) >= m);
f = NaN;
if numel(up) >= 2
    tc = t(up) + (m - x(up))./(x(up + 1) - x(up)).*(t(up + 1) - t(up));
    f = (numel(tc) - 1)/(tc(end) - tc(1));
end
end
