function out = lagrangian_carved_march(k, geometry, background, solid_mesh, amplitude, dt, n_steps, options)
%LAGRANGIAN_CARVED_MARCH  FSI march with the interface following a Lagrangian solid.
%
%   out = lagrangian_carved_march(k, geometry, background, solid_mesh, A, dt, n_steps, ...
%       geometry_iterations = 0, snap = false, verbose = false, ...
%       solid_model = "linear", material = "svk")
%
%   Solid: LinearElasticityProblem, Newmark, on its fixed reference mesh
%   (solid_mesh, fitted to Sigma_0). Fluid: Eulerian, carved from the fixed
%   background grid at every step (SimplicialGrid.split, cut elements
%   discarded), backward Euler on the first step and, after it, BDF2
%   (fluid_scheme = "bdf2", the default: the variable-step BDF2 written as
%   a theta step, L-stable, so the pressure mode excited when elements
%   switch on or off is damped) or Crank-Nicolson (fluid_scheme =
%   "cranknicolson");
%   its history is remapped as in eulerian_carved_march, elements being
%   identified by their background element (a snapped element keeps its
%   polynomial, re-expressed on its new shape; a switched-on one gets the
%   nearest old element's polynomial, projected). Coupling: the curve pair
%   {Sigma_h^{n+1}, Sigma_0} of InterfaceProblem.fsi, data from
%   lagrangian_march_problems on the discrete curve.
%
%   Geometry: before step n+1, Sigma_h^{n+1} = X + d* with the velocity
%   predictor of the discrete solid displacement on Sigma_0,
%     d* = u^n + dt v^n  (predictor = "velocity", the default; "newmark"
%     adds (dt^2/4) a^n and "taylor" (dt^2/2) a^n, see predictor_weight),
%   evaluated by solid_interface_displacement (the TPM extrapolation of
%   the solid polynomials). The fluid is carved once per step, against
%   this predictor (snapped first when snap = true, snap_to_curve).
%   geometry_iterations = m adds m fixed-point corrections: Sigma_h^{n+1}
%   is rebuilt from the solved u^{n+1} and the step re-solved, on the same
%   carved fluid mesh. A Sigma_h built from the degree-k u_h has its
%   normal and Nanson factor only O(h^k) accurate, which caps the coupled
%   errors at O(h^k); interface_geometry = "postprocessed" (the default)
%   builds it instead from the P_k+1 postprocessing u* of u_h with its
%   reconstructed gradient (postprocessed_displacement: from the stress
%   and rotation for the linear solid, from the discrete gradient L_h of
%   (u_h, û_h) for the hyperelastic one; of order k + 2 with a gradient of order k + 1 when
%   that gradient is of order k + 1), the predictor and the geometry
%   iterations applied to the fields it reads alike. Built from the
%   element-wise solid field, Sigma_h jumps at the solid's interface
%   vertices by the jump of that field; a jump can hide a stretch of Sigma
%   from the fluid's closest-point map (get_partitions:unequalRange) or
%   pair a fluid face with a far stretch (a long transfer path).
%   continuous_interface = true (the default) removes the jumps by a
%   linear correction of each piece (solid_interface_displacement,
%   continuous), of the size of the jumps; false keeps them.
%
%   solid_model = "hyperelastic": the solid is HyperelasticProblem
%   (material "svk" | "neo_hookean", lagrangian_march_problems) with the
%   constant stabilization tau = max(30 mu, 2.4 lambda)/h on solid_mesh,
%   lambda and mu the solid's Lame parameters and h = sqrt(2 max|K|) (the
%   leg of a structured right triangle), or tau = solid_tau mu/h when
%   solid_tau is given. A smaller tau makes the tangent of the compressed
%   solid indefinite in mesh-scale modes before the physical surface
%   instability (stretch about 0.5): on a laterally compressed neo-Hookean
%   strip, k = 1 and 2, nu = 0.3 to 0.495, the smallest tau that keeps
%   the smooth mode first is about max(20 mu, 1.6 lambda)/h, independent
%   of the shear modulus; the default is 1.5 times that.
%   Each step is then a coupled
%   Newton loop: Data.iterate starts at the state (u^n, u_hat^n), each
%   FluidStructureCoupling.solve assembles the Newton step of the solid at
%   the iterate (HyperelasticProblem.newton_blocks) within the coupled
%   system, and the solved solid (u, u_hat) is the next iterate. The
%   solid's constitutive and symmetry rows are the only nonlinear rows of
%   the coupled system, so the residual is theirs (newton_residual,
%   local_rows = true), relative to the linear residual of all the solid's
%   local rows at the zero iterate (its source and inertial history; when
%   both vanish, as on the first step of physical data from rest, relative
%   to the constitutive and symmetry rows' terms A1 X of the first
%   iterate, the linear residual of that iterate with u_hat = 0); the loop
%   stops at the solid's newton_tolerance, or on the round-off floor of
%   HyperelasticProblem.solve (a residual below 1e-6 that no longer
%   halves), warning after newton_max_iterations. Geometry: the predictor
%   curve is built before the loop, and with geometry_iterations = m the
%   curve is rebuilt from the iterate before each of the Newton solves 2
%   to m + 1 (lagged by one iterate), on the same carved fluid mesh; the
%   loop runs at least m + 1 solves. out.newton(n) holds iterations,
%   residual, converged, stop_reason and the step lengths alpha of step n.
%   solid_viscosity = mu_v
%   gives the hyperelastic solid the Kelvin-Voigt stress
%   (lagrangian_march_problems).
%
%   line_search = true damps the coupled Newton step. The iterate is the
%   solid's local dofs (sigma, u, gamma) and skeleton u_hat; with delta
%   the solved minus the iterate, the new iterate is iterate + alpha
%   delta for the largest alpha in {1, 1/2, ..., line_search_min_step}
%   with det(I + L_h) > 0 at every volume quadrature point of the solid
%   (HyperelasticProblem.deformation_jacobian) and, after the first solve
%   of the step, a relative residual at most line_search_growth times
%   the iterate's (or below 1e-6). The fluid and the interface unknowns
%   are not damped: they are not part of the linearization point, and
%   each solve recomputes them for the current (damped) solid iterate; a
%   step converges only on an undamped solve, so the returned fluid, solid
%   and state are those of one coupled solve. alpha = 1 is the undamped
%   iterate itself, bit for bit. No admissible alpha stops the loop with
%   stop_reason "failure" (an error unless adaptive_time_step retries).
%
%   adaptive_time_step = true (hyperelastic solid) retries a step whose
%   Newton loop fails (stop_reason "failure" or "max_iterations", or det
%   F <= 0 in the residual without line_search) from its saved state
%   with half the time step, down to min_time_step (0: dt/64), and
%   doubles the step again after time_step_recovery accepted steps, at
%   the first t that is a multiple of the doubled step, never above dt, so
%   that the steps return to the grid of dt; the march ends at n_steps*dt.
%   A step that still fails at min_time_step is an error
%   (lagrangian_carved_march:timeStepFloor, naming the reason). With
%   line_search and without adaptive_time_step, a step that reaches
%   newton_max_iterations on a damped iterate (alpha < 1) is an error
%   (lagrangian_carved_march:dampedNotConverged); the plain loop only
%   warns. The Newmark (solid) and
%   fluid histories store accelerations (Crank-Nicolson) or the previous
%   velocity (BDF2, whose theta follows the step ratio), so a step
%   change needs nothing else; the geometry predictor and the linear beta
%   extrapolation use the actual steps. out.time_steps holds the dt of
%   every accepted step and out.rejected the rejected attempts [t, dt,
%   reason] (1 max_iterations, 2 failure, 3 det F <= 0); snapshots are
%   saved whenever t reaches a multiple of snapshot_every*dt.
%
%   fluid_model = "navier_stokes": the fluid momentum gains the convective
%   term rho (beta.grad) u, in the conservative Oseen form of the fluid's
%   local problems (Data.convective_velocity, convective_stabilization
%   on), and the manufactured fluid source rho (u_F.grad) u_F
%   (lagrangian_march_problems, convection). The theta scheme collocates
%   every operator term at t^{n+1} (its history is the discrete
%   acceleration a^n), so beta approximates u(t^{n+1}); it is always a
%   divergence-free postprocessed velocity u* (divergence_free_velocity,
%   P_{k+1}) of a solve. Each step's u* is carried with the fluid history
%   and remapped to the new carve by remap_fluid_state (at degree k+1: a
%   kept or snapped element keeps its polynomial, a switched-on element
%   takes its source element's, so u* stays divergence-free in every
%   element). The first beta of step n+1 is the extrapolation
%   convective_extrapolation = "linear", 2 u*^n - u*^{n-1} (O(dt^2)), or
%   "constant", u*^n (O(dt)); step 1 convects with beta = 0 and needs a
%   fluid at rest at t = 0. picard_iterations = m then allows up to m
%   Picard updates per step, beta = u* of the latest solve followed by a
%   re-solve, stopping once the relative increment of the fluid velocity
%   dofs between consecutive solves is at most picard_tolerance (0 runs
%   exactly m updates); m = 0 is the semi-implicit scheme, one Oseen solve
%   per step. With the linear solid the Picard solves follow the geometry
%   iterations; with the hyperelastic solid they are merged into the
%   coupled Newton loop: each Newton solve after the first also refreshes
%   beta, and the loop stops when both the Newton and the Picard criteria
%   hold. out.picard(n) holds the updates and increments of step n. A
%   fluid with density 0 is steady Stokes and convects nothing. The
%   fluid's local problems stabilize with S tau, S = mu + rho
%   max|beta.n|/2, while the interface rows use mu tau;
%   consistent_convective_stabilization = true puts tau (1 + rho
%   max|beta.n|/(2 mu)) on the fluid mesh before every solve (beta of that
%   solve) and turns convective_stabilization off, so that both use the
%   same stabilization.
%
%   problems = @(t) ... replaces the manufactured data: a handle returning
%   [fp, sp, pd, reference, sources] at time t in the layout of
%   lagrangian_march_problems (e.g. cavity_problems, physical data); the
%   errors and the geometry error are then measured against its pd.
%
%   out: errors at t = n_steps*dt (eF, eS: compute_l2_error of
%   (sigma, u, gamma)); per step (history): t, max |d_h - u_S| on Sigma_0
%   for the curve of the final solve and for the predictor, the largest
%   |u_h| dof of fluid and solid, the background elements switched on and
%   off, and the lag max |Sigma_h - (X + u_h^{n+1})| of the curve of the
%   final solve behind the solved displacement; the transfer path lengths
%   [min, max]; the final solutions, fluid mesh and curve. out.snap, per
%   accepted step: unsnapped, the interface vertices snap_to_curve left at
%   their carve position (NaN without snap), and path_ratio, the largest
%   fluid transfer path over its element's height normal to the patch
%   face (delta/h_perp); a finite max_path_ratio stops the march with an
%   error once path_ratio exceeds it. Each snapshot also holds the step's
%   dt, newton_iterations and stop_reason (NaN and '' with the linear
%   solid), unsnapped, path_ratio, and solved_sigma_points, the curve
%   rebuilt from the solved displacement (sigma_points is the curve the
%   step used).
%
%See also: lagrangian_march_problems, solid_interface_displacement,
%          eulerian_carved_march, DeformedCurve
arguments
    k, geometry, background SimplicialGrid, solid_mesh SimplicialGrid, amplitude, dt, n_steps
    options.geometry_iterations (1,1) double {mustBeInteger, mustBeNonnegative} = 0
    options.predictor (1,1) string {mustBeMember(options.predictor, ["newmark", "taylor", "velocity"])} = "velocity"
    % the displacement whose image is Sigma_h: the HDG u_h (degree k) or
    % its P_k+1 postprocessing u*
    options.interface_geometry (1,1) string {mustBeMember(options.interface_geometry, ...
        ["displacement", "postprocessed"])} = "postprocessed"
    % Sigma_h made continuous at the solid's interface vertices
    % (solid_interface_displacement, continuous)
    options.continuous_interface (1,1) logical = true
    % hyperelastic solid: tau = solid_tau*mu/h on solid_mesh ([]: max(30 mu, 2.4 lambda)/h)
    options.solid_tau double {mustBeScalarOrEmpty, mustBePositive} = []
    % fluid time scheme after the first (backward Euler) step
    options.fluid_scheme (1,1) string {mustBeMember(options.fluid_scheme, ["cranknicolson", "bdf2"])} = "bdf2"
    options.snap (1,1) logical = false
    options.verbose (1,1) logical = false
    options.n_quadpoints double {mustBeScalarOrEmpty} = []
    % every snapshot_every steps, append the solution to snapshot_file
    % (save_march_snapshots: meshes, dofs, Sigma_h samples) for plots/movies
    options.snapshot_file (1,1) string = ""
    options.snapshot_every (1,1) double {mustBeInteger, mustBePositive} = 1
    options.solid_model (1,1) string {mustBeMember(options.solid_model, ["linear", "hyperelastic"])} = "linear"
    options.material (1,1) string {mustBeMember(options.material, ["svk", "neo_hookean"])} = "svk"
    options.solid_viscosity (1,1) double = 0
    % multiplies the fluid data (lagrangian_march_problems)
    options.fluid_scale (1,1) double = 1
    % the data at time t, [fp, sp, pd, reference, sources] = problems(t)
    options.problems function_handle = function_handle.empty
    options.fluid_model (1,1) string {mustBeMember(options.fluid_model, ["stokes", "navier_stokes"])} = "stokes"
    options.picard_iterations (1,1) double {mustBeNonnegative} = 0
    options.picard_tolerance (1,1) double {mustBeNonnegative} = 1e-8
    options.convective_extrapolation (1,1) string ...
        {mustBeMember(options.convective_extrapolation, ["linear", "constant"])} = "linear"
    % hyperelastic solid: damped coupled Newton (coupled_newton)
    options.line_search (1,1) logical = false
    options.line_search_min_step (1,1) double {mustBePositive} = 1/64
    options.line_search_growth (1,1) double {mustBePositive} = 1
    % hyperelastic solid: retry a failed step with dt/2 (see the help)
    options.adaptive_time_step (1,1) logical = false
    options.min_time_step (1,1) double {mustBeNonnegative} = 0     % 0: dt/64
    options.time_step_recovery (1,1) double {mustBeInteger, mustBePositive} = 4
    % error once a fluid transfer path exceeds this many element heights
    options.max_path_ratio (1,1) double {mustBePositive} = Inf
    % Navier-Stokes: tau(beta) on the fluid mesh (see the help)
    options.consistent_convective_stabilization (1,1) logical = false
end
hyperelastic = options.solid_model == "hyperelastic";
navier_stokes = options.fluid_model == "navier_stokes";
convection = {};
if navier_stokes, convection = {"convection", true}; end
if hyperelastic
    data = @(t) lagrangian_march_problems(k, geometry, t, dt, amplitude, ...
        "solid_model", "hyperelastic", "material", options.material, ...
        "fluid_scale", options.fluid_scale, "solid_viscosity", options.solid_viscosity, convection{:});
    newton = struct("iterations", cell(1, n_steps), "residual", [], "converged", [], "stop_reason", [], ...
        "alpha", []);
else
    data = @(t) lagrangian_march_problems(k, geometry, t, dt, amplitude, ...
        "fluid_scale", options.fluid_scale, convection{:});
end
if ~isempty(options.problems)
    data = options.problems;
end
% [fluid mesh, its background element keys]: the carved fluid (snapped
% when asked) and the centroids of the unsnapped carve, which identify
% each element's background element across steps
carve = @(curve) carve_fluid(background, curve, options.snap);
quad = {"n_quadpoints", k + 2, "n_quadpoints_path", k + 3};

if ~isempty(options.n_quadpoints)
    quad = {"n_quadpoints", options.n_quadpoints, "n_quadpoints_path", options.n_quadpoints + 1};
end

[fp, sp, ~, reference, sources] = data(0);
if hyperelastic
    [lambda, mu] = LinearElasticityProblem.get_lame_parameters(sp.Data.young_modulus, sp.Data.poisson_ratio);
    tau = max(30*mu, 2.4*lambda);
    if ~isempty(options.solid_tau), tau = options.solid_tau*mu; end
    solid_mesh = solid_mesh.add_stabilization(tau/sqrt(2*max(solid_mesh.volume)));
end
make_curve = @(dofs) displaced_curve(reference, solid_mesh, dofs, k, options.continuous_interface);
% Sigma_h of a solid state (displacement, stress, rotation)
solid_curve = @(s) make_curve(s.displacement);
if options.interface_geometry == "postprocessed"
    solid_curve = @(s) displaced_curve(reference, solid_mesh, ...
        sp.postprocessed_displacement(solid_mesh, k, s), k + 1, options.continuous_interface);
end
curve = make_curve(zeros(2*nchoosek(k + 2, 2), size(solid_mesh.elements, 2)));
sp.Data.time_scheme = "newmark";
fp.Data.time_scheme = "euler";
[fm, keys] = carve(curve);
src = sources(curve);
state = FluidStructureCoupling.initial_state(fp, sp, fm, solid_mesh, {curve, reference}, ...
    src{:}, "Order", k, quad{:});
s_samples = linspace(reference.t_range(1), reference.t_range(2), 401);
geometry_error = @(c, pd) max(vecnorm(c.get_2d_points(s_samples) - reference.get_2d_points(s_samples) ...
    - exact_on(pd.u_S, reference.get_2d_points(s_samples))));
paths = [Inf, 0];
history = zeros(n_steps, 8);
snap_history = zeros(n_steps, 2);          % [unsnapped vertices, max path/h_perp]
if navier_stokes
    assert(~any(state.fluid.velocity(:)), "lagrangian_carved_march:initialVelocity", ...
        "fluid_model = ""navier_stokes"" starts from a fluid at rest (beta = 0 on step 1)");
    picard = struct("updates", cell(1, n_steps), "increment", []);
    picard_options = struct("iterations", options.picard_iterations, ...
        "tolerance", options.picard_tolerance);
end
% u* of the last two steps, [u1*; u2*] on the current carve (u*^n, u*^{n-1})
star = {[], []};
v_previous = [];
damping = struct("line_search", options.line_search, "min_step", options.line_search_min_step, ...
    "growth", options.line_search_growth);
adaptive = options.adaptive_time_step;
assert(~adaptive || hyperelastic, "lagrangian_carved_march:adaptiveTimeStep", ...
    "adaptive_time_step needs solid_model = ""hyperelastic"" (its Newton loop decides a failure)");
min_ratio = options.min_time_step/dt;
if min_ratio == 0, min_ratio = 1/64; end
% time in units of dt: q (the time reached), ratio (the next step's dt/dt,
% 1 or a power of 1/2), h_previous (the last accepted step's); all exact
% binary fractions, so t = q dt = n dt while no step is halved; retries,
% the halvings of the step being attempted
q = 0;  ratio = 1;  h_previous = 1;  streak = 0;  retries = 0;
time_steps = zeros(n_steps, 1);
rejected = zeros(0, 3);
n = 0;
while q < n_steps
    n = n + 1;
    h = min(ratio, n_steps - q);
    step_dt = h*dt;
    t = (q + h)*dt;
    [fp, sp, pd, ~, sources] = data(t);
    if adaptive
        fp.Data.time_step = step_dt;
        sp.Data.time_step = step_dt;
    end
    sp.Data.time_scheme = "newmark";
    fp.Data.time_scheme = "cranknicolson";
    if n == 1, fp.Data.time_scheme = "euler"; end
    saved = {state, fm, keys, star, v_previous};
    s = state.solid;
    curve = solid_curve(predicted_solid(s, step_dt, predictor_weight(options.predictor)*step_dt^2));
    predictor_error = geometry_error(curve, pd);
    [fm_old, keys_old] = deal(fm, keys);
    [fm, keys, snap_info] = carve(curve);
    [state.fluid, switched] = remap_fluid_state(state.fluid, fm_old, fm, k, keys_old, keys);
    solve_state = state;
    if options.fluid_scheme == "bdf2" && n > 1
        % BDF2 (variable step, omega = dt_n/dt_{n-1}) as the theta step
        % with theta = (1 + omega)/(1 + 2 omega) from the history
        % w = ((1 + omega)^2 v^n - omega^2 v^{n-1})/(1 + 2 omega) and a zero
        % acceleration: backward Euler from w with the step theta dt
        v_previous = remap_fluid_state(struct("velocity", v_previous), fm_old, fm, k, keys_old, keys).velocity;
        omega = h/h_previous;
        fp.Data.theta = (1 + omega)/(1 + 2*omega);
        solve_state.fluid.velocity = ((1 + omega)^2*state.fluid.velocity - omega^2*v_previous)/(1 + 2*omega);
        solve_state.fluid.acceleration = zeros(size(state.fluid.velocity));
    end
    convects = navier_stokes && fp.Data.density > 0;
    if convects
        for j = find(~cellfun(@isempty, star))
            remapped = remap_fluid_state(struct("velocity", star{j}), fm_old, fm, k + 1, keys_old, keys);
            star{j} = remapped.velocity;
        end
        fp.Data.convective_velocity = extrapolated_beta(star, options.convective_extrapolation, ...
            h/h_previous);
        if options.consistent_convective_stabilization
            fp.Data.convective_stabilization = false;
            tau0 = fm.tau;
            fm = convective_tau(fp, fm, tau0, k);
            picard_options.tau0 = tau0;
        end
    end
    if hyperelastic
        if ~convects, picard_options = struct("iterations", 0, "tolerance", 0); end
        try
            [fs, ss, new_state, gs, curve, newton(n), picard_report] = coupled_newton(fp, sp, fm, ...
                solid_mesh, curve, reference, sources, solid_curve, solve_state, options.geometry_iterations, ...
                k, quad, picard_options, damping);
            failure = newton(n).stop_reason;
        catch err
            if ~adaptive || err.identifier ~= "HyperelasticProblem:inverted", rethrow(err); end
            failure = "inverted";
        end
        failed = failure == "failure" || failure == "inverted" || (adaptive && failure == "max_iterations");
        if failed && adaptive && h/2 >= min_ratio
            % retry from the saved step with half the time step
            rejected(end + 1, :) = [q*dt, step_dt, find(failure == ["max_iterations", "failure", "inverted"])]; %#ok<AGROW>
            if options.verbose
                fprintf("  step %d t %.4f dt %.3e rejected (%s): retry with dt %.3e\n", ...
                    n, t, step_dt, failure, step_dt/2);
            end
            [state, fm, keys, star, v_previous] = saved{:};
            ratio = h/2;  streak = 0;  n = n - 1;  retries = retries + 1;
            continue
        end
        if failed && adaptive
            error("lagrangian_carved_march:timeStepFloor", ...
                "Step %d (t = %g) failed (%s) at dt = %g = dt/2^%d, the floor min_time_step = %g, " + ...
                "after %d halving(s) of this step.", n, t, failure, step_dt, round(log2(1/h)), ...
                min_ratio*dt, retries);
        end
        if failure == "failure" || failure == "inverted"
            error("lagrangian_carved_march:stepFailed", ...
                "The coupled Newton loop of step %d (t = %g, dt = %g) found no admissible iterate.", ...
                n, t, step_dt);
        end
        if failure == "max_iterations" && options.line_search && newton(n).alpha(end) < 1
            % the solved state is the undamped update the line search rejected
            error("lagrangian_carved_march:dampedNotConverged", ...
                "The coupled Newton loop of step %d (t = %g, dt = %g) reached newton_max_iterations " + ...
                "= %d on a damped iterate (alpha = %g, relative residual %.3e).", n, t, step_dt, ...
                newton(n).iterations, newton(n).alpha(end), newton(n).residual(end));
        end
    else
        for iteration = 0:options.geometry_iterations
            if iteration > 0
                curve = solid_curve(new_state.solid);
            end
            src = sources(curve);
            [fs, ss, new_state, gs] = FluidStructureCoupling.solve(fp, sp, fm, solid_mesh, ...
                {curve, reference}, src{:}, "newmark", "quiet", true, "state", solve_state, quad{:});
        end
        if convects
            % Picard on the converged geometry: beta = u* of the last solve
            picard_report = struct("updates", 0, "increment", zeros(1, 0));
            while picard_report.updates < options.picard_iterations
                u_previous = fs.Unknowns.field(2).dof;
                fp.Data.convective_velocity = divergence_free_velocity(fs);
                if options.consistent_convective_stabilization
                    fm = convective_tau(fp, fm, tau0, k);
                end
                [fs, ss, new_state, gs] = FluidStructureCoupling.solve(fp, sp, fm, solid_mesh, ...
                    {curve, reference}, src{:}, "newmark", "quiet", true, "state", solve_state, quad{:});
                picard_report.updates = picard_report.updates + 1;
                picard_report.increment(end + 1) = relative_increment(fs, u_previous);
                if picard_report.increment(end) <= options.picard_tolerance, break; end
            end
        end
    end
    snap_history(n, :) = [snap_info.n_unsnapped, fluid_path_ratio(gs.partition_1, fm)];
    if snap_history(n, 2) > options.max_path_ratio
        error("lagrangian_carved_march:pathRatio", ...
            "At t = %g a fluid transfer path reaches %.3g element heights (max_path_ratio = %g); " + ...
            "%d interface vertices unsnapped.", t, snap_history(n, 2), options.max_path_ratio, ...
            snap_info.n_unsnapped);
    end
    if convects
        picard(n) = picard_report;
        star = {stacked_components(divergence_free_velocity(fs)), star{1}};
    end
    v_previous = state.fluid.velocity;            % v^n on the current carve (bdf2)
    state = new_state;
    q = q + h;
    time_steps(n) = step_dt;
    h_previous = h;
    retries = 0;
    if ratio < 1
        % restore the time step after time_step_recovery accepted steps, at
        % the first t on the doubled step's grid (q a multiple of 2 ratio)
        streak = streak + 1;
        if streak >= options.time_step_recovery && mod(q, 2*ratio) == 0
            ratio = min(1, 2*ratio);  streak = 0;
        end
    end
    % lag of the geometry used behind the solved displacement
    solved_points = solid_curve(state.solid).get_2d_points(s_samples);
    lag = max(vecnorm(curve.get_2d_points(s_samples) - solved_points));
    for part = {gs.partition_1, gs.partition_2}
        l = vecnorm(part{1}.quadrature.interface - part{1}.quadrature.physical, 2, 1);
        paths = [min(paths(1), min(l(:))), max(paths(2), max(l(:)))];
    end
    history(n, :) = [t, geometry_error(curve, pd), predictor_error, ...
        max(abs(fs.Unknowns.field(2).dof), [], "all"), ...
        max(abs(ss.Unknowns.field(2).dof), [], "all"), switched, lag];
    if options.snapshot_file ~= "" && (floor(q/options.snapshot_every) > floor((q - h)/options.snapshot_every) ...
            || q >= n_steps)
        diagnostics = struct("dt", step_dt, "newton_iterations", NaN, "stop_reason", '', ...
            "unsnapped", snap_history(n, 1), "path_ratio", snap_history(n, 2), ...
            "solved_sigma_points", solved_points);
        if hyperelastic
            diagnostics.newton_iterations = newton(n).iterations;
            diagnostics.stop_reason = char(newton(n).stop_reason);
        end
        save_march_snapshots(options.snapshot_file, k, t, fs, ss, state, ...
            curve.get_2d_points(s_samples), diagnostics);
    end
    if options.verbose
        eF = fs.compute_l2_error({pd.sigma_F, pd.u_F, pd.gamma_F});
        eS = ss.compute_l2_error({pd.sigma_S, pd.u_S, pd.gamma_S});
        fprintf("  step %d t %.4f dt %.3e  uF %.3e  uS %.3e  geometry %.2e (predictor %.2e)  on/off %d/%d\n", ...
            n, t, step_dt, eF(2), eS(2), history(n, 2), predictor_error, switched);
        if hyperelastic
            fprintf("    Newton (%s): %s; alpha %s\n", newton(n).stop_reason, ...
                sprintf("%.2e ", newton(n).residual), sprintf("%g ", newton(n).alpha));
        end
        if convects
            fprintf("    Picard %d: %s\n", picard(n).updates, sprintf("%.2e ", picard(n).increment));
        end
    end
end
eF = fs.compute_l2_error({pd.sigma_F, pd.u_F, pd.gamma_F});
eS = ss.compute_l2_error({pd.sigma_S, pd.u_S, pd.gamma_S});
[offset, e_rest] = fluid_stress_offset(fs, pd.sigma_F);
out = struct("eF", eF(:).', "eS", eS(:).', "paths", paths, "fs", fs, "ss", ss, ...
    "fluid_stress_offset", [offset, e_rest], ...
    "fluid_mesh", fm, "curve", curve, "history", history, "time_steps", time_steps(1:n), ...
    "snap", struct("unsnapped", snap_history(1:n, 1), "path_ratio", snap_history(1:n, 2)));
if adaptive
    out.rejected = rejected;
end
if hyperelastic
    out.newton = newton;
end
if navier_stokes
    out.picard = picard;
end
end

function beta = extrapolated_beta(star, extrapolation, rho)
% the first beta of a step from u*^n, u*^{n-1} (stacked, on the current
% carve; [] before they exist): {} (no convection) on step 1. rho is the
% ratio dt^{n+1}/dt^n of the step to the previous one: the linear
% extrapolation is u*^n + rho (u*^n - u*^{n-1})
if isempty(star{1})
    beta = {};
    return
end
b = star{1};
if extrapolation == "linear" && rho ~= 1
    b = (1 + rho)*star{1};
    if ~isempty(star{2}), b = b - rho*star{2}; end
elseif extrapolation == "linear" && ~isempty(star{2})
    b = 2*star{1} - star{2};
elseif extrapolation == "linear"
    b = 2*star{1};                       % u*^0 = 0: the fluid starts at rest
end
n = size(b, 1)/2;
beta = {b(1:n, :); b(n + 1:end, :)};
end

function U = stacked_components(u)
% {u1; u2} -> [u1; u2]
U = cat(1, u{:});
end

function r = relative_increment(fs, u_previous)
% relative change of the fluid velocity dofs between consecutive solves
u = fs.Unknowns.field(2).dof;
r = norm(u - u_previous, "fro")/max(norm(u, "fro"), realmin);
end

function [fs, ss, new_state, gs, curve, report, picard_report] = coupled_newton(fp, sp, fm, ...
    solid_mesh, curve, reference, sources, solid_curve, state, geometry_iterations, k, quad, picard, ...
    damping)
% one step of the coupled Newton loop of the hyperelastic solid (see the
% help): sp is a HyperelasticProblem at t^{n+1}, curve the predictor;
% picard.iterations > 0 refreshes the fluid's beta after each solve until
% the Picard criterion holds (picard_report: updates, increments)
d1 = k + 1;
d2 = nchoosek(k + 2, 2);
n_faces = size(solid_mesh.faces, 2);
indices = struct();
[~, ~, indices.rhs_rows] = get_assemble_indices(sp, solid_mesh, k);
d2peers = 4*d2 + (k > 0)*d1;
rows = [1:d2peers, d2peers + 2*d2 + (1:d2)];    % constitutive and symmetry rows
with_history = sp;                                % the momentum rows (and a viscous
with_history.Data.state = state;                  % solid's rate) read the state
local_rows = @(X, x_hat, varargin) newton_residual(with_history, solid_mesh, k, X, x_hat, ...
    indices, "local_rows", true, varargin{:});
r_ref = max(norm(newton_residual(with_history, solid_mesh, k, ...
    zeros(d2peers + 3*d2, size(solid_mesh.elements, 2)), zeros(2*d1*n_faces, 1), indices, ...
    "local_rows", true, "nonlinear", false), "fro"), realmin);
sp.Data.iterate = struct("displacement", state.solid.displacement, ...
    "displacement_skeleton", state.solid.displacement_skeleton);
sp.Data.state = state;                            % the T2 history of the solid path
history = zeros(1, 0);
alpha = zeros(1, 0);
converged = false;
stop_reason = "max_iterations";
picard_report = struct("updates", 0, "increment", zeros(1, 0));
u_previous = [];
if damping.line_search
    % the accepted iterate: all local dofs [sigma; u; gamma] and the skeleton
    n_elements = size(solid_mesh.elements, 2);
    X_it = [state_block(state.solid, "stress", d2peers, n_elements); state.solid.displacement;
            state_block(state.solid, "rotation", d2, n_elements)];
    skeleton_it = state.solid.displacement_skeleton;
    % the state does not solve the step's linear rows: no residual bound
    % on the first step from it
    r_it = Inf;
end
for iteration = 1:max(sp.Data.newton_max_iterations, geometry_iterations + 1)
    if iteration > 1 && iteration <= geometry_iterations + 1
        curve = solid_curve(sp.Data.iterate);
    end
    src = sources(curve);
    [fs, ss, new_state, gs] = FluidStructureCoupling.solve(fp, sp, fm, solid_mesh, ...
        {curve, reference}, src{:}, "newmark", "quiet", true, "state", state, quad{:});
    dofs = arrayfun(@(f) as_matrix(f.dof), ss.Unknowns.field, "UniformOutput", false);
    X = cat(1, dofs{:});
    skeleton = as_matrix(ss.Unknowns.hybrid(1).dof);
    x_hat = [reshape(skeleton(1:d1, :), [], 1); reshape(skeleton(d1 + 1:end, :), [], 1)];
    if iteration == 1 && r_ref == realmin
        r_X = local_rows(X, 0*x_hat, "nonlinear", false);
        r_ref = max(norm(r_X(rows, :), "fro"), realmin);
    end
    step = 1;
    if damping.line_search
        [step, X, skeleton, relative] = line_search(with_history, solid_mesh, k, local_rows, rows, ...
            r_ref, X_it, skeleton_it, r_it, X, skeleton, d2peers, d2, damping);
        alpha(end + 1) = step; %#ok<AGROW>
        history(end + 1) = relative; %#ok<AGROW>
        if isnan(step)
            stop_reason = "failure";
            break
        end
        [X_it, skeleton_it, r_it] = deal(X, skeleton, relative);
    else
        alpha(end + 1) = 1; %#ok<AGROW>
        r = local_rows(X, x_hat);
        history(end + 1) = norm(r(rows, :), "fro")/r_ref; %#ok<AGROW>
    end
    sp.Data.iterate = struct("displacement", X(d2peers + (1:2*d2), :), ...
        "displacement_skeleton", skeleton);
    picard_done = true;
    if picard.iterations > 0
        if ~isempty(u_previous)
            picard_report.increment(end + 1) = relative_increment(fs, u_previous);
        end
        picard_done = picard_report.updates >= picard.iterations || ...
            (~isempty(picard_report.increment) && picard_report.increment(end) <= picard.tolerance);
        if ~picard_done
            u_previous = fs.Unknowns.field(2).dof;
            fp.Data.convective_velocity = divergence_free_velocity(fs);
            if isfield(picard, "tau0"), fm = convective_tau(fp, fm, picard.tau0, k); end
            picard_report.updates = picard_report.updates + 1;
        end
    end
    if iteration <= geometry_iterations || step < 1
        continue                         % a damped step is never the last
    end
    stagnated = numel(history) > 1 && history(end) <= 1e-6 && history(end) > history(end - 1)/2;
    if picard_done && (history(end) <= sp.Data.newton_tolerance || stagnated)
        converged = true;
        stop_reason = "stagnation";
        if history(end) <= sp.Data.newton_tolerance, stop_reason = "tolerance"; end
        break
    end
end
if ~converged
    warning("HyperelasticProblem:newtonNotConverged", ...
        "Coupled Newton stopped after %d steps at relative residual %.3e.", ...
        iteration, history(end));
end
report = struct("iterations", iteration, "residual", history, "converged", converged, ...
    "stop_reason", stop_reason, "alpha", alpha);
end

function [step, X, skeleton, relative] = line_search(pde, solid_mesh, k, local_rows, rows, r_ref, ...
    X_it, skeleton_it, r_it, X_new, skeleton_new, d2peers, d2, damping)
% the largest step = 1, 1/2, ..., >= damping.min_step from the iterate
% (X_it, skeleton_it) toward the solved (X_new, skeleton_new) whose
% deformation gradient has det F > 0 at every quadrature point
% (deformation_jacobian) and whose relative residual is at most
% damping.growth times the iterate's r_it (or below the round-off floor
% 1e-6); step = NaN when none is. step = 1 returns the solved dofs
% unchanged. relative is the residual of the returned (last tried) point.
step = 1;
while true
    if step == 1
        [X, skeleton] = deal(X_new, skeleton_new);
    else
        X = X_it + step*(X_new - X_it);
        skeleton = skeleton_it + step*(skeleton_new - skeleton_it);
    end
    trial = struct("displacement", X(d2peers + (1:2*d2), :), "displacement_skeleton", skeleton);
    relative = Inf;
    admissible = min(deformation_jacobian(pde, solid_mesh, k, trial), [], "all") > 0;
    if admissible
        r = local_rows(X, trace_vector(skeleton, k + 1));
        relative = norm(r(rows, :), "fro")/r_ref;
        admissible = relative <= damping.growth*r_it || relative <= 1e-6;
    end
    if admissible
        return
    end
    if step/2 < damping.min_step
        step = NaN;
        return
    end
    step = step/2;
end
end

function x_hat = trace_vector(skeleton, d1)
% the global hybrid vector of a 2 d1 x F skeleton array (x block over y block)
x_hat = [reshape(skeleton(1:d1, :), [], 1); reshape(skeleton(d1 + 1:end, :), [], 1)];
end

function A = state_block(s, name, n_rows, n_elements)
% a state.solid array of local dofs, zero when absent
A = zeros(n_rows, n_elements);
if isfield(s, name) && ~isempty(s.(name)) && ~isequal(s.(name), 0)
    A = as_matrix(s.(name));
end
end

function A = as_matrix(A)
if iscell(A), A = cat(1, A{:}); end
end

function [fluid, keys, info] = carve_fluid(background, curve, snap)
% info: snap_to_curve's diagnostics (n_unsnapped NaN without snapping).
% With snapping the normal-deviation prune of split is off: snapping makes
% every interface face a chord of the curve, while the prune, where the
% curve's normal lies between two face directions of the grid, removes
% elements the curve does not cut and can leave a notch deeper than one
% element, whose vertices cannot be snapped without inverting a neighbour
if snap
    fluid = background.split(curve, max_normal_deviation = 90);
else
    fluid = background.split(curve);
end
keys = squeeze(mean(reshape(fluid.coordinates(:, fluid.elements), 2, 3, []), 2)).';
info = struct("n_unsnapped", NaN);
if snap
    [fluid, info] = snap_to_curve(fluid, curve);
end
end

function fm = convective_tau(fp, fm, tau0, k)
% the fluid mesh with tau = tau0 (1 + rho max|beta.n|/(2 mu)) for the
% fluid's beta (tau0 without one); with convective_stabilization off the
% local problem uses S = mu on this tau, and the interface rows read the
% same mesh tau
fm.tau = tau0;
beta = fp.Data.convective_velocity;
if ~isempty(beta)
    [~, ~, ~, max_bn] = convective_matrix(fm, k, beta);
    fm.tau = tau0*(1 + fp.Data.density*max_bn/(2*fp.Data.viscosity));
end
end

function ratio = fluid_path_ratio(partition, mesh)
% the largest transfer path length of the fluid's induced partition over
% the height of its element normal to the patch's face, 2|K|/|F|
l = vecnorm(partition.quadrature.interface - partition.quadrature.physical, 2, 1);
h_perp = 2*mesh.volume(partition.elements(:)) ./ mesh.area(partition.faces(:));
ratio = max(l ./ reshape(h_perp, 1, 1, []), [], "all");
end

function c = predictor_weight(predictor)
% the weight of dt^2 a^n in the geometry predictor. The accelerations of
% the fields without inertia (skeleton, stress, rotation) come from the
% trapezoidal Newmark recursion alone, a^{n+1} = 4/dt^2 (x^{n+1} - x^n -
% dt v^n) - a^n, which no equation pins: any jump of x^{n+1} leaves in
% a an undamped mode alternating in sign from step to step. x and v
% read only a^n + a^{n+1} and do not see it; a predictor reading a^n
% does. "velocity" omits a^n, at a predictor error of order dt^2 like
% the other two
c = struct("newmark", 1/4, "taylor", 1/2, "velocity", 0).(predictor);
end

function p = predicted_solid(s, dt, c)
% the Newmark predictor u + dt v + c a of the solid displacement, and of
% its skeleton trace, stress and rotation where the history holds them
% (absent rates are zero)
p.displacement = s.displacement + dt*s.velocity + c*s.acceleration;
rates = struct("displacement_skeleton", ["velocity_skeleton", "acceleration_skeleton"], ...
    "stress", ["stress_velocity", "stress_acceleration"], ...
    "rotation", ["rotation_velocity", "rotation_acceleration"]);
weights = [dt, c];
for name = string(fieldnames(rates)).'
    if ~isfield(s, name), continue, end
    p.(name) = as_matrix(s.(name));
    for j = 1:2
        rate = rates.(name)(j);
        if isfield(s, rate) && ~isempty(s.(rate)) && ~isequal(s.(rate), 0)
            p.(name) = p.(name) + weights(j)*as_matrix(s.(rate));
        end
    end
end
end

function curve = displaced_curve(reference, solid_mesh, dofs, k, continuous)
[d, d_prime] = solid_interface_displacement(reference, solid_mesh, dofs, k, continuous = continuous);
if ~is_periodic(reference)
    % The open (flat) Sigma ends on the vertical walls x = 0, 1, where the
    % exact horizontal displacement vanishes and the fluid chain ends: the
    % discrete d_1 is corrected by its linear interpolant between the two
    % ends (an O(end error) change), so the ends slide along the walls,
    % and the fluid's wall vertices map to the ends (WallEndedCurve).
    t0 = reference.t_range(1);  L = diff(reference.t_range);
    ends = d([t0, t0 + L]);
    d_wall = @(t) [1; 0] .* (ends(1, 1) + (t - t0)/L*(ends(1, 2) - ends(1, 1)));
    d = @(t) d(t) - d_wall(t);
    d_prime = @(t) d_prime(t) - [1; 0] .* ((ends(1, 2) - ends(1, 1))/L + 0*t);
    curve = WallEndedCurve(reference, d, [0 1], "parametric", true, "derivative", d_prime);
    return
end
curve = DeformedCurve(reference, d, parametric = true, derivative = d_prime);
end

function v = exact_on(u, P)
v = [u{1}(P(1, :), P(2, :)) + 0*P(1, :); u{2}(P(1, :), P(2, :)) + 0*P(1, :)];
end

function [s, switched] = remap_fluid_state(s, old_mesh, new_mesh, k, old_keys, new_keys)
% the volume histories (velocity, acceleration; Dubiner dofs, component
% blocks stacked, one column per element) from old_mesh to new_mesh, both
% carved from the same background grid (eulerian_carved_march's remap,
% with elements identified by their background element, the keys): an
% element present at both steps keeps its polynomial, re-expressed on its
% new shape if snapping moved it (exact); an element that switches on
% gets the L2 projection of the nearest old element's polynomial.
% switched counts the background elements switched on and off.
centroid = @(m) squeeze(mean(reshape(m.coordinates(:, m.elements), 2, 3, []), 2)).';
[present, from] = ismembertol(new_keys, old_keys, 1e-10, "ByRows", true);
switched = [nnz(~present), size(old_mesh.elements, 2) - nnz(present)];
c_old = centroid(old_mesh);
c_new = centroid(new_mesh);
same_shape = present;
same_shape(present) = vecnorm(c_new(present, :) - c_old(from(present), :), 2, 2) < 1e-14;
kept = find(same_shape);
added = find(~same_shape);
source = zeros(size(added));
for j = 1:numel(added)
    if present(added(j))
        source(j) = from(added(j));
    else
        [~, source(j)] = min(vecnorm(c_old - c_new(added(j), :), 2, 2));
    end
end
d2 = nchoosek(k + 2, 2);
[rs, w] = get_2d_reference_quadrature(2*k + 2);
P = dubiner_basis(k, rs);
fit = (P.'*(w(:).*P)) \ (P.'.*w(:).');          % d2 x n_q: values -> dofs
for name = ["velocity", "acceleration"]
    if ~isfield(s, name) || ~isnumeric(s.(name)) || isempty(s.(name)) || isequal(s.(name), 0)
        continue
    end
    D = s.(name);
    E = zeros(size(D, 1), size(new_mesh.elements, 2));
    E(:, kept) = D(:, from(kept));
    for j = 1:numel(added)
        xy = new_mesh.get_2d_physical(rs, added(j));
        rs_old = old_mesh.get_2d_local(xy, source(j));
        Q = dubiner_basis(k, rs_old);
        for c = 0:size(D, 1)/d2 - 1
            rows = c*d2 + (1:d2);
            E(rows, added(j)) = fit*(Q*D(rows, source(j)));
        end
    end
    s.(name) = E;
end
end
