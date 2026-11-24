%-------------------------------------------------------------------------
%Description: tip_mapping = "collapse" of get_partitions on the interior
%             (gapped) Joukowski airfoil pair, and the default "spread".
%
%             On the exact cusp (regularized = 0) the closest point of the
%             vertices behind the trailing edge is the tip itself, so the
%             faces between them have a zero-length image. "spread"
%             stretches them along Sigma and, from h = 0.05 on, leaves two
%             without a patch (orphanFace: the tied run straddles the
%             parameter seam). "collapse" sends each such face whole to
%             the tip. Checked:
%               1. default: fingerprints of the default mapping (free-dof
%                  hash exact; solution, rhs and a matrix probe by norm and
%                  one projection, relative 1e-12)
%               2. regularized airfoil (0.02): nothing collapses and
%                  "collapse" reproduces "spread" bit for bit
%               3. open curve: "collapse" is refused
%                  (map_mesh_to_interface:collapseOpenCurve), "spread"
%                  reports the faces beyond the ends (orphanFace)
%               4. cusp partitions, fluid outside (fluid faces collapse)
%                  and inside (solid faces collapse): no orphaned face,
%                  whole face on the collapsed side, area 0 on the other,
%                  every interface node at the tip, collapsed_side
%               5. the row rule of InterfaceProblem.solve: the kinematic
%                  row must be tested on the collapsed side; the other
%                  orientation and collapse on both sides are refused
%               6. supported configurations: the solution does not depend
%                  on collapsed_normal (bit for bit), and is accurate
%                  (k = 1, bounds from the measured errors, and error
%                  reduction from h = 0.1 to 0.05)
%               7. history transport on the point patches: the minimum-
%                  norm fit of trace_extension_coefficients reproduces the
%                  path integral at the paired point (k = 1, 2)
%               8. a transient step with collapse (newmark solid,
%                  Crank-Nicolson fluid, through initial_state and
%                  interface_datum), fluid outside and inside: finite,
%                  and the t = 0 state and the step independent of
%                  collapsed_normal (bit for bit)
%-------------------------------------------------------------------------
n_failed = 0;
box = [-3 3 -0.4 0.6];
cusp = AirfoilCurve(-1/10, 1/10, 0);

sha = @(v) lower(reshape(dec2hex(typecast(int8( ...
    java.security.MessageDigest.getInstance('SHA-256').digest( ...
    typecast(double(v(:)), 'uint8'))), 'uint8')).', 1, []));
probe = @(n) cos((1:n).');

%% 1. default mapping: fingerprints
% columns: regularized, h, k, free-dof SHA-256, then norm(x), cos'*x,
% norm(rhs), cos'*rhs, norm(A*cos)
reference = { ...
    0.02, 0.1, 1, "e70c06305b4aad25469c58ad1923d181ae62ecd18f9b95a91a447c67d92f10df", ...
        [122.1581328413256, 2.4116219430139907, 447.08111374415188, -61.280915777988128, 25905.29771933037]; ...
    0.02, 0.1, 2, "f033219119c11f0d7ca013f8803171428da749e8b49cb42789f0d8a2ddbf1323", ...
        [117.79361996150885, 9.6140689550076175, 1523.9841403232952, 118.07797563044005, 2085592.588981559]; ...
    0,    0.1, 1, "7fb9e1baf5c4e79113dce91d9b83951d1ff01a1dcffcaca87129c13708801c3b", ...
        [92348.910472171017, -12797.863533766707, 447.09860917562031, 109.86933880764911, 48099.019516413893]};
for c = 1:size(reference, 1)
    [reg, h, k, dof_hash, values] = reference{c, :};
    o = airfoil_solve(reg, h, k, "outer");
    g = o.gs;
    w = probe(size(g.matrix, 2));
    got = [norm(g.solution), w.' * g.solution, norm(g.rhs), w.' * g.rhs, norm(g.matrix * w)];
    label = sprintf("default reg = %g, h = %g, k = %d", reg, h, k);
    n_failed = n_failed + check(label + ": free dofs (hash)", string(sha(g.dof)) == dof_hash);
    n_failed = n_failed + check(label + ": solution, rhs, matrix probe (rel 1e-12)", ...
        all(abs(got - values) <= 1e-12 * max(abs(values), [], "all")));
end

%% 2. regularized airfoil: nothing collapses, collapse == spread
o_s = airfoil_solve(0.02, 0.1, 1, "outer");
o_c = airfoil_solve(0.02, 0.1, 1, "outer", "tip_mapping", "collapse");
n_failed = n_failed + check("reg = 0.02: no collapsed patch", ~any(o_c.gs.partition_1.collapsed));
n_failed = n_failed + check("reg = 0.02: collapse reproduces spread bit for bit", ...
    isequal(o_s.gs.solution, o_c.gs.solution) && isequal(o_s.gs.matrix, o_c.gs.matrix) ...
    && isequal(o_s.gs.rhs, o_c.gs.rhs));

%% 3. open curve: faces beyond the ends of Sigma
segment = FlatCurveHorizontal(0.5, [0.25 0.75]);
[fm_o, sm_o] = interface_meshes(FlatCurveHorizontal(0.5), "interior", 1/16, ...
    box = [0 1 0 1], fluid_side = "lower");
n_failed = n_failed + check("open curve, spread: orphanFace", ...
    error_id(@() get_partitions(segment, fm_o, sm_o)) == "get_partitions:orphanFace");
n_failed = n_failed + check("open curve, collapse: refused", ...
    error_id(@() get_partitions(segment, fm_o, sm_o, tip_mapping = "collapse")) ...
    == "map_mesh_to_interface:collapseOpenCurve");

%% 4. cusp partitions
h = 0.05;
n_failed = n_failed + check("cusp, h = 0.05: spread leaves orphaned faces", ...
    error_id(@() get_partitions(cusp, ...
    interface_meshes(cusp, "interior", h, box = box, fluid_side = "outer"), ...
    second_mesh(cusp, h, box, "outer"))) == "get_partitions:orphanFace");
for fluid_side = ["outer", "inner"]
    [fm, sm] = interface_meshes(cusp, "interior", h, box = box, fluid_side = fluid_side);
    [p1, p2] = get_partitions(cusp, fm, sm, tip_mapping = "collapse", ...
        quadrature_patch = 2, quadrature_path = 3);
    side = 1 + (fluid_side == "inner");
    meshes = {fm, sm};  parts = {p1, p2};
    own = parts{side};  other = parts{3 - side};  own_mesh = meshes{side};
    c = find(p1.collapsed);
    label = "cusp, fluid " + fluid_side;
    n_failed = n_failed + check(sprintf("%s: %d collapsed faces, all on side %d", label, numel(c), side), ...
        ~isempty(c) && all(p1.collapsed_side(c) == side) && isequal(p1.collapsed_side, p2.collapsed_side) ...
        && all(p1.collapsed_side(~p1.collapsed) == 0));
    n_failed = n_failed + check(label + ": no orphaned face", ...
        isempty(setdiff(fm.intfaces(:).', p1.faces(:).')) && isempty(setdiff(sm.intfaces(:).', p2.faces(:).')));
    n_failed = n_failed + check(label + ": collapsed area = face length, paired area = 0", ...
        max(abs(own.area(c) - own_mesh.area(squeeze(own.faces(1, 1, c)).'))) < 1e-14 ...
        && all(other.area(c) == 0));
    t_star = squeeze(p1.vertices_interface(1, :, c));
    xbar = reshape(cusp.get_2d_points(t_star(1)), 2, 1);
    n_failed = n_failed + check(label + ": one t*, every collapsed interface node at the tip", ...
        all(t_star(:) == t_star(1)) && max(abs(p1.quadrature.interface(:, :, c) - xbar), [], "all") < 1e-14 ...
        && isequal(p1.quadrature.interface(:, :, c), p2.quadrature.interface(:, :, c)));
end

%% 5. the row rule
n_failed = n_failed + check("rule: fluid faces collapse, default rows refused", ...
    error_id(@() airfoil_solve(0, 0.1, 1, "outer", "tip_mapping", "collapse")) ...
    == "InterfaceProblem:collapsedRows");
n_failed = n_failed + check("rule: solid faces collapse, swapped rows refused", ...
    error_id(@() airfoil_solve(0, 0.1, 1, "inner", "tip_mapping", "collapse", ...
    "swap_transmission", true)) == "InterfaceProblem:collapsedRows");
[fm, ~] = interface_meshes(cusp, "interior", 0.1, box = box, fluid_side = "outer");
[q1, ~] = get_partitions(cusp, fm, fm, tip_mapping = "collapse");
n_failed = n_failed + check("rule: the same tied mesh on both sides collapses on both", ...
    isequal(unique(q1.collapsed_side(q1.collapsed)), [1, 2]));
n_failed = n_failed + check("rule: collapse on both sides refused", ...
    error_id(@() airfoil_solve(0, 0.1, 1, "outer", "tip_mapping", "collapse", ...
    "swap_transmission", true, meshes = {fm, fm})) == "InterfaceProblem:collapsedBothSides");

%% 6. supported configurations: normal independence and accuracy (k = 1)
% fluid outside: swapped rows; fluid inside: default rows
supported = {"outer", true, [0.1, 0.02]; "inner", false, [0.15, 0.4]};
for s = 1:size(supported, 1)
    [fluid_side, swap, bound] = supported{s, :};
    args = {"tip_mapping", "collapse", "swap_transmission", swap};
    a = airfoil_solve(0, 0.05, 1, fluid_side, args{:}, "collapsed_normal", "path");
    b = airfoil_solve(0, 0.05, 1, fluid_side, args{:}, "collapsed_normal", "face");
    coarse = airfoil_solve(0, 0.1, 1, fluid_side, args{:});
    label = sprintf("fluid %s, swap %d", fluid_side, swap);
    n_failed = n_failed + check(label + ": independent of collapsed_normal (bit for bit)", ...
        isequal(a.gs.solution, b.gs.solution) && isequal(a.gs.matrix, b.gs.matrix) && isequal(a.gs.rhs, b.gs.rhs));
    fprintf("  %s: e(u_F) %.3e -> %.3e, e(u_S) %.3e -> %.3e (h = 0.1 -> 0.05), residue %.1e\n", ...
        label, coarse.eF(2), a.eF(2), coarse.eS(2), a.eS(2), a.gs.residue);
    n_failed = n_failed + check(sprintf("%s: h = 0.05 e(u_F) < %g, e(u_S) < %g", label, bound), ...
        a.eF(2) < bound(1) && a.eS(2) < bound(2) && a.gs.residue < 1e-10);
    n_failed = n_failed + check(label + ": errors drop from h = 0.1 to 0.05 (factor < 0.7)", ...
        a.eF(2) < 0.7 * coarse.eF(2) && a.eS(2) < 0.7 * coarse.eS(2));
end

%% 7. history transport on the point patches
[fm, sm] = interface_meshes(cusp, "interior", 0.05, box = box, fluid_side = "outer");
for k = 1:2
    [~, p2] = get_partitions(cusp, fm, sm, tip_mapping = "collapse", ...
        quadrature_patch = k + 1, quadrature_path = k + 2);
    rng(k);
    n_el = size(sm.elements, 2);
    sigma = randn(PEERSBasis(k).Dimension, n_el);
    gamma = randn(nchoosek(k + 2, 2), n_el);
    solid = LinearElasticityProblem(k, young_modulus = 1, poisson_ratio = 0.3333);
    ext = trace_extension_coefficients(solid, sm, p2, sigma, gamma, Order = k);
    [lambda, mu] = LinearElasticityProblem.get_lame_parameters(1, 0.3333);
    c1 = 1/(2*mu);  c2 = lambda/(2*mu*(2*lambda + 2*mu));
    G = get_affine_map(sm, p2.elements, inverse = true);
    worst = 0;
    for e = find(p2.collapsed)
        K = p2.elements(e);
        reference_value = path_gradient_integral(k, p2, e, G(:, :, e), sigma(:, K), gamma(:, K), c1, c2);
        l = legendre_basis(k, p2.quadrature.local(1, :, e));
        fitted = l * [ext(1:k + 1, e), ext(k + 2:end, e)];
        worst = max(worst, max(abs(fitted - reference_value), [], "all") / max(abs(reference_value), [], "all"));
    end
    n_failed = n_failed + check(sprintf("history k = %d: point patches finite, fit = path integral (rel %.1e)", k, worst), ...
        all(isfinite(ext), "all") && worst < 1e-12);
end

%% 8. transient: one step with collapse, independent of collapsed_normal
% fluid outside (swapped rows) and inside (default rows); the t = 0 state
% (initial_state, through interface_datum) and the step must not depend
% on the arbitrary normal at the collapse point
for fluid_side = ["outer", "inner"]
    try
        a = transient_step(cusp, box, fluid_side, "path");
        b = transient_step(cusp, box, fluid_side, "face");
        fprintf("  transient, fluid %s (newmark / Crank-Nicolson, h = 0.1, dt = 0.01): e(u_F) %.3e, e(u_S) %.3e\n", ...
            fluid_side, a.e_F, a.e_S);
        finite = isfinite(a.e_F) && isfinite(a.e_S);
        same_state = isequal(a.state0, b.state0);
        same_step = isequal(a.x, b.x);
    catch err
        fprintf(2, "  transient, fluid %s: %s\n", fluid_side, err.message);
        [finite, same_state, same_step] = deal(false);
    end
    n_failed = n_failed + check("transient, fluid " + fluid_side + ": step runs and stays finite", finite);
    n_failed = n_failed + check("transient, fluid " + fluid_side + ": t = 0 state independent of collapsed_normal", same_state);
    n_failed = n_failed + check("transient, fluid " + fluid_side + ": step independent of collapsed_normal", same_step);
end

fprintf("\n%d check(s) failed.\n", n_failed);

%=========================================================================%
function out = airfoil_solve(regularized, h, k, fluid_side, varargin)
% Steady manufactured FSI on the interior airfoil pair (smooth fields,
% E = 1, nu = 0.3333, mu = 1, pinned pressure, jacobian patch weight);
% varargin goes to interface_problem, except meshes = {fluid, solid}.
curve = AirfoilCurve(-1/10, 1/10, regularized);
box = [-3 3 -0.4 0.6];
m = find(cellfun(@(a) isequal(a, "meshes"), varargin), 1);
if isempty(m)
    [fm, sm] = interface_meshes(curve, "interior", h, box = box, fluid_side = fluid_side);
else
    [fm, sm] = varargin{m + 1}{:};
    varargin(m:m + 1) = [];
end
E = 1.0; nu = 0.3333; mu = 1.0;
pd = steady_manufactured_data(curve, box, "sin(pi*x)*sin(pi*y)", ...
    "exp(-(x - 0.5)^2 - (y - 0.5)^2)", ["x*sin(y)"; "y*cos(x)"], mu, E, nu, 1.0, 0.0, 0.0);
fp = IncompressibleElasticityProblem(k, time_step = 1.0, source_term = pd.f_F, ...
    datum_dirichlet = pd.u_F, viscosity = mu, density = 0.0, ...
    domain_measure = pd.measure_F, pressure_mean = "pinned");
sp = LinearElasticityProblem(k, time_step = 1.0, time = Inf, source_term = pd.f_S, ...
    datum_dirichlet = pd.u_S, density = 0.0, hooke_modulus = 0, ...
    young_modulus = E, poisson_ratio = nu);
coupled = FluidStructureCoupling.interface_problem(fp, sp, fm, sm, curve, ...
    pd.T_1, pd.T_2, "steady", "patch_weight", "jacobian", varargin{:});
[fs, ss, gs] = solve(coupled, k);
eF = fs.compute_l2_error({pd.sigma_F, pd.u_F, pd.gamma_F});
eS = ss.compute_l2_error({pd.sigma_S, pd.u_S, pd.gamma_S});
out = struct("gs", gs, "eF", eF(:).', "eS", eS(:).');
end

function out = transient_step(curve, box, fluid_side, collapsed_normal)
% A time-independent solution (zero reaction in the sources; T2 = u_F,
% since u_F = du_S/dt + T2 with u_S constant) started from the exact
% fields: newmark solid, Crank-Nicolson fluid, one step, collapse with the
% rows the collapsed side requires. Returns the t = 0 state, the stacked
% field dofs of the step and the errors.
k = 1; dt = 1e-2; h = 0.1;
[fm, sm] = interface_meshes(curve, "interior", h, box = box, fluid_side = fluid_side);
zero = {@(x, y) 0*x; @(x, y) 0*x};
pd = steady_manufactured_data(curve, box, "sin(pi*x)*sin(pi*y)", ...
    "exp(-(x - 0.5)^2 - (y - 0.5)^2)", ["x*sin(y)"; "y*cos(x)"], 1.0, 1.0, 0.3333, dt, 0.0, 0.0);
fp = IncompressibleElasticityProblem(k, time_step = dt, source_term = pd.f_F, ...
    datum_dirichlet = pd.u_F, viscosity = 1.0, density = 1.0, ...
    domain_measure = pd.measure_F, pressure_mean = "pinned", ...
    initial_velocity = pd.u_F, datum_dirichlet_rate = zero);
sp = LinearElasticityProblem(k, time_step = dt, source_term = pd.f_S, ...
    datum_dirichlet = pd.u_S, density = 1.0, hooke_modulus = 0, ...
    young_modulus = 1.0, poisson_ratio = 0.3333, ...
    initial_displacement = pd.u_S, initial_velocity = zero);
sp.Data.time_scheme = "newmark";
fp.Data.time_scheme = "cranknicolson";
% the cusp's coarse initial data trip the traction compatibility warning
% under either mapping; the identity checks do not depend on it
warning_state = warning("off", "FluidStructureCoupling:incompatibleInitialData");
restore = onCleanup(@() warning(warning_state));
mapping = {"tip_mapping", "collapse", "collapsed_normal", collapsed_normal};
state = FluidStructureCoupling.initial_state(fp, sp, fm, sm, curve, pd.T_1, pd.u_F, mapping{:});
[fs, ss] = FluidStructureCoupling.solve(fp, sp, fm, sm, curve, pd.T_1, pd.u_F, "newmark", ...
    "patch_weight", "jacobian", "quiet", true, "state", state, mapping{:}, ...
    "swap_transmission", fluid_side == "outer");
eF = fs.compute_l2_error({pd.sigma_F, pd.u_F, pd.gamma_F});
eS = ss.compute_l2_error({pd.sigma_S, pd.u_S, pd.gamma_S});
dofs = @(sol) [arrayfun(@(f) {f.dof(:)}, sol.Unknowns.field), arrayfun(@(f) {f.dof(:)}, sol.Unknowns.hybrid)];
out = struct("state0", state, "x", {[dofs(fs), dofs(ss)]}, "e_F", eF(2), "e_S", eS(2));
end

function sm = second_mesh(curve, h, box, fluid_side)
[~, sm] = interface_meshes(curve, "interior", h, box = box, fluid_side = fluid_side);
end

function id = error_id(f)
id = "";
try
    f();
catch err
    id = string(err.identifier);
end
end

function failed = check(name, condition)
failed = ~condition;
if failed
    fprintf(2, "  FAIL  %s\n", name);
else
    fprintf("  ok    %s\n", name);
end
end
