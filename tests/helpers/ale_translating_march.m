function out = ale_translating_march(k, fluid_mesh, solid_mesh, amplitude, dt, n_steps, ...
    uu_fluid, p_fluid, u_solid, options)
%ALE_TRANSLATING_MARCH  FSI march with a translating circle on moving meshes.
%
%   out = ale_translating_march(k, fluid_mesh, solid_mesh, A, dt, n_steps, ...
%       uu_fluid, p_fluid, u_solid, delta = [], ...)
%
%   Sigma(t) = Sigma_0 + d(t) e_x, d(t) = A sin(2 pi t), Sigma_0 the circle
%   of radius 0.4 at (0.5, 0.5); data from translating_circle_problems.
%   At t_n the meshes are, from their t = 0 vertices X,
%     solid  X + d(t_n) e_x                (rigid: the solid's d_t at fixed
%                                            X is the material D_t of the data)
%     fluid  X + d(t_n) phi(X_1) e_x       (ALE)
%   with phi = 0 at X_1 = 0, 1 and phi = 1 on [delta, 1 - delta] (quintic
%   C^2 ramps; delta = [] takes the widest ramps that keep phi = 1 on every
%   vertex off the walls x = 0, 1 of the fluid elements touching Sigma_h, so
%   that layer translates rigidly with Sigma (letting it deform, phi < 1
%   at its outer vertices, made the coarsest interpolating march blow up);
%   max phi' = 1.875/delta, so A ~ delta/1.875 and beyond can invert
%   elements, an error), so the fluid interface faces
%   translate with Sigma and the
%   transfer paths keep their t = 0 lengths; the walls x = 0, 1 stay put
%   and y = 0, 1 slide along themselves. The dofs are carried unchanged
%   (fixed reference points), so the fluid step solves the ALE form
%   rho d_t|_X u + rho (beta.grad) u - div sigma = f with beta = -w_h,
%   w_h = d'(t_n) times the P1 interpolant of phi on the current mesh,
%   the mesh velocity of the affine elements (Data.convective_velocity as
%   P1 Dubiner dofs; convective_form = "advective", since div w_h ~= 0 in
%   the ramps; convective_stabilization = false so the fluid keeps the
%   mesh tau the transmission rows use). The fluid is backward Euler
%   on the first step and Crank-Nicolson after (its a0 history is the ALE
%   derivative of the previous step, as solved on the previous mesh); the
%   solid is Newmark with its state carried unchanged. A = 0 skips beta,
%   so the march is the fixed-mesh one.
%
%   out: errors at t = n_steps*dt (eF, eS: compute_l2_error of
%   (sigma, u, gamma) on the moved meshes), the transfer path lengths
%   [min, max] over the run, delta, the last step's solutions, and per
%   step (history) the largest |u_h| dof of the fluid and of the solid and
%   the smallest fluid element area ratio |K(t_n)|/|K(0)|.
%
%See also: translating_circle_problems, FluidStructureCoupling.solve
arguments
    k, fluid_mesh, solid_mesh, amplitude, dt, n_steps, uu_fluid, p_fluid, u_solid
    options.delta double {mustBeScalarOrEmpty} = []
    options.rho_f = 1; options.mu = 1; options.rho_s = 1; options.E = 10; options.nu = 0.3;
    options.patch_weight (1,1) string = "jacobian"
    % k + 1 points per patch under-integrate the interface terms when a
    % patch is a whole face (conforming interpolating pairs): the coarsest
    % square_circle_fitted_none pair then has unstable modes
    options.n_quadpoints double {mustBeScalarOrEmpty} = []
    options.verbose (1,1) logical = false
    % Sigma follows the discrete solid: the meshes and Sigma_h are placed
    % at a predicted translation instead of d(t) (see below)
    options.follow_solid (1,1) logical = false
    % push the fluid interface vertices radially out so every fluid
    % interface chord lies outside Sigma_0 (Omega_F,h inside Omega_F): on
    % interpolating pairs the chords lie inside Sigma, the fluid paths point
    % into their own elements and the kinematic row acts as a wrong-sign
    % Robin condition (an unstable mode on the coarsest pair at k >= 2)
    options.fluid_inside (1,1) logical = false
end
curve_0 = CircularCurve([0.5; 0.5], 0.4);
data = @(t) translating_circle_problems(k, curve_0, amplitude, t, dt, uu_fluid, p_fluid, ...
    u_solid, options.rho_f, options.mu, options.rho_s, options.E, options.nu);
d = @(t) amplitude*sin(2*pi*t);
d_prime = @(t) 2*pi*amplitude*cos(2*pi*t);

if options.fluid_inside
    fluid_mesh = outside_chords(fluid_mesh, curve_0);
end
X = fluid_mesh.coordinates;
interface_vertices = unique(fluid_mesh.faces(1:2, fluid_mesh.intfaces));
delta = options.delta;
if isempty(delta)
    layer = unique(fluid_mesh.elements(:, any(ismember(fluid_mesh.boundaries, fluid_mesh.intfaces), 1)));
    x = X(1, layer);
    x = x(x > 1e-12 & x < 1 - 1e-12);
    delta = min(min(x, 1 - x));
end
phi = ramp(X(1, :), delta);
assert(all(phi(interface_vertices) == 1), "ale_translating_march:delta", ...
    "phi < 1 at a fluid interface vertex: lower delta.");
phi_dofs = p1_dofs(fluid_mesh, phi);

[fp, sp, src, ~, curve] = data(0);
sp.Data.time_scheme = "newmark";
fp.Data.time_scheme = "euler";
quad = {"n_quadpoints", options.n_quadpoints, "n_quadpoints_path", options.n_quadpoints + 1};
if isempty(options.n_quadpoints), quad = {"n_quadpoints", k + 2, "n_quadpoints_path", k + 3}; end
state = FluidStructureCoupling.initial_state(fp, sp, fluid_mesh, solid_mesh, curve, src{:}, "Order", k, quad{:});
paths = [Inf, 0];
% per step: max |u_h| dof of the fluid and the solid, min |K(t)|/|K(0)|
history = zeros(n_steps, 4);
% follow_solid: the frame (d_n, w_n, a_n) is Newmark-integrated from a
% predicted acceleration a_n = 2 b'' (t_{n-1}) - b''(t_{n-2}), b'' the
% mean over Sigma of the solid's acceleration (frame + relative a_S,
% extrapolated to the curve points); the solid source gets the frame
% term rho_S (d'' - a_n), so a frame error is an exact rigid shift of the
% discrete solid and b'' is independent of the frame. Exact at t = 0.
d_n = d(0); w_n = d_prime(0); a_n = -(2*pi)^2*d(0);
b_ddot = a_n; b_ddot_old = a_n;
for n = 1:n_steps
    t = n*dt;
    if options.follow_solid
        a_old = a_n;
        a_n = 2*b_ddot - b_ddot_old;
        d_n = d_n + dt*w_n + dt^2/4*(a_old + a_n);
        w_n = w_n + dt/2*(a_old + a_n);
    else
        d_n = d(t);
        w_n = d_prime(t);
    end
    fm = with_coordinates(fluid_mesh, X + [d_n*phi; 0*phi]);
    assert(all(sign(fm.volume) == sign(fluid_mesh.volume)), "ale_translating_march:inverted", ...
        "The fluid mesh inverts at t = %g: lower A or widen the ramps (delta = %g).", t, delta);
    sm = solid_mesh + [d_n; 0];
    [fp, sp, src, pd, curve] = data(t);
    if options.follow_solid
        % Sigma_h(t_n) = Sigma_0 + d_n; the kinematic datum carries the
        % frame velocity w_n in place of d'(t) (the solid's rigid velocity)
        curve = CircularCurve(curve_0.center + [d_n; 0], curve_0.radius);
        T2 = src{2}{1};
        src{2}{1} = @(x, y) T2(x, y) + (w_n - d_prime(t));
        f = sp.Data.source_term{1};
        sp.Data.source_term{1} = @(x, y) f(x, y) + options.rho_s*(-(2*pi)^2*d(t) - a_n);
    end
    sp.Data.time_scheme = "newmark";
    fp.Data.time_scheme = "cranknicolson";
    if n == 1, fp.Data.time_scheme = "euler"; end
    if amplitude ~= 0 || options.follow_solid
        % the P1 interpolant is invariant under the motion: same dofs
        fp.Data.convective_velocity = {-w_n*phi_dofs; 0*phi_dofs};
        fp.Data.convective_stabilization = false;
        fp.Data.convective_form = "advective";
    end
    [fs, ss, state, gs] = FluidStructureCoupling.solve(fp, sp, fm, sm, curve, src{:}, ...
        "newmark", "quiet", true, "state", state, "patch_weight", options.patch_weight, quad{:});
    for part = {gs.partition_1, gs.partition_2}
        l = vecnorm(part{1}.quadrature.interface - part{1}.quadrature.physical, 2, 1);
        paths = [min(paths(1), min(l(:))), max(paths(2), max(l(:)))];
    end
    history(n, :) = [max(abs(fs.Unknowns.field(2).dof), [], "all"), ...
        max(abs(ss.Unknowns.field(2).dof), [], "all"), min(fm.volume./fluid_mesh.volume), d_n];
    if options.follow_solid
        b_ddot_old = b_ddot;
        b_ddot = a_n + interface_mean(state.solid.acceleration, gs.partition_2, k);
    end
    if options.verbose
        eF = fs.compute_l2_error({pd.sigma_F, pd.u_F, pd.gamma_F});
        eS = ss.compute_l2_error({pd.sigma_S, pd.u_S, pd.gamma_S});
        fprintf("  step %d t %.4f  uF %.3e  uS %.3e\n", n, t, eF(2), eS(2));
    end
end
eF = fs.compute_l2_error({pd.sigma_F, pd.u_F, pd.gamma_F});
eS = ss.compute_l2_error({pd.sigma_S, pd.u_S, pd.gamma_S});
out = struct("eF", eF(:).', "eS", eS(:).', "paths", paths, "fs", fs, "ss", ss, ...
    "h", max([fluid_mesh.area(:); solid_mesh.area(:)]), "delta", delta, ...
    "history", history);
end

function m = interface_mean(dofs, partition, k)
% mean x-component over Sigma of the solid field (Dubiner dofs, the
% component blocks stacked) extrapolated to the partition's curve points
q = partition.quadrature;
[n_q, n_patches] = deal(size(q.interface_2d, 2), size(q.interface_2d, 3));
[~, w] = lgwt(n_q, 0, 1);
weights = w(:) .* reshape(partition.area, 1, []);   % n_q x n_patches
d2 = nchoosek(k + 2, 2);
u = zeros(n_q, n_patches);
for j = 1:n_patches
    u(:, j) = dubiner_basis(k, q.interface_2d(:, :, j)) * dofs(1:d2, partition.elements(j));
end
m = sum(weights .* u, "all") / sum(weights, "all");
end

function phi = ramp(x, delta)
% 1 on [delta, 1 - delta], 0 at 0 and 1, quintic C^2 ramps in between
s = min(1, min(x, 1 - x)/delta);
phi = s.^3.*(10 - 15*s + 6*s.^2);
end

function dofs = p1_dofs(mesh, v)
% Dubiner dofs (3 x N) of the P1 interpolant of the vertex values v;
% element vertex 1, 2, 3 is the reference point (1, 0), (0, 1), (0, 0)
[rs, w] = get_2d_reference_quadrature(1);
p = dubiner_basis(1, rs);
lambda = [rs(1, :); rs(2, :); 1 - rs(1, :) - rs(2, :)].';
values = lambda*reshape(v(mesh.elements), 3, []);
dofs = (p.'*(w.*p)) \ (p.'*(w.*values));
end

function mesh = outside_chords(mesh, circle)
% each interface vertex at radius R/cos(a/2), a its widest adjacent chord
% angle: that chord is then tangent to the circle at its midpoint and the
% others lie outside it
X = mesh.coordinates;
F = mesh.faces(1:2, mesh.intfaces);
for v = unique(F(:)).'
    a = 0;
    for f = find(any(F == v, 1))
        r = X(:, F(:, f)) - circle.center;
        a = max(a, abs(angle(complex(r(1, 1), r(2, 1))/complex(r(1, 2), r(2, 2)))));
    end
    X(:, v) = circle.center + (X(:, v) - circle.center)/norm(X(:, v) - circle.center)*circle.radius/cos(a/2);
end
mesh = with_coordinates(mesh, X);
end

function mesh = with_coordinates(mesh, X)
% the same mesh (connectivity, tags, tau) on the vertices X, with its
% geometry recomputed as the SimplicialGrid constructor does
mesh.coordinates = X;
x_K = reshape(X(1, mesh.elements), 3, []);
y_K = reshape(X(2, mesh.elements), 3, []);
mesh.volume = .5*sum(x_K.*circshift(y_K, -1) - circshift(x_K, -1).*y_K);
dx_e = diff(reshape(X(1, mesh.faces(1:2, :)), 2, []));
dy_e = diff(reshape(X(2, mesh.faces(1:2, :)), 2, []));
mesh.area = sqrt(dx_e.^2 + dy_e.^2);
n_by_face = [dy_e; -dx_e];
mesh.normals = reshape(n_by_face(:, mesh.boundaries), 2*3, []);
end
