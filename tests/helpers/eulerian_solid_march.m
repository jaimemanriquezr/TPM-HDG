function out = eulerian_solid_march(k, background, amplitude, dt, n_steps, uu_fluid, p_fluid, u_solid, options)
%EULERIAN_SOLID_MARCH  FSI march with fluid and solid both carved from one background grid.
%
%   out = eulerian_solid_march(k, background, A, dt, n_steps, uu_fluid, p_fluid, u_solid, ...
%       geometry = "flat", snap = true, fluid_scheme = "bdf2", ...)
%
%   Sigma(t) is prescribed: {y = y0 + d(t)} (geometry = "flat", fluid
%   above, translating_flat_problems) or the circle of the given radius
%   centred at center + d(t) e_x (geometry = "circle", fluid outside,
%   translating_circle_problems), d(t) = A sin(2 pi t). At every t_n both
%   meshes are the two sides of SimplicialGrid.split of the background
%   (cut elements discarded; the solid from solid_background when given),
%   their interface vertices moved onto Sigma(t_n) when snap = true
%   (snap_to_curve, normal_sign = -1 for the solid). Both histories live on
%   the background elements, identified across steps by the centroids
%   (and, for the solid skeleton, the face midpoints) of the unsnapped
%   carve.
%   Fluid: backward Euler on the first step, then BDF2 (fluid_scheme =
%   "bdf2", written as a theta step) or Crank-Nicolson; its velocity and
%   acceleration are remapped as in eulerian_carved_march.
%   Solid: Newmark along the characteristics of the known rigid velocity
%   d'(t) (the displacement measured from the reference translated by
%   d(t)): by default (transport = "translate", history_translate =
%   "direct"; see transport for the defaults under "foot") every volume level of
%   the history is the exact projection onto the new carve of the old field
%   at x - shift, shift = d(t_n) - d(t_n-1); the skeleton levels are
%   translated on the old carve (LinearElasticityProblem.translate_state)
%   and remapped onto the new carve (remap_solid_state), an interface face
%   that is an old face translated by the shift (snapped faces of the flat
%   interface) taking the old face's history unchanged, and the
%   displacement level on the boundary faces of the solid carve is then the
%   trace of the translated u* (skeleton_history = "ustar_displacement").
%   The velocity of the elements that switch on or change shape is the
%   two-step form of u* (history_velocity = "ustar_twostep"), the volume
%   acceleration of the other elements the consistent acceleration
%   (history_acceleration = "consistent_kept"), and the step is the
%   trapezoidal Newmark scheme (solid_bossak_rho = 1).
%   history_velocity: on the elements that switch on or change shape
%   (moved by remap_solid_state), "newmark" keeps the remapped volume
%   velocity and acceleration. The remapped v_h, a_h are best
%   approximations on the old shape, while the displacement history of a
%   moved element is its P_k+1 postprocessing u* (remap_solid_state), and
%   the trapezoidal Newmark carries the difference undamped. "ustar"
%   rebuilds them from the u* of the last displacements, carried
%   along the characteristics (LinearElasticityProblem.translated_field)
%   and re-expressed on the new carve: for dt >= h/4 (h the largest face
%   of the grid the solid is carved from), from the fourth step on, v and a are the one-sided
%   differences of order 2 of the last four u*,
%     v = (3 u*_1 - 4 u*_2 + u*_3)/(2 dt),
%     a = (2 u*_1 - 5 u*_2 + 4 u*_3 - u*_4)/dt^2;
%   otherwise (and in the second and third steps) v is the velocity of the
%   two-step form of the trapezoidal Newmark scheme,
%     v = (u*_1 - u*_2)/dt + dt (a^n + a^(n-1))/4,
%   with a^n, a^(n-1) the last two remapped Newmark accelerations, and a is
%   kept; in the first step the remapped fields are kept. The differences
%   of u* are projected onto P_k. With the trapezoidal update the
%   one-sided differences amplify the modes with omega dt of order one
%   (the largest root of the recursion of one element is 1.07 per step)
%   and their acceleration amplifies the step-to-step noise of u* by
%   dt^(-2): the solid error grows as dt -> 0 at fixed h. The two-step
%   form keeps the scheme's neutral stability and reads the displacement
%   history as the solve does, but is less accurate for dt of order h.
%   "ustar_bdf": the one-sided differences of order 2 at every dt, from the
%   fourth step on; before that the remapped fields are kept.
%   "ustar_bdf_all" (transport = "foot" only): the "ustar_bdf" rule on
%   every solid element.
%   "ustar_twostep" (default), "ustar_twostep_low", "ustar_twostep_mean": the
%   two-step form at every dt, from the second step on,
%     v = Pi_k(u*_1 - u*_2)/dt + dt (F a^n + F a^(n-1))/4,  a = F a^n,
%   with F the element-wise L2 projection (on the new shape) onto P_k
%   (the identity: a is kept), P_k-1 and P_0 respectively; a^(n-1) is
%   carried unprojected.
%   history_acceleration: the volume acceleration history. "translated"
%   keeps the translated and remapped Newmark acceleration. The
%   Newmark step leaves its acceleration in the element momentum identity
%     rho a = f + div sigma - S(u - u_hat)
%   (LinearElasticityProblem.consistent_acceleration), which the
%   separately translated a no longer satisfies for the translated
%   (sigma, u, u_hat); the defect is carried undamped and grows. The
%   consistent variants replace a by the consistent acceleration of the
%   translated and remapped state, with the previous step's source along
%   the characteristic, f(x - shift, t_n-1): "consistent" on every
%   element, before the history_velocity rule (which then reads it as
%   a^n); "consistent_kept" only on the elements that neither switch on
%   nor change shape (default); "consistent_after" on every element, after
%   the rule (the rule's velocity is kept). The skeleton, stress and rotation
%   acceleration levels are left as they are. "balanced" and
%   "balanced_after" (before and after the rule, as the consistent ones)
%   replace a by the globally consistent acceleration M a = f - K u of the
%   translated volume displacement u (LinearElasticityProblem.
%   balanced_acceleration, the same source along the characteristic): its
%   stress, rotation and interior traces are the flux-balanced lift of u
%   with the translated and remapped displacement_skeleton history as the
%   trace on the boundary faces of the solid carve, the trace the previous
%   coupled solve produced; "balanced_levels" is "balanced" with the
%   displacement level's stress, rotation and interior skeleton history
%   replaced by that lift as well.
%   solid_bossak_rho: the spectral radius at infinity rho_inf of the
%   Bossak (WBZ-alpha) solid step
%     M[(1 - am) a^(n+1) + am a^n] + K u^(n+1) = f^(n+1),
%   am = (rho_inf - 1)/(rho_inf + 1), gamma = 1/2 - am, beta = (1 - am)^2/4
%   (1, the default: the trapezoidal Newmark scheme). It is solved as a
%   Newmark step with the density (1 - am) rho and the volume history
%   acceleration a^n - am a^n/((1 - am) c), c = (1/2 - beta)/beta, and the
%   solid levels are then advanced from the true history
%   (get_kinematic_state with beta, gamma). The two-step history rules
%   keep their trapezoidal form.
%   history_translate: "library" carries the history by
%   translate_state on the old carve followed by remap_solid_state, which
%   extrapolates twice: where x - shift leaves the old carve (trailing
%   side) and onto the elements that switch on (leading side). "direct"
%   replaces every volume level (displacement, velocity, acceleration and
%   the three stress and three rotation levels) by the exact projection
%   onto the new carve of the old field at x - shift
%   (direct_translate_state), and carries the u* history of the rules and
%   the previous acceleration the same way
%   (LinearElasticityProblem.translated_field_between); the skeleton
%   levels keep the library path (remap_solid_state still runs in full;
%   its volume levels are then overwritten). history_acceleration and
%   history_velocity then act as with "library". "direct" is the default
%   under transport = "translate" (see transport).
%   initial_perturbation = [seed, amplitude] adds to each volume level
%   (displacement, velocity, acceleration) of the initial solid state a
%   random field (randn with the given seed) of amplitude times that
%   level's dof norm; two marches that differ only by it differ by the
%   homogeneous step map applied to the perturbation. record_solid keeps
%   per step the solid volume dofs after the solve and the moved elements
%   in out.solid_record.
%   history_velocity_scope: the elements the u* rules act on, "moved"
%   (default: every element that switches on or changes shape),
%   "reshaped" (present in both carves, changed shape) or "on" (switched
%   on); the others keep the translated velocity and acceleration.
%   balanced_trace: the boundary trace of the balanced lift, "history"
%   (default: the translated and remapped displacement_skeleton history),
%   "volume" (the face projection of the trace of the translated volume
%   displacement of the owning element) or "ustar" (the same for the u*
%   of the old state carried along the characteristic, P_k+1).
%   record_solid also keeps the history the solve read (volume levels and
%   displacement_skeleton) and the solid mesh of every step.
%   skeleton_history: "library" keeps the translated and
%   remapped skeleton levels; "volume" and "ustar" replace them on the
%   boundary faces of the solid carve by the face projections of the
%   traces of the volume velocity, acceleration and displacement (or of
%   the translated u* for the displacement) of the owning element;
%   "ustar_displacement" (default under transport = "translate"; under
%   "foot" the default is "library", see transport) replaces only the
%   displacement level as
%   "ustar" does, "rates" only the velocity and acceleration levels.
%   slid_subdivisions (default 1), slid_correction (default true) and
%   slid_displacement (default false): passed to remap_solid_state. skeleton_history "rates_velocity" and
%   "rates_acceleration" replace only the velocity, respectively the
%   acceleration, level as "rates" does.
%   swap_transmission (default true): see Coupling. record_solid also
%   keeps the fluid velocity dofs of every step.
%   transport = "foot" replaces the translate and the remap by the
%   semi-Lagrangian transport_solid_state (projections of the history at
%   the feet of the step; the u* history by foot_projection), for the
%   reference moving rigidly with phi_t(X) = c(t) + Q(theta(t)) (X - c0),
%   c(t) = c0 + d(t) e, theta(t) = Theta sin(2 pi t) (rotation = Theta,
%   circle only; Theta = 0 on the flat), whose velocity field is beta(x, t)
%   = d'(t) e + theta'(t) J (x - c(t)), J = [0 -1; 1 0]; the data are
%   translating_circle_problems' with Theta. The foot y of a point x of
%   the step t_n -> t_n+1, the position at t_n of the reference point at
%   x at t_n+1:
%     foot_rule = "exact"  y = c(t_n) + Q(theta_n - theta_n+1)(x - c(t_n+1));
%       "midpoint" (default)  x_m = x - (dt/2) w(x), y = x - dt w(x_m);
%       "arrival"  y = x - dt w(x), locally O(dt^2) only (for a rotation
%       the radius drifts);
%     foot_velocity, the field w of the two last rules: "extrapolated"
%       (default) (3 beta(t_n) - beta(t_n-1))/2, the extrapolation to the
%       midpoint; "extrapolated_end" 2 beta(t_n) - beta(t_n-1), to t_n+1;
%       "prescribed" beta(t_n+1/2). beta(t_n-1) is the prescribed field,
%       also at the first step (t_-1 = -dt).
%   Every foot is affine in x here; its Jacobian is passed as a handle.
%   quadrature_degree and foot_integration ("clipped", default: exact
%   integration on the overlaps of the foot images with the old elements;
%   "quadrature": a rule on each new element) are transport_solid_state's
%   ([] = 2k + 2), also used for the u* history (foot_projection). transport
%   = "translate" (default) is the rule above and needs rotation = 0.
%   Under transport = "foot" the history options keep their meaning with
%   the translate replaced by the feet: the u* history and the previous
%   acceleration of the two-step rules are carried by foot_projection, the
%   consistent accelerations use the previous step's source at the feet,
%   f(foot(x), t_n-1), in place of f(x - shift, t_n-1), and the Bossak
%   step is unchanged (its a^(n-1) is carried by foot_projection as well).
%   history_translate and skeleton_history default to "direct" and
%   "ustar_displacement" under transport = "translate" and to "library"
%   under "foot"; history_translate = "library" is the only value there
%   (the transport already projects every volume level at the feet), and
%   skeleton_history takes its other values with the u* trace carried by
%   foot_projection of the old u* in place of x - shift. The balanced
%   accelerations and interface_feedback are not available under "foot".
%   Of the default history, history_velocity = "ustar_twostep",
%   history_acceleration = "consistent_kept" and solid_bossak_rho = 1 act
%   under "foot" as under "translate"; the direct translate, the u*
%   displacement trace and slid_subdivisions (remap_solid_state is not
%   called) do not. The transport carries the vector components
%   unrotated, (u o foot)(x), and div C sym grad (u o foot) differs from
%   (div C sym grad u) o foot when grad foot is a rotation (the grad div
%   part does not commute), so the consistent acceleration of the
%   transported stress level is not the material acceleration at the
%   feet. consistent_flux (transport = "foot" only) chooses the stress the
%   consistent accelerations read: "stress" ("", the default) the
%   transported stress level, with history_acceleration ignored when
%   rotation ~= 0 and the transported acceleration kept ("translated");
%   "transported" the transported momentum flux of transport_solid_state,
%   the PEERS projection of sigma(foot(x)) grad_foot^-T of the old stress,
%   for which (div sigma)(foot(x)) = div(sigma(foot(x)) grad_foot^-T)
%   (affine foot), and history_acceleration acts at every rotation (the
%   state the solve reads keeps the transported stress level).
%   history_velocity = "ustar_bdf",
%   history_acceleration = "translated" and solid_bossak_rho = 1 give the
%   march with the history rule of order 2 alone.
%   Coupling: the kinematic row is tested with the fluid's interface traces
%   and the traction row with the solid's (swap_transmission), because the
%   two snapped polygons do not match on a curved interface.
%
%   out: errors at t = n_steps dt (eF, eS: compute_l2_error of
%   (sigma, u, gamma)); paths, the transfer path lengths [min, max] over
%   the run; path_ratio, per step the largest path over the element height
%   normal to its face, [fluid, solid]; history, per step [t, max|u_h|
%   fluid, max|u_h| solid, fluid on, fluid off, solid on, solid off];
%   errors, per step the largest of eF and eS (verbose only, else NaN);
%   step_errors, per step the (sigma, u) errors [fluid, solid], i.e.
%   [eF(1), eF(2), eS(1), eS(2)] (verbose only, else NaN);
%   the final fs, ss, fluid_mesh, solid_mesh, state (the coupled state)
%   and keys (of the final carve, see carve_both); h, the largest face.
%
%   measure_interface = true (no effect on the solve) adds out.interface,
%   (n_steps + 1) x 9, one row per t_n, n = 0 (the initial state), ...,
%   n_steps: the interface rebuilt from the solved solid displacement on
%   Sigma_0 (eulerian_interface_displacement, reference = the curve at
%   d = 0, rigid = d(t_n)), measured by interface_curve_metrics against
%   the exact set Sigma_0 + d(t_n) and the exact material points x(s) =
%   X(s) + d(t_n) + u_S(x(s)) (fixed point with pd.u_S). With rotation, X(s)
%   is the reference rotated by theta(t_n) about its centre, so that X(s) +
%   d(t_n) = phi_t(X_0(s)) and the rebuilt interface solves x - u*(x) =
%   phi_t(X_0(s)). interface_reference = "initial" (default) is this
%   construction; "preimage" rebuilds the interface instead on the
%   material pre-image of Sigma(t_n), Y(s) = x(s) - d(t_n) - u_S(x(s)),
%   x(s) = X(s) + d(t_n), whose exact material points are the x(s), on
%   Sigma(t_n): with "initial" the exact points lie off Sigma(t_n) by the
%   normal part of u_S and u* is extrapolated that far outside its
%   elements; with "preimage" the distance (columns 2 and 8) is an error
%   as well. Columns:
%     1 t; from the P_k+1 postprocessing u* (postprocessed_displacement),
%     continuous = true: 2 distance to the exact set, 3 point error,
%     4 measure ratio, 5 monotonicity; from u* uncorrected: 6 monotonicity,
%     7 max jump; from u_h (degree k), continuous = true: 8 distance,
%     9 point error.
%
%   interface_feedback = true carves and snaps both meshes, at every t_n,
%   against a curve built from the solid state instead of the prescribed
%   Sigma(t_n) (the solid history is still translated by the known shift
%   d(t_n) - d(t_n-1)). The initial meshes are carved against the exact
%   Sigma(0), where the initial state is built. Step t_n-1 -> t_n:
%     predicted curve  on the carve of t_n-1, whose interface faces are
%       mapped to the curve it was carved against, Sigma_h^n-1: x_h(s)
%       solves x - d(t_n-1) - u*(x) = X(s), u* the P_k+1 postprocessing
%       (eulerian_interface_displacement, continuous = true); predictor =
%         "none" (default)  x_p = x_h: no lag in the normal direction
%           when the normal motion of the interface is the rigid one, as
%           in the manufactured data of this march; with a normal elastic
%           motion its lag is of order dt, and the choice of predictor
%           must be revisited on data that have one;
%         "pointwise"  x_p(s) = x_h(s) + dt v_h(x_h(s)), v_h
%           the solid's volume velocity on the element of the piece of s,
%           made continuous at the breaks as continuous = true does;
%         "velocity"  x_p(s) from the Newmark predictor p = u + dt v of
%           displacement, stress and rotation and its postprocessing p*
%           in place of u* (the accelerations of stress and rotation come
%           from the Newmark recursion alone, which no equation pins, and
%           the remap feeds an alternating mode of their rates: unstable
%           in feedback, e.g. k = 1, N = 16);
%       carried along the characteristic, x*(s) = x_p(s) + shift: the
%       DeformedCurve over Sigma_0 of D(s) + shift. On the flat (open)
%       curve D_1 is corrected by the linear interpolant of its end values,
%       so that the ends lie on the walls x = 0, 1 (the exact tangential
%       displacement vanishes there), and the curve is a WallEndedCurve;
%     data  the traction datum T1 on the predicted curve, with its normal
%       (eulerian_interface_traction); T2 is pointwise and unchanged;
%     carve, remap, solve on the predicted curve; the solid history's
%       interface faces are carried by remap_solid_state's rule with
%       carry_tolerance (default 0.25 under feedback on the flat, where
%       the snapped faces are near translates of the old ones; 0 on the
%       circle, where they are not even under the prescribed Sigma and a
%       loose match copies traces onto faces a fraction of h away);
%     solved curve  Sigma_h^n rebuilt from the solved solid on the new
%       carve (u*, rigid d(t_n), continuous = true), measured.
%   No geometry iterations. out.feedback, (n_steps + 1) x 14, one row per
%   t_n (row 1, t = 0: the curve rebuilt from the initial state, predicted
%   columns NaN). Each curve is sampled as interface_curve_metrics does,
%   on the pieces of the solved curve for the lag. Columns:
%     1 t; 2 lag max_s |x*(s) - x_h(s)| between the predicted and the
%     solved curve, 3 its normal part (flat: |dy|; circle: the difference
%     of the distances to the exact centre); predicted curve: 4 distance to
%     the exact set, 5 point error (against the exact material points, as
%     in out.interface); solved curve: 6 distance, 7 point error, 8 measure
%     ratio, 9 monotonicity, 10 and 11 (flat; NaN on the circle) x_1 - 0
%     and x_1 - 1 at its ends s = 0, 1; 12 the solid skeleton faces carried
%     by the remap, 13 the largest relative mismatch of a new interface
%     face to the nearest translated old one (remap_solid_state info),
%     14 the wall time of the step (s).
%
%See also: eulerian_carved_march, remap_solid_state, translating_flat_problems
arguments
    k, background SimplicialGrid, amplitude, dt, n_steps, uu_fluid, p_fluid, u_solid
    options.geometry (1,1) string {mustBeMember(options.geometry, ["flat", "circle"])} = "flat"
    options.y0 (1,1) double = 0.45
    options.center (2,1) double = [0.5; 0.5]
    options.radius (1,1) double = 0.2
    % the grid the solid is carved from (default: background)
    options.solid_background = []
    options.snap (1,1) logical = true
    options.fluid_scheme (1,1) string {mustBeMember(options.fluid_scheme, ["bdf2", "cranknicolson"])} = "bdf2"
    options.rho_f = 1; options.mu = 1; options.rho_s = 1; options.E = 10; options.nu = 0.3;
    options.n_quadpoints double {mustBeScalarOrEmpty} = []
    % per step errors in out.errors, printed
    options.verbose (1,1) logical = false
    options.snapshot_file (1,1) string = ""
    options.snapshot_every (1,1) double {mustBeInteger, mustBePositive} = 1
    options.history_velocity (1,1) string {mustBeMember(options.history_velocity, ...
        ["newmark", "ustar", "ustar_bdf", "ustar_bdf_all", "ustar_twostep", "ustar_twostep_low", ...
        "ustar_twostep_mean"])} = "ustar_twostep"
    options.history_acceleration (1,1) string {mustBeMember(options.history_acceleration, ...
        ["translated", "consistent", "consistent_kept", "consistent_after", "balanced", ...
        "balanced_after", "balanced_levels"])} = "consistent_kept"
    options.solid_bossak_rho (1,1) double {mustBeInRange(options.solid_bossak_rho, 0, 1)} = 1
    % "" resolves to "direct" under transport = "translate", "library" under "foot"
    options.history_translate (1,1) string {mustBeMember(options.history_translate, ...
        ["", "library", "direct"])} = ""
    options.initial_perturbation double = []
    options.record_solid (1,1) logical = false
    options.history_velocity_scope (1,1) string {mustBeMember(options.history_velocity_scope, ...
        ["moved", "reshaped", "on"])} = "moved"
    options.balanced_trace (1,1) string {mustBeMember(options.balanced_trace, ...
        ["history", "ustar", "volume"])} = "history"
    options.swap_transmission (1,1) logical = true
    % "" resolves to "ustar_displacement" under transport = "translate", "library" under "foot"
    options.skeleton_history (1,1) string {mustBeMember(options.skeleton_history, ...
        ["", "library", "volume", "ustar", "ustar_displacement", "rates", "rates_velocity", ...
        "rates_acceleration"])} = ""
    options.slid_subdivisions (1,1) double {mustBePositive} = 1
    options.slid_correction (1,1) logical = true
    options.slid_displacement (1,1) logical = false
    % per step metrics of the interface rebuilt from the solid (out.interface)
    options.measure_interface (1,1) logical = false
    % the curve out.interface rebuilds on (see measure_interface in the help)
    options.interface_reference (1,1) string {mustBeMember(options.interface_reference, ...
        ["initial", "preimage"])} = "initial"
    % carve against the interface predicted from the solid state (out.feedback)
    options.interface_feedback (1,1) logical = false
    % the interface predictor of interface_feedback (see the help)
    options.predictor (1,1) string {mustBeMember(options.predictor, ["pointwise", "velocity", "none"])} = "none"
    % remap_solid_state's carry_tolerance; [] is 0.25 under interface_feedback on the flat, else 0
    options.carry_tolerance double {mustBeScalarOrEmpty} = []
    % rotation amplitude Theta of the reference about its centre (circle, transport = "foot")
    options.rotation (1,1) double = 0
    % the transport of the solid history (see the help)
    options.transport (1,1) string {mustBeMember(options.transport, ["translate", "foot"])} = "translate"
    options.foot_rule (1,1) string {mustBeMember(options.foot_rule, ["midpoint", "exact", "arrival"])} = "midpoint"
    options.foot_velocity (1,1) string {mustBeMember(options.foot_velocity, ...
        ["extrapolated", "extrapolated_end", "prescribed"])} = "extrapolated"
    % transport_solid_state's quadrature_degree ([] = its default 2k + 2)
    options.quadrature_degree double {mustBeScalarOrEmpty} = []
    % transport_solid_state's foot_integration
    options.foot_integration (1,1) string {mustBeMember(options.foot_integration, ...
        ["clipped", "quadrature"])} = "clipped"
    % the stress of the consistent accelerations under transport = "foot" ("" = "stress")
    options.consistent_flux (1,1) string {mustBeMember(options.consistent_flux, ...
        ["", "stress", "transported"])} = ""
end
by_foot = options.transport == "foot";
assert(options.rotation == 0 || (by_foot && options.geometry == "circle"), ...
    "eulerian_solid_march:rotation", "rotation ~= 0 needs geometry = ""circle"" and transport = ""foot"".");
assert(~(by_foot && options.interface_feedback), "eulerian_solid_march:feedbackFoot", ...
    "interface_feedback is not available with transport = ""foot"".");
assert(by_foot || options.history_velocity ~= "ustar_bdf_all", "eulerian_solid_march:ustarAll", ...
    "history_velocity = ""ustar_bdf_all"" needs transport = ""foot"".");
assert(~(by_foot && options.history_translate == "direct"), "eulerian_solid_march:directFoot", ...
    "history_translate = ""direct"" is not available with transport = ""foot"".");
assert(by_foot || options.consistent_flux ~= "transported", "eulerian_solid_march:fluxTranslate", ...
    "consistent_flux = ""transported"" needs transport = ""foot"".");
assert(~(by_foot && startsWith(options.history_acceleration, "balanced")), ...
    "eulerian_solid_march:balancedFoot", ...
    "the balanced accelerations are not available with transport = ""foot"".");
% the defaults of the history (see transport in the help)
history_translate = options.history_translate;
skeleton_history = options.skeleton_history;
if history_translate == "", history_translate = "direct"; end
if by_foot
    history_translate = "library";
    if skeleton_history == "", skeleton_history = "library"; end
end
if skeleton_history == "", skeleton_history = "ustar_displacement"; end
transported_flux = options.consistent_flux == "transported";
bossak_am = (options.solid_bossak_rho - 1)/(options.solid_bossak_rho + 1);
acceleration_rule = options.history_acceleration;
if by_foot && options.rotation ~= 0 && ~transported_flux
    % the transported acceleration is kept (see the help)
    acceleration_rule = "translated";
end
% beta = 1/2 (rho_inf = sqrt(2) - 1) leaves no history acceleration in the
% Newmark form of the step (c = 0 below)
assert(abs((1 - bossak_am)^2/4 - 1/2) > 1e-12, "eulerian_solid_march:bossak", ...
    "solid_bossak_rho = sqrt(2) - 1 (beta = 1/2) cannot be written as a Newmark solve.");
rule = options.history_velocity;
by_ustar = rule ~= "newmark";
% the degree F projects the accelerations of the two-step form onto
q = k - (rule == "ustar_twostep_low") - k*(rule == "ustar_twostep_mean");
physics = {options.rho_f, options.mu, options.rho_s, options.E, options.nu};
if options.geometry == "flat"
    data = @(t) translating_flat_problems(k, options.y0, amplitude, t, dt, uu_fluid, p_fluid, ...
        u_solid, physics{:});
    direction = [0; 1];
else
    curve_0 = CircularCurve(options.center, options.radius);
    if options.rotation ~= 0, physics = [physics, {options.rotation}]; end
    data = @(t) translating_circle_problems(k, curve_0, amplitude, t, dt, uu_fluid, p_fluid, ...
        u_solid, physics{:});
    direction = [1; 0];
end
theta = @(t) options.rotation*sin(2*pi*t);
d = @(t) amplitude*sin(2*pi*t);
carve = @(curve) carve_both(background, options.solid_background, curve, options.snap);
quad = {"n_quadpoints", options.n_quadpoints, "n_quadpoints_path", options.n_quadpoints + 1};
if isempty(options.n_quadpoints), quad = {"n_quadpoints", k + 2, "n_quadpoints_path", k + 3}; end

[fp, sp, src, pd, curve] = data(0);
sp.Data.time_scheme = "newmark";
fp.Data.time_scheme = "euler";
[fm, sm, keys] = carve(curve);
state = FluidStructureCoupling.initial_state(fp, sp, fm, sm, curve, src{:}, "Order", k, quad{:});
if ~isempty(options.initial_perturbation)
    stream = RandStream("mt19937ar", Seed = options.initial_perturbation(1));
    for name = ["displacement", "velocity", "acceleration"]
        D = state.solid.(name);
        R = randn(stream, size(D));
        state.solid.(name) = D + options.initial_perturbation(2)*norm(D, "fro")/norm(R, "fro")*R;
    end
end
solid_record = struct("displacement", {}, "velocity", {}, "acceleration", {}, "moved", {}, ...
    "history", {}, "mesh", {}, "fluid_velocity", {}, "velocity_skeleton", {}, ...
    "history_velocity_skeleton", {});
% the previous step's solid source, for the consistent acceleration history
f_previous = sp.Data.source_term;
% u* of the last four displacements, newest first, and the previous
% acceleration history (the remapped Newmark acceleration)
u_stars = {};
a_previous = [];
% the volume acceleration the previous solid solve read as a^(n-1) (Bossak)
a_read = [];
% h of the grid the solid is carved from
solid_grid = background;
if ~isempty(options.solid_background), solid_grid = options.solid_background; end
by_differences = ismember(rule, ["ustar_bdf", "ustar_bdf_all"]) || ...
    (rule == "ustar" && dt >= 0.25*max(solid_grid.area));
two_step = ~ismember(rule, ["ustar_bdf", "ustar_bdf_all"]);
if by_ustar, u_stars = {postprocessed_displacement(sp, sm, k, state.solid)}; end
carry = {};
if ~isempty(options.carry_tolerance)
    carry = {"carry_tolerance", options.carry_tolerance};
elseif options.interface_feedback && options.geometry == "flat"
    carry = {"carry_tolerance", 0.25};
end
if options.measure_interface || options.interface_feedback
    if options.geometry == "flat"
        reference = FlatCurveHorizontal(options.y0);
        exact_set = @(t) struct("type", "flat", "height", options.y0 + d(t));
    else
        reference = curve_0;
        exact_set = @(t) struct("type", "circle", "center", options.center + d(t)*direction, ...
            "radius", options.radius);
    end
    % the reference rotated by theta(t) about its centre (the reference
    % itself without rotation)
    reference_at = @(t) reference;
    if options.rotation ~= 0
        reference_at = @(t) rotated_circle(options.center, options.radius, theta(t));
    end
end
if options.measure_interface
    measure = @(t, sp, sm, solid, curve, pd) interface_measures(t, sp, sm, k, solid, reference_at(t), ...
        curve, d(t)*direction, exact_set(t), pd.u_S, options.interface_reference == "preimage");
    interface = nan(n_steps + 1, 9);
    interface(1, :) = measure(0, sp, sm, state.solid, curve, pd);
end
if options.interface_feedback
    feedback = nan(n_steps + 1, 14);
    [~, feedback(1, :)] = solved_curve(0, [], sp, sm, k, state.solid, reference, curve, [0; 0], ...
        exact_set(0), pd.u_S, options.center);
    feedback(1, 14) = 0;
    curve_carved = curve;
end
paths = [Inf, 0];
path_ratio = zeros(n_steps, 2);
history = zeros(n_steps, 7);
errors = nan(n_steps, 1);
step_errors = nan(n_steps, 4);
v_previous = [];
for n = 1:n_steps
    step_timer = tic;
    t = n*dt;
    [fp, sp, src, pd, curve] = data(t);
    if options.interface_feedback
        % the predicted curve from the state on the carve of t_n-1 (mapped
        % to the curve it was carved against) and T1 on it
        predicted = predicted_curve(sp, sm, k, state.solid, dt, reference, curve_carved, ...
            d(t - dt)*direction, (d(t) - d(t - dt))*direction, options.predictor);
        curve = predicted.curve;
        src{1} = eulerian_interface_traction(curve, pd.sigma_F, pd.sigma_S);
    end
    sp.Data.time_scheme = "newmark";
    if bossak_am ~= 0
        sp.Data.newmark_gamma = 1/2 - bossak_am;
        sp.Data.newmark_beta = (1 - bossak_am)^2/4;
    end
    fp.Data.time_scheme = "cranknicolson";
    if n == 1, fp.Data.time_scheme = "euler"; end
    [fm_old, sm_old, keys_old] = deal(fm, sm, keys);
    [fm, sm, keys] = carve(curve);
    [state.fluid, switched_f] = remap_fluid_state(state.fluid, fm_old, fm, k, keys_old.fluid, keys.fluid);
    direct = history_translate == "direct";
    if by_foot
        [foot, grad_foot] = step_foot(t - dt, dt, amplitude, options.rotation, options.center, ...
            direction, options.foot_rule, options.foot_velocity);
        untranslated = state.solid;
        % the transported momentum flux, when the consistent accelerations read it
        transported = cell(1, 3 + (transported_flux && startsWith(acceleration_rule, "consistent")));
        [transported{:}] = transport_solid_state(sp, state.solid, sm_old, sm, k, ...
            keys_old.solid, keys.solid, keys_old.solid_faces, keys.solid_faces, foot, grad_foot, ...
            quadrature_degree = options.quadrature_degree, foot_integration = options.foot_integration);
        [state.solid, switched_s, moved] = transported{1:3};
        q_star = [];
        if ~isempty(options.quadrature_degree), q_star = options.quadrature_degree + 2; end
        if by_ustar
            U = foot_projection(sm_old, sm, cat(1, u_stars{:}), k + 1, k + 1, foot, ...
                quadrature_degree = q_star, foot_integration = options.foot_integration);
            u_stars = mat2cell(U, 2*nchoosek(k + 3, 2)*ones(1, numel(u_stars))).';
        end
    else
        shift = (d(t) - d(t - dt))*direction;
        untranslated = state.solid;
        state = translate_state(sp, sm_old, state, shift, k);
        [state.solid, switched_s, moved, remap_info] = remap_solid_state(sp, state.solid, sm_old, sm, k, ...
            keys_old.solid, keys.solid, keys_old.solid_faces, keys.solid_faces, ...
            "untranslated", untranslated, "shift", shift, carry{:}, ...
            "slid_subdivisions", options.slid_subdivisions, "slid_correction", options.slid_correction, ...
            "slid_displacement", options.slid_displacement);
        if direct
            state.solid = direct_translate_state(sp, sm_old, sm, k, untranslated, shift, ...
                keys_old.solid, keys.solid, state.solid);
        end
        if by_ustar
            if direct
                U = sp.translated_field_between(sm_old, sm, k + 1, k + 1, cat(1, u_stars{:}), shift);
            else
                U = sp.translated_field(sm_old, k + 1, cat(1, u_stars{:}), shift);
                U = carve_remap(U, k + 1, sm_old, sm, keys_old.solid, keys.solid);
            end
            u_stars = mat2cell(U, 2*nchoosek(k + 3, 2)*ones(1, numel(u_stars))).';
        end
    end
    face_on = accumarray(sm.boundaries(:), 1, [size(sm.faces, 2), 1]) == 1;
    if any(skeleton_history == ["volume", "ustar", "ustar_displacement"])
        % the displacement level of skeleton_history, before the consistent
        % and balanced accelerations, which read it: the trace of u_h or of
        % the translated u*
        if skeleton_history == "volume"
            U_hat = boundary_trace_projection(sm, k, state.solid.displacement, k);
        else
            U_star_t = postprocessed_displacement(sp, sm_old, k, untranslated);
            if by_foot
                U_star_t = foot_projection(sm_old, sm, U_star_t, k + 1, k + 1, foot, ...
                    quadrature_degree = q_star, foot_integration = options.foot_integration);
            elseif direct
                U_star_t = sp.translated_field_between(sm_old, sm, k + 1, k + 1, U_star_t, shift);
            else
                U_star_t = carve_remap(sp.translated_field(sm_old, k + 1, U_star_t, shift), ...
                    k + 1, sm_old, sm, keys_old.solid, keys.solid);
            end
            U_hat = boundary_trace_projection(sm, k + 1, U_star_t, k);
        end
        state.solid.displacement_skeleton(:, face_on) = U_hat(:, face_on);
    end
    consistent = startsWith(acceleration_rule, "consistent");
    balanced = startsWith(acceleration_rule, "balanced");
    if consistent || balanced
        % the acceleration of the element momentum identity of the
        % translated, remapped state, with f(x - shift, t_n-1) (with the
        % feet: f(foot(x), t_n-1))
        sp_previous = sp;
        if by_foot
            sp_previous.Data.source_term = cellfun(@(f) @(x, y) at_feet(f, foot, x, y), ...
                f_previous, UniformOutput = false);
        else
            sp_previous.Data.source_term = cellfun(@(f) @(x, y) f(x - shift(1), y - shift(2)), ...
                f_previous, UniformOutput = false);
        end
        if consistent
            flux_state = state;
            if by_foot && numel(transported) > 3
                % the stress replaced by the transported momentum flux (see
                % consistent_flux in the help)
                flux_state.solid.stress = transported{4};
            end
            a_consistent = consistent_acceleration(sp_previous, sm, k, flux_state);
        else
            % the lift of u with the solid's boundary trace history
            switch options.balanced_trace
                case "history"
                    g_hat = state.solid.displacement_skeleton;
                case "volume"
                    g_hat = boundary_trace_projection(sm, k, state.solid.displacement, k);
                case "ustar"
                    % u* of the old state carried along the characteristic
                    U_star = postprocessed_displacement(sp, sm_old, k, untranslated);
                    if direct
                        U_star = sp.translated_field_between(sm_old, sm, k + 1, k + 1, U_star, shift);
                    else
                        U_star = carve_remap(sp.translated_field(sm_old, k + 1, U_star, shift), ...
                            k + 1, sm_old, sm, keys_old.solid, keys.solid);
                    end
                    g_hat = boundary_trace_projection(sm, k + 1, U_star, k);
            end
            [a_consistent, sigma_b, gamma_b, u_hat_b] = balanced_acceleration(sp_previous, sm, k, ...
                state.solid.displacement, boundary_trace = g_hat);
            if acceleration_rule == "balanced_levels"
                interior = accumarray(sm.boundaries(:), 1, [size(sm.faces, 2), 1]) == 2;
                state.solid.stress = sigma_b;
                state.solid.rotation = gamma_b;
                state.solid.displacement_skeleton(:, interior) = u_hat_b(:, interior);
            end
        end
        if bossak_am ~= 0
            % the Bossak identity (1 - am) a^n + am a^(n-1) = r^n of the
            % rebuilt r^n = M^-1(f - K u), with a^(n-1) the acceleration the
            % previous solve read, carried along the characteristic (in the
            % first step the carried initial acceleration)
            a_hist = state.solid.acceleration;
            if ~isempty(a_read)
                if by_foot
                    a_hist = foot_projection(sm_old, sm, a_read, k, k, foot, ...
                        quadrature_degree = options.quadrature_degree, ...
                        foot_integration = options.foot_integration);
                elseif direct
                    a_hist = sp.translated_field_between(sm_old, sm, k, k, a_read, shift);
                else
                    a_hist = carve_remap(sp.translated_field(sm_old, k, a_read, shift), k, ...
                        sm_old, sm, keys_old.solid, keys.solid);
                end
            end
            a_consistent = (a_consistent - bossak_am*a_hist)/(1 - bossak_am);
        end
        if any(acceleration_rule == ["consistent", "balanced", "balanced_levels"])
            state.solid.acceleration = a_consistent;
        elseif acceleration_rule == "consistent_kept"
            kept = setdiff(1:size(sm.elements, 2), moved);
            state.solid.acceleration(:, kept) = a_consistent(:, kept);
        end
    end
    if rule == "ustar_bdf_all", moved = (1:size(sm.elements, 2)).'; end
    a_newmark = state.solid.acceleration;
    % the elements the velocity rule acts on
    ruled = moved;
    if options.history_velocity_scope ~= "moved"
        present = ismembertol(keys.solid, keys_old.solid, 1e-10, "ByRows", true);
        if options.history_velocity_scope == "reshaped"
            ruled = moved(present(moved));
        else
            ruled = moved(~present(moved));
        end
    end
    if by_ustar && ~isempty(a_previous)
        if by_foot
            a_previous = foot_projection(sm_old, sm, a_previous, k, k, foot, ...
                quadrature_degree = options.quadrature_degree, ...
                foot_integration = options.foot_integration);
        elseif direct
            a_previous = sp.translated_field_between(sm_old, sm, k, k, a_previous, shift);
        else
            a_previous = carve_remap(sp.translated_field(sm_old, k, a_previous, shift), k, ...
                sm_old, sm, keys_old.solid, keys.solid);
        end
        if ~isempty(ruled) && numel(u_stars) == 4 && by_differences
            % the one-sided differences of order 2
            V = (3*u_stars{1} - 4*u_stars{2} + u_stars{3})/(2*dt);
            A = (2*u_stars{1} - 5*u_stars{2} + 4*u_stars{3} - u_stars{4})/dt^2;
            state.solid.velocity(:, ruled) = reexpress(V, k + 1, k, sm, sm, ruled, ruled);
            state.solid.acceleration(:, ruled) = reexpress(A, k + 1, k, sm, sm, ruled, ruled);
        elseif ~isempty(ruled) && two_step
            % the two-step form of the trapezoidal Newmark scheme,
            % u^(n+1) - 2 u^n + u^(n-1) = dt^2 (a^(n+1) + 2 a^n + a^(n-1))/4
            [a_n, a_p] = deal(a_newmark(:, ruled), a_previous(:, ruled));
            if q < k
                a_n = project_down(a_newmark, k, q, sm, ruled);
                a_p = project_down(a_previous, k, q, sm, ruled);
                state.solid.acceleration(:, ruled) = a_n;
            end
            state.solid.velocity(:, ruled) = reexpress((u_stars{1} - u_stars{2})/dt, k + 1, k, ...
                sm, sm, ruled, ruled) + dt/4*(a_n + a_p);
        end
    end
    if any(acceleration_rule == ["consistent_after", "balanced_after"])
        state.solid.acceleration = a_consistent;
        % the acceleration the solve reads: a^(n-1) of the next step's rule
        a_newmark = a_consistent;
    end
    if by_ustar, a_previous = a_newmark; end
    if ~any(skeleton_history == ["library", "ustar_displacement"])
        % the velocity and acceleration levels of skeleton_history: the
        % traces of the volume levels the solve reads, of the owning element
        V_hat = boundary_trace_projection(sm, k, state.solid.velocity, k);
        A_hat = boundary_trace_projection(sm, k, state.solid.acceleration, k);
        if skeleton_history ~= "rates_acceleration"
            state.solid.velocity_skeleton(:, face_on) = V_hat(:, face_on);
        end
        if skeleton_history ~= "rates_velocity"
            state.solid.acceleration_skeleton(:, face_on) = A_hat(:, face_on);
        end
    end
    solve_state = state;
    history_solid = state.solid;
    if options.fluid_scheme == "bdf2" && n > 1
        % BDF2 (constant step) as the theta step with theta = 2/3 from the
        % history w = (4 v^n - v^{n-1})/3 and a zero acceleration
        v_previous = remap_fluid_state(struct("velocity", v_previous), fm_old, fm, k, ...
            keys_old.fluid, keys.fluid).velocity;
        fp.Data.theta = 2/3;
        solve_state.fluid.velocity = (4*state.fluid.velocity - v_previous)/3;
        solve_state.fluid.acceleration = zeros(size(state.fluid.velocity));
    end
    sp_solve = sp;
    if bossak_am ~= 0
        % the Bossak step as a Newmark solve: density (1 - am) rho and the
        % volume history acceleration a0 - am a0/((1 - am) c); the source
        % keeps the true density
        c = (1/2 - sp.Data.newmark_beta)/sp.Data.newmark_beta;
        sp_solve.Data.density = (1 - bossak_am)*sp.Data.density;
        a0 = solve_state.solid.acceleration;
        a_read = a0;
        solve_state.solid.acceleration = a0 - bossak_am*a0/((1 - bossak_am)*c);
    end
    [fs, ss, new_state, gs] = FluidStructureCoupling.solve(fp, sp_solve, fm, sm, curve, src{:}, ...
        "newmark", "quiet", true, "state", solve_state, "patch_weight", "jacobian", quad{:}, ...
        "swap_transmission", options.swap_transmission);
    if bossak_am ~= 0
        % every solid level advanced from the true history
        advanced = get_kinematic_state(sp, ss, state);
        new_state.solid = advanced.solid;
    end
    v_previous = state.fluid.velocity;
    state = new_state;
    f_previous = sp.Data.source_term;
    if by_ustar
        u_stars = [{postprocessed_displacement(sp, sm, k, state.solid)}, u_stars(1:min(3, end))];
    end
    if options.measure_interface
        interface(n + 1, :) = measure(t, sp, sm, state.solid, curve, pd);
    end
    if options.interface_feedback
        [~, feedback(n + 1, :)] = solved_curve(t, predicted, sp, sm, k, state.solid, reference, ...
            curve, d(t)*direction, exact_set(t), pd.u_S, options.center + d(t)*direction);
        feedback(n + 1, 12:14) = [remap_info.n_carried, remap_info.mismatch, toc(step_timer)];
        curve_carved = curve;
    end
    for part = {gs.partition_1, gs.partition_2}
        l = vecnorm(part{1}.quadrature.interface - part{1}.quadrature.physical, 2, 1);
        paths = [min(paths(1), min(l(:))), max(paths(2), max(l(:)))];
    end
    path_ratio(n, :) = [max_path_ratio(gs.partition_1, fm), max_path_ratio(gs.partition_2, sm)];
    history(n, :) = [t, max(abs(fs.Unknowns.field(2).dof), [], "all"), ...
        max(abs(ss.Unknowns.field(2).dof), [], "all"), switched_f, switched_s];
    if options.snapshot_file ~= "" && (mod(n, options.snapshot_every) == 0 || n == n_steps)
        save_march_snapshots(options.snapshot_file, k, t, fs, ss, state, ...
            curve.get_2d_points(linspace(curve.t_range(1), curve.t_range(2), 401)));
    end
    if options.record_solid
        solid_record(n) = struct("displacement", state.solid.displacement, "velocity", ...
            state.solid.velocity, "acceleration", state.solid.acceleration, "moved", moved, ...
            "history", struct("displacement", history_solid.displacement, "velocity", ...
            history_solid.velocity, "acceleration", history_solid.acceleration, ...
            "displacement_skeleton", history_solid.displacement_skeleton), ...
            "mesh", struct("coordinates", sm.coordinates, "elements", sm.elements, ...
            "faces", sm.faces, "boundaries", sm.boundaries), "fluid_velocity", state.fluid.velocity, ...
            "velocity_skeleton", state.solid.velocity_skeleton, ...
            "history_velocity_skeleton", history_solid.velocity_skeleton);
    end
    if options.verbose
        eF = fs.compute_l2_error({pd.sigma_F, pd.u_F, pd.gamma_F});
        eS = ss.compute_l2_error({pd.sigma_S, pd.u_S, pd.gamma_S});
        errors(n) = max([eF(:); eS(:)]);
        step_errors(n, :) = [eF(1), eF(2), eS(1), eS(2)];
        fprintf("  step %d t %.4f  uF %.3e  uS %.3e  sF %.3e  sS %.3e  on/off F %d/%d S %d/%d\n", ...
            n, t, eF(2), eS(2), eF(1), eS(1), switched_f, switched_s);
    end
end
eF = fs.compute_l2_error({pd.sigma_F, pd.u_F, pd.gamma_F});
eS = ss.compute_l2_error({pd.sigma_S, pd.u_S, pd.gamma_S});
out = struct("eF", eF(:).', "eS", eS(:).', "paths", paths, "path_ratio", path_ratio, ...
    "history", history, "errors", errors, "step_errors", step_errors, "fs", fs, "ss", ss, ...
    "fluid_mesh", fm, "solid_mesh", sm, "h", max([fm.area(:); sm.area(:)]), ...
    "state", state, "keys", keys, "solid_record", solid_record);
if options.measure_interface, out.interface = interface; end
if options.interface_feedback, out.feedback = feedback; end
end

function [foot, grad_foot] = step_foot(t_old, dt, amplitude, rotation, center, direction, rule, velocity)
% the foot map of the step t_old -> t_old + dt and its (constant) Jacobian
% (see the help): the rigid field beta(x, t) = a(t) + theta'(t) J x,
% a(t) = d'(t) e - theta'(t) J c(t), is affine in x, and so is every foot
J = [0, -1; 1, 0];
d = @(t) amplitude*sin(2*pi*t);
c = @(t) center + d(t)*direction;
theta = @(t) rotation*sin(2*pi*t);
Q = @(a) [cos(a), -sin(a); sin(a), cos(a)];
t_new = t_old + dt;
if rule == "exact"
    R = Q(theta(t_old) - theta(t_new));
    foot = @(x) c(t_old) + R*(x - c(t_new));
    grad_foot = @(x) repmat(R, 1, 1, size(x, 2));
    return
end
omega = @(t) 2*pi*rotation*cos(2*pi*t);                       % theta'
a = @(t) 2*pi*amplitude*cos(2*pi*t)*direction - omega(t)*J*c(t);
switch velocity
    case "extrapolated"
        a_w = (3*a(t_old) - a(t_old - dt))/2;  w_B = (3*omega(t_old) - omega(t_old - dt))/2*J;
    case "extrapolated_end"
        a_w = 2*a(t_old) - a(t_old - dt);  w_B = (2*omega(t_old) - omega(t_old - dt))*J;
    case "prescribed"
        a_w = a(t_old + dt/2);  w_B = omega(t_old + dt/2)*J;
end
w = @(x) a_w + w_B*x;
if rule == "midpoint"
    foot = @(x) x - dt*w(x - dt/2*w(x));
    R = eye(2) - dt*w_B*(eye(2) - dt/2*w_B);
else
    foot = @(x) x - dt*w(x);
    R = eye(2) - dt*w_B;
end
grad_foot = @(x) repmat(R, 1, 1, size(x, 2));
end

function v = at_feet(f, foot, x, y)
% f(foot(x, y)) for arrays x, y of one size
z = foot([reshape(x, 1, []); reshape(y, 1, [])]);
v = f(z(1, :), z(2, :));
if ~isscalar(v), v = reshape(v, size(x)); end
end

function curve = rotated_circle(center, radius, angle)
% the circle of the given centre and radius, its parameter t at the angle
% 2 pi t + angle: CircularCurve's points rotated by angle about the centre
x_t = sprintf("%.17g + %.17g*cos(2*pi*t + %.17g)", center(1), radius, angle);
y_t = sprintf("%.17g + %.17g*sin(2*pi*t + %.17g)", center(2), radius, angle);
curve = SmoothCurveClosed(x_t, y_t, [0 1]);
end

function predicted = predicted_curve(sp, sm, k, solid, dt, reference, curve_old, rigid, shift, predictor)
% the predicted interface at t_n from the solid state of t_n-1 (see the
% help): x*(s) = X(s) + D(s) + shift
p = solid;
if predictor == "velocity", p = predicted_solid(solid, dt, 0); end
p_star = sp.postprocessed_displacement(sm, k, p);
[D, D_prime, info] = eulerian_interface_displacement(reference, curve_old, sm, p_star, k + 1, ...
    rigid, continuous = true);
if predictor == "pointwise"
    [D, D_prime] = advance_by_velocity(reference, sm, k, solid.velocity, dt, D, D_prime, info);
end
if ~is_periodic(reference)
    % the ends on the walls: D_1 minus the linear interpolant of its end
    % values (as lagrangian_carved_march's displaced_curve)
    t0 = reference.t_range(1);  L = diff(reference.t_range);
    ends = D([t0, t0 + L]);
    D_wall = @(t) [1; 0] .* (ends(1, 1) + (t - t0)/L*(ends(1, 2) - ends(1, 1)));
    D_shifted = @(t) D(t) - D_wall(t) + shift;
    D_prime_wall = @(t) D_prime(t) - [1; 0] .* ((ends(1, 2) - ends(1, 1))/L + 0*t);
    curve = WallEndedCurve(reference, D_shifted, [0 1], "parametric", true, "derivative", D_prime_wall);
else
    D_shifted = @(t) D(t) + shift;
    curve = DeformedCurve(reference, D_shifted, parametric = true, derivative = D_prime);
end
predicted = struct("curve", curve, "position", @(s) reference.get_2d_points(s) + D_shifted(s), ...
    "breaks", info.breaks);
end

function [position, row] = solved_curve(t, predicted, sp, sm, k, solid, reference, curve, rigid, ...
    geometry, u_S, center)
% Sigma_h of the solved state and one row of out.feedback (see the help;
% columns 12-14 are filled by the caller)
u_star = sp.postprocessed_displacement(sm, k, solid);
[D, ~, info] = eulerian_interface_displacement(reference, curve, sm, u_star, k + 1, rigid, ...
    continuous = true);
position = @(s) reference.get_2d_points(s) + D(s);
exact = @(s) exact_material_points(reference.get_2d_points(s), rigid, u_S);
period = {"period", diff(reference.t_range)};
closed = is_periodic(reference);
m = interface_curve_metrics(position, info.breaks, exact, geometry, period{:});
row = nan(1, 14);
row([1, 6:9]) = [t, m.distance, m.point_error, m.measure_ratio, m.monotonicity];
if ~closed
    ends = position(reference.t_range);
    row(10:11) = ends(1, :) - [0, 1];
end
if isempty(predicted), return, end
m_p = interface_curve_metrics(predicted.position, predicted.breaks, exact, geometry, period{:});
s = piece_samples(info.breaks, closed, diff(reference.t_range));
x_h = position(s);
x_p = predicted.position(s);
if closed
    normal = abs(vecnorm(x_p - center, 2, 1) - vecnorm(x_h - center, 2, 1));
else
    normal = abs(x_p(2, :) - x_h(2, :));
end
row(2:5) = [max(vecnorm(x_p - x_h, 2, 1)), max(normal), m_p.distance, m_p.point_error];
end

function [D, D_prime] = advance_by_velocity(reference, sm, k, velocity, dt, D, D_prime, info)
% D + dt V and its derivative, V(s) = v_h(X(s) + D(s)) with v_h the volume
% velocity (Dubiner, degree k) of the element of the piece of s, made
% continuous as eulerian_interface_displacement's continuous option does:
% each piece corrected by the linear function taking its end values to
% the means of the one-sided values (the ends of an open curve left)
periodic = is_periodic(reference);
period = diff(reference.t_range);
breaks = info.breaks;
K = info.elements;
n = numel(K);
if periodic
    s_start = breaks;  s_end = [breaks(2:end), breaks(1) + period];
else
    s_start = breaks(1:n);  s_end = breaks(2:n + 1);
end
d2 = nchoosek(k + 2, 2);
G = reshape(get_affine_map(sm, K, inverse = true), 4, n);
[~, x3] = get_affine_map(sm, K);
C1 = velocity(1:d2, K).';
C2 = velocity(d2 + (1:d2), K).';
D0 = D;  D0_prime = D_prime;
[c_lo, c_hi] = deal(zeros(2, n));
if n > 1
    left = on_piece(s_end, 1:n, D0(s_end));
    right = on_piece(s_start, 1:n, D0(s_start));
    following = [2:n, 1];
    mean_value = (left + right(:, following))/2;
    c_hi = mean_value - left;
    c_lo(:, following) = mean_value - right(:, following);
    if ~periodic
        c_hi(:, end) = 0;
        c_lo(:, 1) = 0;
    end
end
slope = (c_hi - c_lo)./(s_end - s_start);
D = @(s) advanced(s);
D_prime = @(s) advanced_derivative(s);

    function [s, e] = locate(s)
        s = reshape(s, 1, []);
        if periodic, s = breaks(1) + mod(s - breaks(1), period); end
        e = discretize(s, [-Inf, breaks(2:n), Inf]);
    end

    function v = advanced(s)
        [s, e] = locate(s);
        d0 = D0(s);
        v = d0 + dt*(on_piece(s, e, d0) + c_lo(:, e) + slope(:, e).*(s - s_start(e)));
    end

    function v = advanced_derivative(s)
        [s, e] = locate(s);
        d0_prime = D0_prime(s);
        x_prime = reference.tangent(s) + 0*s + d0_prime;
        [~, J] = on_piece(s, e, D0(s));
        v = d0_prime + dt*([J(1, :).*x_prime(1, :) + J(3, :).*x_prime(2, :); ...
             J(2, :).*x_prime(1, :) + J(4, :).*x_prime(2, :)] + slope(:, e));
    end

    function [v, J] = on_piece(s, e, d0)
        % v_h of the element of piece e at x(s) = X(s) + d0, and grad v_h
        % there (rows d1v1 d1v2 d2v1 d2v2)
        e = e + 0*s;
        y = reference.get_2d_points(s) + d0 - x3(:, e);
        r = (G(1, e).*y(1, :) + G(3, e).*y(2, :)).';
        q = (G(2, e).*y(1, :) + G(4, e).*y(2, :)).';
        [B, Br, Bq] = dubiner_basis(k, r, q);
        a = C1(e, :);  b = C2(e, :);
        v = [sum(B.*a, 2), sum(B.*b, 2)].';
        if nargout < 2, return, end
        ar = sum(Br.*a, 2).';  aq = sum(Bq.*a, 2).';
        br = sum(Br.*b, 2).';  bq = sum(Bq.*b, 2).';
        J = [ar.*G(1, e) + aq.*G(2, e); br.*G(1, e) + bq.*G(2, e); ...
             ar.*G(3, e) + aq.*G(4, e); br.*G(3, e) + bq.*G(4, e)];
    end
end

function s = piece_samples(breaks, closed, period)
% the samples of interface_curve_metrics: 9 Chebyshev-Lobatto points of
% each piece [a + delta, b - delta], delta = 1e-12 (b - a), as one row
if closed
    s_lo = breaks;  s_hi = [breaks(2:end), breaks(1) + period];
else
    s_lo = breaks(1:end - 1);  s_hi = breaks(2:end);
end
nodes = -cos(pi*(0:8).'/8);
delta = 1e-12*(s_hi - s_lo);
s = reshape((s_lo + delta) + (nodes + 1)/2 .* ((s_hi - delta) - (s_lo + delta)), 1, []);
end

function p = predicted_solid(s, dt, c)
% the Newmark predictor u + dt v + c a of the solid displacement, and of
% its skeleton trace, stress and rotation where the history holds them
% (absent rates are zero); lagrangian_carved_march's
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

function A = as_matrix(A)
if iscell(A), A = cat(1, A{:}); end
end

function row = interface_measures(t, sp, sm, k, solid, reference, curve, rigid, geometry, u_S, preimage)
% one row of out.interface (see the help)
u_star = sp.postprocessed_displacement(sm, k, solid);
if preimage
    % the material pre-image of Sigma(t), Y(s) = x(s) - rigid - u_S(x(s)),
    % x(s) = X(s) + rigid on Sigma(t): its exact material points lie on
    % Sigma(t), so x_h is solved where u* is defined, not extrapolated by
    % the normal part of u_S
    u_at = @(x) [u_S{1}(x(1, :), x(2, :)) + 0*x(1, :); u_S{2}(x(1, :), x(2, :)) + 0*x(1, :)];
    reference = DeformedCurve(reference, @(s) -u_at(reference.get_2d_points(s) + rigid), parametric = true);
end
[D_c, ~, info] = eulerian_interface_displacement(reference, curve, sm, u_star, k + 1, rigid, ...
    continuous = true);
D_raw = eulerian_interface_displacement(reference, curve, sm, u_star, k + 1, rigid);
[D_h, ~, info_h] = eulerian_interface_displacement(reference, curve, sm, solid.displacement, k, ...
    rigid, continuous = true);
X = @(s) reference.get_2d_points(s);
exact = @(s) exact_material_points(X(s), rigid, u_S);
period = {"period", diff(reference.t_range)};
m_c = interface_curve_metrics(@(s) X(s) + D_c(s), info.breaks, exact, geometry, period{:});
m_raw = interface_curve_metrics(@(s) X(s) + D_raw(s), info.breaks, [], geometry, period{:});
m_h = interface_curve_metrics(@(s) X(s) + D_h(s), info_h.breaks, exact, geometry, period{:});
row = [t, m_c.distance, m_c.point_error, m_c.measure_ratio, m_c.monotonicity, ...
    m_raw.monotonicity, m_raw.max_jump, m_h.distance, m_h.point_error];
end

function x = exact_material_points(X, rigid, u)
% the solution of x - rigid - u(x) = X by the fixed-point iteration (u a
% contraction on the data used here)
x = X + rigid;
for it = 1:500
    x_new = X + rigid + [u{1}(x(1, :), x(2, :)); u{2}(x(1, :), x(2, :))];
    change = max(abs(x_new - x), [], "all");
    x = x_new;
    if change <= 1e-15, break; end
end
end

function [fluid, solid, keys] = carve_both(background, solid_background, curve, snap)
% the two sides of the background split by curve and the keys of the
% unsnapped carve: element centroids (fluid, solid) and solid face
% midpoints. With snapping the normal-deviation prune of split is off
% (see lagrangian_carved_march's carve_fluid)
split = @(grid) grid.split(curve);
if snap
    split = @(grid) grid.split(curve, max_normal_deviation = 90);
end
[fluid, solid] = split(background);
if ~isempty(solid_background)
    [~, solid] = split(solid_background);
end
centroid = @(m) squeeze(mean(reshape(m.coordinates(:, m.elements), 2, 3, []), 2)).';
keys = struct("fluid", centroid(fluid), "solid", centroid(solid), "solid_faces", ...
    ((solid.coordinates(:, solid.faces(1, :)) + solid.coordinates(:, solid.faces(2, :)))/2).');
if snap
    fluid = snap_to_curve(fluid, curve);
    solid = snap_to_curve(solid, curve, normal_sign = -1);
end
end

function ratio = max_path_ratio(partition, mesh)
% the largest transfer path length over the height of its element normal
% to the patch's face, 2|K|/|F|
l = vecnorm(partition.quadrature.interface - partition.quadrature.physical, 2, 1);
h_perp = 2*mesh.volume(partition.elements(:)) ./ mesh.area(partition.faces(:));
ratio = max(l ./ reshape(h_perp, 1, 1, []), [], "all");
end

function [s, switched] = remap_fluid_state(s, old_mesh, new_mesh, k, old_keys, new_keys)
% the fluid volume histories (velocity, acceleration) from old_mesh to
% new_mesh, as eulerian_carved_march's remap_fluid_state: an element
% present at both steps keeps its polynomial, re-expressed on its new
% shape if snapping moved it; an element that switches on gets the L2
% projection of the nearest old element's polynomial. switched counts the
% background elements switched on and off.
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
fit = (P.'*(w(:).*P)) \ (P.'.*w(:).');
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

function U = carve_remap(U, degree, old_mesh, new_mesh, old_keys, new_keys)
% Dubiner fields of the given degree from old_mesh to new_mesh: each
% element's source is itself if present in both carves, else the old
% element with the nearest centroid; the source's polynomial re-expressed
[present, from] = ismembertol(new_keys, old_keys, 1e-10, "ByRows", true);
for K = find(~present).'
    [~, from(K)] = min(vecnorm(old_keys - new_keys(K, :), 2, 2));
end
U = reexpress(U, degree, degree, old_mesh, new_mesh, from, 1:numel(from));
end

function V = project_down(D, k, q, mesh, elements)
% the L2 projections onto P_q of the degree-k Dubiner fields D on the given
% elements, written in the degree-k basis
Dq = zeros(size(D, 1)/nchoosek(k + 2, 2)*nchoosek(q + 2, 2), size(D, 2));
Dq(:, elements) = reexpress(D, k, q, mesh, mesh, elements, elements);
V = reexpress(Dq, q, k, mesh, mesh, elements, elements);
end

function V = reexpress(D, deg_from, deg_to, from_mesh, to_mesh, el_from, el_to)
% the P_deg_to L2 projections on the to_mesh elements el_to of the Dubiner
% fields D (degree deg_from, component blocks) of the from_mesh elements
% el_from, evaluated (extrapolated) there
[rs, w] = get_2d_reference_quadrature(max(deg_from, deg_to) + 1);
P = dubiner_basis(deg_to, rs);
fit = (P.'*(w(:).*P)) \ (P.'.*w(:).');
d_from = nchoosek(deg_from + 2, 2);
d_to = nchoosek(deg_to + 2, 2);
n_c = size(D, 1)/d_from;
V = zeros(n_c*d_to, numel(el_to));
if isempty(el_to), return, end
xy = to_mesh.get_2d_physical(repmat(rs, 1, 1, numel(el_to)), el_to);
r_from = from_mesh.get_2d_local(xy, el_from);
for j = 1:numel(el_to)
    Q = dubiner_basis(deg_from, r_from(:, :, j));
    for c = 0:n_c - 1
        V(c*d_to + (1:d_to), j) = fit*(Q*D(c*d_from + (1:d_from), el_from(j)));
    end
end
end

function G = boundary_trace_projection(mesh, degree, U, k)
% the Legendre P_k dofs (component blocks, 2*(k+1) x faces) on the boundary
% faces of mesh of the traces of the Dubiner fields U (degree `degree`,
% two components) of the owning elements; zero on interior faces
n_faces = size(mesh.faces, 2);
d1 = k + 1;
G = zeros(2*d1, n_faces);
owners = accumarray(mesh.boundaries(:), 1, [n_faces, 1]);
boundary = find(owners == 1);
owner = zeros(n_faces, 1);
owner(mesh.boundaries(:)) = repelem(1:size(mesh.elements, 2), 3);
[xyf, tf, wf] = quadpoints_faces_2d(mesh, k + 2, boundary, output = "cell");
lw = wf .* legendre_basis(k, tf);
fit = (lw.'*legendre_basis(k, tf)) \ lw.';
nq = numel(wf);
rs = mesh.get_2d_local(reshape([xyf{1}(:).'; xyf{2}(:).'], 2, nq, []), owner(boundary));
d = nchoosek(degree + 2, 2);
for j = 1:numel(boundary)
    Q = dubiner_basis(degree, rs(:, :, j));
    for c = 0:1
        G(c*d1 + (1:d1), boundary(j)) = fit*(Q*U(c*d + (1:d), owner(boundary(j))));
    end
end
end
