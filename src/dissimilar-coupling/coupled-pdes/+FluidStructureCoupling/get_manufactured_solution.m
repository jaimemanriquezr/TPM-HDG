function [filename, data] = get_manufactured_solution(uu_fluid, p_fluid, u_structure, interface, viscosity_fluid, young_modulus, poisson_ratio, rho_f, rho_s, beta, time, time_step, options)
%GET_MANUFACTURED_SOLUTION  Symbolic manufactured data for a fluid-structure test.
%
%   [filename, data] = get_manufactured_solution(uu_fluid, p_fluid, ...
%       u_structure, interface, viscosity_fluid, young_modulus, ...
%       poisson_ratio, rho_f, rho_s, beta, time, time_step, ...
%       name = "", swap_domains = false, normalize_pressure = true)
%
%   derives, from the exact solutions on the unit square, the fluid and
%   solid stresses, sources and interface data T_1, T_2 of the coupled
%   problem. The positionals (all optional) are the fluid stream function
%   uu_fluid (u_F = (-d/dy, d/dx) of it) and pressure p_fluid, the 2 x 1
%   solid displacement u_structure (strings in x, y, t), the interface
%   (numeric = height of a FlatCurveHorizontal; FlatCurveHorizontal,
%   SmoothCurveHorizontal, CircularCurve or SquircleCurve), the fluid
%   viscosity, the solid Young modulus and Poisson ratio, the densities
%   rho_f and rho_s, the zeroth-order solid coefficient beta, and the time
%   and time step at which the data are evaluated; time = Inf gives the
%   time-discrete data (rho/time_step terms in place of time derivatives).
%
%   data  struct of symbolic fields (measure_F, viscosity, sigma_F, u_F,
%         gamma_F, f_F, young_modulus, poisson_ratio, sigma_S, u_S,
%         u_t_S, gamma_S, f_S, T_1, T_2).
%   filename  a private scratch .mat file, written under cache_dir() with
%         a unique suffix, holding the same data as function handles. It
%         is not a cache: the data are derived on every call and the
%         caller (load_problems) loads the file and deletes it.
%
%   The options only shape the data and the scratch file name: name
%   replaces the parameter-derived part of the name, swap_domains puts the
%   fluid below or inside and the solid above or outside, and
%   normalize_pressure subtracts the fluid-domain mean of the pressure.
%
%See also: FluidStructureCoupling.load_problems
arguments
    % exact solutions
    uu_fluid (1, 1) = "cos(pi*x)*cos(pi*y)"; %u1 = -duu/dy, u2 = duu/dx
    p_fluid (1, 1) = "cos(pi*x)*cos(pi*y)";

    u_structure (:, 1) = ["cos(pi*x)*cos(pi*y)"; "sin(pi*x)*sin(pi*y)"];

    % flat interface I = {(t, y0)}
    interface = FlatCurveHorizontal(.5);

    viscosity_fluid = 1.;

    young_modulus = 1;
    poisson_ratio = 0.25;

    rho_f = 1.0;
    rho_s = 1.0;

    beta = 0.0;
    time = 0.0;
    time_step = 1e-1;

    options.name string = ""

    % Put the fluid on the lower/inner side and the solid on the
    % upper/outer one. The interface normals do not change: they are
    % geometric, and the solver already flips them when the meshes are
    % exchanged.
    options.swap_domains (1,1) logical = false

    % Subtract the fluid-domain mean from the manufactured pressure;
    % required under pressure_mean = "constrained", unnecessary under
    % "pinned"
    options.normalize_pressure (1,1) logical = true
end

if options.swap_domains && ~(isa(interface, "SmoothCurveHorizontal") || ...
        isa(interface, "FlatCurveHorizontal") || isa(interface, "CircularCurve") || ...
        isnumeric(interface))
    error("get_manufactured_solution:swapNeedsOpenInterface", ...
        "swap_domains is implemented for open graph and circular interfaces; got %s.", ...
        class(interface));
end

% The swapped data differs, so its file name must differ from the
% standard assignment's.
swap_tag = "";
if options.swap_domains
    swap_tag = "_swap";
end
% Likewise the un-normalised data: a different manufactured solution,
% so it gets its own tag.
if ~options.normalize_pressure
    swap_tag = swap_tag + "_rawp";
end

if options.name ~= ""
    filename = cache_dir() + "/fluid_structure_solutions_" + options.name + swap_tag;
else
    sanitize = @(s) regexprep(join(string(s(:)).', "_"), '[^a-zA-Z0-9]', '');
    filename = cache_dir() + "/fluid_structure_solutions_" + ...
        sanitize(uu_fluid)        + "_" + ...
        sanitize(p_fluid)         + "_" + ...
        sanitize(u_structure)     + "_" + ...
        class(interface)          + ...
        "_E"    + young_modulus   + ...
        "_nu"   + poisson_ratio   + ...
        "_mu"   + viscosity_fluid + ...
        "_b"    + beta            + ...
        "_rhof" + rho_f           + ...
        "_rhos" + rho_s           + ...
        "_t"    + time            + ...
        "_dt"   + time_step        + swap_tag;
end
% the descriptive part is capped so that the scratch file name (the token
% and the extension added below) stays within the 255-character limit of
% the file system; the token keeps it unique
prefix = string(cache_dir()) + "/";
if strlength(filename) - strlength(prefix) > 180
    filename = extractBefore(filename, strlength(prefix) + 181);
end
if isnumeric(interface)
    interface = FlatCurveHorizontal(interface);
end
supported = ["FlatCurveHorizontal", "SmoothCurveHorizontal", "CircularCurve", "SquircleCurve"];
if ~any(arrayfun(@(c) isa(interface, c), supported))
    error("get_manufactured_solution:unsupportedInterface", ...
        "No manufactured solution for a %s interface (supported: %s).", ...
        class(interface), strjoin(supported, ", "));
end

mu_F = viscosity_fluid;

E = young_modulus;
nu = poisson_ratio;

mu_S = E / (2 * (1 + nu));
lambda_S = E*nu / ((1 + nu)*(1 - 2*nu));

syms x y t real
strain = @(v) .5 * (jacobian(v, [x y]) + jacobian(v, [x, y]).');
shear = @(v) -.5 * (-diff(v(1), y) + diff(v(2), x));
Div = @(v) [divergence(v(1,:),[x y]); divergence(v(2,:), [x y])];
zeta = @(v) 2*mu_S*v + lambda_S*trace(v)*eye(2);

u_F = [-diff(str2sym(uu_fluid), y); diff(str2sym(uu_fluid), x)];
p = str2sym(p_fluid);
if isa(interface, "FlatCurveHorizontal")
    y0 = interface.y_coord(0);
    if options.swap_domains
        % Fluid below, solid above: only the fluid measure and pressure
        % mean change
        measure_F = y0;
    else
        measure_F = 1 - y0;
    end
    n_F = [0; -1];
    n_S = -n_F;
    p0 = p;
    if options.normalize_pressure
        if options.swap_domains
            p0 = p - 1/y0*int(int(p,x,0,1),y,0,y0);
        else
            p0 = p - 1/(1 - y0)*int(int(p,x,0,1),y,y0,1);
        end
    end
elseif isa(interface, "SmoothCurveHorizontal")
    to_x = @(s) subs(s, symvar(s), x);
    I_sym = to_x(interface.sym_expr);
    measure_F = double(1 - int(I_sym, [0 1]));
    if options.swap_domains
        measure_F = double(int(I_sym, [0 1]));
    end
    % The fluid-domain mean is frozen at t = time; exact, since the mean
    % enters only algebraically
    p0 = p;
    if options.normalize_pressure
        if options.swap_domains
            % Fluid below the graph, solid above.
            mean_p = 1/(measure_F) * int(int(p, y, [0 I_sym]), x, [0 1]);
        else
            mean_p = 1/(measure_F) * int(int(p, y, [I_sym 1]), x, [0 1]);
        end
        p0 = p - double(subs(mean_p, t, time));
    end
    % Not swapped with the domains: the orientation is geometric, and the
    % solver flips the normals when the meshes are exchanged
    n_S = [-diff(I_sym); 1] / sqrt(1 + diff(I_sym)^2);
    n_F = -n_S;
elseif isa(interface, "CircularCurve")
    x0 = interface.center(1);
    y0 = interface.center(2);
    r0 = interface.radius;

    if options.swap_domains
        % Fluid inside the circle; normals not flipped, mirroring the
        % open-graph swap branch
        measure_F = pi*r0^2;
    else
        measure_F = 1 - pi*r0^2;
    end

    p0 = p;
    if options.normalize_pressure
        % The polar integral has no closed form for a general p
        syms r theta
        p_r = subs(p, [x y], [x0 + r*cos(theta) y0 + r*sin(theta)]);
        p_int_square = int(int(p, x, 0, 1), y, 0, 1);
        p_int_circle = int(int(p_r*r, r, 0, r0), theta, 0, 2*pi);
        if options.swap_domains
            p_avg = 1/measure_F * p_int_circle;
        else
            p_avg = 1/measure_F * (p_int_square - p_int_circle);
        end
        p0 = p - double(subs(p_avg, t, time));   % mean frozen at t = time
    end
    n_F = -[(x - x0); (y - y0)] / r0;
    n_S = -n_F;
elseif isa(interface, "SquircleCurve")
    % Rounded square r = R(1 - eps*cos(4*theta)), fluid outside. The
    % normal comes from the implicit form F = rho - R(1 - eps*cos4theta),
    % with cos(4theta) = (X^4 - 6 X^2 Y^2 + Y^4)/rho^4 rational in (X, Y),
    % so it is exact in (x, y).
    x0 = interface.center(1);
    y0 = interface.center(2);
    r0 = interface.radius;
    e0 = interface.eccentricity;

    X = x - x0;  Y = y - y0;
    rho2 = X^2 + Y^2;
    cos4theta = (X^4 - 6*X^2*Y^2 + Y^4) / rho2^2;
    F_level = sqrt(rho2) - r0*(1 - e0*cos4theta);
    grad_F = [diff(F_level, x); diff(F_level, y)];
    n_out = grad_F / sqrt(grad_F(1)^2 + grad_F(2)^2);

    % Fluid is the exterior: its outward normal points into the enclosed
    % region, as in the circle branch
    n_F = -n_out;
    n_S = -n_F;

    % Exact: A = 1/2 int r^2 dtheta = pi R^2 (1 + eps^2/2).
    measure_F = 1 - interface.enclosed_area();

    % Mean pressure over the fluid domain, by quadrature (no elementary
    % form). .*ones: matlabFunction of a constant field returns a scalar
    % regardless of input size, and integral2 requires array outputs.
    p0 = p;
    if options.normalize_pressure
        p_raw = matlabFunction(subs(p, t, time), "Vars", [x y]);
        p_fun = @(xx, yy) p_raw(xx, yy) .* ones(size(xx));
        r_of_theta = @(th) r0*(1 - e0*cos(4*th));
        p_int_square = integral2(p_fun, 0, 1, 0, 1, ...
            AbsTol = 1e-12, RelTol = 1e-12);
        % r runs only out to the boundary, so substitute r = s*r_of_theta(th)
        % with s in [0,1]; the Jacobian r dr dtheta becomes s*r_of_theta^2.
        p_int_squircle = integral2( ...
            @(ss, th) p_fun(x0 + ss.*r_of_theta(th).*cos(th), ...
                            y0 + ss.*r_of_theta(th).*sin(th)) ...
                      .* ss .* r_of_theta(th).^2, ...
            0, 1, 0, 2*pi, AbsTol = 1e-12, RelTol = 1e-12);
        p0 = p - (p_int_square - p_int_squircle) / measure_F;
    end
end

sigma_F = 2*mu_F*strain(u_F) - p0*eye(2);
gamma_F = shear(u_F);


u_S = str2sym(u_structure);

sigma_S = zeta(strain(u_S));
gamma_S = shear(u_S);

if ~isinf(time)
    f_F = rho_f*diff(u_F, t) - Div(sigma_F);
    f_S = rho_s*diff(u_S, t, 2) + beta*u_S - Div(sigma_S);
else
    f_F = rho_f/time_step*u_F - Div(sigma_F);
    f_S = rho_s/time_step^2*u_S + beta*u_S - Div(sigma_S);
end
T_1 = sigma_F * n_F + sigma_S * n_S;

if ~isinf(time) && rho_f*rho_s > 0
    T_2 = u_F - diff(u_S, t);
else
    T_2 = u_F - u_S/time_step;
end

sigma_F = subs(sigma_F, t, time);
gamma_F = subs(gamma_F, t, time);
u_F = subs(u_F, t, time);
f_F = subs(f_F, t, time);

sigma_S = subs(sigma_S, t, time);
gamma_S = subs(gamma_S, t, time);
% u_t_S must be differentiated before t is substituted out of u_S
u_t_S = subs(diff(u_S, t), t, time);
u_S = subs(u_S, t, time);
f_S = subs(f_S, t, time);

T_1 = subs(T_1, t, time);
T_2 = subs(T_2, t, time);


data.measure_F = measure_F;
data.viscosity = mu_F;
data.sigma_F = sigma_F;
data.u_F = u_F;
data.gamma_F = gamma_F;
data.f_F = f_F;

data.young_modulus = E;
data.poisson_ratio = nu;
data.sigma_S = sigma_S;
data.u_S = u_S;
data.u_t_S = u_t_S;
data.gamma_S = gamma_S;
data.f_S = f_S;

data.T_1 = T_1;
data.T_2 = T_2;

% Written to a PRIVATE file (the descriptive name plus a unique token) and
% returned: the data are regenerated on every call, so there is nothing to
% share, and concurrent jobs with the same parameters never touch the same
% file. The caller deletes it after loading (see load_problems).
[~, token] = fileparts(tempname);
scratch = filename + "_" + token;
fID = fopen(scratch + ".txt", 'w');
data_vars = reshape(namedargs2cell(data), 2, []);
for i = 1:size(data_vars, 2)
    print_var(fID, to_var(data_vars{:, i}));
end
fclose(fID);
%===========================================================================================%
txt2mat(scratch + ".txt")
%===========================================================================================%
delete(scratch + ".txt");
filename = scratch + ".mat";
end

function var = to_var(name, v)
if isa(v, "sym")
    value = string(v);
    value_type = "function_handle";
else
    value = v;
    value_type = class(v);
end
var = struct("name", name, ...
    "size", size(v, 1:2), ...
    "value", value, ...
    "type", value_type);
end

function print_var(fID, var)
fprintf(fID, ...
    sprintf("%s %i %i %s\n", var.name, var.size(1), var.size(2), var.type));
switch var.type
    case "double"
        fprintf(fID, ...
            sprintf("%.17g\n", var.value));
    case "function_handle"
        fprintf(fID, ...
            sprintf("%s\n", var.value));
end
end
