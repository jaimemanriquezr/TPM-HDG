function mms = get_manufactured_solution(normal_p, options)
%GET_MANUFACTURED_SOLUTION  Manufactured data for a Stokes-Biot test.
%
%   mms = StokesBiotCoupling.get_manufactured_solution(normal_p, ...
%       stream_function =, pressure_fluid =, displacement =, ...
%       pressure_biot =, viscosity =, young_modulus =, poisson_ratio =, ...
%       biot_alpha =, storage =, diffusion =, hooke_modulus =, ...
%       friction_bjs =, time_step =, rate =, density_fluid =, ...
%       density_biot =, convection =)
%
%   derives, from independent exact fields (strings in x, y, t), the fluid
%   and Biot stresses, rotations, Darcy flux, volume sources and the four
%   interface sources of the coupled problem
%
%     fluid:  sigma_F = 2 mu_f eps(u_F) - p_F I,
%             f_F = rho_F D_t u_F + rho_c (u_F . grad) u_F - div sigma_F,
%             u_F = (-d/dy, d/dx) of the stream function (div u_F = 0)
%     Biot:   sigma_P = 2 mu eps(eta) + lambda (div eta) I - alpha p I,
%             q = -kappa grad p,  f_P = rho_P D_tt eta - div sigma_P + b0 eta,
%             g_P = s0 D_t p + alpha D_t div eta + div q
%     Sigma:  T_1 = sigma_F nu_f + sigma_P nu_p            (I2)
%             T_2 = u_F - D_t eta - q                      (I1, normal part)
%             T_3 = p + nu_f . sigma_F nu_f                (I3)
%             T_4 = sigma_F nu_f + gamma_BJS (u_F - D_t eta)   (I4, tangential part)
%
%   with nu_f = -nu_p and (lambda, mu) from
%   LinearElasticityProblem.get_lame_parameters. The rotations are the
%   scalars gamma = (d_y w_1 - d_x w_2)/2 of the skew part of grad w, as
%   in steady_manufactured_data. D_t = d/dt, D_tt = d^2/dt^2 for rate =
%   "exact", and D_t w = w/time_step, D_tt w = w/time_step^2 for rate =
%   "from_rest" (one Euler step from a zero state at rest). rho_F =
%   density_fluid and rho_P = density_biot (0: quasi-static). The
%   convective term is present only with convection = true, with rho_c =
%   rho_F when rho_F > 0 and 1 otherwise (the convention of
%   IncompressibleElasticityMatrices). The fluid pressure is not normalised (pinned closure). The
%   fields are independent, so all four interface sources are non-zero in
%   general; the coupling rows use only nu_p . T_2 and t . T_4.
%
%   normal_p  2 x 1 string: nu_p, the unit normal of Sigma outward from
%             the poroelastic domain, as expressions in x, y valid in a
%             neighbourhood of Sigma (not normalised here).
%
%   mms.parameters  the numbers: viscosity, young_modulus, poisson_ratio,
%             lame_lambda, lame_mu, biot_alpha, coupling_constant
%             (alpha/(2(lambda + mu))), storage, diffusion, hooke_modulus,
%             friction_bjs, time_step, rate, density_fluid, density_biot,
%             convection.
%   mms.at    handle, data = mms.at(t): struct of handles of (x, y) that
%             broadcast (+ 0*x), cells {2,1} for vectors and {2,2} for
%             tensors: u_F, p_F, sigma_F, gamma_F, f_F, u_P, u_t_P,
%             sigma_P, gamma_P, q_P, p_P, f_P, g_P, T_1, T_2, T_3, T_4.
%   mms.symbolic  the same fields as sym arrays in x, y, t (for exact
%             verification; the handles carry constants such as pi^2 in
%             floating point).
%
%   The symbolic derivation runs once per call of this function; mms.at
%   only binds t.
%
%   Example (flat Sigma y = 1/2, Biot below):
%       mms  = StokesBiotCoupling.get_manufactured_solution(["0"; "1"], storage = 1);
%       data = mms.at(0.0);
%
%See also: BiotProblem.get_manufactured_solution,
%          FluidStructureCoupling.get_manufactured_solution
arguments
    normal_p (2,1) string
    options.stream_function (1,1) string = "sin(pi*x)*sin(pi*y)"
    options.pressure_fluid  (1,1) string = "sin(pi*x)*exp(y - 0.5)"
    options.displacement    (2,1) string = ["x*sin(y)"; "y*cos(x)"]
    options.pressure_biot   (1,1) string = "cos(pi*x)*sin(pi*y)"
    options.viscosity       (1,1) double = 1.0
    options.young_modulus   (1,1) double = 1.0
    options.poisson_ratio   (1,1) double = 0.25
    options.biot_alpha      (1,1) double = 1.0
    options.storage         (1,1) double = 0.0
    options.diffusion       (2,2) double = eye(2)
    options.hooke_modulus   (1,1) double = 0.0
    options.friction_bjs    (1,1) double = 0.0
    options.time_step       (1,1) double = 1.0
    options.rate (1,1) string {mustBeMember(options.rate, ["exact", "from_rest"])} = "exact"
    options.density_fluid (1,1) double {mustBeNonnegative} = 0.0
    options.density_biot  (1,1) double {mustBeNonnegative} = 0.0
    options.convection    (1,1) logical = false
end

[lambda, mu] = LinearElasticityProblem.get_lame_parameters( ...
    options.young_modulus, options.poisson_ratio);
mu_f  = options.viscosity;
alpha = options.biot_alpha;
s0    = options.storage;
kappa = options.diffusion;
b0    = options.hooke_modulus;
gamma_bjs = options.friction_bjs;

syms x y t real
strain = @(v) (jacobian(v, [x y]) + jacobian(v, [x y]).')/2;
shear  = @(v) (diff(v(1), y) - diff(v(2), x))/2;
Div    = @(S) [diff(S(1, 1), x) + diff(S(1, 2), y); diff(S(2, 1), x) + diff(S(2, 2), y)];
div    = @(v) diff(v(1), x) + diff(v(2), y);
if options.rate == "exact"
    D_t = @(w) diff(w, t);
    D_tt = @(w) diff(w, t, 2);
else
    D_t = @(w) w/options.time_step;
    D_tt = @(w) w/options.time_step^2;
end
rho_F = options.density_fluid;
rho_P = options.density_biot;

nu_p = str2sym(normal_p);
nu_f = -nu_p;

% Fluid
psi = str2sym(options.stream_function);
u_F = [-diff(psi, y); diff(psi, x)];
p_F = str2sym(options.pressure_fluid);
sigma_F = 2*mu_f*strain(u_F) - p_F*eye(2);
gamma_F = shear(u_F);
f_F = -Div(sigma_F);
if rho_F > 0
    f_F = f_F + rho_F*D_t(u_F);
end
if options.convection
    rho_c = 1;
    if rho_F > 0, rho_c = rho_F; end
    f_F = f_F + rho_c*jacobian(u_F, [x y])*u_F;
end

% Biot
u_P = str2sym(options.displacement);
p_P = str2sym(options.pressure_biot);
sigma_P = 2*mu*strain(u_P) + lambda*div(u_P)*eye(2) - alpha*p_P*eye(2);
gamma_P = shear(u_P);
q_P = -kappa*[diff(p_P, x); diff(p_P, y)];
f_P = -Div(sigma_P) + b0*u_P;
if rho_P > 0
    f_P = f_P + rho_P*D_tt(u_P);
end
u_t_P = D_t(u_P);
g_P = s0*D_t(p_P) + alpha*D_t(div(u_P)) + div(q_P);

% Interface sources
T_1 = sigma_F*nu_f + sigma_P*nu_p;
T_2 = u_F - u_t_P - q_P;
T_3 = p_P + nu_f.'*sigma_F*nu_f;
T_4 = sigma_F*nu_f + gamma_bjs*(u_F - u_t_P);

symbolic = struct("u_F", u_F, "p_F", p_F, "sigma_F", sigma_F, ...
    "gamma_F", gamma_F, "f_F", f_F, "u_P", u_P, "u_t_P", u_t_P, ...
    "sigma_P", sigma_P, "gamma_P", gamma_P, "q_P", q_P, "p_P", p_P, ...
    "f_P", f_P, "g_P", g_P, "T_1", T_1, "T_2", T_2, "T_3", T_3, "T_4", T_4);

% Handles of (x, y, t), built once; mms.at binds t
names = string(fieldnames(symbolic)).';
handles = struct();
for name = names
    S = symbolic.(name);
    H = cell(size(S));
    for i = 1:numel(S)
        H{i} = matlabFunction(S(i), Vars = [x y t]);
    end
    handles.(name) = H;
end

mms.parameters = struct("viscosity", mu_f, ...
    "young_modulus", options.young_modulus, ...
    "poisson_ratio", options.poisson_ratio, ...
    "lame_lambda", lambda, "lame_mu", mu, ...
    "biot_alpha", alpha, "coupling_constant", alpha/(2*(lambda + mu)), ...
    "storage", s0, "diffusion", kappa, "hooke_modulus", b0, ...
    "friction_bjs", gamma_bjs, "time_step", options.time_step, ...
    "rate", options.rate, "density_fluid", rho_F, "density_biot", rho_P, ...
    "convection", options.convection);
mms.at = @(time) bind_time(handles, time);
mms.symbolic = symbolic;
end

function data = bind_time(handles, time)
% Every field as handles of (x, y) at the given time; scalars are plain
% handles, vectors {2,1} and tensors {2,2} cells of handles.
data = struct();
for name = string(fieldnames(handles)).'
    H = handles.(name);
    B = cell(size(H));
    for i = 1:numel(H)
        B{i} = bind_one(H{i}, time);
    end
    if isscalar(B)
        data.(name) = B{1};
    else
        data.(name) = B;
    end
end
end

function f = bind_one(h, time)
f = @(xx, yy) h(xx, yy, time) + 0*xx;
end
