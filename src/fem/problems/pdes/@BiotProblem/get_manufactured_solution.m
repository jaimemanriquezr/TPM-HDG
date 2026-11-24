function mms = get_manufactured_solution(displacement, pressure, options)
%GET_MANUFACTURED_SOLUTION  Single-domain Biot data from a displacement and a pressure.
%
%   mms = BiotProblem.get_manufactured_solution(displacement, pressure, ...)
%
%   displacement (2x1 strings) and pressure (string) are expressions in
%   x, y, t. With (λ, μ) from (young_modulus, poisson_ratio), κ =
%   diffusion and D_t the time derivative:
%       σ = 2μ ε(η) + λ (∇·η) I − α p I,   γ = (∂_y η_1 − ∂_x η_2)/2,
%       q = −κ ∇p,   f = ρ_P D_tt η − ∇·σ + b₀ η,
%       g = s₀ D_t p + α D_t ∇·η + ∇·q,
%   ρ_P = density. rate = "exact": D_t = ∂_t, D_tt = ∂_tt; rate =
%   "from_rest": D_t w = w/time_step, D_tt w = w/time_step² (one Euler step
%   from a zero state, zero velocity included).
%
%   mms.parameters holds the numbers; data = mms.at(t) returns the fields
%   u_P, u_t_P, v_P (2x1 cells; v_P = ∂_t η for every rate, so mms.at(t0).v_P
%   is the initial velocity of a run started at t0), sigma_P (2x2 cell), gamma_P, q_P (2x1 cell),
%   p_P, f_P (2x1 cell), g_P, each a handle of (x, y) at time t. The
%   symbolic derivation is done once, here.
%
%See also: BiotProblem
arguments
    displacement (2,1) string
    pressure     (1,1) string
    options.young_modulus (1,1) double = 1.0
    options.poisson_ratio (1,1) double = 0.25
    options.biot_alpha    (1,1) double = 1.0
    options.storage       (1,1) double = 0.0
    options.diffusion     (2,2) double = eye(2)
    options.hooke_modulus (1,1) double = 0.0
    options.density       (1,1) double {mustBeNonnegative} = 0.0
    options.time_step     (1,1) double = 1.0
    options.rate (1,1) string {mustBeMember(options.rate, ["exact", "from_rest"])} = "exact"
end
[lambda, mu] = LinearElasticityProblem.get_lame_parameters( ...
    options.young_modulus, options.poisson_ratio);
alpha = options.biot_alpha;
s0 = options.storage;
kappa = options.diffusion;
b0 = options.hooke_modulus;
rho = options.density;

syms x y t real
eta = str2sym(displacement);
p = str2sym(pressure);
grad_eta = jacobian(eta, [x y]);
epsilon = (grad_eta + grad_eta.')/2;
div_eta = trace(grad_eta);
sigma = 2*mu*epsilon + lambda*div_eta*eye(2) - alpha*p*eye(2);
gamma = (diff(eta(1), y) - diff(eta(2), x))/2;
q = -kappa*gradient(p, [x y]);
Div = @(w) [divergence(w(1,:), [x y]); divergence(w(2,:), [x y])];
switch options.rate
    case "exact"
        Dt = @(w) diff(w, t);
        Dtt = @(w) diff(w, t, 2);
    case "from_rest"
        Dt = @(w) w/options.time_step;
        Dtt = @(w) w/options.time_step^2;
end
f = -Div(sigma) + b0*eta;
if rho ~= 0
    f = f + rho*Dtt(eta);
end
eta_t = Dt(eta);
v = diff(eta, t);
g = s0*Dt(p) + alpha*Dt(div_eta) + divergence(q, [x y]);

fields = struct("u_P", eta, "u_t_P", eta_t, "v_P", v, "sigma_P", sigma, "gamma_P", gamma, ...
    "q_P", q, "p_P", p, "f_P", f, "g_P", g);
compiled = struct();
for name = string(fieldnames(fields)).'
    expr = fields.(name);
    handles = cell(size(expr));
    for i = 1:numel(expr)
        handles{i} = matlabFunction(expr(i), Vars = {x, y, t});
    end
    compiled.(name) = handles;
end

mms.parameters = struct("young_modulus", options.young_modulus, ...
    "poisson_ratio", options.poisson_ratio, "lambda", lambda, "mu", mu, ...
    "biot_alpha", alpha, "storage", s0, "diffusion", kappa, ...
    "hooke_modulus", b0, "density", rho, "time_step", options.time_step, "rate", options.rate);
mms.at = @(time) evaluate_at(compiled, time);
end

function data = evaluate_at(compiled, time)
% handles of (x, y) at the given time; scalars are plain handles, vectors
% and tensors cells of handles; constants broadcast through + 0*x
data = struct();
for name = string(fieldnames(compiled)).'
    F = cellfun(@(c) @(x, y) c(x, y, time) + 0*x, compiled.(name), UniformOutput = false);
    if isscalar(F)
        F = F{1};
    end
    data.(name) = F;
end
end
