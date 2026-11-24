function data = get_manufactured_solution(u, young_modulus, poisson_ratio, hooke_modulus, density, time, material, solid_viscosity)
%-------------------------------------------------------------------------
%Description: Manufactured data of the hyperelastic problem from an exact
%             displacement u(x, y, t) (strings in x, y, t):
%                 P = P(I + ∇u) + F S_v,   f = ρ ü − ∇·P + b₀ u,
%             S_v = μ_v(ḞᵀF + FᵀḞ) the Kelvin–Voigt stress (solid_viscosity
%             = μ_v, default 0),
%             with the exact local unknowns of the HDG scheme: the first
%             Piola stress σ = P (field 1) and the rotation
%             γ = (∂_y u_x − ∂_x u_y)/2 (field 3). Returns handles at
%             t = time: sigma = P (2x2), u, u_t, u_tt (2x1), gamma, f (2x1)
%             and the linearized stress linear_stress = Cε(u) (2x2).
%-------------------------------------------------------------------------
arguments (Input)
    u             = ["x^2*y"; "x*y^2"]
    young_modulus = 1
    poisson_ratio = 0.3
    hooke_modulus = 0.0
    density       = 0.0
    time          = 0.0
    material      (1,1) string = "saint_venant_kirchhoff"
    solid_viscosity (1,1) double = 0
end
[lambda, mu] = LinearElasticityProblem.get_lame_parameters(young_modulus, poisson_ratio);
syms x y t real
u = str2sym(u(:));
G = jacobian(u, [x y]);
epsilon = (G + G.')/2;
sigma = 2*mu*epsilon + lambda*trace(epsilon)*eye(2);
gamma = (diff(u(1), y) - diff(u(2), x))/2;
P = HyperelasticProblem.first_piola_symbolic(G, lambda, mu, material);
if solid_viscosity ~= 0
    F = eye(2) + G;
    F_dot = diff(G, t);
    P = P + F*(solid_viscosity*(F_dot.'*F + F.'*F_dot));
end
f = density*diff(u, t, 2) - [divergence(P(1, :), [x y]); divergence(P(2, :), [x y])] ...
    + hooke_modulus*u;

at = @(e) subs(e, t, time);
data = struct( ...
    "young_modulus", young_modulus, ...
    "poisson_ratio", poisson_ratio, ...
    "sigma", {handles(at(P))}, ...
    "u",     {handles(at(u))}, ...
    "u_t",   {handles(at(diff(u, t)))}, ...
    "u_tt",  {handles(at(diff(u, t, 2)))}, ...
    "gamma", handles(at(gamma)), ...
    "f",     {handles(at(f))}, ...
    "linear_stress", {handles(at(sigma))});
end

function h = handles(e)
% vectorized @(x, y) handles of a symbolic array (a single handle for a scalar)
syms x y real
h = cell(size(e));
for i = 1:numel(e)
    g = matlabFunction(e(i), Vars = [x, y]);
    h{i} = @(x, y) g(x, y) + 0*x;
end
if isscalar(h), h = h{1}; end
end
