% Stokes-Biot manufactured data (StokesBiotCoupling.get_manufactured_solution).
%
% Three groups of checks, on a flat Sigma (y = 1/2, Biot below, nu_p = e_y)
% and a curved one (circle centred at (1/2, 1/2), Biot inside, nu_p radial):
%
%   (a) symbolic: the expressions mms.symbolic satisfy every PDE,
%       constitutive law and interface-source definition identically
%       (simplify to 0), for several parameter sets including the limits
%       kappa -> 0 and s0 -> 0 and both rates; the handles of mms.at agree
%       with them to 1e-12 (relative) at sample points;
%   (b) central differences of the handles at 20 random points and 2
%       times (tolerance 1e-5): both momentum equations, div u_F = 0,
%       Darcy's law, the mass equation, D_t eta, symmetry of both stresses;
%   (c) the hydrostatic state (u_F = 0, p_F = P, eta = 0, p = P, alpha = 1,
%       s0 = 0) gives T_1 = T_2 = T_3 = 0 and t . T_4 = 0 to 1e-12 (the
%       normal part of T_4 is -P nu_f, which the coupling never reads).
clc; clear;
n_fail = 0;

circle_normal = ["(x - 1/2)/sqrt((x - 1/2)^2 + (y - 1/2)^2)"; ...
                 "(y - 1/2)/sqrt((x - 1/2)^2 + (y - 1/2)^2)"];
geometries = {"flat", ["0"; "1"]; "circle", circle_normal};

moving = {"stream_function", "sin(pi*x)*sin(pi*y)*cos(t)", ...
    "pressure_fluid", "x^2 - y*t + exp(x*y)", ...
    "displacement", ["x*sin(y)*(1 + t^2)"; "y*cos(x)*exp(-t)"], ...
    "pressure_biot", "cos(pi*x)*sin(pi*y)*(1 + t)"};
sets = {
    "defaults",          {}
    "general",           [moving, {"viscosity", 0.5, "young_modulus", 3, ...
                          "poisson_ratio", 0.3, "biot_alpha", 0.8, "storage", 1, ...
                          "diffusion", [2 0.5; 0.5 1], "hooke_modulus", 0.7, "friction_bjs", 1}]
    "kappa 1e-8, s0 0",  [moving, {"diffusion", 1e-8*eye(2), "storage", 0, ...
                          "friction_bjs", 10, "rate", "from_rest", "time_step", 0.1}]
    "kappa 0, s0 0",     [moving, {"diffusion", zeros(2), "storage", 0, ...
                          "poisson_ratio", 0.4999, "friction_bjs", 2}]
    "s0 0, from rest",   {"storage", 0, "biot_alpha", 1, "rate", "from_rest", "time_step", 0.25}
    };

%% (a) symbolic residuals
fprintf("(a) symbolic residuals\n");
for g = 1:size(geometries, 1)
    for s = 1:size(sets, 1)
        mms = StokesBiotCoupling.get_manufactured_solution(geometries{g, 2}, sets{s, 2}{:});
        [n_exact, n_total] = symbolic_residuals(mms, geometries{g, 2});
        mismatch = handle_mismatch(mms);
        ok = n_exact == n_total && mismatch < 1e-12;
        n_fail = n_fail + ~ok;
        fprintf("  %-6s %-18s %3d/%3d residuals simplify to 0, handles vs symbolic %.1e  %s\n", ...
            geometries{g, 1}, sets{s, 1}, n_exact, n_total, mismatch, pass_str(ok));
    end
end

%% (b) central differences
fprintf("(b) central differences, 20 points x 2 times\n");
rng(7);
X = 0.05 + 0.9*rand(1, 20);  Y = 0.05 + 0.9*rand(1, 20);
for g = 1:size(geometries, 1)
    for s = 1:size(sets, 1)
        mms = StokesBiotCoupling.get_manufactured_solution(geometries{g, 2}, sets{s, 2}{:});
        err = fd_residuals(mms, X, Y, [0.3, 0.7]);
        ok = err < 1e-5;
        n_fail = n_fail + ~ok;
        fprintf("  %-6s %-18s max scaled residual %.1e  %s\n", ...
            geometries{g, 1}, sets{s, 1}, err, pass_str(ok));
    end
end

%% (c) hydrostatic state
fprintf("(c) hydrostatic state, P = 3\n");
for g = 1:size(geometries, 1)
    mms = StokesBiotCoupling.get_manufactured_solution(geometries{g, 2}, ...
        stream_function = "0", pressure_fluid = "3", displacement = ["0"; "0"], ...
        pressure_biot = "3", biot_alpha = 1, storage = 0, friction_bjs = 1, ...
        diffusion = [2 0.5; 0.5 1], young_modulus = 2, poisson_ratio = 0.3);
    nu_p = normal_handles(geometries{g, 2});
    err = 0;
    for time = [0.0, 0.5]
        d = mms.at(time);
        n1 = nu_p{1}(X, Y);  n2 = nu_p{2}(X, Y);
        tangential = -n2.*d.T_4{1}(X, Y) + n1.*d.T_4{2}(X, Y);
        err = max([err, abs(d.T_1{1}(X, Y)), abs(d.T_1{2}(X, Y)), ...
            abs(d.T_2{1}(X, Y)), abs(d.T_2{2}(X, Y)), abs(d.T_3(X, Y)), abs(tangential)]);
    end
    ok = err < 1e-12;
    n_fail = n_fail + ~ok;
    fprintf("  %-6s max |T_1|, |T_2|, |T_3|, |t.T_4| = %.1e  %s\n", geometries{g, 1}, err, pass_str(ok));
end

fprintf("%d check(s) failed.\n", n_fail);

%=========================================================================%
function s = pass_str(ok)
if ok; s = "ok"; else; s = "FAIL"; end
end

function H = normal_handles(normal_p)
syms x y real
n = str2sym(normal_p);
H = {matlabFunction(n(1), Vars = [x y]); matlabFunction(n(2), Vars = [x y])};
H = cellfun(@(h) @(xx, yy) h(xx, yy) + 0*xx, H, UniformOutput = false);
end

function [n_exact, n_total] = symbolic_residuals(mms, normal_p)
% Residuals of the defining equations on mms.symbolic, with the material
% constants as exact rationals of the stored doubles.
syms x y t real
P = mms.parameters;
F = mms.symbolic;
grad = @(v) jacobian(v, [x y]);
Div  = @(S) [diff(S(1, 1), x) + diff(S(1, 2), y); diff(S(2, 1), x) + diff(S(2, 2), y)];
div  = @(v) diff(v(1), x) + diff(v(2), y);
J = [0 1; -1 0];            % skew matrix of the scalar rotation
I = eye(2);
if P.rate == "exact"
    D_t = @(w) diff(w, t);
else
    D_t = @(w) w/P.time_step;
end
lambda = sym(P.lame_lambda);  mu = sym(P.lame_mu);
c1 = 1/(2*mu);  c2 = lambda/(4*mu*(lambda + mu));
A = @(S) c1*S - c2*trace(S)*I;
c_b = P.biot_alpha/(2*(lambda + mu));
nu_p = str2sym(normal_p);  nu_f = -nu_p;
dev = @(S) S - trace(S)/2*I;

R = {
    % fluid
    -Div(F.sigma_F) - F.f_F
    div(F.u_F)
    F.sigma_F(1, 2) - F.sigma_F(2, 1)
    grad(F.u_F) - dev(F.sigma_F)/(2*P.viscosity) - F.gamma_F*J
    F.p_F + trace(F.sigma_F)/2
    % Biot
    -Div(F.sigma_P) + P.hooke_modulus*F.u_P - F.f_P
    F.sigma_P(1, 2) - F.sigma_P(2, 1)
    grad(F.u_P) - A(F.sigma_P + P.biot_alpha*F.p_P*I) - F.gamma_P*J
    A(F.sigma_P + P.biot_alpha*F.p_P*I) - (A(F.sigma_P) + c_b*F.p_P*I)
    F.q_P + P.diffusion*[diff(F.p_P, x); diff(F.p_P, y)]
    F.u_t_P - D_t(F.u_P)
    F.g_P - (P.storage*D_t(F.p_P) + P.biot_alpha*D_t(div(F.u_P)) + div(F.q_P))
    % interface sources
    F.T_1 - (F.sigma_F*nu_f + F.sigma_P*nu_p)
    F.T_2 - (F.u_F - D_t(F.u_P) - F.q_P)
    F.T_3 - (F.p_P + nu_f.'*F.sigma_F*nu_f)
    F.T_4 - (F.sigma_F*nu_f + P.friction_bjs*(F.u_F - D_t(F.u_P)))
    };
r = sym([]);
for i = 1:numel(R)
    r = [r; R{i}(:)]; %#ok<AGROW>
end
z = simplify(r);
n_exact = nnz(arrayfun(@(e) isequal(e, sym(0)), z));
n_total = numel(z);
% The stored numbers themselves
[l0, m0] = LinearElasticityProblem.get_lame_parameters(P.young_modulus, P.poisson_ratio);
n_total = n_total + 3;
n_exact = n_exact + (P.lame_lambda == l0) + (P.lame_mu == m0) + ...
    (abs(P.coupling_constant - double(c_b)) <= 1e-15*abs(double(c_b)));
end

function worst = handle_mismatch(mms)
% Relative difference between the handles of mms.at and mms.symbolic.
syms x y t real
pts = [0.21 0.33 0.4; 0.77 0.61 0.9; 0.13 0.88 1.3];
worst = 0;
for j = 1:size(pts, 1)
    d = mms.at(pts(j, 3));
    for name = string(fieldnames(d)).'
        C = d.(name);
        if ~iscell(C); C = {C}; end
        S = mms.symbolic.(name);
        for i = 1:numel(C)
            v = C{i}(pts(j, 1), pts(j, 2));
            e = double(subs(S(i), [x y t], pts(j, :)));
            worst = max(worst, abs(v - e)/max(1, abs(e)));
        end
    end
end
end

function err = fd_residuals(mms, X, Y, times)
% Max residual, scaled by max(1, |terms|), of the PDEs by central differences.
P = mms.parameters;
h = 1e-4;  ht = 1e-4;
dx = @(f) (f(X + h, Y) - f(X - h, Y))/(2*h);
dy = @(f) (f(X, Y + h) - f(X, Y - h))/(2*h);
err = 0;
scaled = @(r, s) max(abs(r)./max(1, abs(s)));
for time = times
    d = mms.at(time);
    % fluid momentum and incompressibility
    for a = 1:2
        lhs = -dx(d.sigma_F{a, 1}) - dy(d.sigma_F{a, 2});
        err = max(err, scaled(lhs - d.f_F{a}(X, Y), d.f_F{a}(X, Y)));
    end
    err = max(err, scaled(dx(d.u_F{1}) + dy(d.u_F{2}), dx(d.u_F{1})));
    % Biot momentum
    for a = 1:2
        lhs = -dx(d.sigma_P{a, 1}) - dy(d.sigma_P{a, 2}) + P.hooke_modulus*d.u_P{a}(X, Y);
        err = max(err, scaled(lhs - d.f_P{a}(X, Y), d.f_P{a}(X, Y)));
    end
    % stress symmetry
    err = max(err, scaled(d.sigma_F{1, 2}(X, Y) - d.sigma_F{2, 1}(X, Y), d.sigma_F{1, 2}(X, Y)));
    err = max(err, scaled(d.sigma_P{1, 2}(X, Y) - d.sigma_P{2, 1}(X, Y), d.sigma_P{1, 2}(X, Y)));
    % Darcy's law
    gp = [dx(d.p_P); dy(d.p_P)];
    for a = 1:2
        r = d.q_P{a}(X, Y) + P.diffusion(a, :)*gp;
        err = max(err, scaled(r, d.q_P{a}(X, Y)));
    end
    % D_t eta and D_t p
    if P.rate == "exact"
        dp = mms.at(time + ht);  dm = mms.at(time - ht);
        Dt = @(name, a) (value(dp.(name), a, X, Y) - value(dm.(name), a, X, Y))/(2*ht);
    else
        Dt = @(name, a) value(d.(name), a, X, Y)/P.time_step;
    end
    for a = 1:2
        err = max(err, scaled(d.u_t_P{a}(X, Y) - Dt("u_P", a), d.u_t_P{a}(X, Y)));
    end
    % mass equation, D_t div eta = div D_t eta
    div_ut = dx(d.u_t_P{1}) + dy(d.u_t_P{2});
    div_q = dx(d.q_P{1}) + dy(d.q_P{2});
    rhs = P.storage*Dt("p_P", 0) + P.biot_alpha*div_ut + div_q;
    err = max(err, scaled(d.g_P(X, Y) - rhs, d.g_P(X, Y)));
end
end

function v = value(C, a, X, Y)
% Component a of a vector field (a = 0: a scalar field) at the points.
if a == 0; v = C(X, Y); else; v = C{a}(X, Y); end
end
