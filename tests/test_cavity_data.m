% test_cavity_data.m
% Data of the elastic disc in a lid-driven cavity (cavity_problems) and its
% lid-graded background grid (graded_cavity_grid); no solve.
%   1. the ramp r(t) is 0 at t = 0, 1 from ramp_time on, continuous, and
%      the rate datum is its derivative (centred differences);
%   2. lid profiles: "regularized" 16 x^2 (1 - x)^2 vanishes at the top
%      corners and is 1 at x = 1/2; "uniform" is 1 on the lid; at the
%      Gauss points of the boundary faces the lid faces carry the lid
%      value and the side-wall faces 0;
%   3. the default options give the same data as a call without options;
%      viscous_model = "spatial" is stored in the solid's Data;
%   4. solid material: E = 2 G (1 + nu), nu = 0.49;
%   5. graded_cavity_grid: rows y_j = 1 - s (r + (1 - r) s), s = 1 - j/N,
%      heights decreasing toward the lid, vertices in [0, 1]^2 reaching
%      the four sides, element areas positive with total 1.
clc; clear;

n_fail = 0;
k = 1;  dt = 0.01;  T_r = 0.1;
lid_value = @(t, varargin) lid_datum(k, t, dt, 0.5, "datum_dirichlet", varargin{:});
lid_rate = @(t, varargin) lid_datum(k, t, dt, 0.5, "datum_dirichlet_rate", varargin{:});

%% 1. ramp
t = [linspace(0, 1.5*T_r, 31), T_r*(1 - 1e-9), T_r];
r = arrayfun(@(t) lid_value(t), t);
fprintf("ramp: r(0) = %g, r(ramp_time) = %g, max jump on the samples %.2e\n", r(1), r(end), ...
    max(abs(diff(r(1:31)))));
n_fail = n_fail + check(r(1) == 0 && all(r(t >= T_r) == 1), "ramp end values");
n_fail = n_fail + check(abs(r(end - 1) - 1) < 1e-12 && max(abs(diff(r(1:31)))) < 0.1, "ramp continuity");
n_fail = n_fail + check(all(diff(r(1:31)) >= 0), "ramp monotone");
e = 1e-6;
t = [0.013, 0.05, 0.087, 0.15, 0.5];
rate = arrayfun(@(t) lid_rate(t), t);
fd = arrayfun(@(t) (lid_value(t + e) - lid_value(t - e))/(2*e), t);
fprintf("ramp rate: max |rate - centred difference| = %.2e\n", max(abs(rate - fd)));
n_fail = n_fail + check(max(abs(rate - fd)) < 1e-6, "ramp rate");
n_fail = n_fail + check(abs(lid_value(0.05, "ramp_time", 0.2) - (1 - cos(pi/4))/2) < 1e-14 ...
    && lid_value(0, "ramp_time", 0) == 1 && lid_rate(0, "ramp_time", 0) == 0, "ramp_time option");

%% 2. lid profiles
x = linspace(0, 1, 11);
fp = cavity_problems(k, 1, dt, lid_profile = "regularized");
g = fp.Data.datum_dirichlet{1}(x, ones(size(x)));
n_fail = n_fail + check(g(1) == 0 && g(end) == 0 && abs(g(6) - 1) < 1e-14, "regularized lid");
fpu = cavity_problems(k, 1, dt, lid_profile = "uniform");
gu = fpu.Data.datum_dirichlet{1}(x, ones(size(x)));
n_fail = n_fail + check(isequal(size(gu), size(x)) && all(gu == 1), "uniform lid");
n_fail = n_fail + check(all(fpu.Data.datum_dirichlet{1}([0 1 0.5], [0.99 0.5 0]) == 0) ...
    && all(fpu.Data.datum_dirichlet{2}(x, ones(size(x))) == 0), "uniform lid: zero off the lid");
mesh = graded_cavity_grid(4, 0.25);
V = mesh.coordinates;
F = mesh.faces;
lid_faces = find(V(2, F(1, :)) == 1 & V(2, F(2, :)) == 1);
side_faces = find((V(1, F(1, :)) == 0 & V(1, F(2, :)) == 0) | (V(1, F(1, :)) == 1 & V(1, F(2, :)) == 1));
for order = [1, 2]
    xy = quadpoints_faces_2d(mesh, order, lid_faces(:));
    on_lid = fpu.Data.datum_dirichlet{1}(xy{:});
    xy = quadpoints_faces_2d(mesh, order, side_faces(:));
    on_side = fpu.Data.datum_dirichlet{1}(xy{:});
    n_fail = n_fail + check(numel(lid_faces) == 4 && numel(side_faces) == 8 ...
        && all(on_lid == 1, "all") && all(on_side == 0, "all"), ...
        sprintf("uniform lid at the face Gauss points, order %d", order));
end

%% 3. defaults
s = linspace(0, 1, 9);
for t = [0, 0.05, 1]
    a = cell(1, 5);  b = cell(1, 5);
    [a{:}] = cavity_problems(k, t, dt);
    [b{:}] = cavity_problems(k, t, dt, lid_profile = "uniform", ramp_time = 0.1, ...
        viscous_model = "material");
    same = same_value(a{1}.Data, b{1}.Data) && same_value(a{2}.Data, b{2}.Data) ...
        && same_value(a{3}, b{3}) && isequal(a{4}.get_2d_points(s), b{4}.get_2d_points(s)) ...
        && same_value(a{5}(a{4}), b{5}(b{4}));
    n_fail = n_fail + check(same, sprintf("default options change the data at t = %g", t));
end
[~, sp] = cavity_problems(k, 1, dt);
n_fail = n_fail + check(~isfield(sp.Data, "viscous_model") || sp.Data.viscous_model == "material", ...
    "default viscous_model");
[~, sps] = cavity_problems(k, 1, dt, viscous_model = "spatial");
n_fail = n_fail + check(isfield(sps.Data, "viscous_model") && sps.Data.viscous_model == "spatial", ...
    "viscous_model = ""spatial"" not in the solid's Data");

%% 4. material
for G = [0.1, 10]
    [~, sp] = cavity_problems(k, 1, dt, shear_modulus = G);
    n_fail = n_fail + check(sp.Data.poisson_ratio == 0.49 ...
        && abs(sp.Data.young_modulus - 2*G*1.49) < 1e-14*G && sp.Data.material == "neo_hookean", ...
        sprintf("material, G = %g", G));
end

%% 5. graded grid
for N = [8, 16]
    for top_ratio = [0.25, 1]
        mesh = graded_cavity_grid(N, top_ratio);
        V = mesh.coordinates;  T = mesh.elements;
        y = unique(V(2, V(1, :) == 0));
        s = 1 - (0:N)/N;
        rows = diff(y);
        A = ((V(1, T(2, :)) - V(1, T(1, :))).*(V(2, T(3, :)) - V(2, T(1, :))) ...
            - (V(1, T(3, :)) - V(1, T(1, :))).*(V(2, T(2, :)) - V(2, T(1, :))))/2;
        fprintf("graded N = %d, r = %g: top row %.4e (r/N = %.4e), bottom row %.4e, area %.15f\n", ...
            N, top_ratio, rows(end), top_ratio/N, rows(1), sum(abs(A)));
        ok = numel(y) == N + 1 && max(abs(y - (1 - s.*(top_ratio + (1 - top_ratio)*s)))) < 1e-14 ...
            && abs(rows(end) - (top_ratio + (1 - top_ratio)/N)/N) < 1e-14;
        n_fail = n_fail + check(ok, sprintf("row levels, N = %d, r = %g", N, top_ratio));
        if top_ratio < 1
            ok = all(diff(rows) < 0);
        else
            ok = max(abs(rows - 1/N)) < 1e-14;
        end
        n_fail = n_fail + check(ok, sprintf("row heights, N = %d, r = %g", N, top_ratio));
        ok = min(V, [], "all") == 0 && max(V, [], "all") == 1 ...
            && isequal(min(V, [], 2), [0; 0]) && isequal(max(V, [], 2), [1; 1]);
        n_fail = n_fail + check(ok, sprintf("vertices span [0, 1]^2, N = %d, r = %g", N, top_ratio));
        ok = all(sign(A) == sign(A(1))) && min(abs(A)) > 0 && abs(sum(abs(A)) - 1) < 1e-12;
        n_fail = n_fail + check(ok, sprintf("element areas, N = %d, r = %g", N, top_ratio));
    end
end

fprintf("%d check(s) failed.\n", n_fail);
if n_fail > 0
    error("test_cavity_data:failed", "%d check(s) failed", n_fail);
end

function v = lid_datum(k, t, dt, x, name, varargin)
% the x-component of the fluid datum "name" at (x, 1)
fp = cavity_problems(k, t, dt, varargin{:});
v = fp.Data.(name){1}(x, 1);
end

function same = same_value(a, b)
% equality of data, function handles compared by their values at sample points
x = [0, 0.25, 0.5, 0.8, 1, 0.3, 0.6];
y = [1, 1, 1, 1, 0.5, 0, 0.7];
if isa(a, "function_handle") && isa(b, "function_handle")
    if isempty(a) || isempty(b)
        same = isempty(a) && isempty(b);
    else
        same = isequaln(a(x, y), b(x, y));
    end
elseif isstruct(a) && isstruct(b)
    same = isequal(size(a), size(b)) && isequal(sort(fieldnames(a)), sort(fieldnames(b)));
    names = fieldnames(a);
    for i = 1:numel(names)
        if ~same, break; end
        for j = 1:numel(a)
            same = same && same_value(a(j).(names{i}), b(j).(names{i}));
        end
    end
elseif iscell(a) && iscell(b)
    same = isequal(size(a), size(b));
    for i = 1:numel(a)
        if ~same, break; end
        same = same_value(a{i}, b{i});
    end
else
    same = isequaln(a, b);
end
end
