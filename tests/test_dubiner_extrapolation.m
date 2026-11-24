% test_dubiner_extrapolation.m
% The Dubiner functions are polynomials in the natural coordinates (r, s);
% dubiner_eval must return their values and first and second derivatives
% at any (r, s), inside the reference triangle (0, 0), (1, 0), (0, 1) or
% outside it, where points extrapolated from a neighbouring element land.
% For k = 0..4 the six outputs are compared with an independent reference,
% the monomial expansion r^i s^j (i + j <= k) fitted by least squares to
% dubiner_eval on well-interior points and differentiated exactly, at
%   1. interior points (not fitting points);
%   2. the collapsed vertex (0, 1);
%   3. the line s = 1 at r ~= 0, on both sides of the vertex;
%   4. points within a few eps of s = 1 (both sides);
%   5. points at s = 1 - 1e-9;
%   6. points far outside the triangle;
% and 7. a 3-D array of points evaluates as the same points in a column.
% Tolerance: |computed - reference| <= tol*(1 + |reference|), tol = 1e-9.
% The fitting matrix has condition number up to 5.5e3 (k = 4), with fit
% residual 3e-13 and coefficients up to 2.5e3; the reference itself then
% differs from the exact polynomials by up to ~1e-10 at the points farthest
% out (|r|, |s| <= 2, second derivatives, k = 4), well below tol, while the
% defects this test guards against are of size 1e-1 to 1e15.
clc; clear;
n_fail = 0;
tol = 1e-9;

% fitting points, well inside the triangle
[R, S] = meshgrid(linspace(0.02, 0.96, 25));
inside = R + S < 0.98;
r_fit = R(inside); s_fit = S(inside);

groups = { ...
    "interior",       [0.1 0.2; 0.33 0.33; 0.7 0.1; 0.05 0.9; 0.45 0.5]; ...
    "vertex (0, 1)",  [0 1]; ...
    "line s = 1",     [0.3 1; 0.5 1; -0.7 1; 1.5 1]; ...
    "s = 1 +- eps",   [0 1-eps; 0.3 1-eps; -0.7 1+2*eps; 0.3 1+2*eps; 0 1-4*eps; 1.5 1+4*eps]; ...
    "s = 1 - 1e-9",   [0 1-1e-9; 0.3 1-1e-9; -0.7 1-1e-9; 1.5 1-1e-9]; ...
    "far outside",    [1.5 1.3; -0.5 2; 2 -1; -1 -1]};
names = ["P", "dPx", "dPy", "dPxx", "dPyy", "dPxy"];

for k = 0:4
    d2 = nchoosek(k + 2, 2);
    M = monomials(k, r_fit, s_fit);
    C = M \ dubiner_eval(1:d2, r_fit, s_fit);
    for g = 1:size(groups, 1)
        pts = groups{g, 2};
        expected = cell(1, 6);
        [expected{:}] = monomials(k, pts(:, 1), pts(:, 2));
        expected = cellfun(@(A) A*C, expected, UniformOutput = false);
        computed = cell(1, 6);
        [computed{:}] = dubiner_eval(1:d2, pts(:, 1), pts(:, 2));
        d = cellfun(@(a, b) max(abs(a - b)./(1 + abs(b)), [], "all"), computed, expected);
        [d_max, worst] = max(d);
        fprintf("k = %d %-14s max rel. difference %.2e (%s)\n", k, groups{g, 1}, d_max, names(worst));
        if ~(d_max <= tol)
            fprintf("k = %d %s FAILED\n", k, groups{g, 1}); n_fail = n_fail + 1;
        end
    end

    %% 7. 3-D point arrays
    r3 = reshape(linspace(-0.5, 1.5, 12), 3, 2, 2);
    s3 = reshape(linspace(1.2, -0.3, 12), 3, 2, 2);
    computed = cell(1, 6);
    [computed{:}] = dubiner_eval(1:d2, r3, s3);
    d = 0;
    for c = 1:2
        expected = cell(1, 6);
        [expected{:}] = dubiner_eval(1:d2, reshape(r3(:, :, c), [], 1), reshape(s3(:, :, c), [], 1));
        d = max(d, max(cellfun(@(a, b) max(abs(a(:, :, c) - b), [], "all"), computed, expected)));
    end
    fprintf("k = %d 3-D points: size %s, max diff = %.2e\n", k, mat2str(size(computed{1})), d);
    if ~isequal(size(computed{1}), [6, d2, 2]) || d > 0
        fprintf("k = %d 3-D points FAILED\n", k); n_fail = n_fail + 1;
    end
end

fprintf("%d check(s) failed.\n", n_fail);
if n_fail > 0
    error("test_dubiner_extrapolation:failed", "%d check(s) failed", n_fail);
end

function [M, Mr, Ms, Mrr, Mss, Mrs] = monomials(k, r, s)
% r^i s^j, i + j <= k, and their derivatives (dubiner_eval's output order)
p = @(x, e) (e >= 0).*x.^max(e, 0);
[M, Mr, Ms, Mrr, Mss, Mrs] = deal(zeros(numel(r), 0));
for d = 0:k
    for j = 0:d
        i = d - j;
        M(:, end + 1)   = p(r, i).*p(s, j);
        Mr(:, end + 1)  = i*p(r, i - 1).*p(s, j);
        Ms(:, end + 1)  = j*p(r, i).*p(s, j - 1);
        Mrr(:, end + 1) = i*(i - 1)*p(r, i - 2).*p(s, j);
        Mss(:, end + 1) = j*(j - 1)*p(r, i).*p(s, j - 2);
        Mrs(:, end + 1) = i*j*p(r, i - 1).*p(s, j - 1);
    end
end
end
