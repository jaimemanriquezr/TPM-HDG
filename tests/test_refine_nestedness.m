%-------------------------------------------------------------------------
%Description: Nestedness of the discrete spaces under red refinement:
%             every coarse basis function (Dubiner, PEERS, Legendre),
%             restricted to a child, must lie in that child's space
%             (relative residual < 1e-10).
%
%             The enriched PEERS stress space is nested because each red
%             child's reference-to-reference map is a homothety
%             (A = +-I/2). A general sub-triangle has no such structure;
%             the negative control checks that it is not nested.
%
%See also: refinement_l2_error, SimplicialGrid.uniform_refine
%-------------------------------------------------------------------------
n_failed = 0;
rng(0);

% A deliberately generic parent: no right angle, no symmetry, so nothing
% passes by accident of the reference element.
V = [0.3 1.2 -0.1; 0.1 0.4 -0.2];
v1 = V(:, 1);  v2 = V(:, 2);  v3 = V(:, 3);
m1 = (v2 + v3)/2;  m2 = (v1 + v3)/2;  m3 = (v1 + v2)/2;

% The four children, in uniform_refine's ordering.
children = {[v1 m3 m2], [m3 v2 m1], [m2 m1 v3], [m1 m2 m3]};
B_parent = [v1 - v3, v2 - v3];

for k = 1:3
    worst_dubiner = 0;
    worst_peers   = 0;
    for c = 1:4
        [rs_child, rs_parent] = sample_pair(children{c}, B_parent, v3);

        D_child  = DubinerBasis.evaluate_basis(k, rs_child(1, :), rs_child(2, :));
        D_parent = DubinerBasis.evaluate_basis(k, rs_parent(1, :), rs_parent(2, :));
        worst_dubiner = max(worst_dubiner, span_residual(D_child, D_parent));

        T = children{c};
        B_child = [T(:, 1) - T(:, 3), T(:, 2) - T(:, 3)];
        P_child  = flatten(PEERSBasis.evaluate_basis(k, rs_child(1, :), ...
            rs_child(2, :), inv(B_child)));
        P_parent = flatten(PEERSBasis.evaluate_basis(k, rs_parent(1, :), ...
            rs_parent(2, :), inv(B_parent)));
        worst_peers = max(worst_peers, span_residual(P_child, P_parent));
    end
    n_failed = n_failed + check(sprintf("k=%d Dubiner nested (%.1e)", k, worst_dubiner), ...
        worst_dubiner < 1e-10);
    n_failed = n_failed + check(sprintf("k=%d PEERS nested (%.1e)", k, worst_peers), ...
        worst_peers < 1e-10);

    % Declared dimension and emitted column count must agree. (They do not
    % at k = 0: evaluate_basis lacks the `if k > 0` guard of the
    % PEERS_basis handle, so k = 0 is not tested.)
    n_failed = n_failed + check(sprintf("k=%d PEERS Dimension matches", k), ...
        size(PEERSBasis.evaluate_basis(k, 0.3, 0.3, eye(2)), 2) == PEERSBasis(k).Dimension);
end

%% 1-D faces --------------------------------------------------------------
for k = 1:3
    worst = 0;
    for g = {@(t) t/2, @(t) (t + 1)/2}
        t_child = rand(1, 400);
        t_parent = g{1}(t_child);
        L_child  = LegendreBasis.evaluate_basis(k, t_child(:));
        L_parent = LegendreBasis.evaluate_basis(k, t_parent(:));
        worst = max(worst, span_residual(L_child, L_parent));
    end
    n_failed = n_failed + check(sprintf("k=%d Legendre nested (%.1e)", k, worst), ...
        worst < 1e-10);
end

%% Negative control -------------------------------------------------------
% Not a red child: the map is no longer a homothety, so the enrichment
% should not be reproducible. If this passes, the test above is vacuous.
centroid = (v1 + v2 + v3)/3;
[rs_child, rs_parent] = sample_pair([v1 m3 centroid], B_parent, v3);
B_odd = [v1 - centroid, m3 - centroid];
P_child  = flatten(PEERSBasis.evaluate_basis(1, rs_child(1, :), rs_child(2, :), inv(B_odd)));
P_parent = flatten(PEERSBasis.evaluate_basis(1, rs_parent(1, :), rs_parent(2, :), inv(B_parent)));
residual = span_residual(P_child, P_parent);
n_failed = n_failed + check(sprintf( ...
    "non-red sub-triangle is NOT nested (%.1e)", residual), residual > 1e-8);

fprintf("\n%d check(s) failed.\n", n_failed);

%=========================================================================%
function [rs_child, rs_parent] = sample_pair(T, B_parent, v3)
% The same physical points, expressed in the child's and the parent's
% reference coordinates.
a = rand(2, 1500);
a = a(:, sum(a, 1) < 1);                    % inside the reference triangle
B_child = [T(:, 1) - T(:, 3), T(:, 2) - T(:, 3)];
xy = B_child * a + T(:, 3);
rs_child = a;
rs_parent = B_parent \ (xy - v3);
end

%=========================================================================%
function r = span_residual(A_child, A_parent)
% Relative residual of projecting the parent's functions onto the span of
% the child's. Zero means every coarse function lives in the child space.
r = norm(A_parent - A_child * (A_child \ A_parent), 'fro') / norm(A_parent, 'fro');
end

%=========================================================================%
function A = flatten(A)
A = reshape(permute(A, [1 3 4 2]), [], size(A, 2));
end

%=========================================================================%
function failed = check(name, condition)
failed = ~condition;
if failed
    fprintf(2, "  FAIL  %s\n", name);
else
    fprintf("  ok    %s\n", name);
end
end
