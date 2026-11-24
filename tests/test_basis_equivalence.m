%% Basis equivalence: low-level function-handle bases vs DSL AbstractFEBasis
% Every basis the low-level path uses must evaluate identically (same values,
% same column enumeration) to its DSL counterpart, so low-level dof vectors
% can be read as FEField dofs without permutation.
clc; clear;
n_fail = 0;
tol = 0;  % same kernels -> byte-identical is the expectation

% deterministic interior points of the reference triangle
t1 = linspace(0.05, 0.85, 7);
t2 = linspace(0.05, 0.9, 7);
[R, S] = meshgrid(t1, t2);
keep = R + S < 0.95;
r = R(keep).'; s = S(keep).';
pts2 = [r; s];
s1 = linspace(0.03, 0.97, 11);

mesh = StructuredSimplicialGrid(3, 3);
G = get_affine_map(mesh, inverse=true);

for k = 0:4
    %% scalar Dubiner
    P_leg = dubiner_basis(k, r.', s.');
    P_dsl = DubinerBasis(k).evaluate(pts2);
    d = max(abs(P_leg(:) - P_dsl(:)), [], 'all');
    fprintf("k = %d scalar Dubiner:  size %s vs %s, max diff = %.3e\n", ...
        k, mat2str(size(P_leg)), mat2str(size(P_dsl)), d);
    if ~isequal(size(P_leg), size(P_dsl)) || d > tol
        fprintf("k = %d scalar Dubiner FAILED\n", k); n_fail = n_fail + 1;
    end

    %% vector Dubiner
    P_leg = vector_dubiner_basis(k, r.', s.');
    P_dsl = DubinerBasis(k).to_vector(2).evaluate(pts2);
    d = max(abs(P_leg(:) - P_dsl(:)), [], 'all');
    fprintf("k = %d vector Dubiner:  size %s vs %s, max diff = %.3e\n", ...
        k, mat2str(size(P_leg)), mat2str(size(P_dsl)), d);
    if ~isequal(size(squeeze(P_leg)), size(squeeze(P_dsl))) || d > tol
        fprintf("k = %d vector Dubiner FAILED\n", k); n_fail = n_fail + 1;
    end

    %% scalar Legendre (skeleton, 1d)
    P_leg = legendre_basis(k, s1.');
    P_dsl = LegendreBasis(k).evaluate(s1);
    d = max(abs(P_leg(:) - P_dsl(:)), [], 'all');
    fprintf("k = %d scalar Legendre: size %s vs %s, max diff = %.3e\n", ...
        k, mat2str(size(P_leg)), mat2str(size(P_dsl)), d);
    if ~isequal(size(P_leg), size(P_dsl)) || d > tol
        fprintf("k = %d scalar Legendre FAILED\n", k); n_fail = n_fail + 1;
    end

    %% PEERS (mesh-dependent)
    % Enumeration is identical, but the bubble part evaluates through
    % different float paths (DSL bubble_eval vs low-level
    % curl_curlA_bubble_basis), so agreement is machine precision, not
    % byte-identity. A permuted column would show O(1) differences.
    P_leg = PEERS_basis(k, r.', s.', G);
    P_dsl = PEERSBasis(k).evaluate(pts2, mesh);
    d = max(abs(P_leg(:) - P_dsl(:)), [], 'all') / max(abs(P_leg(:)));
    fprintf("k = %d PEERS:           size %s vs %s, max rel diff = %.3e\n", ...
        k, mat2str(size(P_leg)), mat2str(size(P_dsl)), d);
    if ~isequal(size(P_leg, 1:5), size(P_dsl, 1:5)) || d > 1e-12
        fprintf("k = %d PEERS FAILED\n", k); n_fail = n_fail + 1;
    end
end

fprintf("\n%d check(s) failed.\n", n_fail);
