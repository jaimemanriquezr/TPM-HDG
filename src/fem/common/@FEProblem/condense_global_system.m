function S = condense_global_system(pde, mesh)
% Assemble this problem's equations and statically condense them over its
% heterogeneous global-support layout, blockwise per global-space pair.
% Shared by the declarative solve and InterfaceProblem, which replaces the
% interface rows of S.H with transmission rows.
%
% Returns a struct S with:
%   L        global layout (global_dof_layout)
%   H, RHS   condensed sparse system over L.n_dof
%   A11d     local matrix pages          [n_loc, n_loc, 1, 1, n_el]
%   S1f      A11 \ A1f                   [n_loc, 1, 1, 1, n_el]
%   S12      per global space, A11 \ A12_t (cell)
%   d_g, f_g per-space dof dims and face-slot counts
% S12 / d_g / f_g: one entry per global space, plus one for the uncondensed
% modes when any local space declares them; recovery is
% x = S1f - sum_t S12{t} x_t over ALL entries.
%
% Uncondensed modes (FESpace.UncondensedModes): each index set
% d = [i1, ..., im] into a local space's dofs is the 0/1 vector
% v = e_i1 + ... + e_im, pivot i1. With T = I + sum_{j>=2} e_ij e_i1' (unit
% triangular up to a permutation, exact; T(:, i1) = v) the local unknowns
% are x = T y, y = (p; y_w), p the coefficients along the directions. The
% local system A11 x + A12 u = A1f tested with the rows w of T' gives
%   y_w = Aww \ ((T' A1f)_w - (T' A12)_w u - Awp p),
% hence x = S1f - sum_t S12_t u_t - S12p p with S1f = loc_solve(A1f),
% S12_t = loc_solve(A12_t) and S12p = T [-I; Aww \ Awp] (scattered to the
% pivot / w rows). The pivot rows of T' are the local equations tested
% with v, v' A11 x + v' A12 u = v' A1f; substituting x gives the p-rows
%   (v' A12_t - rowp S12_t) u_t + (-rowp S12p) p = v' A1f - rowp S1f,
% rowp = v' A11. p is assembled as one extra element-supported global
% block after the global spaces (d = n_p, f = 1); eliminating it again
% recovers the fully condensed system. Several modes (disjoint index
% sets, in one or several local spaces) are handled by the same formulas.
arguments
    pde  (1,1) FEProblem
    mesh AbstractGrid
end
L = global_dof_layout(pde.GlobalSpaces, mesh);
assert(issorted(L.support ~= "face"), ...
    "FEProblem:globalOrder", ...
    "skeleton-supported global spaces must precede element/scalar ones.");

blocks = assemble(pde.Operator, pde);
[A11, A12, A1f] = deal(blocks.A11, blocks.A12, blocks.A1f);
[A21, A22, A2f] = deal(blocks.A21, blocks.A22, blocks.A2f);

n_el = size(mesh.elements, 2);
n_g  = numel(pde.GlobalSpaces);
d_g  = arrayfun(@(s) s.Dimension, pde.GlobalSpaces);
n_sides = size(mesh.boundaries, 1);
f_g  = 1 + (n_sides - 1)*(L.support == "face");    % face slots per space

A11d = double(A11);

trial_of = @(X) cellfun(@(C) C(1), X.id);
test_of  = @(X) cellfun(@(C) C(2), X.id);

% the linear (b) blocks are not zero-padded by assemble_equation:
% place each at its test space's row offset
dims_loc = arrayfun(@(s) s.Dimension, pde.LocalSpaces);
A1fd = zeros(sum(dims_loc), 1, 1, 1, n_el);
% a form set with no Datum terms assembles to a bare 0, not an
% IndexedMatrix -- treat it as an absent RHS
if ~isa(A1f, 'double')
    for j = 1:numel(pde.LocalSpaces)
        bj = A1f{test_of(A1f) == pde.Id.field(j)};
        if ~isempty(bj.id)
            A1fd(sum(dims_loc(1:j-1)) + (1:dims_loc(j)), 1, 1, 1, :) = double(bj);
        end
    end
end

% the uncondensed local modes, as concatenated-local index sets
loc = pde.LocalSpaces;
off_loc = cumsum([0, dims_loc(1:end-1)]);
dirs = {};
for j = 1:numel(loc)
    for r = 1:numel(loc(j).UncondensedModes)
        ix = reshape(loc(j).UncondensedModes{r}, 1, []);
        assert(all(ix >= 1 & ix <= dims_loc(j) & ix == fix(ix)) && numel(unique(ix)) == numel(ix), ...
            "FEProblem:uncondensedMode", ...
            "an uncondensed mode is a set of distinct dof indices of its space.");
        dirs{end+1} = off_loc(j) + ix; %#ok<AGROW>
    end
end
assert(all(cellfun(@isempty, {pde.GlobalSpaces.UncondensedModes})), ...
    "FEProblem:uncondensedMode", "UncondensedModes needs Condense = true.");
all_ix = [dirs{:}];
assert(numel(unique(all_ix)) == numel(all_ix), "FEProblem:uncondensedMode", ...
    "uncondensed modes must be disjoint index sets.");
n_p = numel(dirs);
n_b = double(n_p > 0);                     % the uncondensed block, if any

if n_p == 0
    loc_solve = @(B) pagemldivide(A11d, B);
else
    % the change of basis x = T y (header); only y_w is eliminated
    n_loc = size(A11d, 1);
    piv = cellfun(@(d) d(1), dirs);        % 1 x n_p
    T = eye(n_loc);
    for r = 1:n_p, T(dirs{r}(2:end), piv(r)) = 1; end
    V = T(:, piv);                         % n_loc x n_p, the directions as columns
    w = setdiff(1:n_loc, piv);
    Ay  = pagemtimes(pagemtimes(T.', A11d), T);
    Aww = Ay(w, w, :, :, :);
    Awp = Ay(w, piv, :, :, :);
    rowp = pagemtimes(V.', A11d);          % v' A11, the kept rows in x coordinates
    lift = @(Yw) scatter_rows(Yw, w, n_loc);
    loc_solve = @(B) pagemtimes(T, lift(pagemldivide(Aww, index_rows(pagemtimes(T.', B), w))));
    Yp = lift(pagemldivide(Aww, Awp));
    for r = 1:n_p, Yp(piv(r), r, :, :, :) = -1; end
    S12p = pagemtimes(T, Yp);              % [n_loc, n_p, 1,1, n_el, 1, 1]
    % the uncondensed block, appended last as an element-supported space
    d_g(end+1) = n_p;  f_g(end+1) = 1;
    L.support(end+1)   = "element";
    L.offset(end+1)    = L.n_dof;
    L.count(end+1)     = n_p*n_el;
    L.footprint{end+1} = L.n_dof + reshape(1:n_p*n_el, n_p, n_el);
    L.n_dof            = L.n_dof + n_p*n_el;
end
S1f = loc_solve(A1fd);                     % [n_loc, 1, 1, 1, n_el]

% blocks t, s <= n_g are the global spaces, n_g + 1 the uncondensed block
[row_list, col_list, val_list] = deal(cell(n_g + n_b, n_g + n_b));
rhs_vals = cell(n_g + n_b, 1);
S12 = cell(1, n_g + n_b);                  % A11 \ A12_t, kept for recovery
for t = 1:n_g + n_b
    if t <= n_g
        gid_t = pde.Id.hybrid(t);
        % a global space with no blocks in the equations (e.g. a coupling-
        % owned scalar filled only by transmission rows) condenses to zeros
        if any(trial_of(A12) == gid_t)
            A12t = double(A12{trial_of(A12) == gid_t});   % [n_loc, d_t, 1,1, n_el, 1, f_t]
        else
            A12t = zeros(size(A11d, 1), d_g(t), 1, 1, n_el, 1, f_g(t));
        end
        S12{t} = loc_solve(A12t);
    else
        S12{t} = S12p;
    end
    for s = 1:n_g + n_b
        if s > n_g
            A21s = rowp;                               % [n_p, n_loc, 1,1, n_el]
        elseif any(test_of(A21) == pde.Id.hybrid(s))
            A21s = double(A21{test_of(A21) == pde.Id.hybrid(s)}); % [d_s, n_loc, 1,1, n_el, f_s, 1]
        else
            A21s = zeros(d_g(s), size(A11d, 1), 1, 1, n_el, f_g(s), 1);
        end
        Cst  = -pagemtimes(A21s, S12{t});          % [d_s, d_t, 1,1, n_el, f_s, f_t]
        % the A22 block: assembled for two global spaces, v' A12_t for the
        % p-rows, absent in the uncondensed block's columns
        D = [];
        if s <= n_g && t <= n_g
            A22st = A22{trial_of(A22) == gid_t & test_of(A22) == pde.Id.hybrid(s)};
            if ~isempty(A22st.id), D = double(A22st); end
        elseif t <= n_g
            D = pagemtimes(V.', A12t);                 % [n_p, d_t, 1,1, n_el, 1, f_t]
        end
        if ~isempty(D)
            if size(D, 7) < size(Cst, 7) && size(D, 6) == size(Cst, 6)
                % trace-mass-like block [.., f, 1]: on the face diagonal
                for i = 1:size(Cst, 6)
                    Cst(:,:,:,:,:,i,i) = Cst(:,:,:,:,:,i,i) + D(:,:,:,:,:,i);
                end
            else
                Cst = Cst + D;
            end
        end
        B = reshape(permute(Cst, [1 6 2 7 5 3 4]), d_g(s)*f_g(s), d_g(t)*f_g(t), n_el);
        row_list{s, t} = repmat(reshape(L.footprint{s}, [], 1, n_el), 1, d_g(t)*f_g(t));
        col_list{s, t} = repmat(reshape(L.footprint{t}, 1, [], n_el), d_g(s)*f_g(s), 1);
        val_list{s, t} = B;
    end
    % RHS block of test space t
    if t > n_g
        A21t = rowp;
    elseif any(test_of(A21) == pde.Id.hybrid(t))
        A21t = double(A21{test_of(A21) == pde.Id.hybrid(t)});
    else
        A21t = zeros(d_g(t), size(A11d, 1), 1, 1, n_el, f_g(t), 1);
    end
    Cf_t = -pagemtimes(A21t, S1f);                 % [d_t, 1, 1,1, n_el, f_t]
    if t > n_g
        Cf_t = Cf_t + pagemtimes(V.', A1fd);       % v' A1f
    elseif ~isa(A2f, 'double')
        A2ft = A2f{test_of(A2f) == pde.Id.hybrid(t)};
        if ~isempty(A2ft.id)
            Cf_t = Cf_t + double(A2ft);
        end
    end
    rhs_vals{t} = reshape(permute(Cf_t, [1 6 5 2 3 4]), [], 1);
end
rows = cellfun(@(x) x(:), row_list(:), 'UniformOutput', false);
cols = cellfun(@(x) x(:), col_list(:), 'UniformOutput', false);
vals = cellfun(@(x) x(:), val_list(:), 'UniformOutput', false);
H    = sparse(vertcat(rows{:}), vertcat(cols{:}), vertcat(vals{:}), L.n_dof, L.n_dof);
rhs_rows = cellfun(@(F) reshape(F, [], 1), L.footprint(:), 'UniformOutput', false);
RHS  = accumarray(vertcat(rhs_rows{:}), vertcat(rhs_vals{:}), [L.n_dof, 1]);

S = struct("L", L, "H", H, "RHS", RHS, "A11d", A11d, "S1f", S1f, ...
    "S12", {S12}, "d_g", d_g, "f_g", f_g);
end

function Y = scatter_rows(Yw, w, n)
% Yw on the rows w of n, zeros on the others (7-D pages)
Y = zeros([n, size(Yw, 2:7)]);
Y(w, :, :, :, :, :, :) = Yw;
end

function Y = index_rows(X, w)
Y = X(w, :, :, :, :, :, :);
end
