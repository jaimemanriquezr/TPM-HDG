function [Hc, Fc] = assemble_transmission(obj, k, S1, S2, mesh_1, mesh_2, ...
    part_1, part_2, off, n_dof, sq, tq)
% Assemble the symbolic transmission forms (obj.transmission, built from
% Extrapolate / PathIntegral / bare trials / trace trials over the
% Measure(mesh_i, "dI") interface measure) into the coupled system rows,
% folding the condensed local solves as the built-in *_rows do: x_loc = S1f - S12 x_hat, so condensed parts enter with -R*S12 and
% the rhs with -R*S1f. Terms route by their trial's own mesh; the
% measure's mesh names the test side.
forms = obj.transmission;
[~, wS]  = lgwt(sq, 0, 1);
[~, wTP] = lgwt(tq, 0, 1);
n_p = size(part_1.faces, 3);

sides = struct( ...
    "mesh", {mesh_1, mesh_2}, "part", {part_1, part_2}, ...
    "S", {S1, S2}, "off", {0, off}, "pde", {obj.problems(1), obj.problems(2)});
for i = 1:2
    pde = sides(i).pde;
    S = sides(i).S;
    sides(i).dims_loc = arrayfun(@(sp) sp.Dimension, pde.LocalSpaces);
    sides(i).n_local = numel(pde.LocalSpaces);
    S12f = cell(1, numel(S.S12));
    for t = 1:numel(S.S12)
        S12f{t} = reshape(permute(S.S12{t}, [1 2 7 5 3 4 6]), ...
            size(S.S12{t}, 1), S.d_g(t)*S.f_g(t), []);
    end
    sides(i).S12f = S12f;
end

side_of = @(mesh) 1 + isequal(mesh, mesh_2);
assert(~isequal(mesh_1, mesh_2), "InterfaceProblem:transmission", ...
    "transmission routing needs distinct side meshes");

[rows_acc, cols_acc, vals_acc] = deal({});
rhs_acc = zeros(n_dof, 1);

for f = 1:numel(forms)
    form = forms(f);
    assert(~isempty(form.test) && ~isempty(form.measure), ...
        "InterfaceProblem:transmission", ...
        "each transmission form needs a test (expr*mu) and a dI measure");
    it = side_of(form.measure.domain);
    ts = sides(it);
    test_basis = form.test.basis;
    [~, ~, ~, test_id] = read_signature(form.test);
    test_space_ix = test_id - ts.n_local;
    tr_test = ts.pde.GlobalSpaces(test_space_ix);
    scalar_test = [test_basis.Support] == "scalar";
    if scalar_test
        lt_fun = @(s_) ones(numel(s_), 1);
        n_ct = 1; d1s = 1;
    else
        [lt_fun, n_ct, d1s] = trace_basis_info(test_basis);
    end
    d1t = n_ct * d1s;

    for e = 1:n_p
        F_t = ts.part.faces(e);
        if scalar_test
            rows_e = ts.S.L.offset(test_space_ix) + 1 + ts.off;
        else
            rows_e = get_global_dof(tr_test, subset=F_t);
            rows_e = rows_e(:) + ts.off;
        end
        s_test = reshape(ts.part.quadrature.local(1, :, e), [], 1);
        l_test = lt_fun(s_test);
        patch_len = ts.part.area(e);
        int_test = @(B) int_blocked(B, l_test, wS, patch_len, n_ct, d1s);

        R = {zeros(d1t, sum(sides(1).dims_loc)), ...
             zeros(d1t, sum(sides(2).dims_loc))};
        touched = [false, false];

        for m = 1:numel(form.terms)
            term = form.terms(m);

            if term.kind == "data"
                % interface source at the test side's physical points,
                % moved to the rhs (sign flip)
                x_if = ts.part.quadrature.interface(:, :, e);
                B = eval_datum(term.trial, x_if, n_ct);
                rhs_acc(rows_e) = rhs_acc(rows_e) ...
                    - reshape(int_test(term.coeff * B), [], 1);
                continue
            end

            is_ = side_of(term.trial.mesh);
            s = sides(is_);
            K = s.part.elements(e);
            coeff = term.coeff;
            if term.stab
                coeff = coeff * s.part.tau(e);
            end

            if term.kind == "trace_data"
                % skeleton data (kinematic history): rhs -= coeff * <h l, mu>
                F_d = s.part.faces(e);
                s_tr = reshape(s.part.quadrature.local(1, :, e), [], 1);
                [ltr_fun, n_ctr, d1tr] = trace_basis_info(term.trial.basis);
                assert(n_ctr == n_ct, "InterfaceProblem:transmission", ...
                    "trace data and test component counts differ");
                l_tr = ltr_fun(s_tr);
                h_dof = term.trial.dof(:, F_d);
                B = zeros(sq, 1, n_ct);
                for c = 1:n_ct
                    B(:, 1, c) = l_tr * h_dof((c-1)*d1tr + (1:d1tr));
                end
                rhs_acc(rows_e) = rhs_acc(rows_e) ...
                    - reshape(int_test(coeff * B), [], 1);
                continue
            end

            if term.kind == "trace" || term.kind == "scalar"
                assert(isempty(term.weight), "InterfaceProblem:transmission", ...
                    "matrix weights on trace/scalar trials not supported");
                [~, ~, ~, id] = read_signature(term.trial);
                if term.kind == "scalar"
                    cols_e = s.S.L.offset(id - s.n_local) + 1 + s.off;
                    l_tr = ones(sq, 1);
                    n_ctr = 1; d1tr = 1;
                else
                    tr_sp = s.pde.GlobalSpaces(id - s.n_local);
                    cols_e = get_global_dof(tr_sp, subset=s.part.faces(e));
                    cols_e = cols_e(:) + s.off;
                    s_tr = reshape(s.part.quadrature.local(1, :, e), [], 1);
                    [ltr_fun, n_ctr, d1tr] = trace_basis_info(term.trial.basis);
                    l_tr = ltr_fun(s_tr);
                end
                % pointwise weights per (test comp, trial comp) pair
                W_ct = trace_pair_weights(term.normal, n_ct, n_ctr, ...
                    outward_normal(s.part, e, is_), sq);
                M = zeros(d1t, n_ctr*d1tr);
                for c = 1:n_ct
                    for cc = 1:n_ctr
                        w_pt = W_ct(:, c, cc);
                        if ~any(w_pt), continue; end
                        M((c-1)*d1s + (1:d1s), (cc-1)*d1tr + (1:d1tr)) = ...
                            patch_len * (l_test' * ((wS .* w_pt) .* (coeff * l_tr)));
                    end
                end
                rows_acc{end+1} = repmat(rows_e, 1, numel(cols_e)); %#ok<AGROW>
                cols_acc{end+1} = repmat(cols_e.', numel(rows_e), 1); %#ok<AGROW>
                vals_acc{end+1} = M; %#ok<AGROW>
                continue
            end

            switch term.kind
                case "onface"
                    rs = s.part.quadrature.local_2d(:, :, e);
                    B = eval_local(term.trial, rs, s.mesh, K);
                case "extrapolate"
                    rs = s.part.quadrature.interface_2d(:, :, e);
                    B = eval_local(term.trial, rs, s.mesh, K);
                case "path"
                    rs = reshape(s.part.quadrature.path_2d(:, :, e, :), 2, []);
                    Bp = apply_op_weight(eval_local(term.trial, rs, s.mesh, K), ...
                        term.op, term.weight);
                    m_vec = s.part.quadrature.interface(:, :, e) ...
                          - s.part.quadrature.physical(:, :, e);      % 2 x sq
                    d = size(Bp, 2);
                    % contract the last value index with m: a tensor
                    % integrand (grad u) yields a vector, a vector (q)
                    % yields a scalar
                    if size(Bp, 4) > 1
                        B = zeros(sq, d, size(Bp, 3));
                        for c = 1:size(Bp, 3)
                            for j = 1:size(Bp, 4)
                                Dp = reshape(Bp(:, :, c, j), sq, tq, d);
                                B(:, :, c) = B(:, :, c) + reshape( ...
                                    sum(Dp .* reshape(wTP, 1, tq), 2), sq, d) ...
                                    .* m_vec(j, :).';
                            end
                        end
                    elseif size(Bp, 3) == 2
                        B = zeros(sq, d);
                        for c = 1:size(Bp, 3)
                            Dp = reshape(Bp(:, :, c), sq, tq, d);
                            B = B + reshape(sum(Dp .* reshape(wTP, 1, tq), 2), ...
                                sq, d) .* m_vec(c, :).';
                        end
                    else
                        % scalar integrand: arc-length integral
                        % int expr(x + t m) |m| dt (the alpha rows)
                        Dp = reshape(Bp(:, :, 1), sq, tq, d);
                        B = reshape(sum(Dp .* reshape(wTP, 1, tq), 2), sq, d) ...
                            .* vecnorm(m_vec).';
                    end
            end
            if term.kind ~= "path"
                B = apply_op_weight(B, term.op, term.weight);
            end

            switch term.normal
                case "dot"      % vector . n -> scalar
                    n_vec = outward_normal(s.part, e, is_);
                    assert(size(B, 3) == 2 && size(B, 4) <= 1, ...
                        "InterfaceProblem:transmission", ...
                        "dot(., InterfaceNormal) expects a vector value");
                    B = B(:, :, 1) .* n_vec(1, :).' + B(:, :, 2) .* n_vec(2, :).';
                case "tangent"  % (t x t) v, t = rotated outward normal
                    n_vec = outward_normal(s.part, e, is_);
                    t_vec = [0, -1; 1, 0] * n_vec;
                    assert(size(B, 3) == 2 && size(B, 4) <= 1, ...
                        "InterfaceProblem:transmission", ...
                        "TangentialProjection expects a vector value");
                    tv = B(:, :, 1) .* t_vec(1, :).' + B(:, :, 2) .* t_vec(2, :).';
                    B = cat(3, tv .* t_vec(1, :).', tv .* t_vec(2, :).');
                case "mtimes"   % sigma*n -> vector; scalar p*n -> vector
                    n_vec = outward_normal(s.part, e, is_);
                    if size(B, 4) > 1          % tensor: contract 2nd index
                        Bn = zeros(sq, size(B, 2), size(B, 3));
                        for c = 1:size(B, 3)
                            for j = 1:size(B, 4)
                                Bn(:, :, c) = Bn(:, :, c) ...
                                    + B(:, :, c, j) .* n_vec(j, :).';
                            end
                        end
                        B = Bn;
                    elseif size(B, 3) == 1     % scalar: promote to p n
                        B = cat(3, B .* n_vec(1, :).', B .* n_vec(2, :).');
                    else
                        error("InterfaceProblem:transmission", ...
                            "vector*InterfaceNormal is ambiguous; use dot()");
                    end
            end

            [~, ~, ~, id] = read_signature(term.trial);
            off_loc = sum(s.dims_loc(1:id-1));
            cols = off_loc + (1:s.dims_loc(id));
            B = reshape(B, size(B, 1), size(B, 2), []);
            R{is_}(:, cols) = R{is_}(:, cols) + int_test(coeff * B);
            touched(is_) = true;
        end

        for is_ = 1:2
            if ~touched(is_), continue; end
            s = sides(is_);
            K = s.part.elements(e);
            S12K = cell2mat(cellfun(@(A) A(:, :, K), s.S12f, ...
                'UniformOutput', false));
            foot = cell2mat(cellfun(@(fp) fp(:, K), ...
                s.S.L.footprint(:).', 'UniformOutput', false).') + s.off;
            V = -R{is_} * S12K;
            rows_acc{end+1} = repmat(rows_e, 1, size(V, 2)); %#ok<AGROW>
            cols_acc{end+1} = repmat(foot.', d1t, 1); %#ok<AGROW>
            vals_acc{end+1} = V; %#ok<AGROW>
            rhs_acc(rows_e) = rhs_acc(rows_e) - R{is_} * s.S.S1f(:, :, 1, 1, K);
        end
    end
end

rows = cellfun(@(v) v(:), rows_acc.', 'UniformOutput', false);
cols = cellfun(@(v) v(:), cols_acc.', 'UniformOutput', false);
vals = cellfun(@(v) v(:), vals_acc.', 'UniformOutput', false);
Hc = sparse(vertcat(rows{:}), vertcat(cols{:}), vertcat(vals{:}), n_dof, n_dof);
Fc = rhs_acc;
end

%% =====================================================================
function M = int_blocked(B, l_test, wS, patch_len, n_ct, d1s)
% component-blocked test integration: B is (sq x d x n_c) with
% n_c == n_ct; rows ordered [comp_1; comp_2; ...]
n_cB = size(B, 3);
assert(n_cB == n_ct, "InterfaceProblem:transmission", ...
    "term has %d value components, test expects %d", n_cB, n_ct);
M = zeros(n_ct * d1s, size(B, 2));
for c = 1:n_ct
    M((c-1)*d1s + (1:d1s), :) = patch_len * (l_test' * (wS .* B(:, :, c)));
end
end

function W = trace_pair_weights(mode, n_ct, n_ctr, n_vec, sq)
% pointwise (test comp, trial comp) weights for trace/scalar columns
t_vec = [0, -1; 1, 0] * n_vec;
W = zeros(sq, n_ct, n_ctr);
switch mode
    case ""        % identity pairing (component-diagonal)
        assert(n_ct == n_ctr, "InterfaceProblem:transmission", ...
            "trace trial and test component counts differ");
        for c = 1:n_ct, W(:, c, c) = 1; end
    case "dot"     % vector trace . n against a scalar test
        assert(n_ct == 1 && n_ctr == 2, "InterfaceProblem:transmission", ...
            "dot(trace, n) expects a vector trace and scalar test");
        W(:, 1, 1) = n_vec(1, :).';
        W(:, 1, 2) = n_vec(2, :).';
    case "mtimes"  % scalar trace/global promoted by n against vector test
        assert(n_ct == 2 && n_ctr == 1, "InterfaceProblem:transmission", ...
            "trace*n expects a scalar trial and vector test");
        W(:, 1, 1) = n_vec(1, :).';
        W(:, 2, 1) = n_vec(2, :).';
    case "tangent" % (t x t) pairing of a vector trace with a vector test
        assert(n_ct == 2 && n_ctr == 2, "InterfaceProblem:transmission", ...
            "TangentialProjection on traces expects vector trace and test");
        for c = 1:2
            for cc = 1:2
                W(:, c, cc) = t_vec(c, :).' .* t_vec(cc, :).';
            end
        end
end
end

function B = eval_datum(datum, x_phys, n_ct)
% source values at physical interface points, (sq x 1 x n_ct)
v = datum.value;
if ~iscell(v), v = {v}; end
assert(numel(v) == n_ct, "InterfaceProblem:transmission", ...
    "source has %d components, test expects %d", numel(v), n_ct);
B = zeros(size(x_phys, 2), 1, n_ct);
for c = 1:n_ct
    B(:, 1, c) = reshape(v{c}(x_phys(1, :), x_phys(2, :)) ...
        + 0*x_phys(1, :), [], 1);
end
end

function n_vec = outward_normal(part, e, side)
n_vec = part.quadrature.normals(:, :, e);
if side == 1, n_vec = -n_vec; end   % side-1 outward (FSI convention)
end

function [l_fun, n_comp, d1s] = trace_basis_info(basis)
% scalar evaluator, component count and scalar dof count of a 1-D trace
% basis (components share the scalar basis; per-face dof layout
% [comp_1; comp_2; ...])
switch class(basis)
    case "LegendreBasis"
        l_fun = @(s) legendre_basis(basis.Order, s);
        n_comp = 1;
        d1s = basis.Dimension;
    case "TensorFEBasis"
        sb = basis.ScalarBasis;
        assert(isa(sb, "LegendreBasis"), "InterfaceProblem:transmission", ...
            "trace basis Tensor(%s) not yet supported symbolically", class(sb));
        l_fun = @(s) legendre_basis(sb.Order, s);
        n_comp = numel(basis.LinearMaps);
        d1s = sb.Dimension;
    otherwise
        error("InterfaceProblem:transmission", ...
            "trace basis %s not yet supported symbolically", class(basis));
end
end

function B = apply_op_weight(B, op, W)
% pointwise operator then left matrix action on the value dims
if op == "trace"
    assert(size(B, 3) == size(B, 4), "InterfaceProblem:transmission", ...
        "trace needs a square tensor value");
    Bt = B(:, :, 1, 1);
    for c = 2:size(B, 3)
        Bt = Bt + B(:, :, c, c);
    end
    B = Bt;
end
if isempty(W)
    return
end
if size(B, 3) == 1 && size(B, 4) <= 1
    % matrix weight on a scalar value promotes to a tensor (W * s)
    Bw = zeros(size(B, 1), size(B, 2), size(W, 1), size(W, 2));
    for c = 1:size(W, 1)
        for j = 1:size(W, 2)
            Bw(:, :, c, j) = W(c, j) * B(:, :, 1);
        end
    end
else
    assert(size(B, 4) <= 1, "InterfaceProblem:transmission", ...
        "matrix weights on tensor values are not supported");
    Bw = zeros(size(B, 1), size(B, 2), size(W, 1));
    for c = 1:size(W, 1)
        for j = 1:size(W, 2)
            Bw(:, :, c) = Bw(:, :, c) + W(c, j) * B(:, :, j);
        end
    end
end
B = Bw;
end

function B = eval_local(trial, rs, mesh, K)
% local trial basis values at reference points of element K,
% value shape kept: (n_q x d x c1 x c2)
V = evaluate(trial, rs, mesh, element_list = K);
B = reshape(V, size(V, 1), size(V, 2), size(V, 3), size(V, 4));
end
