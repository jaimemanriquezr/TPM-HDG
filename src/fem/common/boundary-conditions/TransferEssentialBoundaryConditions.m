classdef TransferEssentialBoundaryConditions < WeakEssentialBoundaryConditions
    % Transfer Path Method (Cockburn & Solano) essential boundary conditions.
    %
    % The mesh boundary Gamma_h does not coincide with the true boundary
    % Gamma. For a point x on a boundary face let x_bar = mapping(x) be the
    % associated point on Gamma (e.g. the closest point). Integrating the
    % gradient along the straight transfer path from x to x_bar,
    %
    %   u(x_bar) = u(x) + int_0^1 grad(u)(x + s*m) m ds,   m = x_bar - x,
    %
    % and replacing grad(u) by its PDE-specific reconstruction G_h from the
    % local field unknowns (extended from the element adjacent to the face
    % by polynomial extrapolation), the essential datum g on Gamma is
    % transferred to the trace unknown u_hat on Gamma_h, per component:
    %
    %   u_hat_a(x) + int_0^1 [G_h(x + s*m) m]_a ds = g_a(x_bar).
    %
    % Tested against the trace space on each boundary face, this couples
    % u_hat of the face and the field unknowns of the adjacent element; the
    % field unknowns are eliminated with the same local solver the HDG
    % scheme uses (Schur complement), so evaluate() returns rows over
    % hybrid unknowns only, in the weak-essential format.
    %
    % The gradient reconstruction is declared through the path_gradient
    % property, in one of two forms.
    %
    % (a) Symbolic: an expression in a trial function of a local space
    %     whose pointwise value is the gradient itself, e.g. for Darcy
    %     (q = -kappa*grad u, q the local trial function of field 1)
    %
    %         path_gradient = -kappa \ q
    %
    %     The expression must be linear in exactly one local trial
    %     function with spatially constant coefficients; its value per
    %     basis function is contracted with the transfer direction m
    %     internally. A column-vector value (2 x 1) declares the gradient
    %     of a scalar unknown; an (n_comp x 2) matrix value the gradient
    %     of a vector unknown.
    %
    % (b) A function handle (see the static factories darcy / elasticity /
    %     velocity_gradient / field_gradient), with signature
    %
    %   Gm = path_gradient(problem, mesh, rs, K_list, m, n_local)
    %
    %   rs      (2 x N x n_b)  local coordinates of the path points in the
    %                          element K_list(b) adjacent to boundary face b
    %   m       (2 x N x n_b)  transfer direction x_bar - x at each point
    %   Gm      (N x n_local x n_b x n_comp)  value of [G_phi(x) m]_a for
    %           each local field basis function phi (columns match the
    %           local field dof layout of the condensed equations)
    properties
        mapping function_handle = @(X) X       % x on Gamma_h -> x_bar on Gamma
        path_gradient                          % expression or handle, see class help
        quadrature_order double = []           % face rule degree (default 2k+2)
    end

    methods
        function obj = TransferEssentialBoundaryConditions(values, tags, mapping, options)
            arguments
                values Datum = Datum.empty;
                tags string = string.empty;
                mapping {mustBeA(mapping, ["function_handle", "LevelSetCurve"])} = @(X) X;
                options.path_gradient {mustBeA(options.path_gradient, ...
                    ["function_handle", "AbstractFEFunction"])} = function_handle.empty;
                options.quadrature_order double = [];
            end
            obj = obj@WeakEssentialBoundaryConditions(values, tags);
            if isa(mapping, "LevelSetCurve")
                mapping = @(X) mapping.project(X);
            end
            obj.mapping = mapping;
            obj.quadrature_order = options.quadrature_order;
            if isempty(options.path_gradient)
                error('TransferEssentialBoundaryConditions:path_gradient', ...
                      ['No gradient reconstruction declared: use one of the ', ...
                       'factories (darcy, elasticity, velocity_gradient) or ', ...
                       'pass path_gradient explicitly.']);
            end
            obj.path_gradient = options.path_gradient;
        end

        function [blocks, index] = get_transfer_blocks(obj, problem)
            % Transfer rows before the local unknowns are eliminated:
            %
            %     A31 * z + A32 * u_hat = A3f
            %
            % over the local field dofs z of the adjacent element and the
            % trace dofs of the tagged face itself. The blocks do not
            % involve the local solver A11, so they stay valid when A11
            % depends on a parameter (e.g. A11 - lambda*B in an
            % eigenproblem). evaluate() is the condensed form of these rows.
            %
            % INDEX carries the bookkeeping needed to place them: .dofs
            % (global trace dofs of the tagged faces), .faces, .K (adjacent
            % element), .sn (which of the element's 3 sides the face is),
            % and .D1 (trace dofs per face).
            arguments
                obj TransferEssentialBoundaryConditions
                problem FEProblem
            end
            mesh = problem.Mesh;
            % the trace space: the skeleton-supported global space (a layout
            % may append element means and scalar multipliers after it)
            L = global_dof_layout(problem.GlobalSpaces, mesh);
            hybrid_basis = problem.GlobalSpaces(L.support == "face");
            assert(isscalar(hybrid_basis), 'TransferEssentialBoundaryConditions:trace', ...
                'expected exactly one skeleton-supported global space');
            hybrid_el = hybrid_basis.Basis;
            if isa(hybrid_el, "TensorFEBasis")
                trace_basis = hybrid_el.ScalarBasis;
                n_comp = numel(hybrid_el.LinearMaps);
            else
                trace_basis = hybrid_el;
                n_comp = 1;
            end
            d1 = trace_basis.Dimension;
            D1 = n_comp * d1;
            k = trace_basis.Order;

            % width of the local field block; taken from the spaces rather
            % than from size(A11, 1) so that this method never assembles the
            % local problem
            n_local = sum(arrayfun(@(b) b.Dimension, problem.LocalSpaces));
            n_el = size(mesh.elements, 2);

            % quadrature rules: the face rule integrates (trace basis) x
            % (path-integrated degree-k reconstruction) products; the path
            % rule integrates the degree-k reconstruction along the path
            if isempty(obj.quadrature_order)
                face_order = 2*k + 2;
            else
                face_order = obj.quadrature_order;
            end
            [s_f, w_f] = get_1d_reference_quadrature(face_order);
            n_q = length(s_f);
            [s_p, w_p] = get_1d_reference_quadrature(k);
            n_p = length(s_p);

            % trace basis in the canonical face parametrization (this is the
            % parametrization the global trace dofs are expressed in)
            l = trace_basis.evaluate(s_f);      % (n_q x d1)
            W_l = w_f(:) .* l;

            % physical face points for all faces at once
            x_face_all = get_physical_points(mesh, s_f);   % (2 x n_q x n_faces)

            n_bc = length(obj.values);
            [bc_dofs, bc_face_list, bc_K, bc_sn] = deal(cell(n_bc, 1));

            % transfer blocks over all elements (zero away from bc faces)
            A31 = zeros(D1, n_local, 1, 1, n_el, 3);
            A32 = zeros(D1, D1, 1, 1, n_el, 3);
            A3f = zeros(D1, 1, 1, 1, n_el, 3);

            for ix = 1:n_bc
                bfaces = get_face_by_tag(mesh, obj.tags(ix));
                n_b = length(bfaces);
                if n_b == 0; continue; end
                [K_b, sn_b] = get_element_by_face(mesh, bfaces);
                if iscell(K_b)
                    K_b = cellfun(@(C) C(1), K_b);
                    sn_b = cellfun(@(C) C(1), sn_b);
                end

                x_face = x_face_all(:, :, bfaces);         % (2 x n_q x n_b)
                x_bar = obj.mapping(reshape(x_face, 2, []));
                x_bar = reshape(x_bar, 2, n_q, n_b);
                m = x_bar - x_face;                        % (2 x n_q x n_b)

                % path points x + s*m, path parameter (dim 1) fastest,
                % and the (constant-along-path) direction at each point
                x_path = reshape(x_face, 2, 1, n_q, n_b) ...
                       + reshape(s_p, 1, n_p) .* reshape(m, 2, 1, n_q, n_b);
                x_path = reshape(x_path, 2, n_p*n_q, n_b);
                m_path = repmat(reshape(m, 2, 1, n_q, n_b), 1, n_p);
                m_path = reshape(m_path, 2, n_p*n_q, n_b);
                rs_path = mesh.get_2d_local(x_path, K_b(:));

                % PDE-specific gradient reconstruction along the paths
                if isa(obj.path_gradient, 'function_handle')
                    Gm = obj.path_gradient(problem, mesh, rs_path, K_b, m_path, n_local);
                else
                    Gm = expression_gradient(obj.path_gradient, problem, mesh, ...
                                             rs_path, K_b, m_path, n_local);
                end
                Gm = reshape(Gm, n_p, n_q, n_local, n_b, size(Gm, 4));

                % path integral: PI(q, :, b, a) = int_0^1 [G_phi m]_a ds
                PI = reshape(sum(w_p(:) .* Gm, 1), n_q, n_local, n_b, n_comp);

                area_b = reshape(mesh.area(bfaces), 1, 1, n_b);

                % transfer coupling, face mass, datum moments
                T = area_b .* pagemtimes(W_l, 'transpose', PI, 'none');
                M_face = area_b .* pagemtimes(W_l, 'transpose', l, 'none');
                M_face = pagekron(eye(n_comp), M_face);

                g_comp = components(obj.values(ix), n_comp);
                G_face = zeros(D1, 1, n_b);
                for a = 1:n_comp
                    g_val = evaluate_datum(g_comp{a}, x_bar);   % (n_q x n_b)
                    G_face((a-1)*d1 + (1:d1), 1, :) = ...
                        area_b .* reshape(W_l.' * g_val, d1, 1, n_b);
                end

                for b = 1:n_b
                    e = K_b(b); i = sn_b(b);
                    for a = 1:n_comp
                        A31((a-1)*d1 + (1:d1), :, 1, 1, e, i) = T(:, :, b, a);
                    end
                    A32(:, :, 1, 1, e, i) = M_face(:, :, b);
                    A3f(:, :, 1, 1, e, i) = G_face(:, :, b);
                end

                bc_face_list{ix} = bfaces(:).';
                bc_K{ix} = K_b(:).';
                bc_sn{ix} = sn_b(:).';
                bc_dofs{ix} = get_global_dof(hybrid_basis, subset=bfaces);
            end

            blocks = struct("A31", A31, "A32", A32, "A3f", A3f);
            index = struct( ...
                "dofs",  [bc_dofs{:}], ...
                "faces", [bc_face_list{:}], ...
                "K",     [bc_K{:}], ...
                "sn",    [bc_sn{:}], ...
                "D1",    D1);
        end

        function [dofs, values_H, values_g, blocks] = evaluate(obj, problem, options)
            % BLOCKS = the get_transfer_blocks output (A31, A32, A3f and
            % index); pass it back as options.blocks to reuse it on the
            % same problem.
            arguments
                obj TransferEssentialBoundaryConditions
                problem FEProblem
                options.S struct = struct.empty       % condensed system; empty: condense here
                options.blocks struct = struct.empty  % get_transfer_blocks output; empty: build here
            end
            mesh = problem.Mesh;
            n_faces = size(mesh.faces, 2);
            if isempty(options.blocks)
                [blocks, index] = obj.get_transfer_blocks(problem);
                blocks.index = index;
            else
                blocks = options.blocks;
            end
            [A31, A32, A3f, index] = deal(blocks.A31, blocks.A32, blocks.A3f, blocks.index);
            D1 = index.D1;
            % the condensed local solver, per global space t (as
            % condense_global_system): x_loc = S1f - sum_t S12{t} x_t
            if isempty(options.S)
                S = condense_global_system(problem, mesh);
            else
                S = options.S;
            end
            S12 = S.S12;  S1f = S.S1f;
            Lay = S.L;
            n_g = numel(problem.GlobalSpaces);
            assert(numel(S12) == n_g, "TransferEssentialBoundaryConditions:uncondensedModes", ...
                "uncondensed local modes are not supported with transfer conditions");
            d_g = S.d_g;  f_g = S.f_g;
            trace_space = find(Lay.support == "face", 1);
            g_el = A3f - pagemtimes(A31, S1f);  % (D1 x 1 x 1 x 1 x n_el x 3)

            dofs = index.dofs;
            all_faces = index.faces;
            all_K = index.K;
            all_sn = index.sn;
            n_all = numel(all_faces);
            N = Lay.n_dof;
            FD = get_global_dof(problem.GlobalSpaces(trace_space), subset=1:n_faces);
            FD = reshape(FD, D1, n_faces);
            [rows_H, cols_H, vals_H] = deal(cell(n_g, 1));
            for t = 1:n_g
                % (D1 x d_t x 1 x 1 x n_el x 3 x f_t): rows over the transfer
                % face i of element e, columns over space t's slots
                H_t = -pagemtimes(A31, S12{t});
                if t == trace_space
                    for i = 1:3
                        H_t(:, :, 1, 1, :, i, i) = H_t(:, :, 1, 1, :, i, i) + A32(:, :, 1, 1, :, i);
                    end
                end
                nt = d_g(t)*f_g(t);
                [rt, ct, vt] = deal(zeros(D1, nt, n_all));
                for b = 1:n_all
                    e = all_K(b); i = all_sn(b);
                    rt(:, :, b) = repmat(FD(:, all_faces(b)), 1, nt);
                    ct(:, :, b) = repmat(reshape(Lay.footprint{t}(:, e), 1, nt), D1, 1);
                    vt(:, :, b) = reshape(H_t(:, :, 1, 1, e, i, :), D1, nt);
                end
                rows_H{t} = rt(:);  cols_H{t} = ct(:);  vals_H{t} = vt(:);
            end
            values_g = zeros(D1, n_all);
            for b = 1:n_all
                values_g(:, b) = g_el(:, 1, 1, 1, all_K(b), all_sn(b));
            end
            values_H = sparse(vertcat(rows_H{:}), vertcat(cols_H{:}), vertcat(vals_H{:}), N, N);
            values_g = values_g(:);
        end
    end
    methods (Static)
        function bc = darcy(values, tags, mapping, kappa, options)
            % Darcy: grad(u) = -kappa^{-1} q, q the (vector) flux field 1
            arguments
                values Datum
                tags string
                mapping {mustBeA(mapping, ["function_handle", "LevelSetCurve"])} = @(X) X;
                kappa double = eye(2);
                options.quadrature_order double = [];
            end
            pg = @(problem, mesh, rs, K_list, m, n_local) ...
                darcy_gradient(problem, mesh, rs, K_list, m, n_local, kappa);
            bc = TransferEssentialBoundaryConditions(values, tags, mapping, ...
                path_gradient = pg, quadrature_order = options.quadrature_order);
        end

        function bc = elasticity(values, tags, mapping, E, nu, options)
            % PEERS elasticity: grad(u) = c1 sigma - c2 tr(sigma) I + gamma,
            % gamma = p [0 1; -1 0] the antisymmetric multiplier (field 3),
            % sigma the PEERS stress (field 1);
            % c1 = 1/(2 mu), c2 = lambda/(4 mu (lambda + mu)).
            % The incompressible (Stokes stress) limit is c2 -> 1/(4 mu).
            arguments
                values Datum
                tags string
                mapping {mustBeA(mapping, ["function_handle", "LevelSetCurve"])} = @(X) X;
                E double = 1.0;
                nu double = 0.3;
                options.quadrature_order double = [];
            end
            lambda = E*nu/((1 + nu)*(1 - 2*nu));
            mu = E/2/(1 + nu);
            c1 = 1/(2*mu);
            c2 = lambda/(4*mu*(lambda + mu));
            pg = @(problem, mesh, rs, K_list, m, n_local) ...
                elasticity_gradient(problem, mesh, rs, K_list, m, n_local, c1, c2);
            bc = TransferEssentialBoundaryConditions(values, tags, mapping, ...
                path_gradient = pg, quadrature_order = options.quadrature_order);
        end

        function bc = incompressible_stress(values, tags, mapping, mu, options)
            % Stress-velocity Stokes / Navier-Stokes (IncompressibleElasticity
            % Matrices): the elasticity reconstruction in the incompressible
            % limit, grad(u) = (1/2mu) sigma - (1/4mu) tr(sigma) I + gamma.
            arguments
                values Datum
                tags string
                mapping {mustBeA(mapping, ["function_handle", "LevelSetCurve"])} = @(X) X;
                mu double = 1.0;
                options.quadrature_order double = [];
            end
            pg = @(problem, mesh, rs, K_list, m, n_local) ...
                elasticity_gradient(problem, mesh, rs, K_list, m, n_local, 1/(2*mu), 1/(4*mu));
            bc = TransferEssentialBoundaryConditions(values, tags, mapping, ...
                path_gradient = pg, quadrature_order = options.quadrature_order);
        end

        function bc = velocity_gradient(values, tags, mapping, options)
            % Velocity-gradient (Stokes L = grad u) formulation: field 1 is
            % a TensorFEBasis whose value is grad(u) itself, so
            % [G_phi m]_a = [A_i m]_a phi_j for linear map A_i.
            arguments
                values Datum
                tags string
                mapping {mustBeA(mapping, ["function_handle", "LevelSetCurve"])} = @(X) X;
                options.quadrature_order double = [];
            end
            pg = @(problem, mesh, rs, K_list, m, n_local) ...
                tensor_field_gradient(problem, mesh, rs, K_list, m, n_local);
            bc = TransferEssentialBoundaryConditions(values, tags, mapping, ...
                path_gradient = pg, quadrature_order = options.quadrature_order);
        end

        function bc = field_gradient(values, tags, mapping, options)
            % The transferred unknown is field 1 itself (a scalar or a
            % to_vector Dubiner field): grad(u) is the gradient of its
            % local polynomial, [G_phi m]_a = grad(phi_j) . m in the dof
            % block of component a.
            arguments
                values Datum
                tags string
                mapping {mustBeA(mapping, ["function_handle", "LevelSetCurve"])} = @(X) X;
                options.quadrature_order double = [];
            end
            pg = @(problem, mesh, rs, K_list, m, n_local) ...
                polynomial_field_gradient(problem, mesh, rs, K_list, m, n_local);
            bc = TransferEssentialBoundaryConditions(values, tags, mapping, ...
                path_gradient = pg, quadrature_order = options.quadrature_order);
        end
    end
end

%% ----------------------- gradient reconstructions ----------------------- %%

function Gm = expression_gradient(expr, problem, mesh, rs, K_b, m, n_local)
    % Lower a symbolic gradient expression to the transfer-path values
    % [G_phi m] in the local dof layout. Sums of terms (LinearCombination)
    % lower term by term into their field blocks; each term must be linear
    % in a single local trial function with constant coefficients,
    % e.g. -kappa \ q, or c1*sigma - c2*trace(sigma)*eye(2) + ...
    if isa(expr, 'LinearCombination')
        Gm = expr.coefficients(1) * ...
            expression_gradient(expr.summands(1), problem, mesh, rs, K_b, m, n_local);
        for i = 2:numel(expr.summands)
            Gm_i = expr.coefficients(i) * ...
                expression_gradient(expr.summands(i), problem, mesh, rs, K_b, m, n_local);
            assert(isequal(size(Gm_i), size(Gm)), ...
                'TransferEssentialBoundaryConditions:path_gradient', ...
                'All terms of the gradient sum must have the same value shape.');
            Gm = Gm + Gm_i;
        end
        return
    end
    [~, ~, ~, id] = read_signature(expr);
    assert(id >= 1 && id <= length(problem.LocalSpaces), ...
        'TransferEssentialBoundaryConditions:path_gradient', ...
        'The gradient expression must be built from a local trial function.');
    dims = arrayfun(@(b) b.Dimension, problem.LocalSpaces);
    offset = sum(dims(1:id-1));

    N = size(rs, 2); n_b = size(rs, 3);
    if problem.LocalSpaces(id).IsPiolaTransformed
        % Element-dependent basis (e.g. PEERS): page b of rs evaluated with
        % element K_b(b)'s map -> (N x d x out1 x out2 x n_b), pages kept
        P = expr.evaluate(rs, mesh, element_list = K_b(:));
        d = size(P, 2);
        assert(d == dims(id), ...
            'TransferEssentialBoundaryConditions:path_gradient', ...
            'The gradient expression must be linear in the trial function.');
        assert(size(P, 5) == n_b, ...
            'TransferEssentialBoundaryConditions:path_gradient', ...
            'Element-dependent term lost its page dimension.');
        P = permute(P, [1 5 2 3 4]);            % (N x n_b x d x out1 x out2)
    else
        P = expr.evaluate(reshape(rs, 2, []));  % (N*n_b x d x out1 x out2)
        d = size(P, 2);
        assert(d == dims(id), ...
            'TransferEssentialBoundaryConditions:path_gradient', ...
            'The gradient expression must be linear in the trial function.');
        P = reshape(P, N, n_b, d, size(P, 3), size(P, 4));
    end
    if size(P, 5) == 1
        % column-vector value (2 x 1): gradient of a scalar unknown;
        % otherwise an (n_comp x 2) matrix, contracted with m over dim 5
        P = permute(P, [1 2 3 5 4]);
    end
    n_comp = size(P, 4);

    mm = permute(m, [2 3 4 5 1]);               % (N x n_b x 1 x 1 x 2)
    G = sum(P .* mm, 5);                        % (N x n_b x d x n_comp)
    Gm = zeros(N, n_local, n_b, n_comp);
    Gm(:, offset + (1:d), :, :) = permute(G, [1 3 2 4]);
end

function Gm = darcy_gradient(problem, mesh, rs, K_list, m, n_local, kappa) %#ok<INUSD>
    % grad(u) m = -(kappa^{-1} q) . m = -(kappa^{-T} m) . q
    flux_basis = problem.LocalSpaces(1).Basis;
    if isa(flux_basis, 'TensorFEBasis')
        scalar_basis = flux_basis.ScalarBasis;
    else
        scalar_basis = flux_basis;
    end
    d = scalar_basis.Dimension;
    N = size(rs, 2); n_b = size(rs, 3);

    Phi = double(scalar_basis.evaluate(reshape(rs, 2, [])));   % (N*n_b x d)
    Phi = reshape(Phi, N, n_b, d);
    a = reshape(kappa.' \ reshape(m, 2, []), 2, N, n_b);

    Gm = zeros(N, n_local, n_b);
    for c = 1:2
        Gm(:, (c-1)*d + (1:d), :) = ...
            -reshape(a(c, :, :), N, 1, n_b) .* permute(Phi, [1 3 2]);
    end
end

function Gm = elasticity_gradient(problem, mesh, rs, K_list, m, n_local, c1, c2)
    % grad(u) m = (c1 sigma - c2 tr(sigma) I) m + gamma m,
    % gamma = p [0 1; -1 0]  =>  gamma m = p [m2; -m1]
    % c1 = 1/(2 mu), c2 = lambda/(4 mu (lambda + mu)) (1/(4 mu) incompressible)

    stress_basis = problem.LocalSpaces(1).Basis;      % PEERSBasis
    k = stress_basis.Order;
    dS = stress_basis.Dimension;
    d2 = nchoosek(k + 2, 2);
    N = size(rs, 2); n_b = size(rs, 3);
    n_u = 2*d2;                                 % displacement block (field 2)

    G_maps = get_affine_map(mesh, K_list(:), inverse = true, output = "array");

    Gm = zeros(N, n_local, n_b, 2);
    for b = 1:n_b
        P = PEERSBasis.evaluate_basis(k, rs(1, :, b), rs(2, :, b), G_maps(:, :, b));
        % P: (N x dS x 2 x 2); apply compliance and contract with m
        trP = P(:, :, 1, 1) + P(:, :, 2, 2);
        AP = c1*P - c2*trP.*reshape(eye(2), 1, 1, 2, 2);
        m_b = reshape(m(:, :, b), 2, N).';      % (N x 2)
        for a = 1:2
            Gm(:, 1:dS, b, a) = AP(:, :, a, 1).*m_b(:, 1) + AP(:, :, a, 2).*m_b(:, 2);
        end

        % antisymmetric multiplier (scalar Dubiner, field 3)
        p = DubinerBasis.evaluate_basis(k, rs(1, :, b).', rs(2, :, b).');
        Gm(:, dS + n_u + (1:d2), b, 1) = p .* m_b(:, 2);
        Gm(:, dS + n_u + (1:d2), b, 2) = -p .* m_b(:, 1);
    end
end

function Gm = tensor_field_gradient(problem, mesh, rs, K_list, m, n_local) %#ok<INUSD>
    % field 1 is a TensorFEBasis representing grad(u) directly:
    % [G_phi m]_a = [A_i m]_a phi_j for dof block i of linear map A_i
    grad_basis = problem.LocalSpaces(1).Basis;
    assert(isa(grad_basis, 'TensorFEBasis'), ...
        'velocity_gradient transfer expects a TensorFEBasis first field.');
    scalar_basis = grad_basis.ScalarBasis;
    maps = grad_basis.LinearMaps;
    d = scalar_basis.Dimension;
    N = size(rs, 2); n_b = size(rs, 3);
    n_comp = size(maps{1}, 1);

    Phi = double(scalar_basis.evaluate(reshape(rs, 2, [])));   % (N*n_b x d)
    Phi = permute(reshape(Phi, N, n_b, d), [1 3 2]);           % (N x d x n_b)

    Gm = zeros(N, n_local, n_b, n_comp);
    for i = 1:numel(maps)
        Am = pagemtimes(maps{i}, reshape(m, 2, 1, N, n_b));    % (n_comp x 1 x N x n_b)
        for a = 1:n_comp
            Gm(:, (i-1)*d + (1:d), :, a) = ...
                reshape(Am(a, 1, :, :), N, 1, n_b) .* Phi;
        end
    end
end

function Gm = polynomial_field_gradient(problem, mesh, rs, K_list, m, n_local)
    % field 1 is a scalar basis or its to_vector copies (component a in
    % dof block a): [G_phi m]_a = (B^-T grad_ref phi) . m on element K
    field_basis = problem.LocalSpaces(1).Basis;
    if isa(field_basis, 'TensorFEBasis')
        scalar_basis = field_basis.ScalarBasis;
        n_comp = numel(field_basis.LinearMaps);
    else
        scalar_basis = field_basis;
        n_comp = 1;
    end
    d = scalar_basis.Dimension;
    N = size(rs, 2); n_b = size(rs, 3);
    iB = get_affine_map(mesh, K_list(:), inverse = true, output = "array");

    Gm = zeros(N, n_local, n_b, n_comp);
    for b = 1:n_b
        [~, dPr, dPs] = DubinerBasis.evaluate_basis(scalar_basis.Order, rs(1, :, b), rs(2, :, b));
        dPx = iB(1, 1, b)*dPr + iB(2, 1, b)*dPs;
        dPy = iB(1, 2, b)*dPr + iB(2, 2, b)*dPs;
        Gphi_m = dPx .* reshape(m(1, :, b), N, 1) + dPy .* reshape(m(2, :, b), N, 1);
        for a = 1:n_comp
            Gm(:, (a-1)*d + (1:d), b, a) = Gphi_m;
        end
    end
end

%% ------------------------------ helpers -------------------------------- %%
function g = evaluate_datum(datum, x)
    % evaluate a scalar Datum at physical points x (2 x n_q x n_b) -> (n_q x n_b)
    v = datum.value;
    x1 = reshape(x(1, :, :), size(x, 2), []);
    x2 = reshape(x(2, :, :), size(x, 2), []);
    switch class(v)
        case 'function_handle'
            g = v(x1, x2);
        case 'double'
            g = v(1) * ones(size(x1));
        otherwise
            error('TransferEssentialBoundaryConditions:datum', ...
                  'Unsupported datum value of class %s.', class(v));
    end
end
