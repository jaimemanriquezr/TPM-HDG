classdef PEERSBasis < AbstractFEBasis
    % Discontinuous 2-by-2 matrix-valued basis: the Dubiner polynomials of
    % degree up to k in each of the four components, enriched with the
    % bubble modes of curl_curl_A_dubiner_bubble_eval (none for k = 0).
    % Piola transformed; Dimension = 4*nchoosek(k+2, 2) + (k+1) for k > 0
    % and 4 for k = 0.
    %
    %   basis = PEERSBasis(k)
    %
    % k is the polynomial order (default empty, giving Dimension 0).
    %
    %   P = evaluate(basis, points, mesh, coordinate_system=, element_list=, output=)
    %
    % points is 2-by-N ("local" reference points by default, or "physical"
    % points, which need mesh). Without element_list, P is
    % N-by-Dimension-by-2-by-2-by-n_el, using the inverse affine map of
    % every element. With element_list, points is 2-by-N-by-n_pages, page b
    % is evaluated on element element_list(b), and only output="array" is
    % supported.
    %
    %   E = compute_l2_error(basis, dof, mesh, exact_solution=)
    %
    % returns the 2-by-2 matrix of componentwise L2 errors; exact_solution
    % is a 2-by-2 cell of @(x, y) handles (default zeros).
    %
    % See also: AbstractFEBasis, DubinerBasis, curl_curl_A_dubiner_bubble_eval
    properties
        Order double
    end

    properties (Constant, Hidden)
        ElementType = "Simplex"
        GeometricInputDimension = 2
        GeometricOutputDimension = [2, 2]
        IsDiscontinuous = true
        IsPiolaTransformed = true
    end

    properties (Dependent)
        Dimension
    end

    methods
        function obj = PEERSBasis(k)
            arguments
                k double = [];
            end
            obj.Order = k;
        end

        function dim_space = get.Dimension(obj)
            if isempty(obj.Order)
                dim_space = 0;
            else
                d2 = nchoosek(obj.Order + obj.GeometricInputDimension, obj.GeometricInputDimension);
                if obj.Order == 0
                    d2star = 0;
                else
                    d2star = nchoosek(obj.Order + 1,1);
                end
                dim_space = 4*d2 + d2star;
            end
        end

        function varargout = evaluate(obj, points, mesh, options)
            % Two modes, mirroring get_2d_physical:
            %   no element_list: points (2 x N), evaluated with the inverse
            %       affine map of every element -> (N x dS x 2 x 2 x n_el);
            %   element_list:    points (2 x N x n_pages), page b evaluated
            %       with element_list(b)'s map -> (N x dS x 2 x 2 x n_pages).
            arguments
                obj         PEERSBasis
                points      (:, :, :) double
                mesh        SimplicialGrid = SimplicialGrid.empty
                options.coordinate_system string = "local"
                options.element_list (:, 1) double = []
                options.output string = "array"
            end
            assert(size(points, 1) == obj.GeometricInputDimension)

            switch options.coordinate_system
                case "physical"
                    assert(~isempty(mesh), "No mesh specified!");
                    switch obj.GeometricInputDimension
                        case 2
                            local_points = mesh.get_2d_local(points, options.element_list);
                    end
                otherwise % "local" or "vertex"
                    local_points = points;
            end

            if isempty(options.element_list)
                assert(size(local_points, 3) == 1)
                args = mat2cell(local_points, ones(obj.GeometricInputDimension, 1));
                args{end+1} = get_affine_map(mesh, inverse=true);

                varargout = cell(1, nargout(@PEERSBasis.evaluate_basis));
                [varargout{:}] = obj.evaluate_basis(obj.Order, args{:}, output=options.output);
            else
                % per-page mode: page b of the points belongs to element_list(b)
                assert(size(local_points, 3) == numel(options.element_list), ...
                    "PEERSBasis:pageMismatch", ...
                    "points must carry one page per element_list entry.");
                assert(options.output == "array", ...
                    "PEERSBasis:cellOutput", ...
                    "element_list mode supports output=""array"" only.");
                G = get_affine_map(mesh, options.element_list, inverse = true, output = "array");
                n_pages = numel(options.element_list);
                P = zeros(size(local_points, 2), obj.Dimension, 2, 2, n_pages);
                for b = 1:n_pages
                    P(:, :, :, :, b) = obj.evaluate_basis(obj.Order, ...
                        local_points(1, :, b), local_points(2, :, b), G(:, :, b));
                end
                varargout = {P};
            end
        end

        function E = compute_l2_error(obj, dof, mesh, options)
            arguments
                obj PEERSBasis
                dof double
                mesh SimplicialGrid
                options.exact_solution cell = repmat({@(x,y) 0*(x.*y)}, [2 2]);
            end
            % componentwise errors (2 x 2); their Frobenius norm is
            % FEField.compute_l2_error
            exact_solution = options.exact_solution;

            [xy, rs, w] = quadpoints_elements_2d(mesh, obj.Order + 1, output="array");

            % (N x dS x 2 x 2 x n_el) basis values times the element dofs
            B = obj.evaluate(rs, mesh);
            P = pagemtimes(B, reshape(dof, obj.Dimension, 1, 1, 1, []));

            E = zeros(size(exact_solution));
            for i = 1:numel(exact_solution)
                [r, c] = ind2sub(size(exact_solution), i);
                P_i = reshape(P(:, 1, r, c, :), size(P, 1), []);
                P_ex = reshape(exact_solution{i}(xy(1,:,:), xy(2,:,:)), size(P_i));
                E(i) = sqrt(sum(mesh.volume .* (w'*(P_i - P_ex).^2)));
            end
        end

    end

    methods (Static)
        function P = evaluate_basis(k, r, s, G, options)
            arguments
                k double
                r double
                s double = r(2, :);
                G (2, 2, :) double = eye(2);
                options.output string = "array";
            end
            if nargin < 3; r = r(1, :); end

            p = DubinerBasis.evaluate_basis(k, r, s);
            %P = [P.*reshape([1, 0; 0, 0], [1 1 2 2]), ...
            %     P.*reshape([0, 0; 1, 0], [1 1 2 2]), ...
            %     P.*reshape([0, 1; 0, 0], [1 1 2 2]), ...
            %     P.*reshape([0, 0; 0, 1], [1 1 2 2])];
            E = mat2cell(reshape(eye(4), [4 1 2 2]), ones(4, 1));
            polynomial_modes = cell2mat(cellfun(@(e) p.*e, E, UniformOutput=false).');
            bubble_modes = curl_curl_A_dubiner_bubble_eval(k, r, s, inverse_affine_map = G);

            P = cat(2, repmat(polynomial_modes, [ones(1, 4) size(G, 3)]), bubble_modes);
            switch options.output
                case "array"
                    return
                case "cell"
                    P = mat2cell(P, size(P, 1), size(P, 2), [1 1], [1 1], size(P, 5));
                    P = squeeze(cellfun(@squeeze, P, 'UniformOutput', 0));
            end
        end

    end
end
