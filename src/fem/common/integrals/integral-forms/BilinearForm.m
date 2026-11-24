classdef BilinearForm < LinearForm
    % Product of one trial and one test operand, integrated into
    % element-wise matrix blocks.
    %
    %   a = BilinearForm(u, v, op)
    %
    % u and v are AbstractFEFunctions, one a trial and the other a test
    % function (otherwise BilinearForm:InputType is raised), in either
    % order. op combines their values: @mtimes (default), @dot or @times.
    % The properties trial and test return the respective operand and
    % factors the operands in the written order.
    %
    %   I = integrate(a, mesh, measure_type=, order=, weight=, faces=)
    %
    % returns the IndexedArray of blocks, indexed by [trial id, test id];
    % measure_type is 'dx' (default), 'dS' or 'ds'. With op = @mtimes the
    % inner value sizes must agree (BilinearForm:mtimes otherwise).
    % Sums of trial or test functions give one block per pair of
    % summands, added together.
    %
    % Multiplying by a Measure gives an IntegralForm; multiplying by data
    % contracts it into the form; transpose and ctranspose reverse the
    % product.
    %
    % See also: LinearForm, Measure, IntegralForm
    properties
        % the two operands in the written order; the roles (trial, test)
        % come from their signatures
        factors (1, :) AbstractFEFunction
        op function_handle
    end

    properties (Dependent)
        trial
        test
    end

    methods
        function obj = BilinearForm(a, b, op)
            arguments
                a AbstractFEFunction
                b AbstractFEFunction
                op function_handle = @mtimes
            end
            [trial, test] = AbstractFEFunction.is_trial_or_test(a, b);
            assert(isscalar(trial) && isscalar(test), "BilinearForm:InputType", ...
                "a bilinear form takes one trial and one test operand (%s, %s)", ...
                a.signature, b.signature);
            % the integrand as written; the test function is the form's argument
            evaluate_form = @(~, varargin) AbstractFEFunction.pagefun(op, ...
                a.evaluate(varargin{:}), b.evaluate(varargin{:}));
            obj = obj@LinearForm(evaluate_form, test);
            obj.factors = [a, b];
            obj.op = op;

            [~, u_dim, u_order, u_id] = read_signature(trial);
            [~, v_dim, v_order, v_id] = read_signature(test);
            b_dim = sprintf("%i-%i", u_dim, v_dim);
            b_order = sprintf("%i-%i", u_order, v_order);
            b_id = sprintf("%i-%i", u_id, v_id);
            obj.signature = "BilinearForm" + sprintf(":%s", [b_dim, b_order, b_id]);
        end

        function u = get.trial(obj)
            u = obj.factors(~is_test(obj.factors));
        end

        function v = get.test(obj)
            v = obj.factors(is_test(obj.factors));
        end

        function I = integrate(obj, mesh, options)
            arguments
                obj BilinearForm
                mesh AbstractGrid = SimplicialGrid(reference_element(), ...
                                                   element_diameter=2);
                options.measure_type string = 'dx';
                options.order = [];
                options.weight = 1;
                options.faces = [];
            end

            op_args = namedargs2cell(options);
            % one block per pair of summands: the summands of a sum of
            % trial (or test) functions may live in different spaces
            f = obj.factors;
            [a_terms, a_coefs] = summands(f(1));
            [b_terms, b_coefs] = summands(f(2));
            I = [];
            for i = 1:numel(a_terms)
                for j = 1:numel(b_terms)
                    I_ij = (a_coefs(i) * b_coefs(j)) * integrate_form( ...
                        a_terms{i}, b_terms{j}, obj.op, mesh, op_args{:});
                    if isempty(I)
                        I = I_ij;
                    else
                        I = I + I_ij;
                    end
                end
            end
        end

        function obj_times_x = mtimes(obj, x)
            if ~isa(obj, "LinearForm")
                obj_times_x = contract(x, obj, "left");
            elseif isa(x, "Measure")
                obj_times_x = mtimes(x, obj);
            else
                obj_times_x = contract(obj, x, "right");
            end
        end

        function obj_t = transpose(obj)
            % (a*b).' = b.'*a.'
            obj_t = transpose_form(obj, @transpose);
        end

        function obj_t = ctranspose(obj)
            % (a*b)' = b'*a'
            obj_t = transpose_form(obj, @ctranspose);
        end

        function obj = contract(obj, x, side)
            % data x times the form on the given side. A scalar scales the
            % trial factor; for a written product (op mtimes) anything else
            % contracts the value dims by associativity:
            % L*(a*b) = (L*a)*b and (a*b)*R = a*(b*R)
            arguments
                obj BilinearForm
                x
                side (1,1) string {mustBeMember(side, ["left", "right"])}
            end
            if isa(x, "double") && isscalar(x) || ~isequal(func2str(obj.op), "mtimes")
                i_trial = find(~is_test(obj.factors));
                obj.factors(i_trial) = x * obj.factors(i_trial);
                return
            end
            if side == "left"
                i = 1;
                pair = {x, obj.factors(1)};
            else
                i = 2;
                pair = {obj.factors(2), x};
            end
            is_double = cellfun(@(c) isa(c, "double"), pair);
            pair(is_double) = cellfun(@ConstantWeight, pair(is_double), ...
                "UniformOutput", false);
            p = Product([pair{:}], @AbstractFEFunction.mtimes_);
            p.signature = obj.factors(i).signature;
            obj.factors(i) = p;
        end
    end
end

function obj_t = transpose_form(obj, op_t)
    switch func2str(obj.op)
        case "mtimes"
            obj_t = BilinearForm(op_t(obj.factors(2)), op_t(obj.factors(1)), @mtimes);
        case "dot"
            obj_t = obj;    % scalar-valued
        case "times"
            obj_t = BilinearForm(op_t(obj.factors(1)), op_t(obj.factors(2)), @times);
    end
end

function tf = is_test(f)
    tf = contains([f.signature], "Test");
end

function [s, c] = summands(f)
    % the summands (a cell) and coefficients of a sum, or {f} and 1
    if isa(f, "LinearCombination")
        s = num2cell(f.summands);
        c = f.coefficients;
    else
        s = {f};
        c = 1;
    end
end

function I = integrate_form(a, b, op, mesh, options)
    arguments
        a AbstractFEFunction
        b AbstractFEFunction
        op function_handle
        mesh AbstractGrid = SimplicialGrid(reference_element(), ...
                                           element_diameter=2);
        options.measure_type string = 'dx';
        options.order = [];
        options.weight = 1;
        options.faces = [];
    end
    op_args = namedargs2cell(options);

    a_is_test = contains(a.signature, "Test");
    if a_is_test
        [u, v] = deal(b, a);
    else
        [u, v] = deal(a, b);
    end
    [~, u_dim, ~, id_u] = read_signature(u);
    [~, v_dim, ~, id_v] = read_signature(v);
    [weights, measure, a_, b_] = mesh.evaluate_quadrature(a, b, op_args{:});

    a = as_pages(a_, a_is_test);
    b = as_pages(b_, ~a_is_test);
    switch char(op)
        case 'mtimes'
            % written order; pagemtimes expands 1x1 pages (a scalar scales)
            if prod(size(a, 1:2)) > 1 && prod(size(b, 1:2)) > 1
                assert(size(a, 2) == size(b, 1), "BilinearForm:mtimes", ...
                    "Inner value dimensions must agree: [%s]*[%s].", ...
                    num2str(size(a, 1:2)), num2str(size(b, 1:2)));
            end
            uv = pagemtimes(a, b);
        case 'dot'
            uv = sum(a .* b, 1:2);
        case 'times'
            uv = a .* b;
    end
    uv = permute(uv, [3 4 5 1 2 6 7]);

    measure = reshape(measure, [1 1 1 1 1 size(measure)]);
    I = measure .* sum(weights .* uv, 1);
    I = reshape(I, size(I, 2:7));

    if isequal(options.measure_type, 'dS')
        if ~any([u_dim, v_dim] == 1)
            I = sum(I, 6);
        elseif u_dim == 1 && v_dim ~= 1
            I = permute(I, [1 2 3 4 5 7 6]);
        elseif u_dim ~= 1 && v_dim == 1
            % I
        else
            n_sides = size(mesh.boundaries, 1);
            I_ = zeros([size(I, 1:5) n_sides n_sides]);
            for i = 1:n_sides
                I_(:, :, :, :, :, i, i) = I(:, :, :, :, :, i);
            end
            I = I_;
        end
    end
    I = IndexedArray(I, [id_u, id_v]);
end

function x = as_pages(x_, is_test)
    % quadrature layout [nq nb p q el (face)] -> [p q nq nb_test nb_trial
    % el (face)]: test dofs on page dim 4, trial dofs on page dim 5, so
    % the two operands broadcast against each other
    if is_test
        x = reshape(x_, [size(x_, 1:2) 1 size(x_, 3:6)]);
    else
        x = reshape(x_, [size(x_, 1) 1 size(x_, 2:6)]);
    end
    x = permute(x, [4 5 1 2 3 6 7]);
end
