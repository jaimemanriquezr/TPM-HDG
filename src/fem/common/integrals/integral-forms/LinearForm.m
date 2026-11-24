classdef LinearForm < AbstractFEFunction
    % A function of a test function, integrated into element-wise vector
    % blocks.
    %
    %   l = LinearForm(input, argument)
    %
    % input is a handle called as input(argument, ...) when the form is
    % evaluated; argument is the test function (or a sum of them). The
    % signature is "LinearForm:id:dim:order" from the argument.
    %
    %   I = integrate(l, mesh, measure_type=, order=, weight=, faces=)
    %
    % returns the IndexedArray of blocks indexed by [-1, test id];
    % measure_type is 'dx' (default), 'dS' or 'ds'. A sum of test
    % functions gives one block per summand, added together.
    %
    % Multiplying by a Measure gives an IntegralForm; multiplying by data
    % scales the form.
    %
    % See also: BilinearForm, Measure, IntegralForm
    properties (Hidden)
        evaluate_form
        argument
    end

    properties (Dependent, Hidden)
        id
    end

    methods
        function obj = LinearForm(input, argument)
            arguments
                input
                argument
            end
            obj.evaluate_form = input;
            obj.argument = argument;

            [~, dim, order, id] = read_signature(argument);
            form_signature = sprintf("LinearForm:%i:%i:%i", id, dim, order);
            obj.signature = form_signature;
        end

        function id = get.id(obj)
            [~, ~, ~, id] = read_signature(obj.argument);
        end

        function l = evaluate(obj, varargin)
            l = obj.evaluate_form(obj.argument, varargin{:});
        end

        function I = integrate(obj, mesh, options)
            arguments
                obj LinearForm
                mesh AbstractGrid = SimplicialGrid(reference_element(), ...
                    element_diameter=2);
                options.measure_type string = 'dx';
                options.order = [];
                options.weight = 1;
                options.faces = [];
            end
            op_args = namedargs2cell(options);
            if isa(obj.argument, "LinearCombination")
                % one block per summand of a sum of test functions: the
                % summands may live in different spaces
                lc = obj.argument;
                I = lc.coefficients(1) * integrate(LinearForm(obj.evaluate_form, ...
                    lc.summands(1)), mesh, op_args{:});
                for j = 2:numel(lc.summands)
                    I = I + lc.coefficients(j) * integrate(LinearForm(obj.evaluate_form, ...
                        lc.summands(j)), mesh, op_args{:});
                end
                return
            end
            [~, l_dim, ~, id_l] = read_signature(obj);
            [weights, measure, l_] = mesh.evaluate_quadrature(obj, op_args{:});

            l = reshape(l_, [size(l_, 1:2) 1 size(l_, 3:6)]);

            measure = reshape(measure, [1 1 1 1 1 size(measure)]);
            I = measure .* sum(weights .* l, 1);
            I = reshape(I, size(I, 2:7));

            if isequal(options.measure_type, 'dS')
                if ~isequal(l_dim, 1)
                    I = sum(I, 6);
                end
            end
            I = IndexedArray(I, [-1, id_l]);
        end

        function obj_times_x = mtimes(obj, x)
            if ~isa(obj, "LinearForm")
                [obj, x] = deal(x, obj);
            end
            switch class(x)
                case "Measure"
                    obj_times_x = mtimes(x, obj);
                otherwise
                    % a linear form has no trial factor: scale the form
                    f = obj.evaluate_form;
                    obj_times_x = obj;
                    obj_times_x.evaluate_form = @(v, varargin) x * f(v, varargin{:});
            end
        end
    end
end
