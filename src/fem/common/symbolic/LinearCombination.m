classdef LinearCombination < AbstractFEFunction
    % Sum of FE expressions with coefficients, sum_i c_i*f_i. Products by
    % scalars, transposition, plus and minus act on the coefficients and
    % summands, so the sum stays visible to the assembly.
    properties
        summands (1, :) AbstractFEFunction
        coefficients (1, :)
    end
    
    methods 
        function obj = LinearCombination(summands, coefficients)
            arguments
                summands = FEFunction.empty;
                coefficients = ones(1, length(summands));
            end
            assert(length(summands) == length(coefficients))
            obj.summands = summands;
            obj.coefficients = coefficients;

            if length(summands) > 1
                signatures = cell(length(summands), 1);
                for i = 1:length(summands)
                    signatures{i} = read_signature(summands(i));
                end
                if isequal(signatures{:})
                    obj.signature = summands(1).signature;
                end
            else
                obj.signature = summands.signature;
            end
        end
        
        function P = evaluate(obj, varargin)
            p = obj.coefficients(1) * obj.summands(1);
            P = evaluate(p, varargin{:});
            for i = 2:length(obj.summands)
                p = obj.coefficients(i) * obj.summands(i);
                P = P + evaluate(p, varargin{:});
            end
        end

        function P = evaluate_boundary(obj, varargin)
            p = obj.coefficients(1) * obj.summands(1);
            P = evaluate_boundary(p, varargin{:});
            for i = 2:length(obj.summands)
                p = obj.coefficients(i) * obj.summands(i);
                P = P + evaluate_boundary(p, varargin{:});
            end
        end
        
        function obj_times_x = mtimes(obj, x)
            if ~isa(obj, 'LinearCombination')
                [obj, x] = deal(x, obj);
            end
            switch class(x)
                case "double"
                    obj_times_x = obj;
                    obj_times_x.coefficients = x * obj.coefficients;
                otherwise
                    if isa(x, "BilinearForm")
                        obj_times_x = contract(x, obj, "left");      % L*(a*b)
                    elseif isa(x, "AbstractFEFunction") && ~isa(x, "WeightFunction")
                        obj_times_x = mtimes@AbstractFEFunction(obj, x);   % written order
                    else
                        obj_times_x = mtimes(x, obj);   % Measure; weights distribute
                    end
            end
        end

        function obj_t = transpose(obj)
            % transposition is linear: keep the sum visible (block split,
            % signature)
            obj_t = obj;
            for i = 1:numel(obj.summands)
                obj_t.summands(i) = transpose(obj.summands(i));
            end
        end

        function obj_t = ctranspose(obj)
            obj_t = obj;
            for i = 1:numel(obj.summands)
                obj_t.summands(i) = ctranspose(obj.summands(i));
            end
        end
        
        function obj = plus(obj, x)
            if isa(x, 'LinearCombination')      % flatten
                obj.summands = [obj.summands x.summands];
                obj.coefficients = [obj.coefficients x.coefficients];
            else
                obj.summands = [obj.summands x];
                obj.coefficients = [obj.coefficients 1];
            end
        end

        function obj = minus(obj, x)
            if isa(x, 'LinearCombination')      % flatten
                obj.summands = [obj.summands x.summands];
                obj.coefficients = [obj.coefficients -x.coefficients];
            else
                obj.summands = [obj.summands x];
                obj.coefficients = [obj.coefficients -1];
            end
        end

        function obj = uplus(obj)
        end

        function m_obj = uminus(obj)
            m_obj = obj;
            m_obj.coefficients = -m_obj.coefficients;
        end
    end
end