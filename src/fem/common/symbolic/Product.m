classdef Product < AbstractFEFunction
    % Product of FE expressions: evaluate folds the binary operation
    % product (@times by default, also @dot, @mtimes) over the factors.
    % Multiplication by a double scales the second factor, by a Measure
    % integrates, by a trial or test operand builds a form, and any other
    % operand multiplies the first factor.
    properties
        factors (1, :) AbstractFEFunction
        product
    end
    
    methods 
        function obj = Product(factors, product)
            arguments
                factors = FEFunction.empty;
                product = @times % @dot, @mtimes
            end
            obj.factors = factors;
            obj.product = product;
        end
        
        function p = evaluate(obj, varargin)
            p = evaluate(obj.factors(1), varargin{:});
            for i = 2:length(obj.factors)
                p = obj.product(p, evaluate(obj.factors(i), varargin{:}));
            end
        end

        function p = evaluate_boundary(obj, varargin)
            p = evaluate_boundary(obj.factors(1), varargin{:});
            for i = 2:length(obj.factors)
                p = obj.product(p, evaluate_boundary(obj.factors(i), varargin{:}));
            end
        end
        
        function obj = mtimes(obj, u)
            if isa(u, "BilinearForm")
                % a product on the left of a form contracts its value
                obj = contract(u, obj, "left");
                return
            end
            if ~isa(obj, "Product")
                u_ = obj;
                obj = u;
                u = u_;
            end
            if isa(u, "double")
                obj.factors(2) = u * obj.factors(2);
            elseif isa(u, "Measure")
                % absorb the measure: integrate this product as integrand
                obj = mtimes(u, obj);
            elseif isa(u, "AbstractFEFunction") ...
                    && contains(u.signature, ["Trial", "Test"])
                % trial/test operand: build the bilinear form exactly as
                % the superclass would (dot(q, n)*v_hat == v_hat*dot(q, n))
                obj = mtimes@AbstractFEFunction(obj, u);
            else
                obj.factors(1) = u * obj.factors(1);
            end
        end
    end
end