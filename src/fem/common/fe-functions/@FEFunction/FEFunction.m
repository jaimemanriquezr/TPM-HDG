classdef FEFunction < AbstractFEFunction
    %FEFUNCTION  A function defined by an evaluation handle.
    %
    %   f = FEFunction(input, signature)
    %
    %   input is either a function_handle, called as input(...) by
    %   evaluate, or an AbstractFEBasis / AbstractFEFunction, whose evaluate
    %   method is used. Any other input raises FEFunction:InvalidInput.
    %   signature (default 'FEFunction') labels the function.
    %
    %   evaluate(f, ...) forwards its arguments to the stored handle.
    %   GeometricOutputDimension is the size (dimensions 3:4) of the value
    %   of f at the point [0; 0].
    %
    %See also: AbstractFEFunction, FEField
    properties (Hidden)
        eval_fun
    end

    properties (Dependent)
        GeometricOutputDimension
    end

    methods
        function obj = FEFunction(input, signature)
            arguments
                input
                signature string = 'FEFunction'
            end
            if isa(input, 'AbstractFEBasis') || isa(input, 'AbstractFEFunction')
                obj.eval_fun = @input.evaluate;
            elseif isa(input, 'function_handle')
                obj.eval_fun = input;
            else
                eid = 'FEFunction:InvalidInput';
                msg = 'Invalid input for FEFunction creation.';
                error(eid, msg);
            end
            obj.signature = signature;
        end

        plot(obj);

        function dim = get.GeometricOutputDimension(obj)
            dim = size(evaluate(obj, [.0; .0]), 3:4);
        end

        function v = evaluate(obj, varargin)
            v = obj.eval_fun(varargin{:});
        end
    end
end
