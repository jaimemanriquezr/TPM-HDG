function obj_n = power(obj, n)
%POWER  Elementwise power obj.^n of the evaluated values.
    arguments
        obj AbstractFEFunction
        n double
    end
    obj_n = FEFunction(@(varargin) power(obj.evaluate(varargin{:}), n), ...
                        obj.signature);
end