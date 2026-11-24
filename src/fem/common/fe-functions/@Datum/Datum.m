classdef Datum < AbstractFEFunction
    %DATUM  A given function used as data or exact solution in forms and errors.
    %
    %   d = Datum(value, name) wraps value, which is a function handle of
    %   (x, y), a cell array of such handles (one per component), a numeric
    %   constant (default 1.0), or a string that is parsed symbolically
    %   with str2sym in the variables x and y. name is optional; for a
    %   string value it defaults to the string itself and is used for plot
    %   titles.
    %
    %   Uses:
    %     Datum(f) * w * Measure(...)     data in a form
    %     compute_l2_error(u, mesh, d)    exact solution of an error
    %     plot(d, mesh)                   plot the datum over a mesh
    %
    %See also: FEField, Measure
    properties
        name string = string.empty
        value
    end
    
    properties (Dependent)
        GeometricInputDimension
    end

    methods
        function obj = Datum(value, name)
            arguments
                value = 1.0
                name = string.empty;
            end

            obj.name = name;
            if isa(value, "string")
                if isempty(name)
                    obj.name = value;
                end
                value = str2sym(value);
            end
            obj.value = value;
            
            switch class(obj.value)
                case 'function_handle'
                    dim = nargin(obj.value);
                    order = -1;
                case 'cell'
                    dim = nargin(obj.value{1});
                    order = -1;
                case 'double'
                    dim = 2;
                    order = 0;
                case 'sym'
                    if isscalar(value)
                        obj.value = sym2MATLAB(value);
                    else
                        obj.value = arrayfun(@sym2MATLAB, value, "UniformOutput", false);
                    end
                    dim = length(symvar(value));
                    order = -1;
            end
            obj.signature = "Datum:" + dim + ":" + order;
        end
        
        function dim = get.GeometricInputDimension(obj)
            [~, dim] = read_signature(obj);
        end
        
        function varargout = value_size(obj, varargin)
            varargout = cell(max(nargout, 1), 1);
            if isempty(varargin)
                [varargout{:}] = size(obj.value);
            else
                [varargout{:}] = size(obj.value, varargin{:});
            end
        end
        
        function f_val = get_value(obj, points)
            f_val = cell(size(obj.value));
            points_cell = mat2cell(points, ones(size(points, 1), 1));
            switch class(obj.value)
                case 'function_handle'
                    f_val{1} = permute(obj.value(points_cell{:}), [2 1 4 3]);
                case 'cell'
                    for i = 1:numel(obj.value)
                        f_val{i} = permute(obj.value{i}(points_cell{:}), [2 1 4 3]);
                    end
                case 'double'
                    for i = 1:numel(obj.value)
                        f_val{i} = permute(obj.value(i) * (points(1, :, :).^0), [2 1 4 3]);
                    end
            end
        end
        
        I = test(obj, mesh, basis, options)
        f = evaluate(obj, points, mesh, options)
        plot(obj, mesh, options)
    end
end

function F = sym2MATLAB(f)
    syms x y
    f = matlabFunction(f, 'Vars', {x, y});
    F = @(r, s) (r.*s).^0 .* feval(f, r, s);
end