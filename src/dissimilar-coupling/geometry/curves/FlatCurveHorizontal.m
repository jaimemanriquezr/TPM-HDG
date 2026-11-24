classdef FlatCurveHorizontal < FlatCurve
    % The horizontal line y = y_coord, parametrized by x over t_range, with
    % the upward unit normal.
    methods
        function obj = FlatCurveHorizontal(y_coord, t_range)
            arguments
                y_coord double = 0.5
                t_range double = [0 1]
            end
            obj = obj@FlatCurve([0; y_coord], [1; 0], t_range);
        end

        function y0 = double(obj)
            y0 = obj.y_coord(0);
        end

        function c = minus(obj, x)
            if ~isa(obj, "FlatCurveHorizontal")
                [obj, x] = deal(x, obj);
                c = minus(x, double(obj));
            else
                c = minus(double(obj), x);
            end
        end

        function c = plus(obj, x)
            if ~isa(obj, "FlatCurveHorizontal")
                [obj, x] = deal(x, obj);
            end
            c = obj.y_coord(0) + x;
        end
    end
end
