classdef FlatCurveVertical < FlatCurve
    % The vertical line x = x_coord, parametrized by y over t_range, with
    % the unit normal pointing to -x (the tangent rotated by +90 degrees).
    methods
        function obj = FlatCurveVertical(x_coord, t_range)
            arguments
                x_coord double = 0.5
                t_range double = [0 1]
            end
            obj = obj@FlatCurve([x_coord; 0], [0; 1], t_range);
        end

        function y0 = double(obj)
            y0 = obj.x_coord(0);
        end

        function c = minus(obj, x)
            if ~isa(obj, "FlatCurveVertical")
                [obj, x] = deal(x, obj);
                c = minus(x, double(obj));
            else
                c = minus(double(obj), x);
            end
        end

        function c = plus(obj, x)
            if ~isa(obj, "FlatCurveVertical")
                [obj, x] = deal(x, obj);
            end
            c = obj.x_coord(0) + x;
        end
    end
end
