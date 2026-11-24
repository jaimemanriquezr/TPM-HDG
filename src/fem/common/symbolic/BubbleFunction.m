classdef BubbleFunction < WeightFunction
    % Element bubble weight (bubble_eval) on the reference coordinates;
    % zero on faces and for points that are not 2-D.
    methods
        function obj = BubbleFunction()
            obj.weight = @bubble_eval;
            obj.signature = "BubbleFunction";
        end
        function b = evaluate(~, x_local, mesh, varargin)
            arguments
                ~
                x_local double
                mesh SimplicialGrid
            end
            
            arguments (Repeating)
                varargin
            end

            if size(x_local, 1) == 2
                b = bubble_eval(x_local(1, :), x_local(2, :));
            else
                b = 0 * x_local;
            end
            b = permute(repmat(b(:), [1 size(mesh.elements, 2)]), [1 3 4 5 2]);
        end
        
        function b = evaluate_boundary(~, points, mesh, varargin)
            b = zeros(size(points, 2), 1, 1, 1, size(mesh.elements, 2), 3);
        end
    end
end