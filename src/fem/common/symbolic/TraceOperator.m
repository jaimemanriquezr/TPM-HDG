classdef TraceOperator < FEOperator
    % FEOperator: trace of a square tensor value.
    methods
        function obj = TraceOperator()
            obj = obj@FEOperator(@TraceOperator.apply_action);
            obj.output_dimension = [1 1];
        end
    end
    
    methods (Static)
        function tr_p = apply_action(p)
            arguments
                p
            end
            assert(size(p, 3) == size(p, 4));
            n_d = size(p, 3);
            tr_p = zeros([size(p, 1:2) 1 1 size(p, 5)]);
            for i = 1:n_d
                tr_p = tr_p + p(:, :, i, i, :);
            end
        end
    end
end