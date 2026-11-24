classdef DeviatorOperator < FEOperator
    % FEOperator: deviatoric part p - tr(p)/n*Id of a square tensor value.
    methods
        function obj = DeviatorOperator()
            obj = obj@FEOperator(@DeviatorOperator.apply_action);
            obj.output_dimension = [1 1];
        end
    end
    
    methods (Static)
        function dev_p = apply_action(p)
            arguments
                p
            end
            assert(size(p, 3) == size(p, 4));
            n_d = size(p, 3);
            tr_p = zeros([size(p, 1:2) 1 1 size(p, 5)]);
            for i = 1:n_d
                tr_p = tr_p + p(:, :, i, i, :);
            end
            
            Id = permute(eye(n_d), [3 4 1 2]);
            dev_p = p - tr_p.*(Id / n_d);
        end
    end
end