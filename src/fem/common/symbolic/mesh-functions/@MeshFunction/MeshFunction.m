classdef MeshFunction < WeightFunction
    % WeightFunction that depends on the mesh: the weight is evaluated as
    % weight(mesh) and broadcast over the points and the elements (or
    % faces), then multiplied by element_weights; evaluate_boundary also
    % multiplies by boundary_weights (one entry per local face).
    properties
        element_weights
        boundary_weights
    end
    
    methods
        function obj = MeshFunction(input, options)
            arguments
                input = @(varargin) 1;
                options.element_weights = 1;
                options.boundary_weights = 1;
            end 
            for field = string(fieldnames(options).')
                obj.(field) = options.(field);
            end
            
            if isequal(size(options.boundary_weights), [3 1])
                obj.boundary_weights = reshape(options.boundary_weights, [1 1 1 1 1 3]);
            end


            switch class(input)
                case "function_handle"
                    obj.weight = input;
                    obj.signature = "MeshFunction:-1";
                case "double"
                    obj.signature = "MeshFunction:0";
                    obj.weight = @(varargin) 1;
                    if size(input, 4) > 1
                        obj.element_weights = 1;
                        obj.boundary_weights = input;
                        obj.signature = obj.signature + ":1";
                        obj.boundary_weights = reshape(obj.boundary_weights, [1 1 size(obj.boundary_weights)]);
                    else
                        obj.element_weights = input;
                        obj.boundary_weights = input(:).' .* obj.boundary_weights;
                        obj.signature = obj.signature + ":2";
                    end
            end
        end
        
        function f_K = evaluate_weight(obj, points, mesh, varargin)
            arguments
                obj         MeshFunction
                points      double 
                mesh        AbstractGrid = SimplicialGrid.empty;
            end
            arguments (Repeating)
                varargin
            end
            one_p = ones(size(points, 2), 1);
            switch size(points, 1)
                case 1
                    one_K = ones(1, 1, 1, 1, size(mesh.faces, 2));
                case 2
                    one_K = ones(1, 1, 1, 1, size(mesh.elements, 2));
            end
            f_K = feval(obj.weight, mesh) .* one_p .* one_K;
        end
        
        function mf_K = evaluate(obj, points, mesh, varargin)
            arguments
                obj         MeshFunction
                points      double 
                mesh        AbstractGrid = SimplicialGrid.empty;
            end
            arguments (Repeating)
                varargin
            end
            f_K = obj.evaluate_weight(points, mesh, varargin{:});
            mf_K = f_K .* obj.element_weights;
        end
        
        function obj = transpose(obj)
            % swap the value dims (3-4) of the carried weight, so e.g.
            % n.' is a row vector and u_hat*n.' spells the outer product
            % (u_hat x n), the trial-side form of <u_hat, sigma n>
            if isa(obj.weight, 'function_handle')
                w0 = obj.weight;
                obj.weight = @(mesh) permute(w0(mesh), [1 2 4 3 5 6]);
            end
            obj.element_weights  = permute(obj.element_weights,  [1 2 4 3 5 6]);
            obj.boundary_weights = permute(obj.boundary_weights, [1 2 4 3 5 6]);
        end

        function obj = ctranspose(obj)
            obj = transpose(obj);
        end

        function mf_dK = evaluate_boundary(obj, points, mesh, varargin)
            arguments
                obj         MeshFunction
                points      double 
                mesh        AbstractGrid = SimplicialGrid.empty;
            end
            arguments (Repeating)
                varargin
            end
            if size(points, 1) == 1
                points = [points; points];
            end
            mf_dK = obj.evaluate(points, mesh, varargin{:}) .* obj.boundary_weights;
        end
    end

end