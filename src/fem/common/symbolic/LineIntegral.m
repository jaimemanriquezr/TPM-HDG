classdef LineIntegral < AbstractFEFunction
    % Transfer-path line integral (TPM, Cockburn & Solano).
    %
    % Given an argument u (an FE function or basis on the volume mesh) and a
    % mapping x -> x_bar from a point on a mesh face to a point on the true
    % boundary, evaluate(x_local) returns, for each face of the mesh,
    %
    %   I(x) = int_0^1 u(x + s*(x_bar - x)) * |x_bar - x| ds
    %        = int_{[x, x_bar]} u dsigma,
    %
    % the line integral of u along the (straight) transfer path from x to
    % x_bar, where x is the physical point on each face with local
    % coordinate x_local. The argument is extended from the first element
    % adjacent to each face (polynomial extrapolation, as in the TPM).
    properties
        argument
        mapping
    end

    methods
        function obj = LineIntegral(argument, mapping)
            arguments
                argument AbstractFEFunction
                mapping function_handle = @(X) X;
            end
            obj.argument = argument;
            obj.mapping = mapping;

            [s, ~, order] = read_signature(argument);
            obj.signature = "LineIntegral-" + s(end) + ":" + 1 + ":" + order;
        end

        function Iu = evaluate(obj, x_local, mesh, options)
            arguments
                obj LineIntegral
                x_local (1, :)
                mesh = obj.argument.mesh;
                options.mapping = obj.mapping;
                options.order = obj.argument.order;
                options.coordinate_system = "local";
                options.output = [];
            end

            % physical points on every face: 2 x n_points x n_faces
            x_physical_list = get_physical_points(mesh, x_local);
            [s, w] = get_1d_reference_quadrature(options.order);
            K_list = cellfun(@(C) C(1), get_element_by_face(mesh));

            n_points = size(x_physical_list, 2);
            n_faces = size(x_physical_list, 3);

            Iu = [];
            for i = 1:n_points
                x_start = reshape(x_physical_list(:, i, :), 2, 1, n_faces);
                x_end = options.mapping(x_start);
                ell = vecnorm(x_end - x_start);              % 1 x 1 x n_faces
                x_path = x_start + s .* (x_end - x_start);   % 2 x n_q x n_faces

                u = obj.argument.evaluate(x_path, mesh, ...
                                          element_list = K_list, ...
                                          coordinate_system = "physical");
                % u carries the quadrature points along dim 1 and the faces
                % along its last dimension; align ell with that last dimension
                ell = reshape(ell, [ones(1, max(ndims(u), 2) - 1), n_faces]);
                Iu_i = ell .* sum(w(:) .* u, 1);
                if isempty(Iu)
                    Iu = zeros([n_points, size(Iu_i, 2:max(ndims(Iu_i), 2))]);
                end
                Iu(i, :) = Iu_i(:);
            end
        end
    end
end
