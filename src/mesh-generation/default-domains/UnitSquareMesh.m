classdef UnitSquareMesh < RectangleMesh
    % Mesh of the unit square [0, 1] x [0, 1] with element diameter h
    % (default 0.1).
    %
    %   mesh = UnitSquareMesh(h, bottom_side = "Dirichlet", right_side = ...,
    %       top_side = ..., left_side = ...)
    %
    % The side options are those of RectangleMesh.
    %
    % See also: RectangleMesh
    methods
        function obj = UnitSquareMesh(h, options)
            arguments
                h = .1
                options.bottom_side = "Dirichlet";
                options.right_side = "Dirichlet";
                options.top_side = "Dirichlet";
                options.left_side = "Dirichlet";
            end
            input_options = namedargs2cell(options);
            obj = obj@RectangleMesh(0, 1, 0, 1, h, input_options{:});
        end
    end
end
