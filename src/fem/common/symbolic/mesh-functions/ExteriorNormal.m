classdef ExteriorNormal < MeshFunction
    % MeshFunction: exterior normal of the faces of each element. With a
    % mesh the normals are built at construction from mesh.normals and
    % mesh.orientation, divided by the face area when normalize is true
    % (the default); without one they are taken from
    % mesh.exterior_unit_normals at evaluation time.
    methods
        function obj = ExteriorNormal(mesh, normalize)
            arguments
                mesh = SimplicialGrid.empty;
                normalize = true
            end
            if isempty(mesh)
                normals = @(mesh) permute(mesh.exterior_unit_normals, [5 6 1 2 4 3]);
            else
                if normalize
                    scaled_orientation = mesh.orientation ./ (mesh.area(mesh.boundaries));
                else
                    scaled_orientation = mesh.orientation;
                end
                normals = reshape(mesh.normals .* kron(scaled_orientation, [1; 1]), 2, 1, size(mesh.elements, 1), []);
                normals = permute(normals, [1 2 4 3]);
            end
            obj = obj@MeshFunction(normals);
        end
    end
end
