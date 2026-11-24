function dx = volume_measure(mesh, weight)
%VOLUME_MEASURE  The 'dx' Measure (elements) on a mesh.
%
%   dx = volume_measure(mesh, weight)
%
%   Returns Measure(mesh, 'dx', weight); weight defaults to 1. If weight is
%   empty, returns instead a function handle that builds the measure on
%   mesh from a weight given later. mesh defaults to an empty
%   SimplicialGrid.
    arguments
        mesh = SimplicialGrid.empty;
        weight double = 1.0;
    end
    if isempty(weight)
        dx = @(varargin) get_volume_measure(mesh, varargin{:});
    else
        dx = Measure(mesh, 'dx', weight);
    end
end

function dx = get_volume_measure(mesh, weight)
    arguments
        mesh = SimplicialGrid(reference_element(), element_diameter=2);
        weight double = 1.0
    end
    dx = Measure(mesh, 'dx', weight);
end