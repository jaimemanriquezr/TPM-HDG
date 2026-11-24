function dS = boundary_measure(mesh, weight)
%BOUNDARY_MEASURE  The 'dS' Measure (boundary of each element) on a mesh.
%
%   dS = boundary_measure(mesh, weight)
%
%   Returns Measure(mesh, 'dS', weight); weight defaults to 1. If weight is
%   empty, returns instead a function handle @(weight) that builds the
%   measure on mesh from a weight given later. mesh defaults to the
%   reference element.
    arguments
        mesh = SimplicialGrid(reference_element(), element_diameter=2);
        weight double = 1.0;
    end
    if isempty(weight)
        dS = @(weight) get_boundary_measure(mesh, weight);
    else
        dS = Measure(mesh, 'dS', weight);
    end
end

function dS = get_boundary_measure(mesh, weight)
    arguments
        mesh = SimplicialGrid(reference_element(), element_diameter=2);
        weight double = 1.0
    end
    dS = Measure(mesh, 'dS', weight);
end