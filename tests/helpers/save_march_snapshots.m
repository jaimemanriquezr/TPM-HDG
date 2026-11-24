function save_march_snapshots(file, k, t, fluid_solution, solid_solution, state, sigma_points, diagnostics)
%SAVE_MARCH_SNAPSHOTS  Append one time level of a coupled march to a .mat file.
%
%   save_march_snapshots(file, k, t, fs, ss, state, sigma_points, diagnostics)
%
%   Appends to file (created on the first call) one struct variable per
%   call, snapshot_000001, snapshot_000002, ..., with plain arrays only,
%   so the file loads without the toolbox: time t, degree k; fluid mesh
%   (coordinates 2 x Nv, elements 3 x N) and its Dubiner dofs (stress,
%   velocity, pressure; one column per element, component blocks stacked);
%   solid mesh and its dofs (stress, displacement) plus the Newmark
%   velocity; sigma_points, the 2 x m samples of the interface curve used
%   at t. diagnostics (optional), a struct of plain arrays, adds its
%   fields to the entry. A call writes its own variable only, so its cost
%   does not grow with the levels already saved, and a killed run keeps
%   them.
%
%See also: lagrangian_carved_march, eulerian_carved_march
entry = struct("t", t, "k", k, ...
    "fluid_coordinates", fluid_solution.Mesh.coordinates, ...
    "fluid_elements", fluid_solution.Mesh.elements, ...
    "fluid_stress", fluid_solution.Unknowns.field(1).dof, ...
    "fluid_velocity", fluid_solution.Unknowns.field(2).dof, ...
    "fluid_pressure", field_or_empty(state.fluid, "pressure"), ...
    "solid_coordinates", solid_solution.Mesh.coordinates, ...
    "solid_elements", solid_solution.Mesh.elements, ...
    "solid_stress", solid_solution.Unknowns.field(1).dof, ...
    "solid_displacement", solid_solution.Unknowns.field(2).dof, ...
    "solid_velocity", field_or_empty(state.solid, "velocity"), ...
    "sigma_points", sigma_points);
if nargin > 7
    for name = string(fieldnames(diagnostics)).'
        entry.(name) = diagnostics.(name);
    end
end
if isfile(file)
    level = numel(who("-file", file, "snapshot_*")) + 1;
    snapshot = struct(sprintf("snapshot_%06d", level), entry);
    save(file, "-struct", "snapshot", "-append");
else
    snapshot = struct("snapshot_000001", entry);
    save(file, "-struct", "snapshot", "-v7.3");
end
end

function v = field_or_empty(s, name)
v = [];
if isstruct(s) && isfield(s, name), v = s.(name); end
end
