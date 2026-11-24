function bc = boundary_condition(tag, options)
% Component-wise boundary-condition specification for one boundary tag.
%
% Each Cartesian component of the trace unknown is either constrained
% (dirichlet: u_hat_c = data) or loaded (neumann: (sigma n)_c = data, the
% c-component of the boundary traction with outward normal n). Mixing the
% two on one boundary expresses slip/symmetry ("u.n = 0, tang(sigma n) = 0")
% and traction-driven inflow ("(sigma n).n = -p_in, tang(u) = 0") conditions
% on axis-aligned boundaries, where components and normal/tangential
% directions coincide.
%
% Example (Meddahi 2025, Eq. 52, inlet at x = 0 with n = (-1, 0)):
%   boundary_condition("Inflow", x = "neumann", x_data = @(x,y) p_in + 0*x)
%   % (sigma n).n = -p_in  <=>  (sigma n)_x = +p_in;  tang(u) = u_y = 0
%
% A three-component trace (the Biot trace (eta_x, eta_y, p)) takes a third
% component through p: "dirichlet" prescribes p_hat = p_data; "neumann"
% tests p_data against the trace basis and adds it to the load of the
% p_hat row. Without p the third component is natural with zero data
% (zero Darcy flux), as before.
%   boundary_condition("Bottom", p = "dirichlet")   % eta = 0, p = 0
%
% Assemble several into a struct array and store it in
% problem.Data.boundary_conditions; the mesh must carry matching face tags
% (see SimplicialGrid.apply_tag). Omitted data defaults to zero.
arguments
    tag             (1,1) string
    options.x       (1,1) string {mustBeMember(options.x, ["dirichlet", "neumann"])} = "dirichlet"
    options.y       (1,1) string {mustBeMember(options.y, ["dirichlet", "neumann"])} = "dirichlet"
    options.x_data  (1,1) function_handle = @(x, y) 0*x
    options.y_data  (1,1) function_handle = @(x, y) 0*x
    options.p       (1,1) string {mustBeMember(options.p, ["none", "dirichlet", "neumann"])} = "none"
    options.p_data  (1,1) function_handle = @(x, y) 0*x
end
bc = struct( ...
    "tag",  tag, ...
    "type", [options.x; options.y], ...
    "data", {{options.x_data; options.y_data}});
if options.p ~= "none"
    bc.type = [bc.type; options.p];
    bc.data = [bc.data; {options.p_data}];
end
end
