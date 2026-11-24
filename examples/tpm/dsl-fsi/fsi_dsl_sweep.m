function result = fsi_dsl_sweep(k, interface, make_meshes, h_levels, material, fields, options)
% One mesh-refinement sweep of the steady FSI problem: for each h in
% h_levels builds [fluid_mesh, solid_mesh] = make_meshes(h), solves, and
% prints the L2 errors (fluid and solid combined) of sigma_h, u_h,
% gamma_h, u_hat and the h-weighted ||P_M u - u_hat||_h, with their EOC
% against 1/sqrt(n_elements).
%
%   fields: struct with uu_fluid (stream function, u_F = curl), p_fluid,
%       u_solid (2x1 string)
%   result: struct with h, errors (one column per quantity, in the order
%       above), eoc (first row NaN), errors_solid and eoc_solid (the solid
%       side alone), n_picard (NaN without navier_stokes)
%   solid_errors: also print the solid-side errors in a table of their own
%   picard_tolerance, picard_max_iterations: stopping rule of the
%       Picard iteration (navier_stokes only), see fsi_dsl_picard
arguments
    k
    interface AbstractCurve
    make_meshes function_handle
    h_levels double
    material struct
    fields struct
    options.navier_stokes (1,1) logical = false
    options.title (1,1) string = ""
    options.solid_errors (1,1) logical = false
    options.picard_tolerance (1,1) double = 1e-10
    options.picard_max_iterations (1,1) double = 30
end
[~, ~, interface_sources, pd] = FluidStructureCoupling.load_problems( ...
    k, interface, 0.0, material.time_step, fields.uu_fluid, fields.p_fluid, ...
    fields.u_solid, 0.0, material.viscosity, 0.0, material.young_modulus, ...
    material.poisson_ratio, 0.0, pressure_mean = "pinned");
if options.navier_stokes
    pd = fsi_dsl_convective_source(pd, fields.uu_fluid);
end

n_levels = numel(h_levels);
[errors, errors_solid, n_picard] = deal(zeros(n_levels, 5), zeros(n_levels, 5), NaN(n_levels, 1));
h_list = zeros(n_levels, 1);
for i = 1:n_levels
    disp(options.title + ": level " + i + "/" + n_levels);
    [fluid_mesh, solid_mesh] = make_meshes(h_levels(i));
    if options.navier_stokes
        [fluid_solution, solid_solution, n_picard(i)] = fsi_dsl_picard(k, ...
            fluid_mesh, solid_mesh, interface, material, pd, interface_sources, ...
            tolerance = options.picard_tolerance, ...
            max_iterations = options.picard_max_iterations);
    else
        coupled = fsi_dsl_coupled_problem(k, fluid_mesh, solid_mesh, interface, ...
            material, pd, interface_sources);
        [fluid_solution, solid_solution] = solve(coupled, k);
    end
    e_F = fsi_dsl_field_errors(fluid_solution, fluid_mesh, pd.sigma_F, pd.u_F, pd.gamma_F);
    e_S = fsi_dsl_field_errors(solid_solution, solid_mesh, pd.sigma_S, pd.u_S, pd.gamma_S);
    errors(i, :) = sqrt(e_F.^2 + e_S.^2).';
    errors_solid(i, :) = e_S.';
    h_list(i) = 1/sqrt(size(fluid_mesh.elements, 2) + size(solid_mesh.elements, 2));
end
eoc = [NaN(1, 5); compute_eoc(errors, h_list)];
eoc_solid = [NaN(1, 5); compute_eoc(errors_solid, h_list)];

print_table(options.title, h_list, errors, eoc, n_picard, options.navier_stokes);
if options.solid_errors
    print_table(options.title + " (solid)", h_list, errors_solid, eoc_solid, ...
        n_picard, options.navier_stokes);
end
result = struct(h = h_list, errors = errors, eoc = eoc, ...
    errors_solid = errors_solid, eoc_solid = eoc_solid, n_picard = n_picard);
end

%=========================================================================%
function print_table(title, h_list, errors, eoc, n_picard, navier_stokes)
fprintf("\n%s\n", title);
hdr = "     h    |  e(σ_h)   | r(σ_h) |   e(u_h)  | r(u_h) |  e(γ_h)   | r(γ_h) |  e(û_h)   | r(û_h) | e(Pu-û_h) | r(Pu-û_h)";
if navier_stokes
    hdr = hdr + " | Picard";
end
fprintf("%s\n%s\n", hdr, repmat('-', 1, strlength(hdr)));
for i = 1:numel(h_list)
    fprintf(" %8.2e |" + join(repmat(" %9.3e | %6.2f", 1, 5), " |"), ...
        h_list(i), reshape([errors(i, :); eoc(i, :)], 1, []));
    if navier_stokes
        fprintf(" | %6d", n_picard(i));
    end
    fprintf("\n");
end
end
