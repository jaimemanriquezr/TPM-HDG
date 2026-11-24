function fsi_dsl_picard_convergence(k_list, n_levels, viscosity)
% Convergence tables of the Navier-Stokes/elasticity FSI path on a flat
% interface y = 1/2 with C^inf manufactured fields:
%   (a) Oseen-FSI with the frozen exact beta = u_F on fitted meshes (no
%       gap), tau = 1 + max|u_F|/(2 mu): EOC of the fluid fields ~ k+1;
%   (b) fsi_dsl_picard on the gap case fitted_split(interface, h, h^2):
%       Picard iterations per level and EOC (fluid and solid combined).
% Prints the tables; nothing is asserted.
arguments
    k_list   = [1, 2]
    n_levels = 3
    viscosity = 1.0
end
y0 = 0.5;  h0 = 1/4;
young_modulus = 1.0;  poisson_ratio = 0.33;  dt = 1.0;
psi = "sin(pi*x)*sin(pi*y)";
interface = FlatCurveHorizontal(y0);
material = struct(viscosity = viscosity, young_modulus = young_modulus, ...
    poisson_ratio = poisson_ratio, time_step = dt);
max_u = pi;     % max |u_F| = max |(-pi sin(pi x) cos(pi y), pi cos(pi x) sin(pi y))|

for k = k_list
    [~, ~, interface_sources, pd] = FluidStructureCoupling.load_problems( ...
        k, interface, 0.0, dt, psi, "0", ["cos(pi*x)"; "sin(pi*y)"], ...
        0.0, viscosity, 0.0, young_modulus, poisson_ratio, 0.0, ...
        pressure_mean = "pinned");
    pd = fsi_dsl_convective_source(pd, psi);

    % (a) frozen exact beta, fitted
    [e_a, h_a] = deal(zeros(n_levels, 2), zeros(n_levels, 1));
    for i = 1:n_levels
        h = h0/2^(i - 1);
        [fluid_mesh, solid_mesh] = fitted_split(interface, h, 0);
        fluid_mesh = add_stabilization(fluid_mesh, 1 + max_u/(2*viscosity));
        coupled = fsi_dsl_coupled_problem(k, fluid_mesh, solid_mesh, interface, ...
            material, pd, interface_sources, Datum(pd.u_F));
        [fs, ~] = solve(coupled, k);
        e_F = fsi_dsl_field_errors(fs, fluid_mesh, pd.sigma_F, pd.u_F, pd.gamma_F);
        e_a(i, :) = e_F(1:2).';
        h_a(i) = h;
    end
    fprintf("\n(a) k = %d, mu = %g, frozen beta = u_F, fitted\n", k, viscosity);
    print_table(h_a, e_a, ["e(sigma_F)", "e(u_F)"], []);

    % (b) Picard, gap h^2
    [e_b, h_b, its] = deal(zeros(n_levels, 3), zeros(n_levels, 1), zeros(n_levels, 1));
    for i = 1:n_levels
        h = h0/2^(i - 1);
        [fluid_mesh, solid_mesh] = fitted_split(interface, h, h^2);
        [fs, ss, its(i), increments] = fsi_dsl_picard(k, fluid_mesh, solid_mesh, ...
            interface, material, pd, interface_sources);
        fprintf("  level %d: increments %s\n", i, sprintf("%.1e ", increments));
        e_F = fsi_dsl_field_errors(fs, fluid_mesh, pd.sigma_F, pd.u_F, pd.gamma_F);
        e_S = fsi_dsl_field_errors(ss, solid_mesh, pd.sigma_S, pd.u_S, pd.gamma_S);
        e = sqrt(e_F.^2 + e_S.^2);
        e_b(i, :) = e([1, 2, 4]).';
        h_b(i) = h;
    end
    fprintf("\n(b) k = %d, mu = %g, Picard, gap h^2\n", k, viscosity);
    print_table(h_b, e_b, ["e(sigma)", "e(u)", "e(u_hat)"], its);
end
end

%=========================================================================%
function print_table(h, e, names, its)
r = [NaN(1, size(e, 2)); log(e(1:end-1, :)./e(2:end, :))./log(h(1:end-1)./h(2:end))];
fprintf("    h    ");
fprintf("| %10s  rate ", names);
if ~isempty(its), fprintf("| its"); end
fprintf("\n");
for i = 1:numel(h)
    fprintf(" %7.4f ", h(i));
    fprintf("| %10.3e %5.2f ", [e(i, :); r(i, :)]);
    if ~isempty(its), fprintf("| %3d", its(i)); end
    fprintf("\n");
end
end
