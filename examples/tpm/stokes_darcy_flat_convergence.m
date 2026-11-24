function stokes_darcy_flat_convergence(k, max_iterations, h0, y0, d, viscosity, diffusion, friction_bjs, uu_stokes, p_stokes, p_darcy)
% Convergence study of steady Stokes-Darcy coupling (Beavers-Joseph-
% Saffman friction) across a flat interface y = y0, fitted meshes with a
% gap d bridged by TPM. Prints combined velocity and pressure L2 errors
% with EOC.
arguments
    k               = 3
    max_iterations  = 8
    h0              = 0.25
    y0              = 0.5
    d               = 1e-3
    viscosity       = 1.0
    diffusion       = eye(2)
    friction_bjs    = 10.0
    uu_stokes       = "x^5*y^5"
    p_stokes        = "x^2*y^2"
    p_darcy         = "x^2 + y^2"
end

data = StokesDarcyCoupling.get_manufactured_solution(...
    uu_stokes, p_stokes, p_darcy, y0, viscosity, diffusion, friction_bjs);

I = FlatCurveHorizontal(y0);
interface = I;

% Outward boundary normal of the unit square, contracted with a vector
% field.
n_x = @(x, y) sign(x - 1/2) .* (abs(x - 1/2) == 1/2);
n_y = @(x, y) sign(y - 1/2) .* (abs(y - 1/2) == 1/2);
dot_n = @(u) @(x, y) u{1}(x, y).*n_x(x, y) + u{2}(x, y).*n_y(x, y);

stokes_problem = StokesProblem(viscosity = data.nu, ...
    source_term = data.f_S, datum_dirichlet = data.u_S, ...
    domain_measure = data.domain_measure);
darcy_problem = DarcyProblem(diffusion = data.kappa, ...
    source = data.f_D, datum_neumann = dot_n(data.u_D));

[error_u, error_p, h_list, n_elements] = deal(zeros(max_iterations, 1));
for i = 1:max_iterations
    disp("Iteration " + i + "/" + max_iterations);
    h_i = h0 * (3/4)^(i + 1);
    % Marker 2 = Neumann on the Darcy side, matching its flux datum.
    [upper_mesh, lower_mesh] = fitted_split(I, [h_i, h_i], d, ...
        boundary_markers = [1, 2]);

    [stokes_solution, darcy_solution] = StokesDarcyCoupling.solve(...
        stokes_problem, darcy_problem, upper_mesh, lower_mesh, ...
        interface, friction_bjs, data.T_1, data.T_2, ...
        Order = k);

    [error_stokes, ~] = stokes_solution.compute_l2_error({data.L_S, data.u_S, data.p_S}, {data.u_S});
    [error_darcy,  ~] = darcy_solution.compute_l2_error({data.u_D, data.p_D}, {data.p_D});

    n_elements(i) = size(upper_mesh.elements, 2) + size(lower_mesh.elements, 2);
    h_list(i)     = 1 / sqrt(n_elements(i));
    error_u(i)    = sqrt(error_stokes(2)^2 + error_darcy(1)^2);
    error_p(i)    = sqrt(error_stokes(3)^2 + error_darcy(2)^2);
end

r_u = [NaN; compute_eoc(error_u, h_list)];
r_p = [NaN; compute_eoc(error_p, h_list)];

hdr = "   i  |    h     |   e(u)    | r(u)   |   e(p)    | r(p)";
fprintf("\n%s\n%s\n", hdr, repmat('-', 1, strlength(hdr)));
for i = 1:max_iterations
    if i == 1
        fprintf(" %4d | %8.2e | %9.3e |   -    | %9.3e |   -\n", ...
            i, h_list(i), error_u(i), error_p(i));
    else
        fprintf(" %4d | %8.2e | %9.3e | %6.2f | %9.3e | %6.2f\n", ...
            i, h_list(i), error_u(i), r_u(i), error_p(i), r_p(i));
    end
end
end
