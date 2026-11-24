function stokes_darcy_flat_example(k, h, y0, d, viscosity, diffusion, friction_bjs, uu_stokes, p_stokes, p_darcy)
% Steady Stokes-Darcy coupling (Beavers-Joseph-Saffman friction) across a
% flat interface y = y0, fitted meshes with a gap d bridged by TPM. Prints
% the L2 errors on each side and plots both solutions.
arguments
    k           = 3
    h           = 1/8
    y0          = 0.5
    d           = 5e-2
    viscosity   = 1.0
    diffusion   = eye(2)
    friction_bjs = 10.0
    uu_stokes   = "x^2*y^2"
    p_stokes    = "x^2 + y^2"
    p_darcy     = "x^3 + y^3"
end

data = StokesDarcyCoupling.get_manufactured_solution(...
    uu_stokes, p_stokes, p_darcy, y0, viscosity, diffusion, friction_bjs);

I = FlatCurveHorizontal(y0);
interface = I;
% Marker 2 = Neumann on the Darcy side: its datum below is the flux
% u_D . n, so the outer boundary must be Neumann. With the default
% all-Dirichlet markers the Darcy trace is pinned to an unset (zero)
% pressure datum and the error saturates at O(1).
[upper_mesh, lower_mesh] = fitted_split(I, [h, h], d, boundary_markers = [1, 2]);

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

[stokes_solution, darcy_solution] = StokesDarcyCoupling.solve(...
    stokes_problem, darcy_problem, upper_mesh, lower_mesh, ...
    interface, friction_bjs, data.T_1, data.T_2, ...
    Order = k);

[error_stokes, ~] = stokes_solution.compute_l2_error({data.L_S, data.u_S, data.p_S}, {data.u_S});
[error_darcy, ~]  = darcy_solution.compute_l2_error({data.u_D, data.p_D}, {data.p_D});
fprintf("Stokes — e(L_h): %g | e(u_h): %g | e(p_h): %g\n", ...
    error_stokes(1), error_stokes(2), error_stokes(3));
fprintf("Darcy  — e(u_h): %g | e(p_h): %g\n", ...
    error_darcy(1), error_darcy(2));

plot(stokes_solution);
plot(darcy_solution);
end
