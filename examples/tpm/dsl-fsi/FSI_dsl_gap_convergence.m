function results = FSI_dsl_gap_convergence(k, s_list, max_iterations, h0, options)
% Convergence of steady FSI on a flat interface y = 1/2 with fitted meshes
% retreating from it, leaving a gap of width h^s that TPM bridges
% (s = Inf: no gap). Smooth manufactured solution with nonzero interface
% data. For each s prints the L2 errors (fluid and solid combined) of
% sigma_h, u_h, gamma_h, u_hat and ||P_M u - u_hat||_h, with their EOC.
%
%   navier_stokes: convective fluid (Picard iteration), with the source
%       of the same exact solution
%   picard_tolerance, picard_max_iterations: stopping rule of the
%       Picard iteration (navier_stokes only), see fsi_dsl_picard
%   results(j): h, errors, eoc, n_picard for s_list(j)
arguments
    k               = 2
    s_list          = [1, 2, Inf]
    max_iterations  = 4
    h0              = 1/4
    options.navier_stokes (1,1) logical = false
    options.picard_tolerance (1,1) double = 1e-10
    options.picard_max_iterations (1,1) double = 30
end
material = struct(viscosity = 1.0, young_modulus = 1.0, poisson_ratio = 0.3, ...
    time_step = 1.0);
fields = struct(uu_fluid = "sin(pi*x)*sin(pi*y)", ...
    p_fluid = "cos(pi*x)*exp(y)", ...
    u_solid = ["sin(pi*x)*exp(y)"; "cos(pi*y)*exp(-x)"]);
interface = FlatCurveHorizontal(0.5);
h_levels = h0*(1/2).^(0:max_iterations - 1);

for j = 1:numel(s_list)
    s = s_list(j);
    results(j) = fsi_dsl_sweep(k, interface, ...
        @(h) fsi_dsl_meshes(interface, "fitted", h, gap = h^s), h_levels, ...
        material, fields, navier_stokes = options.navier_stokes, ...
        picard_tolerance = options.picard_tolerance, ...
        picard_max_iterations = options.picard_max_iterations, ...
        title = sprintf("k = %d, gap = h^%g", k, s));
end
end
