function results = FSI_dsl_lambda_convergence(k, lambda_list, max_iterations, h0, options)
% Convergence of steady FSI towards the incompressible solid limit: the
% solid displacement is divergence-free (the curl of a smooth stream
% function), so the exact stress 2 mu_S eps(u_S) and the source do not
% depend on the Lame lambda, which grows over lambda_list at fixed shear
% modulus mu_S. Fitted meshes, no gap. For each lambda prints the L2
% errors (fluid and solid combined) of sigma_h, u_h, gamma_h, u_hat and
% ||P_M u - u_hat||_h, with their EOC, then the same for the solid alone;
% robust rates do not degrade with lambda.
%
%   geometry: "flat" (y = 1/2), "open" (y = 1/2 + cos(2 pi x)/10) or
%       "closed" (circle of radius 0.2, solid inside)
%   navier_stokes: convective fluid (Picard iteration)
%   picard_tolerance, picard_max_iterations: stopping rule of the
%       Picard iteration (navier_stokes only), see fsi_dsl_picard
%   results(j): h, errors, eoc, errors_solid, eoc_solid, n_picard for
%       lambda_list(j)
arguments
    k               = 3
    lambda_list     = [1, 1e2, 1e4, 1e6]
    max_iterations  = 4
    h0              = []    % 1/4, or 1/8 for "closed"
    options.geometry (1,1) string {mustBeMember(options.geometry, ["flat", "open", "closed"])} = "flat"
    options.shear_modulus (1,1) double {mustBePositive} = 1.0
    options.navier_stokes (1,1) logical = false
    options.picard_tolerance (1,1) double = 1e-10
    options.picard_max_iterations (1,1) double = 30
end
switch options.geometry
    case "flat"
        interface = FlatCurveHorizontal(0.5);
    case "open"
        interface = SmoothCurveHorizontal("0.5 + 0.1*cos(2*pi*t)", mapping = "closest");
    case "closed"
        interface = CircularCurve([0.5; 0.5], 0.2, mapping = "closest");
end
if isempty(h0)
    h0 = 1/4;
    if options.geometry == "closed"
        h0 = 1/8;
    end
end
h_levels = h0*(1/2).^(0:max_iterations - 1);

syms x y
phi_S = sin(pi*x)*cos(pi*y)*exp(x)/pi;
fields = struct(uu_fluid = "sin(pi*x)*sin(pi*y)", ...
    p_fluid = "cos(pi*x)*exp(y)", ...
    u_solid = string([diff(phi_S, y); -diff(phi_S, x)]));

mu_S = options.shear_modulus;
for j = 1:numel(lambda_list)
    lambda = lambda_list(j);
    material = struct(viscosity = 1.0, ...
        young_modulus = mu_S*(3*lambda + 2*mu_S)/(lambda + mu_S), ...
        poisson_ratio = lambda/(2*(lambda + mu_S)), time_step = 1.0);
    results(j) = fsi_dsl_sweep(k, interface, ...
        @(h) fsi_dsl_meshes(interface, "fitted", h), h_levels, material, fields, ...
        navier_stokes = options.navier_stokes, solid_errors = true, ...
        picard_tolerance = options.picard_tolerance, ...
        picard_max_iterations = options.picard_max_iterations, ...
        title = sprintf("k = %d, %s interface, lambda = %g, mu_S = %g", ...
            k, options.geometry, lambda, mu_S));
end
end
