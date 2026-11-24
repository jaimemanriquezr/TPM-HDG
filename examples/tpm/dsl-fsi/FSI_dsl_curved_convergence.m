function results = FSI_dsl_curved_convergence(k, max_iterations, h0, options)
% Convergence of steady FSI on a curved interface, smooth manufactured
% solution with nonzero interface data. Prints the L2 errors (fluid and
% solid combined) of sigma_h, u_h, gamma_h, u_hat and ||P_M u - u_hat||_h,
% with their EOC.
%
%   topology: "open"   y = 1/2 + cos(2 pi x)/10, meeting the walls at right
%                      angles and splitting the square into two simply
%                      connected sets (fluid above); closest-point
%                      interface map
%             "closed" circle of radius 0.2 centred at (1/2, 1/2) (solid
%                      inside, fluid outside)
%   mesh:     "fitted"     both sides meshed to the linear interpolant of
%                          Sigma, no gap (h0 the mesh size)
%             "background" active meshes of two crossed background grids,
%                          h_fluid = h_solid/2 (h0 the solid grid size);
%                          the cut elements are discarded, a gap band
%                          holds Sigma
%   navier_stokes: convective fluid (Picard iteration)
%   picard_tolerance, picard_max_iterations: stopping rule of the
%       Picard iteration (navier_stokes only), see fsi_dsl_picard
arguments
    k               = 2
    max_iterations  = 4
    h0              = []    % 1/4 open fitted, else 1/8
    options.topology (1,1) string {mustBeMember(options.topology, ["open", "closed"])} = "open"
    options.mesh (1,1) string {mustBeMember(options.mesh, ["fitted", "background"])} = "fitted"
    options.navier_stokes (1,1) logical = false
    options.picard_tolerance (1,1) double = 1e-10
    options.picard_max_iterations (1,1) double = 30
end
material = struct(viscosity = 1.0, young_modulus = 1.0, poisson_ratio = 0.3, ...
    time_step = 1.0);
fields = struct(uu_fluid = "sin(pi*x)*sin(pi*y)", ...
    p_fluid = "cos(pi*x)*exp(y)", ...
    u_solid = ["sin(pi*x)*exp(y)"; "cos(pi*y)*exp(-x)"]);
switch options.topology
    case "open"
        interface = SmoothCurveHorizontal("0.5 + 0.1*cos(2*pi*t)", mapping = "closest");
    case "closed"
        interface = CircularCurve([0.5; 0.5], 0.2, mapping = "closest");
end
family = replace(options.mesh, "background", "active");
if isempty(h0)
    % coarser: no solid elements inside the circle (closed background),
    % the circle a 16-gon at both of the first two levels (closed fitted)
    h0 = 1/8;
    if options.topology == "open" && options.mesh == "fitted"
        h0 = 1/4;
    end
end
h_levels = h0*(1/2).^(0:max_iterations - 1);

results = fsi_dsl_sweep(k, interface, ...
    @(h) fsi_dsl_meshes(interface, family, h), h_levels, material, fields, ...
    navier_stokes = options.navier_stokes, ...
    picard_tolerance = options.picard_tolerance, ...
    picard_max_iterations = options.picard_max_iterations, ...
    title = sprintf("k = %d, %s interface, %s meshes", k, ...
        options.topology, options.mesh));
end
