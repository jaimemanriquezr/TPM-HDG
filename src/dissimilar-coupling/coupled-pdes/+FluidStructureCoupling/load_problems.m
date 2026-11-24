function [fluid_problem, solid_problem, interface_sources, problem_data] = load_problems(k, ...
    interface, ...
    time, time_step, ...
    uu_fluid, p_fluid, u_solid, ...
    fluid_density, fluid_viscosity, ...
    solid_density, solid_young_modulus, solid_poisson_ratio, solid_hooke_modulus, options)
%LOAD_PROBLEMS  Manufactured fluid and solid problems for a TPM-coupled FSI test.
%
%   [fluid_problem, solid_problem, interface_sources, problem_data] = ...
%       load_problems(k, interface, time, time_step, uu_fluid, p_fluid, ...
%       u_solid, fluid_density, fluid_viscosity, solid_density, ...
%       solid_young_modulus, solid_poisson_ratio, solid_hooke_modulus, ...
%       name = ..., swap_domains = ..., pressure_mean = ...)
%
%   builds, from symbolic exact solutions, the fluid and solid problems of
%   a coupled test together with the interface data that make the exact
%   solution satisfy the coupling conditions. Positional arguments, with
%   their defaults:
%     k                    polynomial degree (required)
%     interface            the interface curve (required)
%     time                 time at which the data are evaluated (0.0)
%     time_step            time step (1e-3)
%     uu_fluid             string in x, y (t for time dependence): the
%                          fluid stream function, u_F = (-d/dy, d/dx) of
%                          it ("sin(pi*x)*sin(pi*y)")
%     p_fluid              string: the fluid pressure
%                          ("sin(pi*x)*exp(y - 0.5)")
%     u_solid              2 x 1 string array: the solid displacement
%                          (["x*sin(y)"; "y*cos(x)"])
%     fluid_density        (0.0)
%     fluid_viscosity      (1.0)
%     solid_density        (0.0)
%     solid_young_modulus  (1.0)
%     solid_poisson_ratio  (0.25)
%     solid_hooke_modulus  (0.0), the coefficient of the zeroth-order
%                          displacement term of the solid equation
%
%   Name-value options:
%     name            string (default ""): names the scratch file of the
%                     manufactured data; empty derives a name from the
%                     parameters
%     swap_domains    logical (false): exchange the physical roles of the
%                     two subdomains, fluid on the lower or inner side and
%                     solid on the upper or outer side; open interfaces
%                     only, and the caller must swap the meshes to match
%     pressure_mean   "constrained" | "pinned" (default "pinned"): the
%                     pressure-level closure of the fluid problem; the
%                     exact pressure is normalized to zero mean only under
%                     "constrained", which a stand-alone all-Dirichlet
%                     fluid requires
%
%   Outputs:
%     fluid_problem      IncompressibleElasticityProblem with the fluid
%                        source term, Dirichlet datum, domain measure and
%                        the initial data at time: initial_velocity = u_F,
%                        datum_dirichlet_rate = u_t_F
%     solid_problem      LinearElasticityProblem with the solid source
%                        term, Dirichlet datum and the initial data at
%                        time: initial_displacement = u_S,
%                        initial_velocity = u_t_S
%     interface_sources  {T1, T2}, the interface data (traction balance
%                        and kinematic condition) for
%                        FluidStructureCoupling.solve
%     problem_data       struct of the exact fields and sources, as
%                        function handles of (x, y) evaluated at time:
%                        sigma_F, u_F, u_t_F, gamma_F, f_F, measure_F,
%                        viscosity for the fluid; sigma_S, u_S, u_t_S,
%                        u_tt_S, gamma_S, f_S, young_modulus,
%                        poisson_ratio for the solid; and T_1, T_2 and
%                        T_2_t = d/dt T_2 (the source_dirichlet_rate of
%                        FluidStructureCoupling.initial_state)
%
%See also: FluidStructureCoupling.solve, FluidStructureCoupling.get_manufactured_solution
arguments
    k
    interface
    time double = 0.0;
    time_step double = 1e-3;
    uu_fluid string = "sin(pi*x)*sin(pi*y)"
    p_fluid string = "sin(pi*x)*exp(y - 0.5)"
    u_solid string = ["x*sin(y)"; "y*cos(x)"];

    fluid_density double = 0.0
    fluid_viscosity double = 1.0

    solid_density double = 0.0

    solid_young_modulus double = 1.0
    solid_poisson_ratio double = 0.25
    solid_hooke_modulus double = 0.0;

    options.name string = ""
    % Exchange the two subdomains' physical roles: fluid on the lower/inner
    % side, solid on the upper/outer one. Open interfaces only. The caller
    % must swap the meshes to match. See get_manufactured_solution.
    options.swap_domains (1,1) logical = false
    % Pressure-level closure of the fluid problem, passed to
    % IncompressibleElasticityProblem and read back by the coupled
    % solvers. "constrained" is required for a stand-alone all-Dirichlet
    % fluid.
    options.pressure_mean (1,1) string {mustBeMember(options.pressure_mean, ...
        ["constrained", "pinned"])} = "pinned"
end
data = {uu_fluid; p_fluid; u_solid; interface;
    fluid_viscosity;
    solid_young_modulus; solid_poisson_ratio;
    fluid_density; solid_density; solid_hooke_modulus; time; time_step};
% The pressure is normalised only under "constrained", so the exact
% field satisfies the same zero-mean condition as the discrete one
data_filename = FluidStructureCoupling.get_manufactured_solution(data{:}, ...
    name=options.name, swap_domains=options.swap_domains, ...
    normalize_pressure = (options.pressure_mean == "constrained"));
problem_data = load(data_filename);
delete(data_filename);   % a private file, see get_manufactured_solution

% Recompute u_t_S from the u_solid strings rather than trust the file,
% which may hold it as identically zero.
t_sym = sym('t');
u_S_sym = str2sym(u_solid);
% the + 0*x wrapper broadcasts spatially constant components to array shape
vectorize = @(h) @(x, y) h(x, y) + 0*x;
problem_data.u_t_S = arrayfun( ...
    @(c) vectorize(matlabFunction(subs(diff(c, t_sym), t_sym, time), ...
        Vars=[sym('x'), sym('y')])), ...
    u_S_sym(:), UniformOutput=false);

% The time derivatives the initial state needs, from the input strings
% likewise: u_t_F, u_tt_S, and T_2_t from the definition of T_2 in
% get_manufactured_solution (u_F - d/dt u_S, or u_F - u_S/time_step when
% time = Inf or a density vanishes)
x_sym = sym('x');  y_sym = sym('y');
uu_sym = str2sym(uu_fluid);
u_F_sym = [-diff(uu_sym, y_sym); diff(uu_sym, x_sym)];
at_time = @(e) arrayfun(@(c) vectorize(matlabFunction(subs(c, t_sym, time), ...
    Vars=[x_sym, y_sym])), e(:), UniformOutput=false);
problem_data.u_t_F  = at_time(diff(u_F_sym, t_sym));
problem_data.u_tt_S = at_time(diff(u_S_sym, t_sym, 2));
if ~isinf(time) && fluid_density*solid_density > 0
    problem_data.T_2_t = at_time(diff(u_F_sym, t_sym) - diff(u_S_sym, t_sym, 2));
else
    problem_data.T_2_t = at_time(diff(u_F_sym, t_sym) - diff(u_S_sym, t_sym)/time_step);
end

fluid_problem = IncompressibleElasticityProblem(k, time_step = time_step, ...
    source_term = problem_data.f_F, datum_dirichlet = problem_data.u_F, ...
    viscosity = fluid_viscosity, density=fluid_density, ...
    domain_measure = problem_data.measure_F, ...
    pressure_mean = options.pressure_mean, ...
    initial_velocity = problem_data.u_F, datum_dirichlet_rate = problem_data.u_t_F);

solid_problem = LinearElasticityProblem(k, time_step = time_step, time = time, ...
    source_term=problem_data.f_S, datum_dirichlet=problem_data.u_S, ...
    density=solid_density, hooke_modulus = solid_hooke_modulus, ...
    young_modulus=solid_young_modulus, poisson_ratio=solid_poisson_ratio, ...
    initial_displacement = problem_data.u_S, initial_velocity = problem_data.u_t_S);

interface_sources = {problem_data.T_1, problem_data.T_2};
end
