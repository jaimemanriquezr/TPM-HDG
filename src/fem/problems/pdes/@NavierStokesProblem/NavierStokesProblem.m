classdef NavierStokesProblem < IncompressibleElasticityProblem
%-------------------------------------------------------------------------
%Description: Oseen / Navier–Stokes on the stress-velocity (PEERS) fluid:
%               (1/2μ)σ − ε(u) + γ = 0,
%               ρ u_t + ρ(β·∇)u − ∇·σ = f,   ∇·u = 0
%             with β a given convective velocity kept out of the stress:
%             it enters the momentum row and the numerical flux only
%             (Solano & Vargas 2022, eq. 3.4), so the transfer machinery
%             built on σ is unchanged. Steady (density = 0) uses the
%             kinematic form (β·∇)u. β is {bx; by} as function handles or
%             as Dubiner dof arrays (the previous Picard iterate or the
%             previous time step, see solve_navier_stokes). Time
%             stepping, the state and initial_state as in
%             IncompressibleElasticityProblem; the acceleration problem
%             of initial_state takes β = u(0).
%-------------------------------------------------------------------------
    methods
        function obj = NavierStokesProblem(k, label, options)
            arguments (Input)
                k       double          {mustBeScalarOrEmpty} = []
                label   string          {mustBeScalarOrEmpty} = "Navier-Stokes (Stress-Velocity)"
                options.density                 (1,1)   double = 0.0
                options.time                    (1,1)   double = 0.0   % initial time t0: seeds state.time; not read by the solve
                options.time_step               (1,1)   double = 1e-1
                options.initial_velocity        (2,1)   cell   = {@(x,y) 0*(x.*y); @(x,y) 0*(x.*y)}
                options.time_scheme             (1,1)   string {mustBeMember(options.time_scheme, ["euler", "cranknicolson"])} = "euler"
                options.theta                   (1,1)   double = NaN
                options.state                           struct = struct.empty
                options.viscosity               (1,1)   double = 1e-3
                options.domain_measure          (1,1)   double = 1.0
                options.pressure_mean           (1,1)   string {mustBeMember(options.pressure_mean, ["constrained", "pinned"])} = "constrained"
                options.source_term             (2,1)   cell   = {@(x,y) 0*(x.*y); @(x,y) 0*(x.*y)}
                options.datum_dirichlet         (2,1)   cell   = {@(x,y) 0*(x.*y); @(x,y) 0*(x.*y)}
                options.datum_neumann           (2,1)   cell   = {@(x,y) 0*(x.*y); @(x,y) 0*(x.*y)}
                options.datum_dirichlet_rate    (2,1)   cell   = {@(x,y) 0*(x.*y); @(x,y) 0*(x.*y)}
                options.convective_velocity             cell   = {}   % {bx; by}; empty -> Stokes
                options.transfer_mapping {mustBeA(options.transfer_mapping, ["function_handle", "LevelSetCurve"])} = function_handle.empty
                % false: keep S = mu*(mesh tau) under convection; the caller sets tau(beta) on the mesh
                options.convective_stabilization (1,1) logical = true
            end
            beta = options.convective_velocity;
            options = rmfield(options, "convective_velocity");
            args = namedargs2cell(options);
            obj = obj@IncompressibleElasticityProblem(k, label, args{:});
            obj.Data.convective_velocity = beta;
        end
    end

    methods (Static)
        data = get_manufactured_solution(U, p, shear_viscosity, beta, rho, time, filename)
    end
end
