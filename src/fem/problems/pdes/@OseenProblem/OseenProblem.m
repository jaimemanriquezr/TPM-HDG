classdef OseenProblem < StokesProblem
%-------------------------------------------------------------------------
%Description: HDG discretization of the Oseen problem (Solano & Vargas
%             2022, eq. 3.4) in velocity-gradient form:
%               L − ∇u = 0,   −∇·(νL − pI) + (β·∇)u = f,   ∇·u = 0
%             with β a given divergence-free convective velocity. Same
%             unknowns and spaces as StokesProblem; the convective term
%             lives in the momentum row and the numerical flux, and the
%             stabilization is τ = max|β·n|/(2ν) + 1 (their (3.5)).
%             β is {bx; by} as function handles, or as Dubiner dof arrays
%             (d2 x Nelts per component) — the previous iterate of a
%             Picard loop, see solve_navier_stokes.
%-------------------------------------------------------------------------
    methods
        function obj = OseenProblem(k, label, options)
            arguments (Input)
                k       double          {mustBeScalarOrEmpty} = []
                label   string          {mustBeScalarOrEmpty} = "Oseen (Velocity-Pressure)"
                options.viscosity           (1,1)   double = 1
                options.domain_measure      (1,1)   double = 1
                options.source_term         (2,1)   cell   = {@(x,y) 0*(x.*y); @(x,y) 0*(x.*y)}
                options.datum_dirichlet     (2,1)   cell   = {@(x,y) 0*(x.*y); @(x,y) 0*(x.*y)}
                options.datum_neumann       (2,1)   cell   = {@(x,y) 0*(x.*y); @(x,y) 0*(x.*y)}
                options.convective_velocity         cell   = {}   % {bx; by}; empty -> Stokes
                options.pressure_mean       (1,1)   string {mustBeMember(options.pressure_mean, ["constrained", "pinned"])} = "constrained"
                options.transfer_mapping {mustBeA(options.transfer_mapping, ["function_handle", "LevelSetCurve"])} = function_handle.empty
            end
            beta = options.convective_velocity;
            options = rmfield(options, "convective_velocity");
            args = namedargs2cell(options);
            obj = obj@StokesProblem(k, label, args{:});
            obj.Data.convective_velocity = beta;
        end
    end

    methods (Static)
        data = get_manufactured_solution(U, p, kinematic_viscosity, beta, filename)
    end
end
