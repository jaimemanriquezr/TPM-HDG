classdef DarcyProblem < FEProblem
%-------------------------------------------------------------------------
%Description: HDG discretization of the Darcy flow problem:
%             Find (u,p) in V^2 x V such that
%               u + κ*∇p = 0,    −∇·u = f   in Ω
%             Hybrid unknown: pressure trace p̂ on mesh skeleton.
%-------------------------------------------------------------------------
    methods
        function obj = DarcyProblem(options)
            arguments (Input)
                options.Label            string          {mustBeScalarOrEmpty} = "Darcy"  % display label
                options.diffusion (2,2)                                 = eye(2)   % diffusion tensor κ (constant or cell of functions)
                options.source   (1,1)           function_handle = @(x,y) 0*(x.*y)  % source f
                options.datum_dirichlet  (1,1)           function_handle = @(x,y) 0*(x.*y)  % Dirichlet datum g_D
                options.datum_neumann    (1,1)           function_handle = @(x,y) 0*(x.*y)  % Neumann datum q_N
            end
            obj = obj@FEProblem(options.Label, Data=rmfield(options, "Label"));
            obj.FieldsLocal = [struct(name="u_h",       description="velocity");
                                 struct(name="p_h",       description="pressure")];
            obj.FieldsHybrid = [struct(name="\hat{p}_h", description="pressure trace")];
        end
    end

    methods (Static)
        spaces   = get_discrete_spaces(k);

        data = get_manufactured_solution(p, kappa, filename);

        k_I       = get_inverse_tensor(k, x, y);
    end
end
