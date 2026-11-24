function p = oseen_stress(beta, varargin)
% the stress-form Oseen problem: the incompressible-elasticity fluid with a
% prescribed convective velocity (a NavierStokesProblem always iterates)
p = IncompressibleElasticityProblem(varargin{:});
p.Data.convective_velocity = beta;
end
