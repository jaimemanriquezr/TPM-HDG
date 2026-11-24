function theta = get_theta(obj)
% θ of the fluid momentum time discretization (single source of truth,
% shared by get_matrices and get_kinematic_state):
%   data.theta, when set (finite), wins;
%   otherwise data.time_scheme selects a preset:
%     "euler" -> θ = 1 (backward Euler), "cranknicolson" -> θ = 1/2.
arguments
    obj (1,1) IncompressibleElasticityProblem
end
data = obj.Data;
if isfield(data, "theta") && ~isnan(data.theta)
    theta = data.theta;
else
    theta = 1.0;
    if isfield(data, "time_scheme") && data.time_scheme == "cranknicolson"
        theta = 0.5;
    end
end
mustBeInRange(theta, 0, 1, "exclude-lower")   % explicit schemes (θ = 0) unsupported
end
