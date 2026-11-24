classdef Stabilization < MeshFunction
    % HDG stabilization parameter tau as a facet weight function: reads
    % mesh.tau (3 x n_elements, one value per local face) at assemble
    % time, so symbolic forms follow whatever stabilization the mesh
    % carries (add_stabilization) instead of a hard-coded coefficient.
    %
    % Usage:
    %   tau_ = Stabilization();
    %   local = ... + tau_*u*w*dS - tau_*u_hat*w*dS + ...
    methods
        function obj = Stabilization()
            tau = @(mesh) permute(mesh.tau, [3 4 5 6 2 1]);
            obj = obj@MeshFunction(tau);
        end
    end
end
