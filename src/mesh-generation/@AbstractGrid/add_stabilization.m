function obj = add_stabilization(obj, tau)
%ADD_STABILIZATION  Set the HDG stabilization parameter of a mesh.
%
%   mesh = add_stabilization(mesh, tau)
%
%   tau  a scalar c, giving tau = c on every side of every element, or a
%        function handle called as tau(mesh) that returns the
%        3-by-n_elements array (default @S_identity, i.e. 1). The result
%        is stored in mesh.tau.
%
%See also: S_identity
arguments
    obj AbstractGrid
    tau (1,1) {mustBeStabilization} = @S_identity
end
if isa(tau, 'double')
    fun_tau = @(obj) S_identity(obj, tau);
else
    fun_tau = tau;
end
obj.tau = feval(fun_tau, obj);
end

function mustBeStabilization(a)
assert(isa(a,'function_handle') || isa(a,'double'),...
    'mustBeStabilization:nostabilization',...
    'Input must be either a function handle or scalar double.')
end
