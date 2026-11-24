function r = solve_roundoff(global_system)
%SOLVE_ROUNDOFF  Round-off level of a global solve, by one refinement step.
%   r = solve_roundoff(global_system) solves the free-dof system of
%   global_system (the second output of solve: matrix, rhs, dof) again,
%   applies one step of iterative refinement and returns the relative size
%   of the correction, |dx|/|x|. It estimates the error the direct solve
%   leaves in the solution, which grows with the condition number of the
%   system (about 1e9 to 1e10 for the k = 3 transfer systems with a gap of
%   order h), so a patch test can tell round-off from a discretization
%   error.
A = global_system.matrix(global_system.dof, global_system.dof);
b = global_system.rhs(global_system.dof);
x = A\b;
dx = A\(b - A*x);
r = norm(dx)/max(norm(x), realmin);
end
