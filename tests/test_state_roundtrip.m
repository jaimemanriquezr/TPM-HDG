% test_state_roundtrip.m
% The FSI state is plain data: on the polynomial problem of
% test_fsi_transient_polynomial (flat Sigma, fitted meshes retreated by
% 0.02, k = 2), for "newmark" (Crank-Nicolson fluid) and "onestep"
% (backward-Euler fluid):
%  1. one FluidStructureCoupling.solve step from
%     FluidStructureCoupling.initial_state equals, bitwise (dofs and the
%     returned state), the same step from a copy of the initial state
%     written to and read back from a .mat file;
%  2. the second step from the returned state equals, bitwise, the same
%     step from the state advanced by hand (get_kinematic_state of each
%     side, merged: solid, fluid, time and step from the solid);
%  3. after each step the state carries the step's own unknowns:
%     state.fluid.velocity = fluid field(2), state.solid.stress = solid
%     field(1);
%  4. a dof array in a Data.initial_* field, which the coupling would not
%     read, is rejected (FluidStructureCoupling:dofHistory).
clc; clear;
n_fail = 0;
k = 2; h = 1/8; gap = 0.02; dt = 0.1;
E = 1.0; nu = 0.25; visc = 0.4; rho = 1.0;
data = {"x*y*(1 + t)", "(x + y)*(1 + t)", ["x*(1 + t + t^2/2)"; "y*(1 + t + t^2/2)"]};
interface = FlatCurveHorizontal();
[fm, sm] = fitted_split(interface, h, gap);
file = [tempname, '.mat'];
cleanup = onCleanup(@() delete(file));

for scheme = ["newmark", "onestep"]
    [fp, sp, src, d] = problems(k, interface, 0.0, dt, data, rho, visc, E, nu, scheme);
    state0 = FluidStructureCoupling.initial_state(fp, sp, fm, sm, interface, src{:}, ...
        source_dirichlet_rate = d.T_2_t);
    save(file, "state0");
    copy = load(file);

    %% 1. a step from the state and from its .mat copy
    [fp1, sp1, src1] = problems(k, interface, dt, dt, data, rho, visc, E, nu, scheme);
    sp1.Data.time = state0.time;   % the onestep Taylor start is gated by Data.time = 0
    [fs1, ss1, state1] = FluidStructureCoupling.solve(fp1, sp1, fm, sm, interface, src1{:}, ...
        scheme, state = state0);
    [fs1c, ss1c, state1c] = FluidStructureCoupling.solve(fp1, sp1, fm, sm, interface, src1{:}, ...
        scheme, state = copy.state0);
    same_1 = isequal(dofs(fs1, ss1), dofs(fs1c, ss1c)) && isequal(state1, state1c);

    %% 2. the second step from the returned state and from the hand-advanced one
    solid = get_kinematic_state(sp1, ss1, state0);
    fluid = get_kinematic_state(fp1, fs1, state0);
    hand = state0;
    hand.solid = solid.solid;  hand.fluid = fluid.fluid;
    hand.time = solid.time;    hand.step = solid.step;
    [fp2, sp2, src2] = problems(k, interface, 2*dt, dt, data, rho, visc, E, nu, scheme);
    sp2.Data.time = state1.time;
    [fs2, ss2, state2] = FluidStructureCoupling.solve(fp2, sp2, fm, sm, interface, src2{:}, ...
        scheme, state = state1);
    [fs2h, ss2h, state2h] = FluidStructureCoupling.solve(fp2, sp2, fm, sm, interface, src2{:}, ...
        scheme, state = hand);
    same_2 = isequal(hand, state1) && isequal(dofs(fs2, ss2), dofs(fs2h, ss2h)) ...
        && isequal(state2, state2h);

    %% 3. the state carries the step's unknowns
    carried = isequal(state1.fluid.velocity, fs1.Unknowns.field(2).dof) ...
        && isequal(state1.solid.stress, ss1.Unknowns.field(1).dof) ...
        && isequal(state2.fluid.velocity, fs2.Unknowns.field(2).dof) ...
        && isequal(state2.solid.stress, ss2.Unknowns.field(1).dof) ...
        && state2.step == 2 && abs(state2.time - 2*dt) < 1e-15;

    fprintf("%-8s .mat copy, bitwise: %d   hand-advanced state, bitwise: %d   " + ...
        "state carries the step's unknowns: %d\n", scheme, same_1, same_2, carried);
    if ~(same_1 && same_2 && carried)
        fprintf("%s round trip FAILED\n", scheme); n_fail = n_fail + 1;
    end
end

%% 4. dof-array histories on Data are rejected
fp1.Data.initial_velocity = state1.fluid.velocity;
try
    FluidStructureCoupling.solve(fp1, sp1, fm, sm, interface, src1{:}, "onestep", state = state0);
    id = "";
catch err
    id = err.identifier;
end
fprintf("Data.initial_velocity as dofs raises: %s\n", id);
if id ~= "FluidStructureCoupling:dofHistory"
    fprintf("dof-array rejection FAILED\n"); n_fail = n_fail + 1;
end

fprintf("%d check(s) failed.\n", n_fail);
if n_fail > 0
    error("test_state_roundtrip:failed", "%d checks failed", n_fail);
end

%% ------------------------------------------------------------------------
function [fp, sp, src, d] = problems(k, interface, t, dt, data, rho, visc, E, nu, scheme)
% load_problems at time t; newmark + Crank-Nicolson or onestep + backward Euler
[fp, sp, src, d] = FluidStructureCoupling.load_problems(k, interface, t, dt, data{:}, ...
    rho, visc, rho, E, nu, 0);
sp.Data.time_scheme = scheme;
fp.Data.time_scheme = "euler";
if scheme == "newmark", fp.Data.time_scheme = "cranknicolson"; end
end

function D = dofs(fs, ss)
D = {fs.Unknowns.field.dof, fs.Unknowns.hybrid.dof, ss.Unknowns.field.dof, ss.Unknowns.hybrid.dof};
end
