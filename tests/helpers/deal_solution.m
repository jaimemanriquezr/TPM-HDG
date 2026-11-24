function [uF, uS, eF] = deal_solution(cp, k, fluid_mesh, pd)
[fs, ss] = solve(cp, k);
uF = fs.Unknowns.field(2).dof;  uS = ss.Unknowns.field(2).dof;
eF = compute_l2_error(fs.Unknowns.field(2), fluid_mesh, Datum(pd.u_F));
end
