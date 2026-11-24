function pd = fsi_dsl_convective_source(pd, uu_fluid)
% Turns the Stokes fluid source of FluidStructureCoupling.load_problems,
% pd.f_F = -div sigma_F, into the Navier-Stokes one
%     f_F + (grad u_F) u_F,     u_F = [-d_y psi; d_x psi],
% psi = uu_fluid the stream function (the convention of
% NavierStokesProblem.get_manufactured_solution with beta = u). The
% interface data are unchanged: Sigma is a material surface and the
% coupling carries no convective flux across it.
arguments
    pd struct
    uu_fluid (1,1) string
end
syms x y real
psi = str2sym(uu_fluid);
u = [-diff(psi, y); diff(psi, x)];
convection = jacobian(u, [x y])*u;
f_F = pd.f_F;
for i = 1:numel(f_F)
    c_i = matlabFunction(convection(i), Vars = [x y]);
    f_i = f_F{i};
    % + 0*x broadcasts a spatially constant term to array shape
    f_F{i} = @(x, y) f_i(x, y) + c_i(x, y) + 0*x;
end
pd.f_F = f_F;
end
