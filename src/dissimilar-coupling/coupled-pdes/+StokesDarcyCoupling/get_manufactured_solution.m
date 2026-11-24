function data = get_manufactured_solution(uu_stokes, p_stokes, p_darcy, interface, nu, kappa, wbjs)
    %GET_MANUFACTURED_SOLUTION  Manufactured data for a Stokes-Darcy test.
    %
    %   data = get_manufactured_solution(uu_stokes, p_stokes, p_darcy, ...
    %       interface, nu, kappa, wbjs)
    %
    %   derives, on the unit square with the Stokes domain above the
    %   interface and the Darcy domain below, the exact fields, sources and
    %   interface data of the coupled problem. All arguments are optional:
    %   uu_stokes the stream function (u_S = (-d/dy, d/dx) of it, strings
    %   in x, y), p_stokes and p_darcy the pressures, interface a
    %   FlatCurveHorizontal, a SmoothCurveHorizontal or a numeric height
    %   (a flat interface), nu the viscosity, kappa the 2 x 2 diffusion
    %   tensor, and wbjs the Beavers-Joseph-Saffman friction coefficient.
    %   The Stokes pressure is shifted to zero mean over its domain.
    %
    %   data  struct of the exact fields and sources: domain_measure, nu,
    %         p_S, u_S, L_S, f_S, kappa, p_D, u_D, f_D, wbjs, T_1, T_2.
    arguments
        % exact solutions
        uu_stokes = "cos(pi*x)*cos(pi*y)";  % stream function: u1 = -duu/dy, u2 = duu/dx
        p_stokes  = "cos(pi*x)*cos(pi*y)";
        p_darcy   = "cos(pi*x)*cos(pi*y)";
        % interface (flat or curved)
        interface = FlatCurveHorizontal(0.5);
        % viscosity ν
        nu    = 1.0;
        % diffusion tensor κ
        kappa = eye(2);
        % Beavers-Joseph-Saffman friction
        wbjs  = 1;
    end

    if isnumeric(interface)
        interface = FlatCurveHorizontal(interface);
    end

    syms x y real

    p_S = str2sym(p_stokes);
    u   = str2sym(uu_stokes);
    p_D = str2sym(p_darcy);

    u_S = [-diff(u, y); diff(u, x)];
    u_D = -kappa * gradient(p_D, [x y]);

    if isa(interface, "FlatCurveHorizontal")
        y0 = interface.y_coord(0);
        domain_measure = 1 - y0;
        p_S = p_S - 1/domain_measure * int(int(p_S, x, 0, 1), y, y0, 1);
        n_S = [0; -1];   % outward from Stokes (upper domain), pointing down
        n_D = [0;  1];   % outward from Darcy  (lower domain), pointing up
    elseif isa(interface, "SmoothCurveHorizontal")
        to_x = @(s) subs(s, symvar(s), x);
        I_sym = to_x(interface.sym_expr);
        domain_measure = double(1 - int(I_sym, x, 0, 1));
        p_S = p_S - double(1/domain_measure * int(int(p_S, y, I_sym, 1), x, 0, 1));
        n_D = [-diff(I_sym); 1] / sqrt(1 + diff(I_sym)^2);  % outward from Darcy
        n_S = -n_D;  % outward from Stokes
    else
        error("StokesDarcyCoupling.get_manufactured_solution: unsupported interface type %s", class(interface));
    end

    t_S = [0 -1; 1 0] * n_S;

    assert(divergence(u_S, [x, y]) == 0, "Velocity field is not divergence free!")

    L_S   = jacobian(u_S, [x, y]);
    sigma = nu * L_S - p_S * eye(2);

    T_1 = sigma * n_S + wbjs * (u_S' * t_S) * t_S - p_D * n_D;
    T_2 = u_S' * n_S + u_D' * n_D;

    f_S = -nu * laplacian(u_S) + gradient(p_S, [x y]);
    f_D = divergence(u_D, [x y]);

    filename = [tempname, '.txt'];
    fID = fopen(filename, 'w');
    write_field(fID, "domain_measure", [1 1], "double", domain_measure);
    write_field(fID, "nu",    [1 1], "double",          nu);
    write_field(fID, "p_S",   [1 1], "function_handle", p_S);
    write_field(fID, "u_S",   [2 1], "function_handle", u_S);
    write_field(fID, "L_S",   [2 2], "function_handle", L_S);
    write_field(fID, "f_S",   [2 1], "function_handle", f_S);
    write_field(fID, "kappa", [2 2], "double",          kappa);
    write_field(fID, "p_D",   [1 1], "function_handle", p_D);
    write_field(fID, "u_D",   [2 1], "function_handle", u_D);
    write_field(fID, "f_D",   [1 1], "function_handle", f_D);
    write_field(fID, "wbjs",  [1 1], "double",          wbjs);
    write_field(fID, "T_1",   [2 1], "function_handle", T_1);
    write_field(fID, "T_2",   [1 1], "function_handle", T_2);
    fclose(fID);

    data = rmfield(txt2mat(filename), "filename");
    delete(filename);
    delete(regexprep(filename, '\.txt$', '.mat'));
end
