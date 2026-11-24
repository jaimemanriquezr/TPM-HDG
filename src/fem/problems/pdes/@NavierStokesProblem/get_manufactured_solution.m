function data = get_manufactured_solution(U, p, shear_viscosity, beta, rho, time, filename)
%-------------------------------------------------------------------------
%Description: Manufactured data for the stress-velocity Oseen/Navier–
%             Stokes problem from a stream function U, pressure p and
%             convective velocity beta:
%                 u = (-∂U/∂y, ∂U/∂x),  σ = 2μ ε(u) - p I,
%                 f = ρ u_t + ρ_c (β·∇)u - ∇·σ,   ρ_c = ρ (ρ > 0) or 1
%             beta = [] takes β = u (manufactured Navier–Stokes).
%             Returns data with fields viscosity, sigma, u, gamma, f, beta.
%             p is shifted to zero mean over the unit square [0,1]^2;
%             on another domain, subtract its mean there.
%-------------------------------------------------------------------------
    arguments (Input)
        U               = "x^2*y^2";
        p               = "x";
        shear_viscosity = 0.5;
        beta            = [];       % ["bx", "by"] symbolic strings; [] -> β = u
        rho             = 0.0;
        time            = 0.0;
        filename        = cache_dir() + "/stress_navier_stokes_solutions.txt";
    end
    arguments (Output)
        data    struct
    end
    mu = shear_viscosity;
    rho_c = 1;
    if rho > 0; rho_c = rho; end

    syms x y t real
    u = [-diff(str2sym(U), y); diff(str2sym(U), x)];
    % the mean over the unit square; kept symbolic when p depends on t
    p = str2sym(p);
    p_mean = int(int(p, x, 0, 1), y, 0, 1);
    if isempty(symvar(p_mean)), p_mean = double(p_mean); end
    p = p - p_mean;
    if isempty(beta)
        b = u;
    else
        b = [str2sym(beta(1)); str2sym(beta(2))];
    end

    epsilon     = .5*(jacobian(u, [x y]) + jacobian(u, [x y]).');
    sigma       = 2*mu*epsilon - p*eye(2);
    gamma       = -.5*(-diff(u(1), y) + diff(u(2), x));
    minus_div   = @(w) [-divergence(w(1,:), [x y]); -divergence(w(2,:), [x y])];
    f           = rho*diff(u, t) + rho_c*jacobian(u, [x y])*b + minus_div(sigma);

    sigma = subs(sigma, t, time);
    gamma = subs(gamma, t, time);
    u_t   = subs(diff(u, t), t, time);
    u     = subs(u,     t, time);
    f     = subs(f,     t, time);
    p     = subs(p,     t, time);
    b     = subs(b,     t, time);

    if nargout > 0
        % private scratch file, removed below: concurrent calls do not collide
        filename = string(tempname) + ".txt";
    end
    fID = fopen(filename, 'w');
    write_field(fID, "viscosity", [1 1], "double",          mu);
    write_field(fID, "sigma",     [2 2], "function_handle", sigma);
    write_field(fID, "u",         [2 1], "function_handle", u);
    write_field(fID, "gamma",     [1 1], "function_handle", gamma);
    write_field(fID, "f",         [2 1], "function_handle", f);
    write_field(fID, "beta",      [2 1], "function_handle", b);
    write_field(fID, "u_t",       [2 1], "function_handle", u_t);
    write_field(fID, "p",         [1 1], "function_handle", p);
    fclose(fID);

    data = rmfield(txt2mat(filename), "filename");
    delete(filename);
    if nargout > 0
        delete(regexprep(filename, '\.txt$', '.mat'));
    end
end
