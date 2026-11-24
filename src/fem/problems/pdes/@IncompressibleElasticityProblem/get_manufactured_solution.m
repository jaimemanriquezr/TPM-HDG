function data = get_manufactured_solution(U, p, shear_viscosity, rho, time, filename)
%-------------------------------------------------------------------------
%Description: Derives the manufactured source and boundary data for the
%             incompressible elasticity (stress-velocity) problem from a
%             stream function U and pressure p:
%                 u = (-∂U/∂y, ∂U/∂x)    (divergence-free by construction)
%                 ε = ½(∇u + ∇u^T),   σ = 2μ*ε - p*I
%                 f = ρ*∂u/∂t - ∇·σ
%             Writes results to a text file; returns data struct if requested.
%             p is shifted to zero mean over the unit square [0,1]^2;
%             on another domain, subtract its mean there.
%-------------------------------------------------------------------------
    arguments (Input)
        U               = "x^2*y^2";               % stream function; u = (-∂U/∂y, ∂U/∂x)
        p               = "x";                      % pressure; mean-zero normalisation applied
        shear_viscosity = 0.5;                      % μ
        rho             = 0.0;                       % ρ  (0 → steady; >0 → transient)
        time            = 0.0;                       % t  at which symbolic t is evaluated
        filename        = cache_dir() + "/stress_stokes_solutions.txt";
    end
    arguments (Output)
        data    struct  % solution data struct with fields: viscosity, sigma, u, gamma, f
    end
    mu = shear_viscosity;

    syms x y t real
    u = [-diff(str2sym(U), y); diff(str2sym(U), x)];
    % the mean over the unit square; kept symbolic when p depends on t
    p = str2sym(p);
    p_mean = int(int(p, x, 0, 1), y, 0, 1);
    if isempty(symvar(p_mean)), p_mean = double(p_mean); end
    p = p - p_mean;

    epsilon     = .5*(jacobian(u, [x y]) + jacobian(u, [x y]).');
    sigma       = 2*mu*epsilon - p*eye(2);
    gamma       = -.5*(-diff(u(1), y) + diff(u(2), x));
    minus_div   = @(w) [-divergence(w(1,:), [x y]); -divergence(w(2,:), [x y])];
    f           = rho*diff(u, t) + minus_div(sigma);

    sigma = subs(sigma, t, time);
    gamma = subs(gamma, t, time);
    u_t   = subs(diff(u, t), t, time);
    u     = subs(u,     t, time);
    f     = subs(f,     t, time);
    p     = subs(p,     t, time);

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
    write_field(fID, "u_t",       [2 1], "function_handle", u_t);
    write_field(fID, "p",         [1 1], "function_handle", p);
    fclose(fID);

    if nargout > 0
        data = rmfield(txt2mat(filename), "filename");
        delete(filename);
        delete(regexprep(filename, '\.txt$', '.mat'));
    else
        txt2mat(filename);
    end
end
