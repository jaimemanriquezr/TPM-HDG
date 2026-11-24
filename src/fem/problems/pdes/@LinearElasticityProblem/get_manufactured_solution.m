function data = get_manufactured_solution(u, young_modulus, poisson_ratio, hooke_modulus, density, time, filename)
%-------------------------------------------------------------------------
%Description: Derives the manufactured source and boundary data for the
%             linear elasticity problem from a given exact displacement u:
%                 ε = ½(∇u + ∇u^T),   σ = 2μ*ε + λ*tr(ε)*I
%                 f = ρ*ü - ∇·σ + b₀*u
%             Supports static (density=0) and dynamic (density>0) problems.
%             Writes results to a text file; returns data struct if requested.
%-------------------------------------------------------------------------
    arguments (Input)
        u             = ["x^2*y"; "x*y^2"];         % (2×1) symbolic string for displacement
        young_modulus = 1;                           % E
        poisson_ratio = 0.4;                         % ν  (→ λ, μ via standard conversion)
        hooke_modulus = 0.0;                         % b₀ (Winkler / Robin stiffness)
        density       = 0.0;                         % ρ  (0 → static; >0 → transient)
        time          = 0.0;                         % t  at which symbolic t is evaluated
        filename      = cache_dir() + "/elasticity_solutions.txt";
    end
    arguments (Output)
        data    struct  % solution data struct with fields: young_modulus, poisson_ratio, sigma, u, u_t, gamma, f
    end
    E   = young_modulus;
    nu  = poisson_ratio;
    b0  = hooke_modulus;
    rho = density;

    mu     = E / (2*(1 + nu));
    lambda = E*nu / ((1 + nu)*(1 - 2*nu));

    syms x y t real
    u       = str2sym(u);
    epsilon = .5*(jacobian(u, [x y]) + jacobian(u, [x y]).');
    gamma   = -.5*(-diff(u(1), y) + diff(u(2), x));
    sigma   = 2*mu*epsilon + lambda*trace(epsilon)*eye(2);

    Div = @(w) [divergence(w(1,:), [x y]); divergence(w(2,:), [x y])];
    f   = rho*diff(u, t, 2) - Div(sigma) + b0*u;
    u_t = diff(u, t);

    sigma = subs(sigma, t, time);
    gamma = subs(gamma, t, time);
    u     = subs(u,     t, time);
    u_t   = subs(u_t,   t, time);
    f     = subs(f,     t, time);

    if nargout > 0
        % private scratch file, removed below: concurrent calls do not collide
        filename = string(tempname) + ".txt";
    end
    fID = fopen(filename, 'w');
    write_field(fID, "young_modulus", [1 1], "double",          E);
    write_field(fID, "poisson_ratio", [1 1], "double",          nu);
    write_field(fID, "sigma",         [2 2], "function_handle", sigma);
    write_field(fID, "u",             [2 1], "function_handle", u);
    write_field(fID, "u_t",           [2 1], "function_handle", u_t);
    write_field(fID, "gamma",         [1 1], "function_handle", gamma);
    write_field(fID, "f",             [2 1], "function_handle", f);
    fclose(fID);

    if nargout > 0
        data = rmfield(txt2mat(filename), "filename");
        delete(filename);
        delete(regexprep(filename, '\.txt$', '.mat'));
    else
        txt2mat(filename);
    end
end
