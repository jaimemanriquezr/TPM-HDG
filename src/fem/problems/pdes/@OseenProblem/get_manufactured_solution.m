function data = get_manufactured_solution(U, p, kinematic_viscosity, beta, filename)
%-------------------------------------------------------------------------
%Description: Manufactured source and boundary data for the Oseen problem
%             from a stream function U, a pressure p and a convective
%             velocity beta:
%                 u = (-∂U/∂y, ∂U/∂x),  L = ∇u,  σ = ν L - p I
%                 f = -∇·σ + (β·∇)u
%             beta = [] takes β = u, i.e. the manufactured Navier–Stokes
%             solution. Returns data with fields nu, L, u, p, f, beta.
%             p is shifted to zero mean over the unit square [0,1]^2;
%             on another domain, subtract its mean there.
%-------------------------------------------------------------------------
    arguments (Input)
        U                   = "x^2*y^2";
        p                   = "x";
        kinematic_viscosity = 1.0;
        beta                = [];        % ["bx", "by"] symbolic strings; [] -> β = u
        filename            = cache_dir() + "/oseen_solutions.txt";
    end
    arguments (Output)
        data    struct
    end
    nu = kinematic_viscosity;

    syms x y real
    u = [-diff(str2sym(U), y); diff(str2sym(U), x)];
    p = str2sym(p) - double(int(int(str2sym(p), x, 0, 1), y, 0, 1));
    if isempty(beta)
        b = u;
    else
        b = [str2sym(beta(1)); str2sym(beta(2))];
    end

    L     = jacobian(u, [x y]);
    sigma = nu*L - p*eye(2);
    f     = [-divergence(sigma(1,:), [x y]); -divergence(sigma(2,:), [x y])] + L*b;

    if nargout > 0
        % private scratch file, removed below: concurrent calls do not collide
        filename = string(tempname) + ".txt";
    end
    fID = fopen(filename, 'w');
    write_field(fID, "nu",   [1 1], "double",          nu);
    write_field(fID, "L",    [2 2], "function_handle", L);
    write_field(fID, "u",    [2 1], "function_handle", u);
    write_field(fID, "p",    [1 1], "function_handle", p);
    write_field(fID, "f",    [2 1], "function_handle", f);
    write_field(fID, "beta", [2 1], "function_handle", b);
    fclose(fID);

    data = rmfield(txt2mat(filename), "filename");
    delete(filename);
    if nargout > 0
        delete(regexprep(filename, '\.txt$', '.mat'));
    end
end
