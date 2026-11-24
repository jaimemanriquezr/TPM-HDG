function data = get_manufactured_solution(U, p, kinematic_viscosity, filename)
%-------------------------------------------------------------------------
%Description: Derives the manufactured source and boundary data for the
%             Stokes problem from a stream function U and pressure p:
%                 u = (-∂U/∂y, ∂U/∂x)    (divergence-free by construction)
%                 L = ∇u,   σ = ν*L - p*I
%                 f = -∇·σ
%             Writes results to a text file; returns data struct if requested.
%             p is shifted to zero mean over the unit square [0,1]^2;
%             on another domain, subtract its mean there.
%-------------------------------------------------------------------------
    arguments (Input)
        U                   = "x^2*y^2";           % stream function; u = (-∂U/∂y, ∂U/∂x)
        p                   = "x";                  % pressure; mean-zero normalisation applied
        kinematic_viscosity = 1.0;                  % ν
        filename            = cache_dir() + "/stokes_solutions.txt";
    end
    arguments (Output)
        data    struct  % solution data struct with fields: nu, L, u, p, f
    end
    nu = kinematic_viscosity;

    syms x y real
    u = [-diff(str2sym(U), y); diff(str2sym(U), x)];
    p = str2sym(p) - double(int(int(str2sym(p), x, 0, 1), y, 0, 1));

    L     = jacobian(u, [x y]);
    sigma = nu*L - p*eye(2);
    f     = [-divergence(sigma(1,:), [x y]); -divergence(sigma(2,:), [x y])];

    if nargout > 0
        % private scratch file, removed below: concurrent calls do not collide
        filename = string(tempname) + ".txt";
    end
    fID = fopen(filename, 'w');
    write_field(fID, "nu", [1 1], "double",          nu);
    write_field(fID, "L",  [2 2], "function_handle", L);
    write_field(fID, "u",  [2 1], "function_handle", u);
    write_field(fID, "p",  [1 1], "function_handle", p);
    write_field(fID, "f",  [2 1], "function_handle", f);
    fclose(fID);

    if nargout > 0
        data = rmfield(txt2mat(filename), "filename");
        delete(filename);
        delete(regexprep(filename, '\.txt$', '.mat'));
    else
        txt2mat(filename);
    end
end
