function data = get_manufactured_solution(p, kappa, filename)
%-------------------------------------------------------------------------
%Description: Derives the manufactured source and boundary data for the
%             Darcy problem from a given exact pressure p and tensor κ:
%                 u = -κ * ∇p,    f = ∇·u
%             Writes results to a text file; returns data struct if requested.
%-------------------------------------------------------------------------
    arguments (Input)
        p        (1,1) string   = "x^2*y^2";           % symbolic expression for exact pressure
        kappa    (2,2) double   = eye(2);               % diffusion tensor κ (constant)
        filename (1,1) string   = cache_dir() + "/darcy_solutions.txt";  % output file path
    end
    arguments (Output)
        data    struct  % solution data struct with fields: kappa, u, p, f (only when nargout > 0)
    end
    syms x y real
    p = str2sym(p);
    u = -kappa * gradient(p, [x y]);
    f = divergence(u, [x y]);

    if nargout > 0
        % private scratch file, removed below: concurrent calls do not collide
        filename = string(tempname) + ".txt";
    end
    fID = fopen(filename, 'w');
    write_field(fID, "kappa", [2 2], "double", kappa);
    write_field(fID, "u", [2 1], "function_handle", u);
    write_field(fID, "p", [1 1], "function_handle", p);
    write_field(fID, "f", [1 1], "function_handle", f);
    fclose(fID);

    if nargout > 0
        data = rmfield(txt2mat(filename), "filename");
        delete(filename);
        delete(regexprep(filename, '\.txt$', '.mat'));
    else
        txt2mat(filename);
    end
end