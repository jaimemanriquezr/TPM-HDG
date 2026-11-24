function u = uminus(u)
%UMINUS  Negation of a field: -u has dof -u.dof.
    arguments
        u FEField
    end
    u.dof = -u.dof;
end