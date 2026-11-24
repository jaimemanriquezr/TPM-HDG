function u_plus_x = plus(u, x)
%PLUS  Sum of a field and another field or function.
%
%   u + x adds the dof when both are FEFields (x with dof); any other
%   combination is delegated to AbstractFEFunction/plus.
    arguments
        u FEField
        x
    end
    if ~isa(u, 'FEField')
        [u, x] = deal(x, u);
    end
    
    if isa(x, 'FEField') && ~isempty(x.dof)
        u_plus_x = u;
        u_plus_x.dof = u.dof + x.dof;
    else
        u_plus_x = plus@AbstractFEFunction(u, x);
    end
end