function u_minus_x = minus(u, x)
%MINUS  Difference of a field and another field or function.
%
%   u - x subtracts the dof when both are FEFields (x with dof); any
%   other combination is delegated to AbstractFEFunction/minus.
    arguments
        u FEField
        x
    end
    if ~isa(u, 'FEField')
        [u, x] = deal(x, u);
    end
    
    if isa(x, 'FEField') && ~isempty(x.dof)
        u_minus_x = u;
        u_minus_x.dof = u.dof - x.dof;
    else
        u_minus_x = minus@AbstractFEFunction(u, x);
    end
end