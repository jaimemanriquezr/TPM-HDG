function mustBeNumericOrCell(x)
%MUSTBENUMERICORCELL  Validate that the input is numeric or a cell array.
    if ~isnumeric(x) && ~isa(x, 'cell')
        eid = 'InputValidator:IsNumericOrCell';
        msg = 'Input must be numeric or cell array.';
        error(eid, msg)
    end
end