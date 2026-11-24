function mustBeBoundaryWeight(x)
%MUSTBEBOUNDARYWEIGHT  Validate that a weight has a dimension of size 3 (one entry per local face).
    if ~(ismember(3, size(x)))
        eid = 'InputValidation:notBoundaryFunction';
        msg = 'Weights must be either of size 3-by-* or *-by-*-by-3-by-*';
        error(eid, msg);
    end
end