function A_minus_B = minus(A, B)
%MINUS  Difference A - B, computed as plus(A, -B).
    arguments
        A IndexedMatrix
        B IndexedMatrix
    end
    A_minus_B = plus(A, -B);
end