function minus_A = uminus(A)
%UMINUS  Negation of every block.
    minus_A = IndexedMatrix(cellfun(@uminus, A.Arrays, "UniformOutput", false));
end