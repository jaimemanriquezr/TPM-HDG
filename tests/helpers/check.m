function n = check(ok, what)
%CHECK  Count a failed check: print FAILED to stderr when ~ok, return 1 or 0.
n = ~ok;
if ~ok, fprintf(2, "FAILED: %s\n", what); end
end
