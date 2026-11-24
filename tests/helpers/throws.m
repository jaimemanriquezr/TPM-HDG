function tf = throws(f, id)
%THROWS  True when f() raises an error with identifier id.
% match the identifier: catching any exception would also pass if the call
% failed for an unrelated reason, e.g. a renamed option
try
    f();  tf = false;
catch err
    tf = err.identifier == id;
end
end
