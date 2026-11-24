function initpath()
filename = mfilename("fullpath");
call_path = fileparts(filename);
addsubpaths = @(path) addpath(genpath(call_path + "/" + path));

% tests/ is on the path so the verification scripts run by name;
% experiments/ exists only in the development repository
folders = ["src", "examples", "experiments", "tests"];
for path = folders
    if isfolder(call_path + "/" + path)
        addsubpaths(path)
    end
end
end
