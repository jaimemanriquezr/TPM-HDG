function plan = buildfile
% BUILDFILE  Build tasks of the TPM-HDG toolbox (run with buildtool).
%
%   buildtool folders   update the member folders of resources/mpackage.json
%                       to the current src/, examples/ and tests/ trees
%   buildtool check     fail if the member folders or the Contents.m version
%                       disagree with resources/mpackage.json
%   buildtool package   check, then write release/TPM-HDG.mltbx
%
% resources/mpackage.json is the single source of the name, ID, version,
% summary and author: the same repository installs as an MPM package
% (mpminstall(folder)) and as the .mltbx add-on.
plan = buildplan(localfunctions);
plan("package").Dependencies = "check";
plan.DefaultTasks = "package";
end

function foldersTask(context)
% Update the member folders of resources/mpackage.json
% Only the "folders" array is rewritten: a full jsondecode/jsonencode round
% trip, and saving through matlab.mpm.Package, both turn one-element
% arrays (supportedPlatforms) into objects, which ToolboxOptions rejects.
root = context.Plan.RootFolder;
entries = compose('    {"path": "%s", "languages": ["matlab"]}', member_folders(root));
block = sprintf('"folders": [\n%s\n  ],\n  "dependencies"', ...
    char(strjoin(entries, sprintf(',\n'))));
text = fileread(spec_file(root));
[i, j] = regexp(text, '"folders"\s*:\s*\[.*?\]\s*,\s*"dependencies"', 'once');
assert(~isempty(i), "buildfile:folders", "No folders array in mpackage.json.");
writelines(string([text(1:i-1), block, text(j+1:end)]), spec_file(root), WriteMode="overwrite");
end

function checkTask(context)
% Check resources/mpackage.json against the tree and Contents.m
root = context.Plan.RootFolder;
spec = read_spec(root);
listed = string({spec.folders.path});
expected = member_folders(root);
missing = setdiff(expected, listed);
stale = setdiff(listed, expected);
assert(isempty(missing) && isempty(stale), "buildfile:folders", ...
    "mpackage.json folders are out of date (missing: %s; stale: %s). " + ...
    "Run buildtool folders.", join(missing, ", "), join(stale, ", "));
contents = fileread(fullfile(root, "Contents.m"));
assert(contains(contents, "% Version " + spec.version + " "), "buildfile:version", ...
    "Contents.m does not state Version %s.", spec.version);
end

function packageTask(context)
% Write release/TPM-HDG.mltbx
root = context.Plan.RootFolder;
spec = read_spec(root);
opts = matlab.addons.toolbox.ToolboxOptions(root, spec.id);
opts.ToolboxName = spec.displayName;
opts.ToolboxVersion = spec.version;
opts.Summary = spec.summary;
opts.AuthorName = spec.provider.name;
opts.MinimumMatlabRelease = "R2024b";
% resources/project (the MATLAB project's files) is development-only
opts.ToolboxFiles = fullfile(root, ["src", "examples", "tests", "resources/mpackage.json", ...
    "Contents.m", "initpath.m", "README.md", "LICENSE"]);   % TODO.md stays in the repository
opts.ToolboxMatlabPath = fullfile(root, member_folders(root));
opts.ProductDependencies = ["Symbolic Math Toolbox", "Parallel Computing Toolbox"];
if ~isfolder(fullfile(root, "release"))
    mkdir(fullfile(root, "release"));
end
opts.OutputFile = fullfile(root, "release", spec.displayName + ".mltbx");
matlab.addons.toolbox.packageToolbox(opts);
end

function folders = member_folders(root)
% The folders initpath puts on the path: src/, examples/ and tests/ and
% their subfolders, without class (@), namespace (+) and private folders,
% relative to root with forward slashes
folders = strings(1, 0);
for top = ["src", "examples", "tests"]
    paths = split(string(genpath(fullfile(root, top))), pathsep).';
    paths = paths(paths ~= "");
    folders = [folders, replace(extractAfter(paths, strlength(root) + 1), "\", "/")]; %#ok<AGROW>
end
folders = sort(folders);
end

function spec = read_spec(root)
spec = jsondecode(fileread(spec_file(root)));
end

function f = spec_file(root)
f = fullfile(root, "resources", "mpackage.json");
end
