function clearElastixWorkspaceForNewRun(elastixDir)
%CLEARELASTIXWORKSPACEFORNEWRUN Remove stale elastix outputs so inversion matches the log.
%
%   invertElastixTransformCP expects the number of TransformParameters*.txt files to
%   match the number of -p parameter files recorded in elastix.log. Leftover files from
%   a previous run break that invariant.

if nargin < 1 || isempty(elastixDir) || exist(elastixDir, 'dir') ~= 7
    return
end
patterns = {'TransformParameters.*.txt', 'elastix.log'};
for pi = 1:numel(patterns)
    L = dir(fullfile(elastixDir, patterns{pi}));
    for k = 1:numel(L)
        if ~L(k).isdir
            delete(fullfile(elastixDir, L(k).name));
        end
    end
end
end
