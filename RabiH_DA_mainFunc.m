function results = RabiH_DA_mainFunc(cfg)
% Top-level organizer for the Rabi--Hofstadter data-availability calculation.

if nargin < 1 || isempty(cfg)
    error('A cfg structure is required. Use run_four_representative_points.m as the entry script.');
end

if ~isfield(cfg,'outputDir') || isempty(cfg.outputDir)
    cfg.outputDir = fullfile(pwd, 'results_data_availability_four_points');
end
if exist(cfg.outputDir,'dir') ~= 7
    mkdir(cfg.outputDir);
end

results = struct();
results.outputDir = cfg.outputDir;
results.cfg = cfg;

if isfield(cfg,'runSemiclassical') && cfg.runSemiclassical
    fprintf('\n=== Semiclassical block ===\n');
    semiOutDir = fullfile(cfg.outputDir, 'semiclassical');
    semi = RabiH_DA_semiclassicalBlock(cfg.common, cfg.points, cfg.semiclassical, semiOutDir);
    results.semiclassical = semi;
else
    semi = [];
end

if isfield(cfg,'runIPEPS') && cfg.runIPEPS
    fprintf('\n=== iPEPS block ===\n');
    ipepsOutDir = fullfile(cfg.outputDir, 'ipeps');
    ipeps = RabiH_DA_iPEPSBlock(cfg.common, cfg.points, cfg.iPEPS, ipepsOutDir);
    results.iPEPS = ipeps;
else
    ipeps = [];
end

combined = struct([]);
if ~isempty(semi) && isfield(semi,'diagnostics')
    combined = local_append_with_source(combined, semi.diagnostics, 'semiclassical');
end
if ~isempty(ipeps) && isfield(ipeps,'diagnostics')
    combined = local_append_with_source(combined, ipeps.diagnostics, 'iPEPS');
end

results.combinedDiagnostics = combined;

save(fullfile(cfg.outputDir,'combined_diagnostics.mat'), 'results', 'combined', 'cfg');
if ~isempty(combined)
    T = struct2table(combined);
    writetable(T, fullfile(cfg.outputDir,'combined_diagnostics.csv'));
end
end

function combined = local_append_with_source(combined, arr, sourceName)
if isempty(arr)
    return;
end
for k = 1:numel(arr)
    x = arr(k);
    x.source = sourceName;
    if isempty(combined)
        combined = x;
    else
        combined(end+1) = x; %#ok<AGROW>
    end
end
end
