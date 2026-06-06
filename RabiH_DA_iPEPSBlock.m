function result = RabiH_DA_iPEPSBlock(common, points, ipepsCfg, outDir)
% Integrated iPEPS finite-resource block.
% The compact RHDC backend is included as backend_RHDC_compact.zip and is unpacked automatically.

if exist(outDir,'dir') ~= 7
    mkdir(outDir);
end

local_prepare_backend();

cfg = RHDC_DefaultConfig();

% Physical parameters.
cfg.g0  = common.g;
cfg.J1  = common.J1;
cfg.Phi = common.Phi;

% Resource controls.
cfg.DList   = ipepsCfg.DList;
cfg.chiForD = ipepsCfg.chiForD;
cfg.dph     = ipepsCfg.dph;

% Manuscript safe route.
cfg.observableMode = 'splitFactoredCTM_coherentNNN_SAFE';
cfg.diagonalEstimator = 'coherent';
cfg.useExactTorus = false;
cfg.allowPlaquetteCTM = false;

cfg.tauList = ipepsCfg.tauList;
cfg.nSweeps = ipepsCfg.nSweeps;
cfg.maxIter = ipepsCfg.maxIter;
cfg.alphaSeedList = ipepsCfg.alphaSeedList;
cfg.alphaPatternList = ipepsCfg.alphaPatternList;

cfg.useCache = ipepsCfg.useCache;
cfg.forceRecompute = ipepsCfg.forceRecompute;
cfg.saveUnitCell = false;
cfg.outputDir = fullfile(outDir, 'raw_RHDC');

if exist(cfg.outputDir,'dir') ~= 7
    mkdir(cfg.outputDir);
end

fprintf('[iPEPS] DList=%s chiForD=%s dph=%d\n', mat2str(cfg.DList), mat2str(cfg.chiForD), cfg.dph);

recordsRaw = RHDC_RunJobs(cfg, points);
diagnostics = RabiH_DA_diagnosticsBlock(recordsRaw);

save(fullfile(outDir,'iPEPS_results.mat'), 'common', 'points', 'ipepsCfg', 'cfg', 'recordsRaw', 'diagnostics');
writetable(struct2table(diagnostics), fullfile(outDir,'iPEPS_results.csv'));

result = struct();
result.recordsRaw = recordsRaw;
result.diagnostics = diagnostics;
result.outputDir = outDir;
end

%% ========================================================================
function local_prepare_backend()
thisFile = mfilename('fullpath');
thisDir = fileparts(thisFile);
backendDir = fullfile(thisDir, '.rhdc_backend');
codeBaseDir = fullfile(backendDir, 'RabiH_FiniteResourceConvergence_Fig', 'code_base');

if exist(codeBaseDir,'dir') ~= 7
    zipFile = fullfile(thisDir, 'backend_RHDC_compact.zip');
    if exist(zipFile,'file') ~= 2
        error('backend_RHDC_compact.zip not found in the package folder.');
    end
    if exist(backendDir,'dir') ~= 7
        mkdir(backendDir);
    end
    fprintf('[iPEPS] unpacking compact RHDC backend...\n');
    unzip(zipFile, backendDir);
end

addpath(genpath(codeBaseDir));

if exist('RHDC_DefaultConfig','file') ~= 2 || exist('RHDC_RunJobs','file') ~= 2
    error('Failed to add the compact RHDC backend to the MATLAB path.');
end
end
