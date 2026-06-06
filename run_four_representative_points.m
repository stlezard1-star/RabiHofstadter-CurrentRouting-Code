%% run_four_representative_points.m
% Entry script for the data-availability package.
% It calls RabiH_DA_mainFunc to compute four representative points.

clear; clc; close all;

cfg = struct();

% Run switches.
cfg.runSemiclassical = true;
cfg.runIPEPS = true;

% Common physical parameters.
cfg.common.omega = 1.0;
cfg.common.Delta = 1.0;
cfg.common.g = 1.40;
cfg.common.J1 = 0.16;
cfg.common.Phi = pi/2;

% Four representative points: pair shown in figures is (thetaOverPi, lambda).
cfg.points = [ ...
    struct('label','square_dominant',      'description','square-dominant representative point',      'lambda',+0.60,'thetaOverPi',+0.00), ...
    struct('label','triangular_enhanced',  'description','triangular-enhanced representative point',  'lambda',-0.60,'thetaOverPi',+0.00), ...
    struct('label','diagonal_selected',    'description','diagonal-selected representative point',    'lambda',+0.60,'thetaOverPi',+0.25), ...
    struct('label','mixed',                'description','mixed routing representative point',        'lambda',-0.60,'thetaOverPi',+0.25) ...
    ];

% Semiclassical controls.
cfg.semiclassical.nStarts = 40;
cfg.semiclassical.randomSeed = 1;
cfg.semiclassical.maxIter = 200000;
cfg.semiclassical.maxFunEvals = 200000;

% iPEPS finite-resource controls. Change these to test convergence settings.
cfg.iPEPS.DList   = 2;
cfg.iPEPS.chiForD = 8;
cfg.iPEPS.dph     = 3;

cfg.iPEPS.tauList = [0.05 0.025 0.0125];
cfg.iPEPS.nSweeps = 2;
cfg.iPEPS.maxIter = 8;
cfg.iPEPS.alphaSeedList = [0 0.8 1.2];
cfg.iPEPS.alphaPatternList = {'uniform'};
cfg.iPEPS.useCache = true;
cfg.iPEPS.forceRecompute = false;

cfg.outputDir = fullfile(pwd, 'results_data_availability_four_points');

results = RabiH_DA_mainFunc(cfg);

disp('Data-availability run finished.');
disp(['Output directory: ' results.outputDir]);
