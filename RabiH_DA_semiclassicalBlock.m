function result = RabiH_DA_semiclassicalBlock(common, points, semiCfg, outDir)
% Integrated semiclassical four-site coherent-state solver.

if exist(outDir,'dir') ~= 7
    mkdir(outDir);
end

if ~isfield(semiCfg,'nStarts'), semiCfg.nStarts = 32; end
if ~isfield(semiCfg,'randomSeed'), semiCfg.randomSeed = 1; end
if ~isfield(semiCfg,'maxIter'), semiCfg.maxIter = 200000; end
if ~isfield(semiCfg,'maxFunEvals'), semiCfg.maxFunEvals = 200000; end

rng(semiCfg.randomSeed);

raw = repmat(local_empty_raw(), numel(points), 1);
for k = 1:numel(points)
    fprintf('[semiclassical] %s: lambda=%+.3f theta/pi=%+.3f\n', ...
        points(k).label, points(k).lambda, points(k).thetaOverPi);
    raw(k) = local_solve_one(common, points(k), semiCfg);
end

diagnostics = RabiH_DA_diagnosticsBlock(raw);

save(fullfile(outDir,'semiclassical_results.mat'), 'common', 'points', 'semiCfg', 'raw', 'diagnostics');
writetable(struct2table(diagnostics), fullfile(outDir,'semiclassical_results.csv'));

result = struct();
result.raw = raw;
result.diagnostics = diagnostics;
result.outputDir = outDir;
end

%% ========================================================================
function r = local_solve_one(common, point, semiCfg)

ph = local_phases(common.Phi, point.thetaOverPi*pi);
J2 = point.lambda*common.J1;

starts = local_initial_conditions(semiCfg.nStarts);
bestE = inf;
bestX = [];

opts = optimset('Display','off', ...
    'MaxFunEvals', semiCfg.maxFunEvals, ...
    'MaxIter', semiCfg.maxIter, ...
    'TolX',1e-11, ...
    'TolFun',1e-11);

for s = 1:size(starts,1)
    x0 = starts(s,:);
    [xopt, eopt] = fminsearch(@(x)local_energy(x,common,ph,J2), x0, opts);
    if eopt < bestE
        bestE = eopt;
        bestX = xopt;
    end
end

a = bestX(1:4) + 1i*bestX(5:8);

I12 = local_current(a,1,2,common.J1,ph.A12);
I23 = local_current(a,2,3,common.J1,ph.A23);
I34 = local_current(a,3,4,common.J1,ph.A34);
I41 = local_current(a,4,1,common.J1,ph.A41);
I13 = local_current(a,1,3,J2,ph.B13);
I24 = local_current(a,2,4,J2,ph.B24);

r = local_empty_raw();
r.pointLabel = point.label;
r.pointDescription = point.description;
r.lambda = point.lambda;
r.thetaOverPi = point.thetaOverPi;
r.validObs = true;
r.energy = bestE;
r.energyPerSite = bestE/4;
r.meanPhotonNumber = mean(abs(a).^2);
r.rhoSR = mean(abs(a).^2);

r.alpha1 = a(1); r.alpha2 = a(2); r.alpha3 = a(3); r.alpha4 = a(4);
r.j12 = I12; r.j23 = I23; r.j34 = I34; r.j41 = I41; r.j13 = I13; r.j24 = I24;
r.chiSquare = I12 + I23 + I34 + I41;
r.chiTri123 = I12 + I23 - I13;
r.chiTri134 = I13 + I34 + I41;
r.chiTri124 = I12 + I24 + I41;
r.chiTri234 = I23 + I34 - I24;
end

function E = local_energy(x,common,ph,J2)
a = x(1:4) + 1i*x(5:8);
E = 0;
for i = 1:4
    xi = real(a(i));
    E = E + common.omega*abs(a(i))^2 - 0.5*sqrt(common.Delta^2 + 16*common.g^2*xi^2);
end
E = E + local_bond_energy(a,1,2,common.J1,ph.A12);
E = E + local_bond_energy(a,2,3,common.J1,ph.A23);
E = E + local_bond_energy(a,3,4,common.J1,ph.A34);
E = E + local_bond_energy(a,4,1,common.J1,ph.A41);
E = E + local_bond_energy(a,1,3,J2,ph.B13);
E = E + local_bond_energy(a,2,4,J2,ph.B24);
end

function Eb = local_bond_energy(a,i,j,J,phi)
Eb = -2*J*real(exp(1i*phi)*conj(a(i))*a(j));
end

function I = local_current(a,i,j,J,phi)
I = 2*J*imag(exp(1i*phi)*conj(a(i))*a(j));
end

function ph = local_phases(Phi, theta)
% Four-site convention used for the semiclassical current-network guide.
ph = struct();
ph.A12 = 0;
ph.A23 = Phi;
ph.A34 = 0;
ph.A41 = 0;
ph.B13 = Phi/2 + theta;
ph.B24 = Phi/2 - theta;
end

function X = local_initial_conditions(nStarts)
baseAmp = [0 0.4 0.8 1.2];
patterns = cell(0,1);
for k = 1:numel(baseAmp)
    a = baseAmp(k);
    patterns{end+1} = [ a a a a 0 0 0 0]; %#ok<AGROW>
    patterns{end+1} = [-a a -a a 0 0 0 0]; %#ok<AGROW>
    patterns{end+1} = [ a -a a -a 0 0 0 0]; %#ok<AGROW>
    patterns{end+1} = [ a a -a -a 0 0 0 0]; %#ok<AGROW>
end
n = max(nStarts, numel(patterns));
X = zeros(n,8);
for k = 1:numel(patterns)
    X(k,:) = patterns{k};
end
for k = numel(patterns)+1:n
    X(k,:) = 0.5*randn(1,8);
end
end

function r = local_empty_raw()
r = struct();
r.pointLabel = '';
r.pointDescription = '';
r.lambda = NaN;
r.thetaOverPi = NaN;
r.validObs = false;
r.energy = NaN;
r.energyPerSite = NaN;
r.meanPhotonNumber = NaN;
r.rhoSR = NaN;
r.alpha1 = NaN; r.alpha2 = NaN; r.alpha3 = NaN; r.alpha4 = NaN;
r.j12 = NaN; r.j23 = NaN; r.j34 = NaN; r.j41 = NaN; r.j13 = NaN; r.j24 = NaN;
r.chiSquare = NaN;
r.chiTri123 = NaN; r.chiTri134 = NaN; r.chiTri124 = NaN; r.chiTri234 = NaN;
end
