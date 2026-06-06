function out = RabiH_DA_diagnosticsBlock(in)
% Compute current diagnostics from records containing currents/chiralities.

if isempty(in)
    out = struct([]);
    return;
end

tmpl = local_empty_record();
out = repmat(tmpl, numel(in), 1);

for k = 1:numel(in)
    r = in(k);
    e = tmpl;

    e.pointLabel = local_getstr_any(r, {'pointLabel','label'});
    e.description = local_getstr_any(r, {'pointDescription','description'});
    e.lambda = local_getnum(r,'lambda');
    e.thetaOverPi = local_getnum(r,'thetaOverPi');
    e.D = local_getnum(r,'D');
    e.chi = local_getnum(r,'chi');
    e.dph = local_getnum(r,'dph');
    e.validObs = local_getbool_any(r, {'validObs','ok'});
    if ~isfield(r,'validObs') && ~isfield(r,'ok')
        e.validObs = true;
    end

    e.energy = local_getnum(r,'energy');
    e.energyPerSite = local_getnum(r,'energyPerSite');
    e.meanPhotonNumber = local_getnum_any(r, {'meanPhotonNumber','photonNumber','nbar'});
    e.rhoSR = local_getnum_any(r, {'rhoSR','superradiantDensity'});

    e.j12 = local_getnum_any(r, {'j12','I12'});
    e.j23 = local_getnum_any(r, {'j23','I23'});
    e.j34 = local_getnum_any(r, {'j34','I34'});
    e.j41 = local_getnum_any(r, {'j41','I41'});
    e.j13 = local_getnum_any(r, {'j13','I13'});
    e.j24 = local_getnum_any(r, {'j24','I24'});

    e.chiSquare = local_getnum_any(r, {'chiSquare','chiSq'});
    e.chiTri123 = local_getnum(r,'chiTri123');
    e.chiTri134 = local_getnum(r,'chiTri134');
    e.chiTri124 = local_getnum(r,'chiTri124');
    e.chiTri234 = local_getnum(r,'chiTri234');

    if ~isfinite(e.chiSquare) && all(isfinite([e.j12 e.j23 e.j34 e.j41]))
        e.chiSquare = e.j12 + e.j23 + e.j34 + e.j41;
    end

    if ~all(isfinite([e.chiTri123 e.chiTri134 e.chiTri124 e.chiTri234]))
        if all(isfinite([e.j12 e.j23 e.j34 e.j41 e.j13 e.j24]))
            e.chiTri123 = e.j12 + e.j23 - e.j13;
            e.chiTri134 = e.j13 + e.j34 + e.j41;
            e.chiTri124 = e.j12 + e.j24 + e.j41;
            e.chiTri234 = e.j23 + e.j34 - e.j24;
        end
    end

    if all(isfinite([e.chiTri123 e.chiTri134 e.chiTri124 e.chiTri234]))
        e.meanAbsChiTri = (abs(e.chiTri123)+abs(e.chiTri134)+abs(e.chiTri124)+abs(e.chiTri234))/4;
    end

    e.Ttri = 2*e.meanAbsChiTri - abs(e.chiSquare);
    e.Ctri = 2*e.meanAbsChiTri + abs(e.chiSquare);
    if isfinite(e.Ctri) && e.Ctri > 0
        e.Qtri = e.Ttri/e.Ctri;
        e.QtriPlus = max(e.Qtri, 0);
    end

    e.Cdiag = abs(e.j13) + abs(e.j24);
    e.SdiagNumerator = abs(e.j13) - abs(e.j24);
    e.AdiagNumerator = e.j13 - e.j24;
    if isfinite(e.Cdiag) && e.Cdiag > 0
        e.Sdiag = e.SdiagNumerator/e.Cdiag;
        e.Adiag = e.AdiagNumerator/e.Cdiag;
    end

    e.PtriNumerator = abs(e.chiTri123-e.chiTri134) + abs(e.chiTri124-e.chiTri234);
    e.PtriDenominator = abs(e.chiTri123)+abs(e.chiTri134)+abs(e.chiTri124)+abs(e.chiTri234);
    if isfinite(e.PtriDenominator) && e.PtriDenominator > 0
        e.Ptri = e.PtriNumerator/e.PtriDenominator;
    end

    e.loopError = local_loop_error(e);
    e.sectorLabel = local_sector_label(e.lambda, e.thetaOverPi);

    out(k) = e;
end
end

function e = local_empty_record()
e = struct();
e.pointLabel = '';
e.description = '';
e.source = '';
e.lambda = NaN; e.thetaOverPi = NaN;
e.D = NaN; e.chi = NaN; e.dph = NaN;
e.validObs = false;
e.energy = NaN; e.energyPerSite = NaN;
e.meanPhotonNumber = NaN; e.rhoSR = NaN;
e.j12 = NaN; e.j23 = NaN; e.j34 = NaN; e.j41 = NaN; e.j13 = NaN; e.j24 = NaN;
e.chiSquare = NaN;
e.chiTri123 = NaN; e.chiTri134 = NaN; e.chiTri124 = NaN; e.chiTri234 = NaN;
e.meanAbsChiTri = NaN;
e.Ttri = NaN; e.Ctri = NaN; e.Qtri = NaN; e.QtriPlus = NaN;
e.Cdiag = NaN; e.SdiagNumerator = NaN; e.Sdiag = NaN; e.AdiagNumerator = NaN; e.Adiag = NaN;
e.PtriNumerator = NaN; e.PtriDenominator = NaN; e.Ptri = NaN;
e.loopError = NaN;
e.sectorLabel = '';
end

function err = local_loop_error(e)
err = NaN;
if all(isfinite([e.chiTri123 e.chiTri134 e.chiTri124 e.chiTri234 e.chiSquare]))
    floorVal = 1e-12;
    denom = max(abs(e.chiSquare), floorVal);
    err1 = abs(e.chiTri123 + e.chiTri134 - e.chiSquare)/denom;
    err2 = abs(e.chiTri124 + e.chiTri234 - e.chiSquare)/denom;
    err = max(err1, err2);
end
end

function label = local_sector_label(lambda, thetaOverPi)
if ~isfinite(lambda) || ~isfinite(thetaOverPi)
    label = '';
elseif lambda < 0 && abs(thetaOverPi) > 1e-12
    label = 'mixed';
elseif lambda < 0
    label = 'triangular-enhanced';
elseif abs(thetaOverPi) > 1e-12
    label = 'diagonal-selected';
else
    label = 'square-dominant';
end
end

function x = local_getnum(s,name)
x = NaN;
if isstruct(s) && isfield(s,name)
    v = s.(name);
    if isempty(v), return; end
    if isnumeric(v) || islogical(v)
        v = v(:);
        if ~isempty(v), x = double(v(1)); end
    elseif ischar(v)
        x = str2double(v);
    end
end
end

function x = local_getnum_any(s,names)
x = NaN;
for i = 1:numel(names)
    x = local_getnum(s,names{i});
    if isfinite(x), return; end
end
end

function b = local_getbool_any(s,names)
b = false;
for i = 1:numel(names)
    name = names{i};
    if isstruct(s) && isfield(s,name)
        v = s.(name);
        if isempty(v), continue; end
        if isnumeric(v) || islogical(v)
            v = v(:); b = logical(v(1)); return;
        elseif ischar(v)
            b = strcmpi(v,'true') || strcmpi(v,'1') || strcmpi(v,'ok'); return;
        end
    end
end
end

function txt = local_getstr_any(s,names)
txt = '';
for i = 1:numel(names)
    name = names{i};
    if isstruct(s) && isfield(s,name)
        v = s.(name);
        if isempty(v), continue; end
        if ischar(v), txt = v; return; end
        if iscell(v) && ~isempty(v), txt = char(v{1}); return; end
    end
end
end
