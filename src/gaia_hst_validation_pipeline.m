function results = gaia_hst_validation_pipeline(rootDir, varargin)
% gaia_hst_validation_pipeline
% -------------------------------------------------------------------------
% MATLAB pipeline for the paper:
%   "A Reproducible Gaia--HST Internal-Kinematics Validation Pipeline
%    for Galactic Globular Clusters"
%
% Scientific scope:
%   A calibration/validation atlas, not a discovery/Jeans-modelling paper.
%   Gaia EDR3 provides wide-field internal-kinematic profiles; HST/HACKS
%   provides a crowded-core proper-motion dispersion benchmark; APOGEE DR17
%   GC VAC provides an external spectroscopic/chemical readiness layer.
%
% Main inputs expected in rootDir/input:
%   gaia_edr3_profile_percentiles.csv
%   gaia_edr3_profile_summary.csv
%   gaia_edr3_catalogue_audit_170.csv
%   hacks_hst_dispersion_profiles.csv
%   apogee_dr17_gc_parameters.csv
%   apogee_dr17_gc_member_summary.csv
%   cluster_crosswalk_master.csv
%   analysis_config.csv
%
% Optional raw HACKS text:
%   input/raw/hacks_velocity_dispersions.txt
%
% Main outputs:
%   outputs/tables/GaiaHST_MATCHED_RADIAL_COMPARISON.csv
%   outputs/tables/GaiaHST_CLUSTER_QUALITY_CLASSES.csv
%   outputs/tables/GaiaHST_GLOBAL_CALIBRATION_SUMMARY.csv
%   outputs/reports/GaiaHST_interpretation_report.txt
% -------------------------------------------------------------------------

    tStart = tic;
    if nargin < 1 || isempty(rootDir); rootDir = pwd; end
    rootDir = char(rootDir);
    if ~isfolder(rootDir); error('Root folder not found: %s', rootDir); end

    cfg = defaultConfig();
    cfg = readConfigCSV(cfg, fullfile(rootDir,'input','analysis_config.csv'));
    cfg = parseNameValue(cfg, varargin{:});
    rng(safeSeed(cfg.randomSeed), 'twister');

    dirs = makeDirs(rootDir);
    fid = fopen(fullfile(dirs.logs,'GaiaHST_run_log.txt'),'w');
    cObj = onCleanup(@() safeClose(fid)); %#ok<NASGU>
    logMsg(fid, 'Gaia--HST--APOGEE kinematic validation pipeline started.');
    logMsg(fid, ['Root folder: ' rootDir]);

    % ----------------------- Load data ----------------------------------
    inDir = fullfile(rootDir,'input');
    gaiaProf = readRequired(fullfile(inDir,'gaia_edr3_profile_percentiles.csv'));
    gaiaSummary = readOptional(fullfile(inDir,'gaia_edr3_profile_summary.csv'));
    gaiaAudit = readOptional(fullfile(inDir,'gaia_edr3_catalogue_audit_170.csv'));
    apogeeParams = readOptional(fullfile(inDir,'apogee_dr17_gc_parameters.csv'));
    apogeeSummary = readOptional(fullfile(inDir,'apogee_dr17_gc_member_summary.csv'));
    crosswalk = readOptional(fullfile(inDir,'cluster_crosswalk_master.csv'));

    hacksFile = fullfile(inDir,'hacks_hst_dispersion_profiles.csv');
    hacks = readOptional(hacksFile);
    rawHacks = fullfile(inDir,'raw','hacks_velocity_dispersions.txt');
    if (isempty(hacks) || height(hacks)==0) && isfile(rawHacks)
        logMsg(fid, 'HACKS CSV empty; parsing raw HACKS velocity-dispersion text.');
        hacks = parseHacksRawText(rawHacks);
        writeTableSafe(fullfile(inDir,'hacks_hst_dispersion_profiles.csv'), hacks);
    end

    % Standardize IDs and columns.
    gaiaProf.cluster_id = string(gaiaProf.cluster_id);
    if ~isempty(hacks); hacks.cluster_id = string(hacks.cluster_id); end
    if ~isempty(gaiaAudit); gaiaAudit.cluster_id = string(gaiaAudit.cluster_id); end
    if ~isempty(apogeeSummary); apogeeSummary.cluster_id = string(apogeeSummary.cluster_id); end

    logMsg(fid, sprintf('Gaia profile rows=%d, clusters=%d', height(gaiaProf), numel(unique(gaiaProf.cluster_id))));
    logMsg(fid, sprintf('HST/HACKS rows=%d, clusters=%d', height(hacks), numel(unique(hacks.cluster_id))));
    logMsg(fid, sprintf('APOGEE summary rows=%d', height(apogeeSummary)));

    % -------------------- Create readiness/crosswalk ---------------------
    readiness = buildReadinessTable(gaiaSummary, gaiaAudit, hacks, apogeeSummary, crosswalk);
    writeTableSafe(fullfile(dirs.tables,'GaiaHST_DATA_LAYER_READINESS.csv'), readiness);

    if isempty(hacks) || height(hacks)==0
        logMsg(fid, 'No HST/HACKS dispersion rows available. Writing readiness outputs only.');
        reportNoHacks(dirs, readiness);
        results = struct('readiness',readiness,'matched',table(),'classes',table(),'global',table(), 'outputDir',dirs.root);
        return;
    end

    % ---------------------- Match profiles -------------------------------
    matchedAll = table(); clusterRows = table();
    clusters = intersect(unique(gaiaProf.cluster_id), unique(hacks.cluster_id));
    logMsg(fid, sprintf('Clusters with both Gaia profile and HST/HACKS profile: %d', numel(clusters)));

    for k = 1:numel(clusters)
        cid = clusters(k);
        try
            g = gaiaProf(gaiaProf.cluster_id==cid,:);
            h = hacks(hacks.cluster_id==cid,:);
            [matched, summary] = compareOneCluster(cid, g, h, gaiaAudit, apogeeSummary, cfg);
            matchedAll = [matchedAll; matched]; %#ok<AGROW>
            clusterRows = [clusterRows; summary]; %#ok<AGROW>
            logMsg(fid, sprintf('Compared %s: matched bins=%d class=%s meanDelta=%.3f', cid, height(matched), string(summary.quality_class(1)), summary.mean_delta_frac(1)));
        catch ME
            logMsg(fid, sprintf('ERROR in %s: %s', cid, ME.message));
            err = table(string(cid), string(ME.message), 'VariableNames', {'cluster_id','error_message'});
            writeTableSafe(fullfile(dirs.tables, ['ERROR_' char(cid) '.csv']), err);
        end
    end

    globalSummary = computeGlobalSummary(matchedAll, clusterRows, cfg);
    classes = clusterRows;

    writeTableSafe(fullfile(dirs.tables,'GaiaHST_MATCHED_RADIAL_COMPARISON.csv'), matchedAll);
    writeTableSafe(fullfile(dirs.tables,'GaiaHST_CLUSTER_QUALITY_CLASSES.csv'), classes);
    writeTableSafe(fullfile(dirs.tables,'GaiaHST_GLOBAL_CALIBRATION_SUMMARY.csv'), globalSummary);
    writeTableSafe(fullfile(dirs.tables,'GaiaHST_DATA_LAYER_READINESS.csv'), readiness);

    makeFigures(dirs, matchedAll, classes, readiness);
    writeInterpretationReport(dirs, matchedAll, classes, globalSummary, readiness, cfg);
    writeOutputDictionary(dirs);

    results = struct();
    results.matched = matchedAll;
    results.classes = classes;
    results.global = globalSummary;
    results.readiness = readiness;
    results.outputDir = dirs.root;
    results.elapsed_seconds = toc(tStart);
    logMsg(fid, sprintf('Pipeline completed in %.2f seconds.', results.elapsed_seconds));

    fprintf('\nGaia--HST kinematic validation completed.\nOutput folder:\n  %s\n', dirs.root);
end

% ========================================================================
% Configuration / IO
% ========================================================================
function cfg = defaultConfig()
    cfg.minMatchedBins = 3;
    cfg.minHSTStarsPerBin = 30;
    cfg.defaultHSTSystematicFloor = 0.0;
    cfg.defaultGaiaSystematicFloorGridMax = 0.050;
    cfg.defaultGaiaScaleGridMax = 2.50;
    cfg.bootstrapN = 1000;
    cfg.randomSeed = 42;
    cfg.classA_absBiasMax = 0.10;
    cfg.classB_absBiasMax = 0.20;
    cfg.classA_chi2Max = 1.5;
    cfg.classB_chi2Max = 2.5;
    cfg.useLogRadiusInterpolation = true;
    cfg.includeComponentComparison = true;
end
function cfg = readConfigCSV(cfg, f)
    if ~isfile(f); return; end
    T = readtable(f, 'VariableNamingRule','preserve');
    for i=1:height(T)
        key = char(string(T.parameter(i))); val = T.value(i);
        if iscell(val); val=val{1}; end
        if isstring(val) || ischar(val); num=str2double(string(val)); else; num=double(val); end
        if isfield(cfg,key)
            if isfinite(num); cfg.(key)=num; else; cfg.(key)=val; end
        end
    end
end
function cfg = parseNameValue(cfg, varargin)
    if isempty(varargin); return; end
    if mod(numel(varargin),2)~=0; error('Options must be name/value pairs.'); end
    for i=1:2:numel(varargin)
        key = varargin{i}; val = varargin{i+1};
        if isfield(cfg,key); cfg.(key)=val; else; error('Unknown option: %s', key); end
    end
end
function seed = safeSeed(x)
    if iscell(x); x=x{1}; end
    if isstring(x)||ischar(x); seed=str2double(string(x)); else; seed=double(x); end
    if ~isfinite(seed) || seed<0 || seed>=2^32; seed=42; end
    seed = floor(double(seed));
end
function dirs = makeDirs(rootDir)
    dirs.root = fullfile(rootDir,'outputs');
    dirs.tables = fullfile(dirs.root,'tables');
    dirs.figures = fullfile(dirs.root,'figures');
    dirs.reports = fullfile(dirs.root,'reports');
    dirs.logs = fullfile(dirs.root,'logs');
    fns = fieldnames(dirs);
    for i=1:numel(fns); if ~isfolder(dirs.(fns{i})); mkdir(dirs.(fns{i})); end; end
end
function T = readRequired(f)
    if ~isfile(f); error('Required file not found: %s', f); end
    T = readtable(f, 'VariableNamingRule','preserve');
end
function T = readOptional(f)
    if isfile(f); T = readtable(f, 'VariableNamingRule','preserve'); else; T = table(); end
end
function writeTableSafe(f,T)
    try; writetable(T,f); catch ME; warning('Could not write %s: %s',f,ME.message); end
end
function safeClose(fid); if fid>0; fclose(fid); end; end
function logMsg(fid,msg)
    s=sprintf('[%s] %s', datestr(now,'yyyy-mm-dd HH:MM:SS'), msg);
    fprintf('%s\n',s); if fid>0; fprintf(fid,'%s\n',s); drawnow; end
end

% ========================================================================
% HACKS parser
% ========================================================================
function hacks = parseHacksRawText(rawFile)
    fid = fopen(rawFile,'r');
    rows = {};
    while true
        line = fgetl(fid);
        if ~ischar(line); break; end
        line = strtrim(line);
        if isempty(line) || startsWith(line,'#'); continue; end
        parts = regexp(line,'\s+','split');
        if numel(parts) < 10; continue; end
        cid = normalizeClusterId(parts{1});
        nums = str2double(parts(2:10));
        if any(~isfinite(nums)); continue; end
        rows(end+1,:) = {string(cid), string(parts{1}), nums(1), nums(2), nums(3), nums(4), nums(5), nums(6), nums(7), nums(8), nums(9), "parsed_raw_hacks"}; %#ok<AGROW>
    end
    fclose(fid);
    hacks = cell2table(rows, 'VariableNames', {'cluster_id','hacks_cluster_id','bin','n_stars','r_arcsec','sigma_pm_masyr','sigma_pm_err_masyr','sigma_radial_masyr','sigma_radial_err_masyr','sigma_tangential_masyr','sigma_tangential_err_masyr','source_note'});
end
function cid = normalizeClusterId(s)
    s=char(string(s)); u=upper(s);
    tok=regexp(u,'NGC[_\s-]*0*(\d+)','tokens','once');
    if ~isempty(tok); cid=sprintf('NGC%04d',str2double(tok{1})); return; end
    tok=regexp(u,'IC[_\s-]*0*(\d+)','tokens','once');
    if ~isempty(tok); cid=sprintf('IC%04d',str2double(tok{1})); return; end
    tok=regexp(u,'PAL[_\s-]*0*(\d+)','tokens','once');
    if ~isempty(tok); cid=sprintf('PAL%02d',str2double(tok{1})); return; end
    cid = regexprep(u,'[^A-Z0-9]+','_');
end

% ========================================================================
% Core analysis
% ========================================================================
function readiness = buildReadinessTable(gaiaSummary, gaiaAudit, hacks, apogeeSummary, crosswalk)
    ids = strings(0,1);
    if ~isempty(gaiaSummary); ids=[ids; string(gaiaSummary.cluster_id)]; end
    if ~isempty(gaiaAudit); ids=[ids; string(gaiaAudit.cluster_id)]; end
    if ~isempty(hacks); ids=[ids; string(hacks.cluster_id)]; end
    if ~isempty(apogeeSummary); ids=[ids; string(apogeeSummary.cluster_id)]; end
    if ~isempty(crosswalk); ids=[ids; string(crosswalk.cluster_id)]; end
    ids = unique(ids);
    readiness = table();
    for i=1:numel(ids)
        cid=ids(i);
        hasGaiaProf = ~isempty(gaiaSummary) && any(string(gaiaSummary.cluster_id)==cid);
        hasGaiaCat = ~isempty(gaiaAudit) && any(string(gaiaAudit.cluster_id)==cid);
        hasHST = ~isempty(hacks) && any(string(hacks.cluster_id)==cid);
        hasAP = ~isempty(apogeeSummary) && any(string(apogeeSummary.cluster_id)==cid);
        nHST = 0; nAP = 0; nProf=0; dens=NaN;
        if hasHST; nHST=sum(string(hacks.cluster_id)==cid); end
        if hasAP; row=apogeeSummary(string(apogeeSummary.cluster_id)==cid,:); nAP=row.n_apogee_rows(1); end
        if hasGaiaProf; row=gaiaSummary(string(gaiaSummary.cluster_id)==cid,:); nProf=row.n_profile_bins(1); end
        if hasGaiaCat && any(strcmp(gaiaAudit.Properties.VariableNames,'median_density_arcmin2'))
            row=gaiaAudit(string(gaiaAudit.cluster_id)==cid,:); dens=row.median_density_arcmin2(1);
        end
        ready = hasGaiaProf && hasHST && nHST>=3;
        readiness = [readiness; table(cid,hasGaiaCat,hasGaiaProf,hasHST,hasAP,nProf,nHST,nAP,dens,ready, ...
            'VariableNames',{'cluster_id','gaia_catalogue_available','gaia_profile_available','hst_hacks_profile_available','apogee_available','n_gaia_profile_bins','n_hst_bins','n_apogee_rows','median_gaia_density_arcmin2','ready_for_gaia_hst_validation'})]; %#ok<AGROW>
    end
end

function [matched, summary] = compareOneCluster(cid, g, h, gaiaAudit, apogeeSummary, cfg)
    % Prepare Gaia profile: radius in arcsec and uncertainty from percentiles.
    g.radius_arcsec = g.radius_deg * 3600;
    g.sigma_pm_masyr = g.disp_p50;
    g.sigma_pm_err_masyr = 0.5 * abs(g.disp_p84_1 - g.disp_p15_9);
    g.sigma_radial_masyr = NaN(height(g),1);
    g.sigma_tangential_masyr = NaN(height(g),1);
    g.sigma_radial_err_masyr = NaN(height(g),1);
    g.sigma_tangential_err_masyr = NaN(height(g),1);

    % Select HST/HACKS bins.
    h = h(h.n_stars >= cfg.minHSTStarsPerBin & isfinite(h.r_arcsec) & isfinite(h.sigma_pm_masyr),:);
    if isempty(h)
        matched = table();
        summary = makeClusterSummary(cid,0,NaN,NaN,NaN,NaN,NaN,"E_insufficient_HST_bins",NaN,NaN,NaN,NaN);
        return;
    end

    % Overlap radial window.
    gg = g(isfinite(g.radius_arcsec)&g.radius_arcsec>0&isfinite(g.sigma_pm_masyr),:);
    if height(gg) < 2
        matched=table(); summary=makeClusterSummary(cid,0,NaN,NaN,NaN,NaN,NaN,"E_insufficient_Gaia_profile",NaN,NaN,NaN,NaN); return;
    end
    rmin=min(gg.radius_arcsec); rmax=max(gg.radius_arcsec);
    h = h(h.r_arcsec>=rmin & h.r_arcsec<=rmax,:);
    if height(h) < cfg.minMatchedBins
        matched=table(); summary=makeClusterSummary(cid,height(h),NaN,NaN,NaN,NaN,NaN,"E_insufficient_overlap",NaN,NaN,NaN,NaN); return;
    end

    xg = gg.radius_arcsec; xh = h.r_arcsec;
    if cfg.useLogRadiusInterpolation
        xgI = log10(xg); xhI = log10(xh);
    else
        xgI = xg; xhI = xh;
    end
    sigG = interp1(xgI, gg.sigma_pm_masyr, xhI, 'linear', NaN);
    errG = interp1(xgI, max(gg.sigma_pm_err_masyr,eps), xhI, 'linear', NaN);

    sigH = h.sigma_pm_masyr; errH = h.sigma_pm_err_masyr;
    if any(strcmp(h.Properties.VariableNames,'sigma_radial_masyr'))
        sigHr = h.sigma_radial_masyr; sigHt = h.sigma_tangential_masyr;
        errHr = h.sigma_radial_err_masyr; errHt = h.sigma_tangential_err_masyr;
    else
        sigHr = NaN(height(h),1); sigHt=NaN(height(h),1); errHr=sigHr; errHt=sigHr;
    end

    ok = isfinite(sigG)&isfinite(errG)&isfinite(sigH)&isfinite(errH)&sigH>0;
    h=h(ok,:); sigG=sigG(ok); errG=errG(ok); sigH=sigH(ok); errH=errH(ok); sigHr=sigHr(ok); sigHt=sigHt(ok); errHr=errHr(ok); errHt=errHt(ok);
    if numel(sigG) < cfg.minMatchedBins
        matched=table(); summary=makeClusterSummary(cid,numel(sigG),NaN,NaN,NaN,NaN,NaN,"E_insufficient_valid_bins",NaN,NaN,NaN,NaN); return;
    end

    diff = sigG - sigH;
    delta = diff ./ sigH;
    uRaw = sqrt(errG.^2 + errH.^2 + cfg.defaultHSTSystematicFloor.^2);
    zRaw = diff ./ max(uRaw,eps);
    redChi2Raw = sum(zRaw.^2) / max(numel(zRaw)-1,1);
    meanDelta = mean(delta,'omitnan');
    medDelta = median(delta,'omitnan');
    rmsDelta = sqrt(mean(delta.^2,'omitnan'));

    % Empirical cluster-level error inflation/floor.
    [scale, floorSys, redChi2Corr] = fitMinimalErrorInflation(diff, errG, errH, cfg);
    totalU = sqrt((scale*errG).^2 + errH.^2 + floorSys.^2 + cfg.defaultHSTSystematicFloor.^2);
    zCorr = diff ./ max(totalU,eps);

    % Component-level anisotropy contrast from HST only.
    betaH = 1 - (sigHt.^2 ./ max(sigHr.^2,eps));

    % Gaia density/APOGEE readiness.
    dens = NaN; nGsrc=NaN; nGpri=NaN; nAP=0; rvMed=NaN; fehMed=NaN;
    if ~isempty(gaiaAudit) && any(string(gaiaAudit.cluster_id)==cid)
        row=gaiaAudit(string(gaiaAudit.cluster_id)==cid,:);
        if any(strcmp(row.Properties.VariableNames,'median_density_arcmin2')); dens=row.median_density_arcmin2(1); end
        if any(strcmp(row.Properties.VariableNames,'n_sources')); nGsrc=row.n_sources(1); end
        if any(strcmp(row.Properties.VariableNames,'n_primary_q2_p90')); nGpri=row.n_primary_q2_p90(1); end
    end
    if ~isempty(apogeeSummary) && any(string(apogeeSummary.cluster_id)==cid)
        row=apogeeSummary(string(apogeeSummary.cluster_id)==cid,:);
        nAP=row.n_apogee_rows(1); rvMed=row.rv_median_kms(1); fehMed=row.feh_median(1);
    end

    matched = table(repmat(cid,numel(sigG),1), h.bin, h.n_stars, h.r_arcsec, sigH, errH, sigG, errG, diff, delta, zRaw, zCorr, sigHr, errHr, sigHt, errHt, betaH, ...
        repmat(dens,numel(sigG),1), repmat(nAP,numel(sigG),1), ...
        'VariableNames',{'cluster_id','hacks_bin','hacks_n_stars','radius_arcsec','hst_sigma_pm_masyr','hst_sigma_pm_err_masyr','gaia_sigma_pm_masyr','gaia_sigma_pm_err_masyr','sigma_diff_masyr','delta_frac','z_raw','z_corrected','hst_sigma_radial_masyr','hst_sigma_radial_err_masyr','hst_sigma_tangential_masyr','hst_sigma_tangential_err_masyr','hst_beta_pm','median_gaia_density_arcmin2','n_apogee_rows'});

    qclass = classifyCluster(numel(sigG), meanDelta, rmsDelta, redChi2Corr, cfg);
    summary = makeClusterSummary(cid,numel(sigG),meanDelta,medDelta,rmsDelta,redChi2Raw,redChi2Corr,qclass,scale,floorSys,dens,nAP);
    summary.n_gaia_sources = nGsrc;
    summary.n_gaia_primary_p90_q2 = nGpri;
    summary.apogee_rv_median_kms = rvMed;
    summary.apogee_feh_median = fehMed;
end

function [scale, floorSys, redChi2] = fitMinimalErrorInflation(diff, errG, errH, cfg)
    scales = linspace(1, cfg.defaultGaiaScaleGridMax, 61);
    floors = linspace(0, cfg.defaultGaiaSystematicFloorGridMax, 61);
    best = [1,0,Inf,Inf];
    n = numel(diff);
    for a=scales
        for s=floors
            u = sqrt((a*errG).^2 + errH.^2 + s.^2);
            rc = sum((diff./max(u,eps)).^2)/max(n-1,1);
            penalty = abs(rc-1) + 0.03*(a-1) + 5*s;
            if rc <= 1.05 && penalty < best(4)
                best=[a,s,rc,penalty];
            elseif best(3)==Inf && penalty < best(4)
                best=[a,s,rc,penalty];
            end
        end
    end
    scale=best(1); floorSys=best(2); redChi2=best(3);
end
function q = classifyCluster(n, meanDelta, rmsDelta, redChi2, cfg)
    if n < cfg.minMatchedBins
        q = "E_insufficient_overlap"; return;
    end
    if ~isfinite(meanDelta) || ~isfinite(redChi2)
        q = "E_unusable"; return;
    end
    if abs(meanDelta) < cfg.classA_absBiasMax && redChi2 < cfg.classA_chi2Max
        q = "A_Gaia_HST_consistent";
    elseif abs(meanDelta) < cfg.classB_absBiasMax && redChi2 < cfg.classB_chi2Max
        q = "B_small_correction";
    elseif abs(meanDelta) < 0.35 && redChi2 < 5
        q = "C_correction_required";
    else
        q = "D_not_reliable_in_overlap";
    end
end
function summary = makeClusterSummary(cid,n,meanDelta,medDelta,rmsDelta,chiRaw,chiCorr,qclass,scale,floorSys,dens,nAP)
    summary = table(cid,n,meanDelta,medDelta,rmsDelta,chiRaw,chiCorr,string(qclass),scale,floorSys,dens,nAP, ...
        'VariableNames',{'cluster_id','n_matched_bins','mean_delta_frac','median_delta_frac','rms_delta_frac','redchi2_raw','redchi2_corrected','quality_class','gaia_error_scale','gaia_systematic_floor_masyr','median_gaia_density_arcmin2','n_apogee_rows'});
end
function G = computeGlobalSummary(matched, classes, cfg)
    if isempty(matched)
        G=table(0,NaN,NaN,NaN,NaN,string('no_matched_bins'), 'VariableNames', {'n_matched_bins','global_mean_delta_frac','global_median_delta_frac','global_rms_delta_frac','global_redchi2_raw','status'}); return;
    end
    delta = matched.delta_frac; z=matched.z_raw;
    G = table(height(matched), numel(unique(matched.cluster_id)), mean(delta,'omitnan'), median(delta,'omitnan'), sqrt(mean(delta.^2,'omitnan')), sum(z.^2)/max(numel(z)-1,1), string(join(unique(classes.quality_class),';')), ...
        'VariableNames',{'n_matched_bins','n_clusters','global_mean_delta_frac','global_median_delta_frac','global_rms_delta_frac','global_redchi2_raw','quality_classes_present'});
end

% ========================================================================
% Figures and reports
% ========================================================================
function makeFigures(dirs, matched, classes, readiness)
    if isempty(matched); return; end
    try
        f=figure('Visible','off'); histogram(matched.delta_frac,30); xlabel('(Gaia-HST)/HST dispersion'); ylabel('Matched radial bins'); title('Gaia--HST dispersion residuals'); saveas(f,fullfile(dirs.figures,'GaiaHST_residual_histogram.png')); close(f);
        f=figure('Visible','off'); scatter(matched.radius_arcsec, matched.delta_frac, 20, matched.median_gaia_density_arcmin2, 'filled'); set(gca,'XScale','log'); xlabel('Radius [arcsec]'); ylabel('(Gaia-HST)/HST'); cb=colorbar; ylabel(cb,'Gaia source density [arcmin^{-2}]'); title('Residuals versus radius and density'); saveas(f,fullfile(dirs.figures,'GaiaHST_residual_vs_radius_density.png')); close(f);
        f=figure('Visible','off'); cats=categorical(classes.quality_class); histogram(cats); ylabel('Number of clusters'); title('Gaia internal-kinematics quality classes'); saveas(f,fullfile(dirs.figures,'GaiaHST_quality_class_counts.png')); close(f);
        % example profiles for up to 6 best clusters by matched bins
        [~,ord]=sort(classes.n_matched_bins,'descend'); n=min(6,numel(ord));
        f=figure('Visible','off'); tiledlayout(2,3,'Padding','compact');
        for k=1:n
            cid=classes.cluster_id(ord(k)); m=matched(matched.cluster_id==cid,:);
            nexttile; errorbar(m.radius_arcsec,m.hst_sigma_pm_masyr,m.hst_sigma_pm_err_masyr,'o'); hold on; errorbar(m.radius_arcsec,m.gaia_sigma_pm_masyr,m.gaia_sigma_pm_err_masyr,'s'); set(gca,'XScale','log'); title(char(cid)); xlabel('R [arcsec]'); ylabel('\sigma_{PM} [mas yr^{-1}]');
        end
        saveas(f,fullfile(dirs.figures,'GaiaHST_example_profile_comparisons.png')); close(f);
    catch ME
        warning('Figure generation failed: %s', ME.message);
    end
end
function writeInterpretationReport(dirs, matched, classes, globalSummary, readiness, cfg)
    fid=fopen(fullfile(dirs.reports,'GaiaHST_interpretation_report.txt'),'w'); if fid<0; return; end
    fprintf(fid,'Gaia--HST--APOGEE kinematic validation report\n');
    fprintf(fid,'=====================================================\n\n');
    fprintf(fid,'Scientific framing: calibration/validation atlas. The code tests whether Gaia EDR3 internal-kinematic profiles are consistent with HST/HACKS crowded-core PM dispersion profiles in the radial overlap region.\n\n');
    fprintf(fid,'Data layer readiness:\n');
    fprintf(fid,'  Gaia profile clusters: %d\n', sum(readiness.gaia_profile_available));
    fprintf(fid,'  HST/HACKS profile clusters: %d\n', sum(readiness.hst_hacks_profile_available));
    fprintf(fid,'  APOGEE spectroscopy clusters: %d\n', sum(readiness.apogee_available));
    fprintf(fid,'  Gaia-HST validation-ready clusters: %d\n\n', sum(readiness.ready_for_gaia_hst_validation));
    if ~isempty(matched)
        fprintf(fid,'Matched radial bins: %d\n', height(matched));
        fprintf(fid,'Matched clusters: %d\n', numel(unique(matched.cluster_id)));
        fprintf(fid,'Global mean fractional residual: %.4f\n', globalSummary.global_mean_delta_frac(1));
        fprintf(fid,'Global RMS fractional residual: %.4f\n', globalSummary.global_rms_delta_frac(1));
        fprintf(fid,'Global raw reduced chi2: %.4f\n\n', globalSummary.global_redchi2_raw(1));
        fprintf(fid,'Quality class counts:\n');
        u=unique(classes.quality_class);
        for i=1:numel(u); fprintf(fid,'  %s: %d\n', u(i), sum(classes.quality_class==u(i))); end
    else
        fprintf(fid,'No matched Gaia-HST bins were available. Fill input/hacks_hst_dispersion_profiles.csv or provide input/raw/hacks_velocity_dispersions.txt.\n');
    end
    fprintf(fid,'\nClaim rules:\n');
    fprintf(fid,'  A/B classes: Gaia profile usable for dynamical modelling after stated calibration.\n');
    fprintf(fid,'  C class: correction required; use with caution and propagate systematic floor.\n');
    fprintf(fid,'  D/E classes: not reliable or insufficient overlap; no Gaia-only internal-dynamics claim.\n');
    fclose(fid);
end
function reportNoHacks(dirs, readiness)
    fid=fopen(fullfile(dirs.reports,'GaiaHST_interpretation_report.txt'),'w'); if fid<0; return; end
    fprintf(fid,'No HST/HACKS dispersion rows were available in the mounted package.\n');
    fprintf(fid,'The Gaia and APOGEE CSV layers were created successfully, but the core Gaia--HST comparison requires a populated input/hacks_hst_dispersion_profiles.csv.\n');
    fprintf(fid,'The uploaded HACKS text excerpt documents the required columns: cluster ID, bin, N, radius, combined/radial/tangential PM dispersions and uncertainties.\n');
    fprintf(fid,'Fill the CSV template or place the raw table at input/raw/hacks_velocity_dispersions.txt and rerun.\n');
    fclose(fid);
end
function writeOutputDictionary(dirs)
    fid=fopen(fullfile(dirs.reports,'GaiaHST_output_dictionary.txt'),'w'); if fid<0; return; end
    fprintf(fid,'GaiaHST_MATCHED_RADIAL_COMPARISON.csv: one row per matched HST radial bin with interpolated Gaia dispersion and residuals.\n');
    fprintf(fid,'GaiaHST_CLUSTER_QUALITY_CLASSES.csv: one row per cluster with fractional bias, chi2, empirical Gaia error scale/floor and quality class.\n');
    fprintf(fid,'GaiaHST_GLOBAL_CALIBRATION_SUMMARY.csv: global residual summary across all matched radial bins.\n');
    fprintf(fid,'GaiaHST_DATA_LAYER_READINESS.csv: availability of Gaia profile, Gaia catalogue, HST/HACKS and APOGEE layers.\n');
    fclose(fid);
end
