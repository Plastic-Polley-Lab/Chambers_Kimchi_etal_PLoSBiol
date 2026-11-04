function plot_ach_aligned_to_fme_continuous_bt(save_path, varargin)
% PLOT_ACH_ALIGNED_TO_FME_CONTINUOUS_BT
%
% Manuscript-quality figure: ACh (470 nm) aligned to FME bouts, with 405 nm
% reference overlaid as a dashed same-color line on all ACh panels.
%
% Usage:
%   plot_ach_aligned_to_fme_continuous_bt(save_path)
%   plot_ach_aligned_to_fme_continuous_bt(save_path, 'pooling', 'per_session')
%   plot_ach_aligned_to_fme_continuous_bt(save_path, 'pooling', 'pooled')
%   plot_ach_aligned_to_fme_continuous_bt(save_path, 'pooling', 'per_mouse')
%
%   'per_mouse' averages sessions within each mouse before statistics, so
%   each data point is one mouse. Mouse ID is parsed from the filename
%   prefix (e.g. 'sub1034_ses27_run3_traces.mat' -> mouse 'sub1034').
%   plot_ach_aligned_to_fme_continuous_bt(save_path, 'good_sessions', gs)
%
% Outputs: two PDF figures saved to save_path
%   ach_aligned_dff.pdf   — ΔF/F traces + AUC
%   ach_aligned_zscore.pdf — z-scored traces + AUC

%% -------------------------------------------------------------------------
%  Parse inputs
%% -------------------------------------------------------------------------
p = inputParser;
addParameter(p, 'pooling',       'per_session', ...
    @(x) ismember(x, {'per_session','pooled','per_mouse'}));
addParameter(p, 'good_sessions', [], @(x) isempty(x) || isstruct(x));
addParameter(p, 'gfp_path',      '', @(x) ischar(x) || isstring(x));
parse(p, varargin{:});
poolingMode   = p.Results.pooling;
good_sessions = p.Results.good_sessions;
gfpPath       = char(p.Results.gfp_path);

fprintf('Pooling mode: %s\n', poolingMode);

%% -------------------------------------------------------------------------
%  Parameters
%% -------------------------------------------------------------------------
fsAch          = 1017;
fsFme          = 30;
sampFactor     = fsAch / fsFme;

minBoutSamples = round(10  * sampFactor);
minIsolSamples = round(30  * sampFactor);
medFiltWidth   = 500;
boutThreshPct  = 85;
shortMaxSec    = 1.0;
longMinSec     = 2.0;

preEvt   = 1.0;
postEvt  = 4.0;
preS     = round(preEvt  * fsAch);
postS    = round(postEvt * fsAch);
tEvt     = (-preS:postS) / fsAch;
nEvtSamp = preS + postS + 1;
blStart  = 1;
blEnd    = preS;

aucWinShort = [0, 1];
aucWinLong  = [0, 4.0];
aucWinMid   = [0, 2.0];   % intermediate (1-2 s) bouts

%% -------------------------------------------------------------------------
%  Colors / style
%% -------------------------------------------------------------------------
shortColor   = [0.00, 0.45, 0.70];
longColor    = [0.90, 0.62, 0.00];
regionColors = [
    0.00, 0.45, 0.70;   % AC  — blue
    0.90, 0.62, 0.00;   % VC  — orange
    0.00, 0.62, 0.45;   % PFC — teal
];
regionNames  = {'AC','VC','PFC'};
nR           = numel(regionNames);
% 3 bout categories are ACCUMULATED: short, long, and the intermediate
% 1-2 s band. The two main figures show only the first two ("published")
% categories; the mid band is reported separately in its own summary.
nBouts           = 3;   % categories accumulated
nPublishedBouts  = 2;   % categories shown in the main trace/AUC figures
nSigs        = 2 + nR;   % FME, Pupil, AC, VC, PFC

boutColors   = {shortColor, longColor, [0.55 0.35 0.65]};   % 3rd = mid band
boutTitleStr = {sprintf('Short bouts (<%.0f s)',   shortMaxSec), ...
                sprintf('Long bouts (\\geq%.0f s)', longMinSec), ...
                sprintf('Mid bouts (%.0f-%.0f s)', shortMaxSec, longMinSec)};
bLabels      = {sprintf('Short (<%.0f s)',  shortMaxSec), ...
                sprintf('Long (\\geq%.0f s)', longMinSec), ...
                sprintf('Mid (%.0f-%.0f s)', shortMaxSec, longMinSec)};

% Figure / font sizes
FIG_W   = 7.0;   % inches
FIG_H   = 6.0;
FONT_AX = 8;     % tick labels
FONT_LB = 9;     % axis labels
FONT_TI = 9;     % panel titles
LW_SIG  = 1.5;   % signal line width
LW_REF  = 0.8;   % reference line width
CI_ALPHA = 0.15; % CI shading alpha

%% -------------------------------------------------------------------------
%  Locate signal files
%% -------------------------------------------------------------------------
if ~isempty(good_sessions) && isfield(good_sessions, 'files')
    allFiles = good_sessions.files;
    % Keep only signal files (no 'ref' in name)
    isSig  = arrayfun(@(f) isempty(strfind(f.name, 'ref')), allFiles);
    sigFiles = allFiles(isSig);
    fprintf('Using QC file list — %d signal files\n', numel(sigFiles));
else
    sigFiles = dir(fullfile(save_path, '*_traces.mat'));
    % Exclude ref files
    isRef    = arrayfun(@(f) ~isempty(strfind(f.name, 'ref')), sigFiles);
    sigFiles = sigFiles(~isRef);
end
if isempty(sigFiles)
    error('No signal trace files found in %s', save_path);
end
nSessions = numel(sigFiles);
fprintf('Found %d signal session files\n', nSessions);

%% -------------------------------------------------------------------------
%  Pair each signal file with its ref counterpart
%% -------------------------------------------------------------------------
refFiles = cell(nSessions, 1);   % path strings, empty if not found
for i = 1:nSessions
    sigName = sigFiles(i).name;                          % e.g. sub1_ses1_run1_traces.mat
    refName = strrep(sigName, '_traces.mat', '_ref_traces.mat');
    refPath = fullfile(sigFiles(i).folder, refName);
    if exist(refPath, 'file')
        refFiles{i} = refPath;
    else
        warning('No ref file found for %s — signal only', sigName);
    end
end
nRefFound = sum(~cellfun(@isempty, refFiles));
fprintf('Paired ref files: %d / %d sessions\n', nRefFound, nSessions);

%% -------------------------------------------------------------------------
%  Validate QC mask
%% -------------------------------------------------------------------------
useQC = false;
if ~isempty(good_sessions)
    maskLen = numel(good_sessions.mask{1});
    if maskLen == nSessions
        useQC = true;
        fprintf('QC mask applied:\n');
        for r = 1:nR
            fprintf('  %s: %d / %d sessions pass\n', ...
                regionNames{r}, sum(good_sessions.mask{r}), nSessions);
        end
    else
        warning('QC mask length (%d) != session count (%d). Ignoring.', ...
            maskLen, nSessions);
    end
end

%% -------------------------------------------------------------------------
%  Storage — signal and ref accumulate separately
%% -------------------------------------------------------------------------
% Dimensions: {bout, signal_index}
%   signal index: 1=FME, 2=Pupil, 3=AC, 4=VC, 5=PFC
traceData_raw    = cell(nBouts, nSigs);
traceData_zscore = cell(nBouts, nSigs);
refData_raw      = cell(nBouts, nR);    % ref only for ACh regions
refData_zscore   = cell(nBouts, nR);
aucData_raw      = cell(nBouts, nR);
aucData_zscore   = cell(nBouts, nR);

% Mouse-ID tracking (one entry per accumulated row) — used by 'per_mouse'
mouseID_trace    = cell(nBouts, nSigs);
mouseID_auc      = cell(nBouts, nR);
mouseID_ref      = cell(nBouts, nR);

% Which regions each accumulated FME row came from ([nRows x nR] logical).
% Lets the FME traces be grouped by cohort for the Fig 1F-I comparability
% check (R1 minor): different regions are recorded in different mice, so the
% reviewer asked whether the movements themselves are comparable across them.
regionMask_fme   = cell(nBouts, 1);
for b = 1:nBouts
    regionMask_fme{b} = false(0, nR);
end
fmeGFP           = cell(nBouts, 1);   % GFP-cohort FME traces, if gfp_path given

% GFP-cohort FLUORESCENCE aligned to that cohort's own FME bouts, per region.
% Previously only GFP *FME* was collected (for the cohort-comparability panel);
% the GFP signal itself was never aligned to bouts, so there were no control
% statistics for Fig 1. One row per GFP mouse, same AUC windows as ACh.
aucGFP_raw    = cell(nBouts, nR);
aucGFP_zscore = cell(nBouts, nR);
for b = 1:nBouts
    for r = 1:nR
        aucGFP_raw{b,r}    = [];
        aucGFP_zscore{b,r} = [];
    end
end
for b = 1:nBouts
    fmeGFP{b} = [];
end

% All detected bout durations BEFORE the short/long filter, so the
% distribution (and the size of the excluded 1-2 s band) can be inspected.
allBoutDurs      = [];
allBoutDurMouse  = {};

for b = 1:nBouts
    for s = 1:nSigs
        traceData_raw{b,s}    = [];
        traceData_zscore{b,s} = [];
    end
    for r = 1:nR
        refData_raw{b,r}    = [];
        refData_zscore{b,r} = [];
        aucData_raw{b,r}    = [];
        aucData_zscore{b,r} = [];
    end
end

boutCounts = zeros(nSessions, nBouts);

%% =========================================================================
%  Main loop
%% =========================================================================
for i = 1:nSessions
    fprintf('Processing session %d / %d  (%s)\n', i, nSessions, sigFiles(i).name);

    % --- Parse mouse ID from filename (e.g. 'sub1034_ses27_run3_traces.mat') ---
    tokMouse = regexp(sigFiles(i).name, '^(sub[0-9]+)', 'tokens', 'once');
    if isempty(tokMouse)
        thisMouseID = sigFiles(i).name;   % fallback: treat file as its own mouse
        warning('Could not parse mouse ID from %s — treating as unique mouse.', ...
            sigFiles(i).name);
    else
        thisMouseID = tokMouse{1};
    end

    % --- Load signal ---
    tmp = load(fullfile(sigFiles(i).folder, sigFiles(i).name), 'traces');
    traces = tmp.traces;

    fmeRes      = traces.fme;
    pupilRes    = traces.pupil;
    achCh       = {traces.achCh1, traces.achCh2};
    nAchSamp    = traces.nAchSamp;
    chForRegion = traces.chForRegion;

    % --- Load ref (if available) ---
    hasRef   = ~isempty(refFiles{i});
    refChAll = {};
    if hasRef
        tmpR     = load(refFiles{i}, 'traces');
        refChAll = {tmpR.traces.achCh1, tmpR.traces.achCh2};
        clear tmpR;
    end

    % --- QC mask ---
    regionPassQC = true(nR, 1);
    if useQC
        for r = 1:nR
            regionPassQC(r) = good_sessions.mask{r}(i);
        end
    end
    skipped = regionNames(~regionPassQC);
    if ~isempty(skipped)
        fprintf('  QC fail — skipping: %s\n', strjoin(skipped, ', '));
    end

    %----------------------------------------------------------------------
    %  Bout detection  (unchanged from original)
    %----------------------------------------------------------------------
    cheek     = medfilt1(fmeRes, medFiltWidth);
    validMask = ~isnan(cheek);
    thresh    = prctile(cheek(validMask), boutThreshPct);

    movement  = cheek > thresh;
    stillness = cheek < thresh;
    movement(~validMask)  = false;
    stillness(~validMask) = false;

    runStarts = find(diff(movement) > 0) + 1;
    runEnds   = find(diff(movement) < 0) + 1;
    if length(runStarts) > length(runEnds); runStarts(end) = []; end
    if length(runStarts) < length(runEnds); runEnds(1)    = []; end

    eligMask_ = (runEnds - runStarts) >= minBoutSamples;
    boutStarts = runStarts(eligMask_);
    boutEnds   = runEnds(eligMask_);

    movElig = false(size(movement));
    for k = 1:length(boutStarts)
        winEnd = boutStarts(k) + round(29 * sampFactor);
        if winEnd <= length(cheek) && ...
           sum(~isnan(cheek(boutStarts(k):winEnd))) > round(28 * sampFactor)
            movElig(boutStarts(k):boutEnds(k)) = true;
        end
    end

    eligStarts = find(diff(movElig) > 0);
    eligEnds   = find(diff(movElig) < 0);

    tooEarly = eligStarts < minIsolSamples + 1;
    eligEnds(tooEarly)   = [];
    eligStarts(tooEarly) = [];
    if isempty(eligStarts)
        fprintf('  No eligible bouts\n');
        continue;
    end

    discardBout = false(size(eligStarts));
    for k = 1:length(eligStarts)
        precStill = stillness(eligStarts(k)-minIsolSamples : ...
            eligStarts(k) - round(11 * sampFactor));
        if sum(precStill) <= minIsolSamples - 20 * sampFactor
            discardBout(k) = true;
        end
    end
    eligStarts(discardBout) = [];
    eligEnds(discardBout)   = [];

    tooLate = eligStarts > length(movement) - round(45 * sampFactor);
    eligStarts(tooLate) = [];
    eligEnds(tooLate)   = [];

    noRoom = eligStarts <= preS | (eligStarts + postS) > nAchSamp;
    eligStarts(noRoom) = [];
    eligEnds(noRoom)   = [];
    if isempty(eligStarts)
        fprintf('  No bouts with sufficient window\n');
        continue;
    end

    boutDurs   = (eligEnds - eligStarts) / fsAch;

    % Record every eligible bout duration before short/long selection
    allBoutDurs     = [allBoutDurs; boutDurs(:)];
    allBoutDurMouse = [allBoutDurMouse; repmat({thisMouseID}, numel(boutDurs), 1)];

    isShort    = boutDurs <  shortMaxSec;
    isLong     = boutDurs >= longMinSec;
    isMid      = boutDurs >= shortMaxSec & boutDurs < longMinSec;
    keepBouts  = isShort | isLong | isMid;
    eligStarts = eligStarts(keepBouts);
    eligEnds   = eligEnds(keepBouts);
    isShort    = isShort(keepBouts);
    isLong     = isLong(keepBouts);
    isMid      = isMid(keepBouts);
    if isempty(eligStarts); continue; end

    boutCounts(i,1) = sum(isShort);
    boutCounts(i,2) = sum(isLong);
    boutCounts(i,3) = sum(isMid);
    fprintf('  Short: %d  Long: %d  Mid(%.0f-%.0f s): %d\n', ...
        sum(isShort), sum(isLong), shortMaxSec, longMinSec, sum(isMid));

    %----------------------------------------------------------------------
    %  Snippet extraction
    %----------------------------------------------------------------------
    nElig  = numel(eligStarts);
    idxMtx = eligStarts(:) + (-preS:postS);   % [nElig x nEvtSamp]

    fmeSnipAll   = fmeRes(idxMtx);
    pupilSnipAll = pupilRes(idxMtx);

    % Signal snippets per region
    achSnip   = cell(nR, 1);
    achSnip_z = cell(nR, 1);
    refSnip   = cell(nR, 1);
    refSnip_z = cell(nR, 1);

    for r = 1:nR
        if chForRegion(r) == 0 || ~regionPassQC(r)
            achSnip{r}   = nan(nElig, nEvtSamp);
            achSnip_z{r} = nan(nElig, nEvtSamp);
            refSnip{r}   = nan(nElig, nEvtSamp);
            refSnip_z{r} = nan(nElig, nEvtSamp);
        else
            ch = chForRegion(r);

            % Signal
            M      = achCh{ch}(idxMtx);
            mu_r   = mean(M, 2, 'omitnan');
            sd_r   = std(M,  0, 2, 'omitnan');
            sd_r(sd_r == 0) = 1;
            achSnip{r}   = M;
            achSnip_z{r} = (M - mu_r) ./ sd_r;

            % Reference
            if hasRef
                Mref     = refChAll{ch}(idxMtx);
                mu_rr    = mean(Mref, 2, 'omitnan');
                sd_rr    = std(Mref,  0, 2, 'omitnan');
                sd_rr(sd_rr == 0) = 1;
                refSnip{r}   = Mref;
                refSnip_z{r} = (Mref - mu_rr) ./ sd_rr;
            else
                refSnip{r}   = nan(nElig, nEvtSamp);
                refSnip_z{r} = nan(nElig, nEvtSamp);
            end
        end
    end

    % AUC window indices
    aucIdx = {tEvt >= aucWinShort(1) & tEvt <= aucWinShort(2), ...
              tEvt >= aucWinLong(1)  & tEvt <= aucWinLong(2), ...
              tEvt >= aucWinMid(1)   & tEvt <= aucWinMid(2)};

    %----------------------------------------------------------------------
    %  Accumulate by bout type
    %----------------------------------------------------------------------
    boutTypeIdx = {find(isShort), find(isLong), find(isMid)};

    for b = 1:nBouts
        bIdx = boutTypeIdx{b};
        if isempty(bIdx); continue; end

        % FME
        Mfme = fmeSnipAll(bIdx, :);
        Mfme = Mfme(all(isfinite(Mfme), 2), :);
        if ~isempty(Mfme)
            switch poolingMode
                case {'per_session','per_mouse'}
                    traceData_raw{b,1}    = [traceData_raw{b,1};    mean(Mfme, 1, 'omitnan')];
                    traceData_zscore{b,1} = [traceData_zscore{b,1}; mean(Mfme, 1, 'omitnan')];
                    mouseID_trace{b,1}    = [mouseID_trace{b,1}; {thisMouseID}];
                    regionMask_fme{b}     = [regionMask_fme{b}; chForRegion(:)' > 0];
                case 'pooled'
                    traceData_raw{b,1}    = [traceData_raw{b,1};    Mfme];
                    traceData_zscore{b,1} = [traceData_zscore{b,1}; Mfme];
            end
        end

        % Pupil
        Mpup = pupilSnipAll(bIdx, :);
        Mpup = Mpup(all(isfinite(Mpup), 2), :);
        if ~isempty(Mpup)
            switch poolingMode
                case {'per_session','per_mouse'}
                    traceData_raw{b,2}    = [traceData_raw{b,2};    mean(Mpup, 1, 'omitnan')];
                    traceData_zscore{b,2} = [traceData_zscore{b,2}; mean(Mpup, 1, 'omitnan')];
                    mouseID_trace{b,2}    = [mouseID_trace{b,2}; {thisMouseID}];
                case 'pooled'
                    traceData_raw{b,2}    = [traceData_raw{b,2};    Mpup];
                    traceData_zscore{b,2} = [traceData_zscore{b,2}; Mpup];
            end
        end

        % ACh + ref per region
        for r = 1:nR

            % -- Signal --
            Mraw = achSnip{r}(bIdx, :);
            Mraw = Mraw(all(isfinite(Mraw), 2), :);
            if ~isempty(Mraw)
                switch poolingMode
                    case {'per_session','per_mouse'}
                        muRaw = mean(Mraw, 1, 'omitnan');
                        muRaw = muRaw - mean(muRaw(blStart:blEnd), 'omitnan');
                        traceData_raw{b, 2+r} = [traceData_raw{b, 2+r}; muRaw];
                        aucData_raw{b,r} = [aucData_raw{b,r}; ...
                            trapz(tEvt(aucIdx{b}), muRaw(aucIdx{b}))];
                        mouseID_trace{b,2+r} = [mouseID_trace{b,2+r}; {thisMouseID}];
                        mouseID_auc{b,r}     = [mouseID_auc{b,r};     {thisMouseID}];
                    case 'pooled'
                        blM  = mean(Mraw(:, blStart:blEnd), 2, 'omitnan');
                        Mraw = Mraw - blM;
                        traceData_raw{b, 2+r} = [traceData_raw{b, 2+r}; Mraw];
                        for ei = 1:size(Mraw, 1)
                            aucData_raw{b,r} = [aucData_raw{b,r}; ...
                                trapz(tEvt(aucIdx{b}), Mraw(ei, aucIdx{b}))];
                        end
                end
            end

            Mz = achSnip_z{r}(bIdx, :);
            Mz = Mz(all(isfinite(Mz), 2), :);
            if ~isempty(Mz)
                switch poolingMode
                    case {'per_session','per_mouse'}
                        muZ = mean(Mz, 1, 'omitnan');
                        muZ = muZ - mean(muZ(blStart:blEnd), 'omitnan');
                        traceData_zscore{b, 2+r} = [traceData_zscore{b, 2+r}; muZ];
                        aucData_zscore{b,r} = [aucData_zscore{b,r}; ...
                            trapz(tEvt(aucIdx{b}), muZ(aucIdx{b}))];
                    case 'pooled'
                        blM = mean(Mz(:, blStart:blEnd), 2, 'omitnan');
                        Mz  = Mz - blM;
                        traceData_zscore{b, 2+r} = [traceData_zscore{b, 2+r}; Mz];
                        for ei = 1:size(Mz, 1)
                            aucData_zscore{b,r} = [aucData_zscore{b,r}; ...
                                trapz(tEvt(aucIdx{b}), Mz(ei, aucIdx{b}))];
                        end
                end
            end

            % -- Reference --
            Rraw = refSnip{r}(bIdx, :);
            Rraw = Rraw(all(isfinite(Rraw), 2), :);
            if ~isempty(Rraw)
                switch poolingMode
                    case {'per_session','per_mouse'}
                        muRef = mean(Rraw, 1, 'omitnan');
                        muRef = muRef - mean(muRef(blStart:blEnd), 'omitnan');
                        refData_raw{b,r} = [refData_raw{b,r}; muRef];
                        mouseID_ref{b,r} = [mouseID_ref{b,r}; {thisMouseID}];
                    case 'pooled'
                        blM  = mean(Rraw(:, blStart:blEnd), 2, 'omitnan');
                        Rraw = Rraw - blM;
                        refData_raw{b,r} = [refData_raw{b,r}; Rraw];
                end
            end

            Rz = refSnip_z{r}(bIdx, :);
            Rz = Rz(all(isfinite(Rz), 2), :);
            if ~isempty(Rz)
                switch poolingMode
                    case {'per_session','per_mouse'}
                        muRz = mean(Rz, 1, 'omitnan');
                        muRz = muRz - mean(muRz(blStart:blEnd), 'omitnan');
                        refData_zscore{b,r} = [refData_zscore{b,r}; muRz];
                    case 'pooled'
                        blM = mean(Rz(:, blStart:blEnd), 2, 'omitnan');
                        Rz  = Rz - blM;
                        refData_zscore{b,r} = [refData_zscore{b,r}; Rz];
                end
            end
        end
    end
end


%% =========================================================================
%  per_mouse: collapse session-level rows to one row per mouse
%  Sessions from the same mouse are averaged before any statistics or plots.
%% =========================================================================
if strcmp(poolingMode, 'per_mouse')
    fprintf('\nCollapsing sessions within mouse (per_mouse mode)...\n');

    for b = 1:nBouts
        % --- Trace data (FME, Pupil, and each ACh region) ---
        for sg = 1:nSigs
            ids = mouseID_trace{b,sg};
            if isempty(ids); continue; end
            uIDs = unique(ids, 'stable');
            Traw = traceData_raw{b,sg};
            Tz   = traceData_zscore{b,sg};
            newRaw = []; newZ = [];
            for u = 1:numel(uIDs)
                sel = strcmp(ids, uIDs{u});
                if ~isempty(Traw) && size(Traw,1) == numel(ids)
                    newRaw = [newRaw; mean(Traw(sel,:), 1, 'omitnan')]; %#ok
                end
                if ~isempty(Tz) && size(Tz,1) == numel(ids)
                    newZ   = [newZ;   mean(Tz(sel,:),   1, 'omitnan')]; %#ok
                end
            end
            if ~isempty(newRaw); traceData_raw{b,sg}    = newRaw; end
            if ~isempty(newZ);   traceData_zscore{b,sg} = newZ;   end
        end

        % --- AUC data per region ---
        for r = 1:nR
            ids = mouseID_auc{b,r};
            if isempty(ids); continue; end
            uIDs = unique(ids, 'stable');
            Araw = aucData_raw{b,r};
            Az   = aucData_zscore{b,r};
            newRaw = []; newZ = [];
            for u = 1:numel(uIDs)
                sel = strcmp(ids, uIDs{u});
                if ~isempty(Araw) && numel(Araw) == numel(ids)
                    newRaw(end+1,1) = mean(Araw(sel), 'omitnan'); %#ok
                end
                if ~isempty(Az) && numel(Az) == numel(ids)
                    newZ(end+1,1)   = mean(Az(sel),   'omitnan'); %#ok
                end
            end
            if ~isempty(newRaw); aucData_raw{b,r}    = newRaw; end
            if ~isempty(newZ);   aucData_zscore{b,r} = newZ;   end
        end

        % --- Reference (405 nm) traces per region ---
        for r = 1:nR
            ids = mouseID_ref{b,r};
            if isempty(ids); continue; end
            uIDs = unique(ids, 'stable');
            Rraw = refData_raw{b,r};
            Rz   = refData_zscore{b,r};
            newRaw = []; newZ = [];
            for u = 1:numel(uIDs)
                sel = strcmp(ids, uIDs{u});
                if ~isempty(Rraw) && size(Rraw,1) == numel(ids)
                    newRaw = [newRaw; mean(Rraw(sel,:), 1, 'omitnan')]; %#ok
                end
                if ~isempty(Rz) && size(Rz,1) == numel(ids)
                    newZ   = [newZ;   mean(Rz(sel,:),   1, 'omitnan')]; %#ok
                end
            end
            if ~isempty(newRaw); refData_raw{b,r}    = newRaw; end
            if ~isempty(newZ);   refData_zscore{b,r} = newZ;   end
        end
    end

    % Report how many mice / sessions contributed
    allIDs = {};
    for b = 1:nBouts
        for r = 1:nR
            if ~isempty(mouseID_auc{b,r})
                allIDs = [allIDs; mouseID_auc{b,r}]; %#ok
            end
        end
    end
    if ~isempty(allIDs)
        fprintf('  %d sessions collapsed to %d unique mice\n', ...
            numel(allIDs), numel(unique(allIDs)));
    end
end

% Output filename tag (defined here because the duration block below is the
% first thing that writes a file)
qcTag = '';
if useQC; qcTag = '_qc'; end

%% =========================================================================
%  Optional: FME from the GFP control cohort
%  Same bout detection, run over a second directory, so the GFP mice's
%  movements can be compared against the GRAB-ACh mice's (R1 minor).
%% =========================================================================
if ~isempty(gfpPath)
    gfpFiles = dir(fullfile(gfpPath, '*_traces.mat'));
    isRefG   = arrayfun(@(f) ~isempty(strfind(f.name, 'ref')), gfpFiles);
    gfpFiles = gfpFiles(~isRefG);
    fprintf('\nGFP cohort: %d sessions in %s\n', numel(gfpFiles), gfpPath);

    for ig = 1:numel(gfpFiles)
        tmpG = load(fullfile(gfpFiles(ig).folder, gfpFiles(ig).name), 'traces');
        tg   = tmpG.traces;
        if ~isfield(tg,'fme'); continue; end
        fmeG = tg.fme(:)';

        cheekG = medfilt1(fmeG, medFiltWidth);
        vmG    = ~isnan(cheekG);
        if ~any(vmG); continue; end
        thG    = prctile(cheekG(vmG), boutThreshPct);
        movG   = cheekG > thG;  movG(~vmG) = false;

        rsG = find(diff(movG) > 0) + 1;
        reG = find(diff(movG) < 0) + 1;
        if numel(rsG) > numel(reG); rsG(end) = []; end
        if numel(rsG) < numel(reG); reG(1)   = []; end
        if isempty(rsG); continue; end

        okG = (reG - rsG) >= minBoutSamples;
        rsG = rsG(okG); reG = reG(okG);
        keepG = rsG > preS & reG + postS <= numel(fmeG);
        rsG = rsG(keepG); reG = reG(keepG);
        if isempty(rsG); continue; end

        durG = (reG - rsG) / fsAch;
        selG = {durG < shortMaxSec, durG >= longMinSec, ...
                durG >= shortMaxSec & durG < longMinSec};

        for b = 1:nBouts
            idxG = find(selG{b});
            if isempty(idxG); continue; end
            % NOTE: no baseline subtraction here — the ACh path stores raw
            % FME snippets (fmeSnipAll = fmeRes(idxMtx)), so subtracting the
            % pre-bout baseline for GFP only would remove a large positive
            % offset and make the cohorts incomparable. Any baseline handling
            % is applied identically to both further below.
            MG = nan(numel(idxG), nEvtSamp);
            for k = 1:numel(idxG)
                MG(k,:) = fmeG(rsG(idxG(k))-preS : rsG(idxG(k))+postS);
            end
            MG = MG(all(isfinite(MG),2),:);
            if isempty(MG); continue; end
            fmeGFP{b} = [fmeGFP{b}; mean(MG, 1, 'omitnan')];
        end

        % --- GFP FLUORESCENCE aligned to this cohort's own FME bouts -----
        % Same alignment, same baseline subtraction and same AUC windows as
        % the ACh path, so the two cohorts are treated identically. GFP mice
        % have their own bouts and their own fibre placement; the comparison
        % is between cohorts, not within mouse.
        if isfield(tg, 'achCh1') && isfield(tg, 'chForRegion')
            achChG    = {tg.achCh1, tg.achCh2};
            cfrG      = tg.chForRegion;
            for r = 1:nR
                if numel(cfrG) < r || cfrG(r) == 0; continue; end
                sigG = achChG{cfrG(r)};
                if isempty(sigG); continue; end
                sigG = sigG(:)';
                for b = 1:nBouts
                    idxG = find(selG{b});
                    if isempty(idxG); continue; end
                    okIdx = rsG(idxG) > preS & rsG(idxG) + postS <= numel(sigG);
                    idxG  = idxG(okIdx);
                    if isempty(idxG); continue; end

                    SG = nan(numel(idxG), nEvtSamp);
                    for k = 1:numel(idxG)
                        SG(k,:) = sigG(rsG(idxG(k))-preS : rsG(idxG(k))+postS);
                    end
                    SG = SG(all(isfinite(SG),2),:);
                    if isempty(SG); continue; end

                    % raw: mean across bouts, baseline-subtracted (as ACh)
                    muG = mean(SG, 1, 'omitnan');
                    muG = muG - mean(muG(blStart:blEnd), 'omitnan');
                    aucGFP_raw{b,r} = [aucGFP_raw{b,r}; ...
                        trapz(tEvt(aucIdx{b}), muG(aucIdx{b}))];

                    % z-scored per snippet, then baseline re-zeroed (as ACh)
                    sdG = std(SG, 0, 2, 'omitnan'); sdG(sdG==0) = 1;
                    SGz = (SG - mean(SG, 2, 'omitnan')) ./ sdG;
                    muGz = mean(SGz, 1, 'omitnan');
                    muGz = muGz - mean(muGz(blStart:blEnd), 'omitnan');
                    aucGFP_zscore{b,r} = [aucGFP_zscore{b,r}; ...
                        trapz(tEvt(aucIdx{b}), muGz(aucIdx{b}))];
                end
            end
        end
    end
end

%% =========================================================================
%  FME comparability across cohorts (R1 minor, Fig 1F-I)
%  "Since different regions are recorded in different mice, it would ensure
%   that the dynamics of FME for these events are comparable in the data used
%   for different regions and in ACh vs. control."
%
%  The reviewer's question is about DYNAMICS - whether a brief facial twitch
%  is the same kind of event across mice - not about how large the movements
%  are. Absolute FME amplitude is not comparable across animals anyway: it
%  depends on camera position, ROI placement and lighting, none of which are
%  standardised between mice. Group differences in amplitude are therefore
%  expected and uninformative, and they do not track the regional ACh
%  differences.
%
%  The primary panel shows Z-SCORED traces, using the same per-snippet
%  convention applied to the ACh data ((M - rowMean) ./ rowSD over the full
%  -1 to +4 s window, then the pre-bout baseline re-zeroed), so FME and ACh
%  can be read in comparable units and both sit at zero before bout onset.
%  Raw-amplitude traces are shown alongside.
%
%  NOTE ON UNITS: the shape metrics (time to peak, FWHM) and the FME AUC are
%  all computed from the RAW baseline-subtracted traces, not from the z-scored
%  display. Time to peak and FWHM are invariant to the linear rescaling that
%  z-scoring applies, so this does not affect them; AUC is reported in raw
%  units so it remains interpretable as integrated motion energy.
%
%  Both cohorts are baseline-subtracted here in exactly the same way; the
%  stored snippets are raw in both paths.
%% =========================================================================
grpNames  = [regionNames, {'GFP'}];
grpColors = [regionColors; 0.45 0.45 0.45];
nGrp      = numel(grpNames);

hasFmeGrp = false;
for b = 1:2
    if ~isempty(regionMask_fme{b}); hasFmeGrp = true; end
end

if hasFmeGrp
    figure('Units','normalized','Position',[0.05 0.14 0.80 0.72]);
    sgtitle('FME dynamics for the events used in Fig 1F-I', ...
        'FontSize', 14, 'FontWeight', 'bold');

    fmeStatRows = {};

    for b = 1:2   % short (twitch) and long (whisk)
        F  = traceData_zscore{b,1};
        RM = regionMask_fme{b};

        % ---- gather each cohort's baseline-subtracted mean traces ----
        grpTrace = cell(nGrp,1);
        for g = 1:nGrp
            if g <= nR
                if isempty(F) || isempty(RM); continue; end
                nUse = min(size(F,1), size(RM,1));
                M = F(1:nUse,:);
                M = M(RM(1:nUse, g), :);
            else
                M = fmeGFP{b};
            end
            if isempty(M); continue; end
            M = M(all(isfinite(M),2),:);
            if size(M,1) < 2; continue; end
            M = M - mean(M(:, blStart:blEnd), 2, 'omitnan');   % identical for all
            grpTrace{g} = M;
        end

        % ---- Row 1: z-scored, matching how ACh is plotted ----
        % Same convention as the ACh snippets: each trace is z-scored using
        % its own mean and SD across the full -1 to +4 s window
        % (achSnip_z = (M - mu_r) ./ sd_r), so FME and ACh are directly
        % comparable in units.
        ax = subplot(2, 2, b);
        hold(ax, 'on');
        for g = 1:nGrp
            M = grpTrace{g};
            if isempty(M); continue; end
            muRow = mean(M, 2, 'omitnan');
            sdRow = std(M,  0, 2, 'omitnan');
            sdRow(sdRow == 0) = NaN;
            Mn = (M - muRow) ./ sdRow;
            % Re-zero the pre-bout baseline AFTER z-scoring, as the ACh path
            % does (muZ = muZ - mean(muZ(blStart:blEnd))). Z-scoring subtracts
            % the whole-window mean, which includes the post-bout rise, so
            % without this step the baseline sits below zero.
            Mn = Mn - mean(Mn(:, blStart:blEnd), 2, 'omitnan');
            Mn = Mn(all(isfinite(Mn),2),:);
            if size(Mn,1) < 2; continue; end
            mu  = mean(Mn,1,'omitnan');
            sem = std(Mn,0,1,'omitnan')/sqrt(size(Mn,1));
            fill(ax, [tEvt, fliplr(tEvt)], [mu+sem, fliplr(mu-sem)], ...
                grpColors(g,:), 'FaceAlpha', 0.16, 'EdgeColor','none');
            plot(ax, tEvt, mu, '-', 'Color', grpColors(g,:), 'LineWidth', 2.0);
        end
        xline(ax, 0, '--k', 'Alpha', 0.4, 'LineWidth', 1.0);
        yline(ax, 0, '-', 'Color', [0.85 0.85 0.85], 'LineWidth', 0.8);
        xlim(ax, [tEvt(1) tEvt(end)]);
        if b == 1; ylabel(ax, 'FME (z-score)', 'FontSize', 12); end
        title(ax, sprintf('%s — z-scored', boutTitleStr{b}), ...
            'FontSize', 12, 'FontWeight','bold');
        box off; set(ax, 'TickDir','out', 'FontSize', 11);
        if b == 2
            lh = gobjects(nGrp,1);
            for g = 1:nGrp
                lh(g) = plot(ax, nan, nan, '-', 'Color', grpColors(g,:), 'LineWidth', 2);
            end
            legend(ax, lh, grpNames, 'FontSize', 10, 'Box','off', 'Location','northeast');
        end

        % ---- Row 2: raw amplitude (secondary) ----
        ax2 = subplot(2, 2, 2+b);
        hold(ax2, 'on');
        for g = 1:nGrp
            M = grpTrace{g};
            if isempty(M); continue; end
            mu  = mean(M,1,'omitnan');
            sem = std(M,0,1,'omitnan')/sqrt(size(M,1));
            fill(ax2, [tEvt, fliplr(tEvt)], [mu+sem, fliplr(mu-sem)], ...
                grpColors(g,:), 'FaceAlpha', 0.16, 'EdgeColor','none');
            plot(ax2, tEvt, mu, '-', 'Color', grpColors(g,:), 'LineWidth', 2.0);
        end
        xline(ax2, 0, '--k', 'Alpha', 0.4, 'LineWidth', 1.0);
        yline(ax2, 0, '-', 'Color', [0.85 0.85 0.85], 'LineWidth', 0.8);
        xlim(ax2, [tEvt(1) tEvt(end)]);
        xlabel(ax2, 'Time re bout onset (s)', 'FontSize', 12);
        if b == 1; ylabel(ax2, 'FME (a.u., raw)', 'FontSize', 12); end
        title(ax2, sprintf('%s — amplitude', boutTitleStr{b}), ...
            'FontSize', 11, 'FontWeight','normal');
        box off; set(ax2, 'TickDir','out', 'FontSize', 11);

        % ---- shape metrics + amplitude, per cohort ----
        switch b
            case 1; wRange = aucWinShort;
            case 2; wRange = aucWinLong;
            otherwise; wRange = aucWinMid;
        end
        aucW  = tEvt >= wRange(1) & tEvt <= wRange(2);
        postM = tEvt >= 0;
        tPost = tEvt(postM);

        for g = 1:nGrp
            M = grpTrace{g};
            if isempty(M); continue; end
            nG = size(M,1);
            ttp = nan(nG,1); fwhm = nan(nG,1);
            for k = 1:nG
                seg = M(k, postM);
                [pkV, pkI] = max(seg);
                if ~isfinite(pkV) || pkV <= 0; continue; end
                ttp(k) = tPost(pkI);
                above  = seg >= pkV/2;
                iL = find(above(1:pkI), 1, 'first');
                iR = pkI - 1 + find(~above(pkI:end), 1, 'first');
                if isempty(iR); iR = numel(seg); end
                if ~isempty(iL); fwhm(k) = tPost(iR) - tPost(iL); end
            end
            aucV = trapz(tEvt(aucW), M(:, aucW), 2);
            fmeStatRows{end+1} = {boutTitleStr{b}, grpNames{g}, nG, ...
                aucV, ttp, fwhm}; %#ok
        end
    end

    % ---- Report: shape first, amplitude second ----
    fprintf('\n============================================================\n');
    fprintf('FME comparability across cohorts (events used in Fig 1F-I)\n');
    fprintf('PRIMARY: shape metrics. Amplitude reported but not comparable\n');
    fprintf('across animals (camera/ROI placement is not standardised).\n');
    fprintf('Kruskal-Wallis across cohorts.\n');
    fprintf('============================================================\n');
    for b = 1:2
        sel = find(cellfun(@(x) strcmp(x{1}, boutTitleStr{b}), fmeStatRows));
        if numel(sel) < 2; continue; end
        fprintf('\n%s:\n', boutTitleStr{b});
        fprintf('  %-5s %5s  %-22s %-22s %s\n', ...
            'Grp','n','time-to-peak (s)','FWHM (s)','FME AUC (a.u.*s)');
        vT = []; vF = []; vA = []; gI = [];
        for k = 1:numel(sel)
            rw = fmeStatRows{sel(k)};
            t = rw{5}(isfinite(rw{5}));  f = rw{6}(isfinite(rw{6}));  a = rw{4}(isfinite(rw{4}));
            fprintf('  %-5s %5d  %8.3f +/- %-9.3f %8.3f +/- %-9.3f %8.2f +/- %.2f\n', ...
                rw{2}, rw{3}, ...
                mean(t,'omitnan'), std(t,'omitnan')/max(sqrt(numel(t)),1), ...
                mean(f,'omitnan'), std(f,'omitnan')/max(sqrt(numel(f)),1), ...
                mean(a,'omitnan'), std(a,'omitnan')/max(sqrt(numel(a)),1));
            vT = [vT; t(:)]; vF = [vF; f(:)]; vA = [vA; a(:)]; %#ok
            gI = [gI; repmat(k, numel(a), 1)]; %#ok
        end
        if numel(unique(gI)) >= 2
            pA = kruskalwallis(vA, gI, 'off');
            nT = min([numel(vT), numel(gI)]);
            nF = min([numel(vF), numel(gI)]);
            pT = NaN; pF = NaN;
            if nT >= 6; pT = kruskalwallis(vT(1:nT), gI(1:nT), 'off'); end
            if nF >= 6; pF = kruskalwallis(vF(1:nF), gI(1:nF), 'off'); end
            if isfinite(pT) && pT > 0.05 && isfinite(pF) && pF > 0.05
                shapeNote = '   <- dynamics comparable across cohorts';
            else
                shapeNote = '';
            end
            fprintf('  Kruskal-Wallis  time-to-peak p = %.4f   FWHM p = %.4f%s\n', ...
                pT, pF, shapeNote);
            fprintf('                  amplitude    p = %.4f   (see caveat above)\n', pA);
        end
    end
    fprintf('\n============================================================\n');
end

%% =========================================================================
%  Bout duration distribution (addresses R3 minor: "what about 1-2 s events?")
%  Shows every detected bout duration before the short/long selection, with
%  the two cutoffs marked, so the size and position of the excluded
%  intermediate band is explicit. A clearly bimodal distribution would
%  justify the split on its own; a smooth unimodal one means the cutoffs are
%  a deliberate conservative choice rather than a natural boundary, and that
%  should be stated as such.
%% =========================================================================
if ~isempty(allBoutDurs)
    dAll = allBoutDurs(isfinite(allBoutDurs));
    nTot = numel(dAll);

    nShortD = sum(dAll <  shortMaxSec);
    nMidD   = sum(dAll >= shortMaxSec & dAll < longMinSec);
    nLongD  = sum(dAll >= longMinSec);

    figure('Units','normalized','Position',[0.03 0.22 0.94 0.56]);
    sgtitle('FME Bout Duration Distribution', 'FontSize', 20, 'FontWeight','bold');

    % --- Linear-scale histogram ---
    ax1 = subplot(1,3,1); hold(ax1,'on');
    edgesD = 0 : 0.1 : min(max(dAll), 8);
    histogram(ax1, dAll, edgesD, 'FaceColor', [0.45 0.45 0.45], ...
        'EdgeColor','none', 'FaceAlpha', 0.85);
    yl = ylim(ax1);
    patch(ax1, [shortMaxSec longMinSec longMinSec shortMaxSec], ...
        [yl(1) yl(1) yl(2) yl(2)], [0.90 0.35 0.35], ...
        'FaceAlpha', 0.16, 'EdgeColor','none');
    xline(ax1, shortMaxSec, '--k', 'LineWidth', 1.8);
    xline(ax1, longMinSec,  '--k', 'LineWidth', 1.8);
    xlabel(ax1, 'Bout duration (s)', 'FontSize', 16);
    ylabel(ax1, 'Count', 'FontSize', 16);
    title(ax1, 'Linear scale', 'FontSize', 17, 'FontWeight','bold');
    box off; set(ax1,'TickDir','out','FontSize',14,'LineWidth',1.1);

    % --- Log-spaced histogram ---
    % Movement bout durations are typically heavy-tailed; log binning shows
    % whether the distribution is genuinely bimodal (which would justify the
    % short/long split) or smoothly unimodal (cutoffs are a conservative
    % choice, not a natural boundary).
    axL = subplot(1,3,2); hold(axL,'on');
    dPos = dAll(dAll > 0);
    if numel(dPos) > 10
        edgesLog = logspace(log10(min(dPos)), log10(max(dPos)), 40);
        histogram(axL, dPos, edgesLog, 'FaceColor', [0.45 0.45 0.45], ...
            'EdgeColor','none', 'FaceAlpha', 0.85);
        set(axL, 'XScale', 'log');
        ylL = ylim(axL);
        patch(axL, [shortMaxSec longMinSec longMinSec shortMaxSec], ...
            [ylL(1) ylL(1) ylL(2) ylL(2)], [0.90 0.35 0.35], ...
            'FaceAlpha', 0.16, 'EdgeColor','none');
        xline(axL, shortMaxSec, '--k', 'LineWidth', 1.8);
        xline(axL, longMinSec,  '--k', 'LineWidth', 1.8);
        xlim(axL, [min(dPos) max(dPos)]);
    end
    xlabel(axL, 'Bout duration (s), log scale', 'FontSize', 16);
    ylabel(axL, 'Count', 'FontSize', 16);
    title(axL, 'Log scale', 'FontSize', 17, 'FontWeight','bold');
    box off; set(axL,'TickDir','out','FontSize',14,'LineWidth',1.1);

    % --- Cumulative distribution ---
    ax2 = subplot(1,3,3); hold(ax2,'on');
    ds = sort(dAll);
    plot(ax2, ds, (1:numel(ds))/numel(ds), '-k', 'LineWidth', 2.8);
    yl2 = [0 1];
    patch(ax2, [shortMaxSec longMinSec longMinSec shortMaxSec], ...
        [yl2(1) yl2(1) yl2(2) yl2(2)], [0.90 0.35 0.35], ...
        'FaceAlpha', 0.16, 'EdgeColor','none');
    xline(ax2, shortMaxSec, '--k', 'LineWidth', 1.8);
    xline(ax2, longMinSec,  '--k', 'LineWidth', 1.8);
    xlim(ax2, [0 min(max(dAll), 8)]); ylim(ax2, yl2);
    xlabel(ax2, 'Bout duration (s)', 'FontSize', 16);
    ylabel(ax2, 'Cumulative probability', 'FontSize', 16);
    title(ax2, 'Cumulative', 'FontSize', 17, 'FontWeight','bold');
    box off; set(ax2,'TickDir','out','FontSize',14,'LineWidth',1.1);

    % No footer annotation, so panels can occupy the full figure height
    for axH = [ax1, axL, ax2]
        pos = get(axH,'Position');
        set(axH,'Position',[pos(1), pos(2)+0.04, pos(3), pos(4)-0.02]);
    end

    % (Counts and percentages are printed to the command window instead of
    %  being annotated on the figure, which crowded the panels.)


    % --- Command-window summary ---
    fprintf('\n============================================================\n');
    fprintf('FME bout duration distribution (all detected bouts)\n');
    fprintf('============================================================\n');
    fprintf('  Total bouts detected:      %d\n', nTot);
    fprintf('  Short  (< %.1f s):          %5d  (%.1f%%)\n', shortMaxSec, nShortD, 100*nShortD/nTot);
    fprintf('  EXCLUDED (%.1f-%.1f s):      %5d  (%.1f%%)\n', shortMaxSec, longMinSec, nMidD, 100*nMidD/nTot);
    fprintf('  Long   (>= %.1f s):         %5d  (%.1f%%)\n', longMinSec, nLongD, 100*nLongD/nTot);
    fprintf('\n  Duration percentiles (s):\n');
    for pc = [5 10 25 50 75 90 95 99]
        fprintf('    p%-3d = %.2f\n', pc, prctile(dAll, pc));
    end
    fprintf('\n  Median short bout: %.2f s   Median long bout: %.2f s\n', ...
        median(dAll(dAll <  shortMaxSec)), median(dAll(dAll >= longMinSec)));
    fprintf('============================================================\n');

    % --- CSV of raw durations for further inspection ---
    csvD = fullfile(save_path, sprintf('fme_bout_durations%s.csv', qcTag));
    fidD = fopen(csvD, 'w');
    if fidD == -1
        warning('Could not write bout duration CSV.');
    else
        fprintf(fidD, 'mouse_id,duration_s,category\n');
        for k = 1:nTot
            if dAll(k) < shortMaxSec
                cat = 'short';
            elseif dAll(k) < longMinSec
                cat = 'excluded';
            else
                cat = 'long';
            end
            if k <= numel(allBoutDurMouse)
                mid = allBoutDurMouse{k};
            else
                mid = 'unknown';
            end
            fprintf(fidD, '%s,%.4f,%s\n', mid, dAll(k), cat);
        end
        fclose(fidD);
        fprintf('Bout durations saved to: %s\n', csvD);
    end
end

%% =========================================================================
%  Intermediate (1-2 s) bout band — supplementary summary
%  ------------------------------------------------------------------------
%  Bout durations form a smooth heavy-tailed continuum, so the 1 s and 2 s
%  criteria isolate unambiguous exemplars rather than dividing natural
%  classes. This section reports the intermediate band so its behaviour is
%  documented, without asserting that it differs in kind from either
%  neighbour. AUC window for this band is 0-2 s (aucWinMid).
%% =========================================================================
bMid = 3;   % index of the mid band in the accumulated cells

hasMidData = ~isempty(traceData_zscore{bMid,1}) || ...
             any(~cellfun(@isempty, aucData_zscore(bMid,:)));

if hasMidData
    figure('Units','normalized','Position',[0.06 0.18 0.88 0.60]);
    sgtitle(sprintf('Intermediate FME bouts (%.0f-%.0f s): traces and AUC', ...
        shortMaxSec, longMinSec), 'FontSize', 13, 'FontWeight','bold');

    tEvtPlot = (-preS:postS) / fsAch;
    sigNamesMid = [{'FME','Pupil'}, regionNames];

    % --- Traces: FME, pupil, then each ACh region ---
    nPanels = 2 + nR;
    for sg = 1:nPanels
        ax = subplot(2, nPanels, sg);
        hold(ax,'on');
        M = traceData_zscore{bMid, sg};
        if isempty(M)
            title(ax, sigNamesMid{sg}, 'FontSize', 10); axis(ax,'off'); continue;
        end
        M = M(all(isfinite(M),2),:);
        if size(M,1) < 2
            title(ax, sigNamesMid{sg}, 'FontSize', 10); continue;
        end
        mu  = mean(M,1,'omitnan');
        sem = std(M,0,1,'omitnan') / sqrt(size(M,1));
        if sg <= 2
            col = [0.35 0.35 0.35];
        else
            col = regionColors(sg-2,:);
        end
        fill(ax, [tEvtPlot, fliplr(tEvtPlot)], [mu+sem, fliplr(mu-sem)], ...
            col, 'FaceAlpha', 0.20, 'EdgeColor','none');
        plot(ax, tEvtPlot, mu, '-', 'Color', col, 'LineWidth', 1.8);
        xline(ax, 0, '--k','Alpha',0.4,'LineWidth',0.8);
        yline(ax, 0, '-','Color',[0.85 0.85 0.85],'LineWidth',0.7);
        xlim(ax, [tEvtPlot(1) tEvtPlot(end)]);
        title(ax, sprintf('%s (n=%d)', sigNamesMid{sg}, size(M,1)), 'FontSize', 10);
        if sg == 1; ylabel(ax, 'z-score', 'FontSize', 9); end
        xlabel(ax, 'Time re bout onset (s)', 'FontSize', 8);
        box off; set(ax,'TickDir','out','FontSize',8);
    end

    % --- AUC per region, with one-sided Wilcoxon vs zero ---
    axA = subplot(2, nPanels, [nPanels+1, 2*nPanels]);
    hold(axA,'on');
    midRows = {};
    for r = 1:nR
        v = aucData_zscore{bMid, r};
        v = v(isfinite(v));
        if numel(v) < 3; continue; end
        muV  = mean(v,'omitnan');
        semV = std(v,0,'omitnan')/sqrt(numel(v));
        pV   = signrank(v, 0, 'tail', 'right');
        bar(axA, r, muV, 0.55, 'FaceColor', regionColors(r,:), ...
            'EdgeColor','none','FaceAlpha',0.75);
        errorbar(axA, r, muV, semV, 'k','LineWidth',1.4,'CapSize',5,'LineStyle','none');
        midRows{end+1} = {regionNames{r}, numel(v), muV, semV, ...
            median(v,'omitnan'), pV}; %#ok
    end
    yline(axA, 0, '-','Color',[0.75 0.75 0.75],'LineWidth',0.8);
    set(axA,'XTick',1:nR,'XTickLabel',regionNames,'TickDir','out','FontSize',9);
    xlim(axA,[0.4 nR+0.6]);
    ylabel(axA, sprintf('ACh AUC %.0f-%.0f s (z*s)', aucWinMid(1), aucWinMid(2)), ...
        'FontSize', 9);
    title(axA, 'ACh AUC, intermediate bouts', 'FontSize', 10, 'FontWeight','bold');
    box off;

    % FDR across regions, then annotate
    if ~isempty(midRows)
        pRawMid = cellfun(@(x) x{6}, midRows);
        nMidT   = numel(pRawMid);
        [psM, idM] = sort(pRawMid(:));
        adjM = psM .* nMidT ./ (1:nMidT)';
        adjM = flipud(cummin(flipud(adjM)));
        adjM = min(adjM, 1);
        pAdjMid = nan(nMidT,1); pAdjMid(idM) = adjM;

        fprintf('\n============================================================\n');
        fprintf('Intermediate bouts (%.0f-%.0f s): ACh AUC vs zero\n', ...
            shortMaxSec, longMinSec);
        fprintf('One-sided Wilcoxon, FDR-corrected across regions.\n');
        fprintf('Pooling mode: %s   AUC window: %.0f-%.0f s\n', ...
            poolingMode, aucWinMid(1), aucWinMid(2));
        fprintf('============================================================\n');
        for k = 1:numel(midRows)
            rw = midRows{k};
            fprintf('  %-5s n=%3d  mean=%+.4f +/- %.4f  median=%+.4f  p_raw=%.4f  p_FDR=%.4f%s\n', ...
                rw{1}, rw{2}, rw{3}, rw{4}, rw{5}, rw{6}, pAdjMid(k), ...
                char(' ' + (pAdjMid(k)<0.05)*('*'-' ')));
        end
        fprintf('============================================================\n');

        csvM = fullfile(save_path, sprintf('ach_midbouts_%s%s.csv', poolingMode, qcTag));
        fidM = fopen(csvM, 'w');
        if fidM == -1
            warning('Could not write intermediate-bout CSV.');
        else
            fprintf(fidM, ['bout_band_s,region,n_obs,mean_AUC,sem_AUC,median_AUC,' ...
                           'auc_win_start_s,auc_win_end_s,p_raw_onesided,p_FDR,sig_FDR05\n']);
            for k = 1:numel(midRows)
                rw = midRows{k};
                fprintf(fidM, '%.1f-%.1f,%s,%d,%.6f,%.6f,%.6f,%.1f,%.1f,%.6f,%.6f,%d\n', ...
                    shortMaxSec, longMinSec, rw{1}, rw{2}, rw{3}, rw{4}, rw{5}, ...
                    aucWinMid(1), aucWinMid(2), rw{6}, pAdjMid(k), pAdjMid(k)<0.05);
            end
            fclose(fidM);
            fprintf('Intermediate-bout stats saved to: %s\n', csvM);
        end
    end

    annotation('textbox',[0.01 0 0.99 0.04], ...
        'String', sprintf(['Intermediate-duration bouts (%.0f-%.0f s), excluded from the ' ...
            'main short/long comparison.  Reported for completeness: bout durations form a ' ...
            'continuum, and this band is shown to document its behaviour rather than to ' ...
            'claim a distinct category.  AUC window %.0f-%.0f s.'], ...
            shortMaxSec, longMinSec, aucWinMid(1), aucWinMid(2)), ...
        'EdgeColor','none','FontSize',8,'Color',[0.4 0.4 0.4], ...
        'HorizontalAlignment','center');
end

%% -------------------------------------------------------------------------
%  Summary
%% -------------------------------------------------------------------------
fprintf('\n--- Summary (%s) ---\n', poolingMode);
sigNames = [{'FME','Pupil'}, regionNames];
fprintf('%-10s  %8s  %8s\n', 'Signal', 'Short', 'Long');
fprintf('%s\n', repmat('-', 30, 1));
for s = 1:nSigs
    fprintf('%-10s  %8d  %8d\n', sigNames{s}, ...
        size(traceData_raw{1,s}, 1), size(traceData_raw{2,s}, 1));
end

%% =========================================================================
%  Generate and save figures
%% =========================================================================

% Trace figures
fig1 = make_trace_figure(traceData_raw,    refData_raw,    '\DeltaF/F');
fig2 = make_trace_figure(traceData_zscore, refData_zscore, 'z-score');

% AUC figures
[fig3, statRows_dff]    = make_auc_figure(aucData_raw,    '\DeltaF/F \times s', '\DeltaF/F');
[fig4, statRows_zscore] = make_auc_figure(aucData_zscore, 'z \times s',         'z-score');

% ------------------------------------------------------------------
%  Above-zero test: one-sided Wilcoxon signed-rank per region,
%  per bout type.  Unit = per-session mean AUC (one value per
%  session, as stored in aucData when poolingMode = 'per_session').
%  FDR correction (Benjamini-Hochberg) across all region × bout
%  combinations within each signal type.
% ------------------------------------------------------------------
boutTypeLabels = {sprintf('Short bouts (<%g s)', shortMaxSec), ...
                  sprintf('Long bouts (>=%g s)', longMinSec), ...
                  sprintf('Mid bouts (%g-%g s)', shortMaxSec, longMinSec)};

for sigType = 1:2   % 1 = raw DF/F, 2 = z-score
    if sigType == 1
        aucD    = aucData_raw;
        unitStr = 'dFF';
        csvSuffix = 'dff';
    else
        aucD    = aucData_zscore;
        unitStr = 'zscore';
        csvSuffix = 'zscore';
    end

    % Collect raw p-values across all bout × region combinations
    aboveZeroRows = {};
    pRaws = [];
    for b = 1:nBouts
        for r = 1:nR
            v = aucD{b,r};
            v = v(isfinite(v));
            if numel(v) < 3
                p = NaN;
            else
                p = signrank(v, 0, 'tail', 'right');
            end
            pRaws(end+1) = p; %#ok
            aboveZeroRows{end+1} = {boutTypeLabels{b}, unitStr, ...
                regionNames{r}, numel(v), ...
                mean(v,'omitnan'), median(v,'omitnan'), ...
                std(v,0,'omitnan'), p}; %#ok
        end
    end

    % FDR correction (Benjamini-Hochberg)
    nT        = numel(pRaws);
    finIdx    = find(isfinite(pRaws));
    pFin      = pRaws(finIdx);
    [ps, idx] = sort(pFin(:));
    adj       = ps .* numel(pFin) ./ (1:numel(pFin))';
    adj       = flipud(cummin(flipud(adj)));
    adj       = min(adj, 1);
    pAdjFin   = nan(size(pFin)); pAdjFin(idx) = adj;
    pAdj      = nan(nT,1); pAdj(finIdx) = pAdjFin;

    % ---------------------------------------------------------------
    %  CONTROLS: same above-zero test for 405 nm and GFP
    %  The experimental cohort is tested per brain-region group; the
    %  controls must be shown the same way or a reader cannot tell whether
    %  the regional pattern is specific to the ACh signal.
    % ---------------------------------------------------------------
    if strcmpi(unitStr,'dFF')
        refD = refData_raw;    gfpD = aucGFP_raw;
    else
        refD = refData_zscore; gfpD = aucGFP_zscore;
    end

    ctrlRows = {}; pCtrl = [];
    for b = 1:nBouts
        for r = 1:nR
            % 405 nm — AUC from the stored per-mouse reference traces
            vR = [];
            if ~isempty(refD{b,r})
                M = refD{b,r};
                M = M(all(isfinite(M),2),:);
                if ~isempty(M)
                    vR = trapz(tEvt(aucIdx{b}), M(:, aucIdx{b}), 2);
                end
            end
            % GFP — AUC computed during the GFP pass
            vG = gfpD{b,r}; vG = vG(isfinite(vG));

            for cIdx = 1:2
                if cIdx == 1; v = vR; cName = '405nm';
                else;         v = vG; cName = 'GFP';   end
                v = v(isfinite(v));
                if numel(v) < 3
                    pC = NaN;
                else
                    pC = signrank(v, 0, 'tail', 'right');
                end
                pCtrl(end+1) = pC; %#ok
                ctrlRows{end+1} = {boutTypeLabels{b}, unitStr, cName, ...
                    regionNames{r}, numel(v), ...
                    mean(v,'omitnan'), median(v,'omitnan'), ...
                    std(v,0,'omitnan'), pC}; %#ok
            end
        end
    end

    % FDR within the control family (kept separate from the ACh family so
    % neither correction is diluted by the other)
    pAdjCtrl = nan(numel(pCtrl),1);
    fin = isfinite(pCtrl);
    if any(fin)
        pf = pCtrl(fin); nf = numel(pf);
        [ps2, ix2] = sort(pf(:));
        a2 = ps2 .* nf ./ (1:nf)';
        a2 = flipud(cummin(flipud(a2)));
        tmp = nan(nf,1); tmp(ix2) = min(a2,1);
        pAdjCtrl(fin) = tmp;
    end

    csvC = fullfile(save_path, sprintf('control_above_zero_stats_%s_%s%s.csv', ...
        csvSuffix, poolingMode, qcTag));
    fidC = fopen(csvC, 'w');
    if fidC ~= -1
        fprintf(fidC, ['cohort,test,test_tail,comparison,' ...
                       'multiple_comparison_correction,correction_family,' ...
                       'bout_type,signal_unit,region,pooling_mode,' ...
                       'n_obs,mean,median,sd,p_raw,p_FDR,' ...
                       'significant_FDR05\n']);
        for k = 1:numel(ctrlRows)
            rw = ctrlRows{k};
            fprintf(fidC, ['%s,Wilcoxon signed-rank,one-sided (greater),vs zero,' ...
                           'Benjamini-Hochberg FDR,control above-zero (%d tests),' ...
                           '%s,%s,%s,%s,%d,%.6f,%.6f,%.6f,%.6f,%.6f,%d\n'], ...
                rw{3}, numel(ctrlRows), ...
                rw{1}, rw{2}, rw{4}, poolingMode, rw{5}, ...
                rw{6}, rw{7}, rw{8}, rw{9}, pAdjCtrl(k), pAdjCtrl(k)<0.05);
        end
        fclose(fidC);
        fprintf('Control above-zero stats saved: %s\n', csvC);
    end

    % ---------------------------------------------------------------
    %  ACh vs CONTROL: the comparison that shows the signal is ACh-specific
    %  Unpaired Wilcoxon rank-sum — different cohorts of mice, so unpaired.
    %  Two separate null results do not by themselves establish a difference;
    %  this tests it directly.
    % ---------------------------------------------------------------
    cmpRows = {}; pCmp = [];
    for b = 1:nBouts
        for r = 1:nR
            vA = aucD{b,r}; vA = vA(isfinite(vA));

            vR = [];
            if ~isempty(refD{b,r})
                M = refD{b,r}; M = M(all(isfinite(M),2),:);
                if ~isempty(M); vR = trapz(tEvt(aucIdx{b}), M(:, aucIdx{b}), 2); end
            end
            vG = gfpD{b,r}; vG = vG(isfinite(vG));

            for cIdx = 1:2
                if cIdx == 1; vC = vR(isfinite(vR)); cName = '405nm';
                else;         vC = vG;               cName = 'GFP';   end
                if numel(vA) < 3 || numel(vC) < 3
                    pV = NaN;
                else
                    pV = ranksum(vA, vC);
                end
                pCmp(end+1) = pV; %#ok
                cmpRows{end+1} = {boutTypeLabels{b}, unitStr, regionNames{r}, ...
                    cName, numel(vA), numel(vC), ...
                    mean(vA,'omitnan'), mean(vC,'omitnan'), ...
                    mean(vA,'omitnan') - mean(vC,'omitnan'), pV}; %#ok
            end
        end
    end

    pAdjCmp = nan(numel(pCmp),1);
    fin = isfinite(pCmp);
    if any(fin)
        pf = pCmp(fin); nf = numel(pf);
        [ps3, ix3] = sort(pf(:));
        a3 = ps3 .* nf ./ (1:nf)';
        a3 = flipud(cummin(flipud(a3)));
        tmp = nan(nf,1); tmp(ix3) = min(a3,1);
        pAdjCmp(fin) = tmp;
    end

    csvX = fullfile(save_path, sprintf('ach_vs_control_stats_%s_%s%s.csv', ...
        csvSuffix, poolingMode, qcTag));
    fidX = fopen(csvX, 'w');
    if fidX ~= -1
        fprintf(fidX, ['test,test_tail,comparison,' ...
                       'multiple_comparison_correction,correction_family,' ...
                       'bout_type,signal_unit,region,control_cohort,' ...
                       'n_ACh,n_control,mean_ACh,mean_control,mean_diff,' ...
                       'p_raw,p_FDR,significant_FDR05\n']);
        for k = 1:numel(cmpRows)
            rw = cmpRows{k};
            fprintf(fidX, ['Wilcoxon rank-sum (unpaired),two-sided,ACh vs %s,' ...
                           'Benjamini-Hochberg FDR,ACh vs control (%d tests),' ...
                           '%s,%s,%s,%s,%d,%d,%.6f,%.6f,%.6f,%.6f,%.6f,%d\n'], ...
                rw{4}, numel(cmpRows), ...
                rw{1}, rw{2}, rw{3}, rw{4}, rw{5}, rw{6}, ...
                rw{7}, rw{8}, rw{9}, rw{10}, pAdjCmp(k), pAdjCmp(k)<0.05);
        end
        fclose(fidX);
        fprintf('ACh vs control stats saved: %s\n', csvX);
    end

    % Write CSV
    csvPath = fullfile(save_path, sprintf('ach_above_zero_stats_%s_%s%s.csv', csvSuffix, poolingMode, qcTag));
    fid = fopen(csvPath, 'w');
    fprintf(fid, ['cohort,test,test_tail,comparison,multiple_comparison_correction,' ...
                  'correction_family,bout_type,signal_unit,region,pooling_mode,' ...
                  'n_obs,mean,median,sd,p_raw,p_FDR,significant_FDR05\n']);
    for k = 1:numel(aboveZeroRows)
        row = aboveZeroRows{k};
        fprintf(fid, ['ACh,Wilcoxon signed-rank,one-sided (greater),vs zero,' ...
                      'Benjamini-Hochberg FDR,ACh above-zero (%d tests),' ...
                      '%s,%s,%s,%s,%d,%.6f,%.6f,%.6f,%.6f,%.6f,%d\n'], ...
            numel(aboveZeroRows), ...
            row{1}, row{2}, row{3}, poolingMode, row{4}, row{5}, row{6}, row{7}, ...
            row{8}, pAdj(k), pAdj(k)<0.05);
    end
    fclose(fid);
    fprintf('Above-zero stats saved: %s\n', csvPath);
end

% Save figures
save_fig(fig1, fullfile(save_path, sprintf('ach_aligned_dff%s.pdf',    qcTag)));
save_fig(fig2, fullfile(save_path, sprintf('ach_aligned_zscore%s.pdf', qcTag)));
save_fig(fig3, fullfile(save_path, sprintf('ach_auc_dff%s.pdf',        qcTag)));
save_fig(fig4, fullfile(save_path, sprintf('ach_auc_zscore%s.pdf',     qcTag)));

% Save statistics tables
save_stats_table(statRows_dff, ...
    fullfile(save_path, sprintf('ach_auc_stats_dff%s.csv',    qcTag)));
save_stats_table(statRows_zscore, ...
    fullfile(save_path, sprintf('ach_auc_stats_zscore%s.csv', qcTag)));

fprintf('\nFigures and stats saved to %s\n', save_path);

%% =========================================================================
%  Nested helpers
%% =========================================================================

    % -----------------------------------------------------------------
    function fig = make_trace_figure(traceData, refData, yUnitStr)
    % 2×2 layout: FME | Pupil / ACh short | ACh long
    % ACh panels: solid = 470 nm signal, dashed = 405 nm ref

        nBoot = 1000;
        rng(42);    % permutation p-values must not change between runs

        fig = figure('Units','inches','Position',[1 1 FIG_W FIG_H], ...
            'Color','w','PaperUnits','inches','PaperSize',[FIG_W FIG_H]);

        panelLabels = {'a','b','c','d'};
        axH = gobjects(2,2);

        % --- FME (a) ---
        axH(1,1) = subplot(2,2,1);
        ax = axH(1,1);
        hold(ax,'on');
        legH = []; legStr = {};
        for b = 1:nPublishedBouts
            M = traceData{b,1};
            if isempty(M) || size(M,1) < 2; continue; end
            M  = M(all(isfinite(M),2),:);
            mu = mean(M, 1, 'omitnan');
            ci = bootstrap_ci(M, nBoot);
            plot_with_ci(ax, tEvt, mu, ci, boutColors{b}, LW_SIG, CI_ALPHA);
            legH(end+1)   = plot(ax, nan, nan, '-', 'Color', boutColors{b}, 'LineWidth', LW_SIG);
            legStr{end+1} = sprintf('%s (n=%d)', bLabels{b}, size(M,1));
        end
        format_trace_ax(ax, tEvt, 'FME (a.u.)', false, FONT_AX, FONT_LB);
        if ~isempty(legH)
            legend(ax, legH, legStr, 'Location','northwest', ...
                'FontSize',FONT_AX,'Box','off');
        end
        title(ax,'Facial motion energy','FontSize',FONT_TI,'FontWeight','bold');

        % --- Pupil (b) ---
        axH(1,2) = subplot(2,2,2);
        ax = axH(1,2);
        hold(ax,'on');
        allPupVals = [];
        for b = 1:nPublishedBouts
            M = traceData{b,2};
            if isempty(M) || size(M,1) < 2; continue; end
            M  = M(all(isfinite(M),2),:);
            ci = bootstrap_ci(M, nBoot);
            plot_with_ci(ax, tEvt, mean(M,1,'omitnan'), ci, boutColors{b}, LW_SIG, CI_ALPHA);
            allPupVals = [allPupVals; ci(2,:)'; ci(1,:)']; %#ok
        end
        format_trace_ax(ax, tEvt, 'Pupil diameter (norm.)', false, FONT_AX, FONT_LB);
        if ~isempty(allPupVals) && range(allPupVals) > 0
            yPad = 0.15 * range(allPupVals);
            ylim(ax, [min(allPupVals)-yPad, max(allPupVals)+yPad]);
        end
        title(ax,'Pupil diameter','FontSize',FONT_TI,'FontWeight','bold');

        % --- ACh short (c) ---
        axH(2,1) = subplot(2,2,3);
        ax = axH(2,1);
        hold(ax,'on');
        legH = []; legStr = {};
        for r = 1:nR
            M = traceData{1, 2+r};
            if isempty(M) || size(M,1) < 2; continue; end
            M   = M(all(isfinite(M),2),:);
            mu  = mean(M, 1, 'omitnan');
            ci  = bootstrap_ci(M, nBoot);
            col = regionColors(r,:);
            plot_with_ci(ax, tEvt, mu, ci, col, LW_SIG, CI_ALPHA);
            h = plot(ax, nan, nan, '-', 'Color', col, 'LineWidth', LW_SIG);

            % Ref overlay (dashed, same color, no CI)
            R = refData{1,r};
            if ~isempty(R) && size(R,1) >= 2
                R   = R(all(isfinite(R),2),:);
                muR = mean(R, 1, 'omitnan');
                plot(ax, tEvt, muR, '--', 'Color', col, 'LineWidth', LW_REF);
            end

            legH(end+1)   = h; %#ok
            legStr{end+1} = sprintf('%s (n=%d)', regionNames{r}, size(M,1)); %#ok
        end
        % Add ref legend entry once
        if nRefFound > 0
            legH(end+1)   = plot(ax, nan, nan, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', LW_REF);
            legStr{end+1} = '405 nm ref';
        end
        format_trace_ax(ax, tEvt, ['ACh (' yUnitStr ')'], true, FONT_AX, FONT_LB);
        if ~isempty(legH)
            legend(ax, legH, legStr, 'Location','best', ...
                'FontSize',FONT_AX,'Box','off');
        end
        title(ax, boutTitleStr{1}, 'FontSize',FONT_TI,'FontWeight','bold', ...
            'Color',boutColors{1});

        % --- ACh long (d) ---
        axH(2,2) = subplot(2,2,4);
        ax = axH(2,2);
        hold(ax,'on');
        legH = []; legStr = {};
        for r = 1:nR
            M = traceData{2, 2+r};
            if isempty(M) || size(M,1) < 2; continue; end
            M   = M(all(isfinite(M),2),:);
            mu  = mean(M, 1, 'omitnan');
            ci  = bootstrap_ci(M, nBoot);
            col = regionColors(r,:);
            plot_with_ci(ax, tEvt, mu, ci, col, LW_SIG, CI_ALPHA);
            h = plot(ax, nan, nan, '-', 'Color', col, 'LineWidth', LW_SIG);

            % Ref overlay
            R = refData{2,r};
            if ~isempty(R) && size(R,1) >= 2
                R   = R(all(isfinite(R),2),:);
                muR = mean(R, 1, 'omitnan');
                plot(ax, tEvt, muR, '--', 'Color', col, 'LineWidth', LW_REF);
            end

            legH(end+1)   = h; %#ok
            legStr{end+1} = sprintf('%s (n=%d)', regionNames{r}, size(M,1)); %#ok
        end
        if nRefFound > 0
            legH(end+1)   = plot(ax, nan, nan, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', LW_REF);
            legStr{end+1} = '405 nm ref';
        end
        format_trace_ax(ax, tEvt, ['ACh (' yUnitStr ')'], true, FONT_AX, FONT_LB);
        if ~isempty(legH)
            legend(ax, legH, legStr, 'Location','best', ...
                'FontSize',FONT_AX,'Box','off');
        end
        title(ax, boutTitleStr{2}, 'FontSize',FONT_TI,'FontWeight','bold', ...
            'Color',boutColors{2});

        % Panel labels a–d
        axOrder = [axH(1,1), axH(1,2), axH(2,1), axH(2,2)];
        for k = 1:4
            add_panel_label(axOrder(k), panelLabels{k});
        end

        % Tighten layout
        set(fig,'Units','inches');
        tightfig_simple(fig);
    end

    % -----------------------------------------------------------------
    function [fig, statRows] = make_auc_figure(aucData, yLabel, unitStr)
    % 1×2 box plot: short | long bouts
    % Standard box (median, IQR, whiskers = 1.5×IQR)
    % Permutation test (n=5000), Bonferroni corrected across 3 pairs
    % Significant pairs annotated with asterisks (* ** ***)
    % Returns statRows: cell array of table rows for CSV export

        nBoot       = 5000;
        rng(42);    % permutation p-values must not change between runs
        alphaNom    = 0.05;
        pairs_local = [1,2; 1,3; 2,3];
        nPairsLocal = size(pairs_local, 1);
        aucWinLabels = {sprintf('0\x20130.5 s'), sprintf('0\x20134.0 s')};

        figW = FIG_W * 0.75;
        figH = FIG_H * 0.60;
        fig  = figure('Units','inches','Position',[1 1 figW figH], ...
            'Color','w','PaperUnits','inches','PaperSize',[figW figH]);

        panelLabels = {'a','b'};
        statRows    = {};

        for b = 1:2
            ax = subplot(1,2,b);
            hold(ax,'on');

            % ----------------------------------------------------------
            % Draw box plots manually so we can colour each box
            % individually and avoid the boxplot() legend clutter
            % ----------------------------------------------------------
            whiskerTops = nan(nR,1);   % track upper whisker for bracket placement

            for r = 1:nR
                vals = aucData{b,r};
                vals = vals(isfinite(vals));
                if isempty(vals); continue; end

                col = regionColors(r,:);
                n   = numel(vals);

                q1  = prctile(vals, 25);
                med = median(vals);
                q3  = prctile(vals, 75);
                iqr_ = q3 - q1;
                wLow = max(min(vals), q1 - 1.5*iqr_);
                wHigh= min(max(vals), q3 + 1.5*iqr_);
                outliers = vals(vals < wLow | vals > wHigh);

                bw = 0.40;   % box half-width

                % Box body
                patch(ax, [r-bw r+bw r+bw r-bw], [q1 q1 q3 q3], col, ...
                    'FaceAlpha',0.25,'EdgeColor',col,'LineWidth',1.2);

                % Median line
                line(ax, [r-bw, r+bw], [med, med], ...
                    'Color',col,'LineWidth',2.0);

                % Whiskers
                line(ax, [r, r],       [wLow,  q1],  'Color',col,'LineWidth',1.0);
                line(ax, [r, r],       [q3,    wHigh],'Color',col,'LineWidth',1.0);
                line(ax, [r-bw*0.4, r+bw*0.4], [wLow,  wLow], 'Color',col,'LineWidth',1.0);
                line(ax, [r-bw*0.4, r+bw*0.4], [wHigh, wHigh],'Color',col,'LineWidth',1.0);

                % Outliers
                if ~isempty(outliers)
                    scatter(ax, repmat(r, numel(outliers),1), outliers, ...
                        18, col, 'o', 'LineWidth',0.8, ...
                        'MarkerEdgeAlpha',0.7,'MarkerFaceColor','none');
                end

                whiskerTops(r) = wHigh;

                % n label below lower whisker
                yLab = wLow - 0.06 * abs(iqr_ + eps);
                text(ax, r, yLab, sprintf('n=%d', n), ...
                    'HorizontalAlignment','center', ...
                    'FontSize',FONT_AX-1,'Color',col);
            end

            % ----------------------------------------------------------
            % Permutation tests — run all pairs first, then apply
            % Benjamini-Hochberg FDR correction
            % ----------------------------------------------------------

            % Step 1: collect raw p-values and metadata for all valid pairs
            pairResults = struct(...
                'r1',{}, 'r2',{}, 'v1',{}, 'v2',{}, ...
                'obsDiff',{}, 'pVal',{}, 'pPerm',{}, 'pRank',{});

            for pr = 1:nPairsLocal
                r1 = pairs_local(pr,1);  r2 = pairs_local(pr,2);
                v1 = aucData{b,r1}; v1 = v1(isfinite(v1));
                v2 = aucData{b,r2}; v2 = v2(isfinite(v2));
                if numel(v1) < 2 || numel(v2) < 2; continue; end

                obsDiff  = mean(v1,'omitnan') - mean(v2,'omitnan');
                pooled   = [v1(:); v2(:)];
                n1       = numel(v1);
                permDiff = zeros(nBoot,1);
                for bt = 1:nBoot
                    perm         = pooled(randperm(numel(pooled)));
                    permDiff(bt) = mean(perm(1:n1)) - mean(perm(n1+1:end));
                end
                pPerm = mean(abs(permDiff) >= abs(obsDiff));
                pPerm = max(pPerm, 1/nBoot);   % floor at resolution limit

                % PRIMARY TEST: Wilcoxon rank-sum.
                % The permutation test above compares DIFFERENCES OF MEANS,
                % which at n = 11 vs 5 is sensitive to individual large
                % values. Rank-sum depends only on ordering and is the test
                % used elsewhere in the paper (Fig 4H region comparisons, the
                % ACh-vs-control tests above), so using it here keeps the
                % manuscript internally consistent. The permutation result is
                % retained alongside as confirmation.
                if numel(v1) >= 3 && numel(v2) >= 3
                    pRank = ranksum(v1, v2);
                else
                    pRank = NaN;
                end
                pVal = pRank;          % drives the FDR correction and the stars

                idx = numel(pairResults) + 1;
                pairResults(idx).r1      = r1;
                pairResults(idx).r2      = r2;
                pairResults(idx).v1      = v1;
                pairResults(idx).v2      = v2;
                pairResults(idx).obsDiff = obsDiff;
                pairResults(idx).pVal    = pVal;
                pairResults(idx).pPerm   = pPerm;
                pairResults(idx).pRank   = pRank;
            end

            % Step 2: Benjamini-Hochberg FDR across the 3 regional pairs.
            % Previously Holm-Bonferroni. Every other family in the paper is
            % BH-corrected (above-zero tests in this same file, Fig 2, Fig 3,
            % Fig 4, coherence), so a single Holm exception was an
            % inconsistency a reader would have to query. BH is now used
            % throughout; no per-analysis rule needs stating in the Methods.
            nValid   = numel(pairResults);
            pAdjPair = nan(1, nValid);
            isSigAll = false(1, nValid);
            if nValid > 0
                pVals = [pairResults.pVal];
                fin   = isfinite(pVals);
                if any(fin)
                    pf = pVals(fin); nf = numel(pf);
                    [pSorted, sortIdx] = sort(pf(:));
                    adj = pSorted .* nf ./ (1:nf)';
                    adj = flipud(cummin(flipud(adj)));   % enforce monotonicity
                    adj = min(adj, 1);
                    tmp = nan(nf,1); tmp(sortIdx) = adj;
                    pAdjPair(fin) = tmp;
                end
                isSigAll = pAdjPair < alphaNom;
            end

            % Step 3: draw brackets in original display order
            yTop = max(whiskerTops, [], 'omitnan');
            if isnan(yTop) || yTop == 0; yTop = 1; end
            bracketStep = max(0.14 * abs(yTop), eps);

            for pr = 1:nValid
                res   = pairResults(pr);
                r1    = res.r1;  r2 = res.r2;
                pVal  = res.pVal;
                isSig = isSigAll(pr);

                % BH-adjusted p for this pair. Both the raw and the adjusted
                % value go to the CSV, explicitly labelled — storing only the
                % raw p made it easy to transcribe the wrong number into the
                % manuscript.
                pAdjThis = pAdjPair(pr);

                % Asterisk label — graded on the HOLM-ADJUSTED p-value.
                % Previously the gate was Holm-corrected but the number of
                % asterisks was taken from the RAW p, so a pair with raw
                % p = 0.0137 (adjusted 0.0273) was drawn as '**' when the
                % corrected value only supports '*'. Grading on one number
                % while reporting another overstates the result.
                % Asterisks graded on the ADJUSTED p, matching the gate.
                % Grading on the raw p while gating on the corrected one
                % overstates the result (it drew ** where * was warranted).
                if ~isSig
                    sigStr = 'ns';
                elseif pAdjThis < 0.001
                    sigStr = '***';
                elseif pAdjThis < 0.01
                    sigStr = '**';
                else
                    sigStr = '*';
                end

                bCol = [0.15 0.15 0.15];
                if ~isSig; bCol = [0.65 0.65 0.65]; end

                yB = yTop + pr * bracketStep;
                line(ax, [r1 r1 r2 r2], ...
                    [yB-bracketStep*0.25, yB, yB, yB-bracketStep*0.25], ...
                    'Color',bCol,'LineWidth',0.8,'Clipping','off');
                text(ax, (r1+r2)/2, yB + bracketStep*0.05, sigStr, ...
                    'HorizontalAlignment','center', ...
                    'FontSize', FONT_AX + (1*isSig), ...
                    'FontWeight','bold', ...
                    'Color',bCol,'Clipping','off');

                % Accumulate stats row
                v1 = res.v1; v2 = res.v2;
                statRows{end+1} = { ...
                    boutTitleStr{b}, unitStr, ...
                    regionNames{r1}, regionNames{r2}, ...
                    numel(v1), numel(v2), ...
                    mean(v1,'omitnan'), median(v1), std(v1,'omitnan'), ...
                    mean(v2,'omitnan'), median(v2), std(v2,'omitnan'), ...
                    res.obsDiff, pVal, pAdjThis, res.pPerm, isSig, sigStr, alphaNom}; %#ok
            end

            % ----------------------------------------------------------
            % Axes formatting
            % ----------------------------------------------------------
            yline(ax, 0, '-', 'Color',[0.78 0.78 0.78], 'LineWidth',0.7);
            set(ax, ...
                'XTick',1:nR,'XTickLabel',regionNames, ...
                'TickDir','out','FontSize',FONT_AX,'Box','off', ...
                'TickLength',[0.02 0.02]);
            ylabel(ax, sprintf('AUC (%s)', yLabel), 'FontSize',FONT_LB);
            title(ax, sprintf('%s  |  %s', boutTitleStr{b}, aucWinLabels{b}), ...
                'FontSize',FONT_TI,'FontWeight','bold','Color',boutColors{b});
            xlim(ax, [0.4, nR+0.6]);

            add_panel_label(ax, panelLabels{b});
        end

        tightfig_simple(fig);
    end

    % -----------------------------------------------------------------
    function save_stats_table(statRows, fpath)
    % Write one CSV row per pairwise comparison across both bout types
    % Columns: bout_type, signal_unit, region_1, region_2,
    %          n_1, n_2, mean_1, median_1, sd_1,
    %          mean_2, median_2, sd_2,
    %          mean_diff, p_uncorrected, p_bonferroni, significance,
    %          alpha_corrected

        if isempty(statRows)
            warning('No statistics rows to save.');
            return;
        end

        % Test identity is recorded in the table itself so the stats
        % spreadsheet does not depend on remembering which test produced
        % which column.
        hdr = ['cohort,test,test_tail,comparison,' ...
               'multiple_comparison_correction,correction_family,' ...
               'bout_type,signal_unit,region_1,region_2,' ...
               'n_1,n_2,mean_1,median_1,sd_1,' ...
               'mean_2,median_2,sd_2,' ...
               'mean_diff,p_raw_ranksum,p_FDR_adjusted,p_permutation_raw,' ...
               'significant_FDR05,significance,alpha'];

        fid = fopen(fpath, 'w');
        if fid == -1
            error('Cannot open file for writing: %s', fpath);
        end
        fprintf(fid, '%s\n', hdr);

        for k = 1:numel(statRows)
            row = statRows{k};
            fprintf(fid, ['ACh,Wilcoxon rank-sum (unpaired),' ...
                          'two-sided,between regions,Benjamini-Hochberg FDR,' ...
                          'ACh between-region (3 pairs per bout type),' ...
                          '%s,%s,%s,%s,%d,%d,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,' ...
                          '%.6f,%.6f,%.6f,%.6f,%d,%s,%.6f\n'], ...
                row{1},  row{2},  row{3},  row{4}, ...   % labels
                row{5},  row{6},  ...                     % n
                row{7},  row{8},  row{9},  ...            % group 1 stats
                row{10}, row{11}, row{12}, ...            % group 2 stats
                row{13}, ...                              % mean_diff
                row{14}, row{15}, row{16}, ...            % p_raw, p_holm, p_perm
                row{17}, row{18}, row{19});               % sig, significance, alpha
        end
        fclose(fid);
        fprintf('Stats table saved: %s\n', fpath);
    end

    % -----------------------------------------------------------------
    function format_trace_ax(ax, t, ylabelStr, isBottom, fontSz, labelSz)
        yL = ylim(ax);
        if range(yL) < 1e-10; yL = [-1 1]; ylim(ax, yL); end

        % Pre-stimulus shading
        patch(ax, [t(1) 0 0 t(1)], [yL(1) yL(1) yL(2) yL(2)], ...
            [0.92 0.92 0.92], 'FaceAlpha',0.30,'EdgeColor','none');

        % Reference lines
        xline(ax, 0, '--', 'Color',[0.4 0.4 0.4], 'LineWidth',0.8, 'Alpha',0.6);
        yline(ax, 0, '-',  'Color',[0.80 0.80 0.80], 'LineWidth',0.7);

        xlim(ax, [t(1) t(end)]);
        ylabel(ax, ylabelStr, 'FontSize',labelSz);
        set(ax, 'TickDir','out','FontSize',fontSz,'Box','off', ...
            'TickLength',[0.02 0.02]);
        if isBottom
            xlabel(ax,'Time from bout onset (s)','FontSize',labelSz);
        else
            set(ax,'XTickLabel',[]);
        end
    end

    % -----------------------------------------------------------------
    function add_panel_label(ax, lbl)
        pos = get(ax,'Position');
        annotation('textbox', ...
            [pos(1)-0.045, pos(2)+pos(4)-0.01, 0.04, 0.04], ...
            'String', lbl, ...
            'FontSize', 11, 'FontWeight','bold', ...
            'EdgeColor','none','BackgroundColor','none', ...
            'HorizontalAlignment','center');
    end

    % -----------------------------------------------------------------
    function ci = bootstrap_ci(M, nBoot_)
        nS = size(M,1);  nT = size(M,2);
        bm = zeros(nBoot_, nT);
        for b_ = 1:nBoot_
            bm(b_,:) = mean(M(randi(nS,nS,1),:), 1, 'omitnan');
        end
        ci = prctile(bm, [2.5 97.5], 1);
    end

    % -----------------------------------------------------------------
    function plot_with_ci(ax_, t_, mu_, ci_, color_, lw_, alpha_)
        tF = t_(:)';
        fill([tF, fliplr(tF)], [ci_(2,:), fliplr(ci_(1,:))], color_, ...
            'FaceAlpha',alpha_,'EdgeColor','none','Parent',ax_);
        plot(ax_, t_, mu_, '-', 'Color',color_, 'LineWidth',lw_);
    end

    % -----------------------------------------------------------------
    function save_fig(fig, fpath)
        set(fig,'PaperPositionMode','auto');
        print(fig, fpath, '-dpdf', '-painters', '-r0');
        fprintf('Saved: %s\n', fpath);
    end

    % -----------------------------------------------------------------
    function tightfig_simple(fig)
    % Reduce excess whitespace around subplots
        set(fig,'Units','normalized');
        axArr = findall(fig,'type','axes');
        for k = 1:numel(axArr)
            axArr(k).TightInset;   % force update
        end
        % Use built-in if available (R2020a+)
        if ~verLessThan('matlab','9.8')
            set(fig,'Units','inches');
        end
    end

end
