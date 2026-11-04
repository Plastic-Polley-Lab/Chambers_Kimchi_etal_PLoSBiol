function plot_ach_fme_stimulus_vs_spontaneous(ach_stim, all_fiber_names_ach, fme_stim, save_path, modOnly, poolTrials, pupil_stim, risingEdgeOnly, subjectIDs)
% PLOT_ACH_FME_STIMULUS_VS_SPONTANEOUS
%
% Usage:
%   plot_ach_fme_stimulus_vs_spontaneous(ach_stim, all_fiber_names_ach, ...
%       fme_stim, save_path)
%   plot_ach_fme_stimulus_vs_spontaneous(ach_stim, all_fiber_names_ach, ...
%       fme_stim, save_path, modOnly)
%   plot_ach_fme_stimulus_vs_spontaneous(ach_stim, all_fiber_names_ach, ...
%       fme_stim, save_path, modOnly, poolTrials)
%   plot_ach_fme_stimulus_vs_spontaneous(ach_stim, all_fiber_names_ach, ...
%       fme_stim, save_path, modOnly, poolTrials, pupil_stim)
%
%   pupil_stim (optional): trial-structured pupil data. If supplied, the
%   pre-bout arousal check also reports pupil in addition to FME.
%
%   risingEdgeOnly (optional, default false): if true, evoked FME onsets must
%   be RISING EDGES — the sample immediately before onset must be below
%   threshold — matching how spontaneous bout onsets are defined. This removes
%   a definitional asymmetry: spontaneous onsets are rising edges by
%   construction (diff of the above-threshold mask), whereas the default
%   evoked rule takes the first post-stimulus supra-threshold sample, which
%   can land mid-movement if the animal was already moving at stimulus onset.
%   That asymmetry inflates pre-bout FME for evoked bouts. Set true to make
%   the two conditions symmetric.
%
% Inputs:
%   ach_stim            — trial-structured ACh data (stimulus-evoked)
%   all_fiber_names_ach — fiber region names per mouse, cell{m}
%   fme_stim            — trial-structured FME data (stimulus-evoked)
%   save_path           — directory containing *_traces.mat files
%                         (same files used by plot_ach_aligned_to_fme_continuous_bt)
%                         Each file contains a struct 'traces' with fields:
%                           traces.fme          — continuous FME [1 x nSamp]
%                           traces.pupil        — continuous pupil [1 x nSamp] (optional)
%                           traces.achCh1/2     — continuous ACh per channel
%                           traces.chForRegion  — [1 x nR] channel index per region
%                           traces.fsAch        — ACh sampling rate
%   modOnly             — (optional, default false) if true, restrict evoked
%                         trials to moderate intensity only:
%                           AudRipple / Tones : nearest intensity to 50 dB
%                           VisGrating        : nearest contrast to 0.33
%   poolTrials          — (optional, default false) if true, boxplot figure
%                         uses all pooled matched trial pairs with Wilcoxon
%                         signed-rank on pairs; if false (default), uses
%                         per-animal means with Wilcoxon vs zero on differences
%
% For each stimulus type and brain region, compares the ACh response on
% stimulus-evoked trials where FME occurred to spontaneous FME bouts of
% matched magnitude extracted from the continuous silent recordings.
%
% Bout detection from continuous FME matches plot_example_session exactly:
%   1. Median filter the raw FME trace
%   2. Percentile threshold (75th pct x multiplier)
%   3. Merge short sub-threshold gaps (< 0.5 s)
%   4. Extract bout onsets/offsets
%
% Figure 1 — Scatter: FME AUC vs ACh AUC, evoked vs spontaneous
% Figure 2 — Boxplots: ACh AUC evoked vs matched spontaneous, with stats

%% -------------------------------------------------------------------------
%  Parse optional arguments
%% -------------------------------------------------------------------------
if nargin < 5 || isempty(modOnly)
    modOnly = false;
end
if nargin < 6 || isempty(poolTrials)
    poolTrials = false;
end

%% -------------------------------------------------------------------------
%  Moderate-intensity targets (used when modOnly == true)
%  Same nearest-neighbour logic as plot_ach_selectivity_magnitude.m
%% -------------------------------------------------------------------------
targetInt = containers.Map(...
    {'AudRipple', 'Tones', 'VisGrating'}, ...
    {50,          50,      0.33});

%% -------------------------------------------------------------------------
%  Reproducibility: fix random seed for stochastic matching
%% -------------------------------------------------------------------------
rng(42);

%% -------------------------------------------------------------------------
%  Parameters
%% -------------------------------------------------------------------------
nMice       = numel(ach_stim);
nTimepoints = 8139;
t           = linspace(-2, 6, nTimepoints);
baselineIdx = t >= -2 & t < 0;
aucIdx      = t >= 0  & t <= 1;

nNative     = 241;
tNative     = linspace(-2, 6, nNative);
blNative    = tNative >= -2 & tNative < 0;
aucNative   = tNative >= 0  & tNative <= 1;

regions     = {'AC', 'VC', 'PFC'};
nR          = numel(regions);

stimTypes   = {'AudRipple', 'Tones', 'VisGrating'};
intFields   = {'int',       'int',   'contrast'};
stimLabels  = {'Aud. Ripple', 'Tones', 'Vis. Grating'};
nS          = numel(stimTypes);

% FME response threshold: evoked trials must exceed this FME AUC
fmeThreshold = 0.1;

% Matching tolerance: spontaneous FME AUC must be within this fraction
% of the evoked FME AUC to be considered a valid match
matchTol = 0.15;   % 15%

regionColors = [
    0.00, 0.45, 0.70;   % blue   -> AC
    0.90, 0.62, 0.00;   % orange -> VC
    0.00, 0.62, 0.45;   % teal   -> PFC
];
stimColors = [
    0.00, 0.45, 0.70;
    0.90, 0.62, 0.00;
    0.80, 0.47, 0.65;
];
markerSz = 30;

% Subject IDs by trial-structured index, from the extraction table.
% Row m of that table IS entry m of ach_stim/fme_stim/pupil_stim.
if nargin < 9 || isempty(subjectIDs)
    subjectIDs = {'1034','1035','1036','1037','1045','1048','1053','1055', ...
                  '1108','1109','1110','1111','1117','1118','1123','1125', ...
                  '1127','1128','1129'};
end

if nargin < 8 || isempty(risingEdgeOnly)
    risingEdgeOnly = false;
end
if risingEdgeOnly
    fprintf('Evoked onset rule: RISING EDGE (symmetric with spontaneous).\n');
else
    fprintf('Evoked onset rule: first supra-threshold sample (default).\n');
end

% Build figure title suffix based on modOnly flag and onset rule
if modOnly
    titleSuffix = ' — Moderate Intensity Only';
else
    titleSuffix = '';
end
if risingEdgeOnly
    titleSuffix = [titleSuffix ' | rising-edge onsets'];
end

%% -------------------------------------------------------------------------
%  Data collection — pooled across animals
%  For each stimulus type and region, one row per matched trial:
%    fmeAUC_evoked{s}        = [nTrials x 1]
%    achAUC_evoked{s, r}     = [nTrials x 1]   (pooled, for scatter only)
%    fmeAUC_matched{s}       = [nTrials x 1]
%    achAUC_matched{s, r}    = [nTrials x 1]   (pooled, for scatter only)
%  Per-animal means (one row per animal) for boxplots and stats:
%    achAUC_evoked_mn{s, r}  = [nMice x 1]
%    achAUC_matched_mn{s, r} = [nMice x 1]
%  Traces aligned to onset:
%    fmeTrace_evoked{s}      = [nTrials x nNative]
%    fmeTrace_matched{s}     = [nTrials x nNative]
%    achTrace_evoked{s, r}   = [nTrials x nNative]
%    achTrace_matched{s, r}  = [nTrials x nNative]
%    achTrace_diff{s, r}     = [nTrials x nNative]  evoked minus matched spont
%% -------------------------------------------------------------------------
% Diagnostic storage: trial counts and intensity distributions
% diagData{s} has fields: intAll, intAboveThresh, intMatched, fmeAll, fmeAboveThresh, fmeMatched
diagData = cell(nS, 1);
for s = 1:nS
    diagData{s}.intAll          = [];   % intensity of all mod-intensity trials
    diagData{s}.intAboveThresh  = [];   % intensity of trials with FME > threshold
    diagData{s}.intMatched      = [];   % intensity of successfully matched trials
    diagData{s}.fmeAll          = [];   % FME AUC of all mod-intensity trials
    diagData{s}.fmeAboveThresh  = [];   % FME AUC of trials above threshold
    diagData{s}.fmeMatched      = [];   % FME AUC of matched trials
end
achAUC_evoked      = cell(nS, nR);
fmeAUC_matched     = cell(nS, 1);
achAUC_matched     = cell(nS, nR);
achAUC_evoked_mn   = cell(nS, nR);   % per-animal means
achTrace_evoked_mn  = cell(nS, nR);  % [nMice x nNative] per-animal mean traces
% Inputs for the matching bootstrap (R3 comment 2c). One entry per
% mouse x stim x region: the FME AUCs that drive the match, and the ACh AUCs
% that the match selects among. Captured during the main loop so the
% bootstrap can re-run the matching afterwards without touching the main path.
bootIn = {};
achTrace_matched_mn = cell(nS, nR);
if risingEdgeOnly
    onsetTag = '_risingEdge';
else
    onsetTag = '';
end

hasPupil = (nargin >= 7) && ~isempty(pupil_stim);
if hasPupil
    fprintf('Pupil data supplied — pre-bout state check will include pupil.\n');
else
    fprintf('No pupil_stim supplied — pre-bout state check uses FME only.\n');
end

preFme_evoked_mn   = cell(nS, 1);    % per-animal mean pre-bout FME (evoked)
preFme_matched_mn  = cell(nS, 1);    % per-animal mean pre-bout FME (spontaneous)
prePup_evoked_mn   = cell(nS, 1);    % per-animal mean pre-bout pupil (evoked)
prePup_matched_mn  = cell(nS, 1);    % per-animal mean pre-bout pupil (spontaneous)
achAUC_matched_mn  = cell(nS, nR);   % per-animal means
fmeTrace_evoked    = cell(nS, 1);
fmeTrace_matched   = cell(nS, 1);
achTrace_evoked    = cell(nS, nR);
achTrace_matched   = cell(nS, nR);
achTrace_diff      = cell(nS, nR);   % evoked minus matched spontaneous (per trial)
achTrace_diff_mn   = cell(nS, nR);   % [nMice x nNative] per-animal difference

for s = 1:nS
    fmeAUC_evoked{s}    = [];
    fmeAUC_matched{s}   = [];
    fmeTrace_evoked{s}  = [];
    fmeTrace_matched{s} = [];
    preFme_evoked_mn{s}  = [];
    preFme_matched_mn{s} = [];
    prePup_evoked_mn{s}  = [];
    prePup_matched_mn{s} = [];
    for r = 1:nR
        achAUC_evoked{s,r}     = [];
        achAUC_matched{s,r}    = [];
        achAUC_evoked_mn{s,r}  = [];
        achAUC_matched_mn{s,r} = [];
        achTrace_evoked_mn{s,r}  = [];
        achTrace_matched_mn{s,r} = [];
        achTrace_evoked{s,r}   = [];
        achTrace_matched{s,r}  = [];
        achTrace_diff{s,r}     = [];
        achTrace_diff_mn{s,r}  = [];
    end
end

%% -------------------------------------------------------------------------
%  Load continuous traces from save_path
%  Matches plot_ach_aligned_to_fme_continuous_bt loading pattern exactly
%% -------------------------------------------------------------------------
sigFiles = dir(fullfile(save_path, '*_traces.mat'));
isRef    = arrayfun(@(f) ~isempty(strfind(f.name, 'ref')), sigFiles);
sigFiles = sigFiles(~isRef);
if isempty(sigFiles)
    error('No signal trace files found in %s', save_path);
end
nSessions = numel(sigFiles);
fprintf('Found %d signal session files\n', nSessions);

% -------------------------------------------------------------------------
%  Bout detection parameters — identical to plot_example_session
% -------------------------------------------------------------------------
boutThreshPct  = 75;    % percentile threshold
boutThreshMult = 1.0;   % multiplier on threshold
boutMedFilt    = 0.5;   % median filter width in SECONDS (applied at actual fs)
minGapSec      = 0.5;   % merge sub-threshold gaps shorter than this (s)
boutAUCWinSec  = 1.0;   % fixed AUC window length for each bout (s)

for m = 1:nMice
    if isempty(ach_stim{m}) || isempty(fme_stim{m}); continue; end

    % Match this mouse to its continuous session file BY SUBJECT ID.
    % Positional matching (sigFiles(m)) is WRONG for this dataset: sub1055
    % has trial-structured data but no continuous file, so every index from
    % that point on is offset by one. Eight of thirteen mice were previously
    % paired with a different animal's continuous recording.
    if numel(subjectIDs) < m || isempty(subjectIDs{m})
        fprintf('  Mouse %d: no subject ID supplied.\n', m);
        continue;
    end
    thisSubj = subjectIDs{m};
    fileIdx  = find(~cellfun(@isempty, strfind({sigFiles.name}, thisSubj)), 1);
    if isempty(fileIdx)
        fprintf('  Mouse %2d (sub%s): no continuous file — excluded.\n', m, thisSubj);
        continue;
    end

    % --- Load continuous traces (same pattern as plot_ach_aligned_to_fme_continuous_bt) ---
    tmp    = load(fullfile(sigFiles(fileIdx).folder, sigFiles(fileIdx).name), 'traces');
    traces = tmp.traces;

    % ---- LINKING CHECK -------------------------------------------------
    % The trial-structured index m is matched to a continuous session file
    % BY POSITION (sigFiles(m)). That assumption is not guaranteed: for this
    % dataset the ACh cohort mapping turned out to be m -> sorted position m
    % for m <= 6 and m-1 thereafter, because one trial-structured slot has no
    % continuous file. If the same offset applies here, some mice would be
    % matched against ANOTHER animal's continuous recording, and the matching
    % analysis would be pairing evoked and spontaneous bouts from different
    % animals.
    %
    % Verification: traces.chForRegion gives the regions present in the
    % continuous file; all_fiber_names_ach{m} gives them for the
    % trial-structured entry. They must agree.
    linkRegions = {'AC','VC','PFC'};
    contRegs = {};
    if isfield(traces, 'chForRegion') && ~isempty(traces.chForRegion)
        cfr = traces.chForRegion;
        contRegs = linkRegions(cfr(:)' > 0);
    end
    trialRegs = {};
    if numel(all_fiber_names_ach) >= m && ~isempty(all_fiber_names_ach{m})
        trialRegs = all_fiber_names_ach{m}(:)';
    end
    if isempty(contRegs)
        linkFlag = '(no chForRegion field)';
    elseif isempty(trialRegs)
        linkFlag = '(no fibre names)';
    elseif isequal(sort(lower(contRegs)), sort(lower(trialRegs)))
        linkFlag = 'MATCH';
    else
        linkFlag = '*** MISMATCH ***';
    end
    fprintf('  LINK  mouse %2d (sub%s) -> %-26s trial:[%-12s] cont:[%-12s] %s\n', ...
        m, thisSubj, sigFiles(fileIdx).name, strjoin(trialRegs, ','), ...
        strjoin(contRegs, ','), linkFlag);
    % --------------------------------------------------------------------

    fmeCont     = traces.fme(:)';
    if isfield(traces, 'pupil') && ~isempty(traces.pupil)
        pupCont = traces.pupil(:)';
    else
        pupCont = [];
    end
    achCh       = {traces.achCh1, traces.achCh2};
    chForRegion = traces.chForRegion;
    fsAch_cont  = traces.fsAch;
    fsFme_cont  = fsAch_cont;   % FME is stored at ACh rate in traces struct

    nContFme  = numel(fmeCont);
    nContAch  = numel(achCh{1});
    tContFme  = (0:nContFme-1) / fsFme_cont;
    tContAch  = (0:nContAch-1) / fsAch_cont;

    % -----------------------------------------------------------------------
    %  Spontaneous FME bout detection from continuous trace
    %  Mirrors plot_example_session FME processing exactly
    % -----------------------------------------------------------------------

    % Step 1: median filter — width in seconds × actual sample rate
    medFiltSamp = round(boutMedFilt * fsFme_cont);   % boutMedFilt is now in seconds
    fmeFilt = medfilt1(fmeCont, medFiltSamp);
    fmeFilt(isnan(fmeCont)) = NaN;

    % Step 1b: global z-score so AUC values are in comparable units
    % across mice regardless of raw FME pixel scale
    fmeGlobalMu = mean(fmeFilt(~isnan(fmeFilt)), 'omitnan');
    fmeGlobalSd = std( fmeFilt(~isnan(fmeFilt)), 'omitnan');
    if fmeGlobalSd == 0; fmeGlobalSd = 1; end
    fmeFiltZ = (fmeFilt - fmeGlobalMu) ./ fmeGlobalSd;

    % Step 2: percentile threshold on z-scored filtered FME
    validSamps = fmeFiltZ(~isnan(fmeFiltZ));
    if numel(validSamps) < 100
        fprintf('  Mouse %d: insufficient continuous FME data.\n', m);
        continue;
    end
    thresh = prctile(validSamps, boutThreshPct) * boutThreshMult;
    fprintf('    FME z-scored range: %.2f–%.2f  thresh (p%d): %.2f\n', ...
        min(validSamps), max(validSamps), boutThreshPct, thresh);

    % Step 3: above-threshold mask on z-scored trace
    above = fmeFiltZ > thresh;
    above(isnan(fmeFiltZ)) = false;
    above = above(:)';

    % Step 4: merge short sub-threshold gaps
    minGapSamp = round(minGapSec * fsFme_cont);
    d_above    = diff([1, above, 1]);
    falls      = find(d_above == -1);
    rises      = find(d_above ==  1);
    for gi = 1:numel(falls)
        if gi <= numel(rises) && (rises(gi) - falls(gi)) < minGapSamp
            above(falls(gi):rises(gi)-1) = true;
        end
    end

    % Step 5: bout onset/offset
    d    = diff([0, above, 0]);
    bOns = find(d ==  1);
    bOffs= find(d == -1) - 1;
    if isempty(bOns)
        fprintf('  Mouse %d: no spontaneous bouts detected.\n', m);
        continue;
    end
    fprintf('  Mouse %d: %d spontaneous bouts detected.\n', m, numel(bOns));

    % Step 6: for each bout compute FME AUC and concurrent ACh AUC
    % All windows are onset-aligned: t=0 = bout onset.
    % FME/ACh display window: [−preFmeSec, +boutAUCWinSec] around onset.
    % AUC window             : onset to onset + boutAUCWinSec.
    % ACh baseline           : 2 s ending at onset.
    preFmeSec   = 0.5;
    preFmeSamp  = round(preFmeSec  * fsFme_cont);
    aucWinSamp  = round(boutAUCWinSec * fsFme_cont);
    tResamp     = linspace(-preFmeSec, boutAUCWinSec, nNative);  % t=0 = onset

    fmeAUC_sil  = [];
    achAUC_sil  = cell(nR, 1);
    fmeTr_sil   = [];            % [nBouts x nNative]
    pupTr_sil  = [];   % pre-bout pupil traces for spontaneous bouts
    achTr_sil   = cell(nR, 1);  % [nBouts x nNative]
    for r = 1:nR; achAUC_sil{r} = []; achTr_sil{r} = []; end

    for ib = 1:numel(bOns)
        onsFme = bOns(ib);   % onset sample in continuous FME (already threshold-crossing)

        % Skip bouts longer than 1 s
        boutDurSec = (bOffs(ib) - bOns(ib)) / fsFme_cont;
        if boutDurSec > boutAUCWinSec; continue; end

        % FME display window: preFmeSec before onset to boutAUCWinSec after onset
        fOn  = onsFme - preFmeSamp;
        fOff = onsFme + aucWinSamp;
        if fOn < 1 || fOff > nContFme; continue; end

        tFmeWin  = tContFme(fOn:fOff);
        fmeSegZ  = fmeFiltZ(fOn:fOff);

        % AUC from onset to onset + boutAUCWinSec
        boutAUC = trapz(tContFme(onsFme:onsFme+aucWinSamp), ...
                        fmeFiltZ(onsFme:onsFme+aucWinSamp));
        if boutAUC <= fmeThreshold; continue; end

        % Resample FME display trace to nNative, t=0 at onset
        tFmeCent  = tFmeWin - tContFme(onsFme);
        fmeResamp = interp1(tFmeCent, fmeSegZ, tResamp, 'linear', 'extrap');

        % Pupil display trace over the same window (if available)
        pupResamp = nan(1, nNative);
        if ~isempty(pupCont) && numel(pupCont) >= fOff
            pupSeg    = pupCont(fOn:fOff);
            pupResamp = interp1(tFmeCent, pupSeg, tResamp, 'linear', 'extrap');
        end

        % ACh windows anchored to onset time
        tOnset  = tContFme(onsFme);
        tAchOn  = tOnset;
        tAchOff = tOnset + boutAUCWinSec;
        if tAchOn < tContAch(1) || tAchOff > tContAch(end); continue; end

        achIdxOn  = find(tContAch >= tAchOn,  1, 'first');
        achIdxOff = find(tContAch <= tAchOff, 1, 'last');
        if isempty(achIdxOn) || isempty(achIdxOff) || achIdxOff <= achIdxOn; continue; end

        % ACh baseline: 2 s ending at onset
        blAchOn  = find(tContAch >= tOnset - 2, 1, 'first');
        blAchOff = find(tContAch <= tOnset,     1, 'last');

        % ACh display window mirrors FME display window
        achDispOn  = find(tContAch >= tOnset - preFmeSec,     1, 'first');
        achDispOff = find(tContAch <= tOnset + boutAUCWinSec, 1, 'last');

        achBoutAUC = nan(nR, 1);
        achBoutTr  = nan(nR, nNative);
        for r = 1:nR
            ch = chForRegion(r);
            if ch == 0; continue; end
            achRaw = achCh{ch};
            if isempty(blAchOn) || isempty(blAchOff) || blAchOff <= blAchOn; continue; end
            blMu = mean(achRaw(blAchOn:blAchOff), 'omitnan');
            blSd = std( achRaw(blAchOn:blAchOff), 'omitnan');
            if blSd == 0; blSd = 1; end
            % AUC: onset to onset+boutAUCWinSec
            achBoutAUC(r) = trapz(tContAch(achIdxOn:achIdxOff), ...
                                  (achRaw(achIdxOn:achIdxOff) - blMu) / blSd);
            % Display trace: t=0 at onset
            achDispSeg    = (achRaw(achDispOn:achDispOff) - blMu) / blSd;
            tDispSeg      = tContAch(achDispOn:achDispOff) - tOnset;
            achBoutTr(r,:)= interp1(tDispSeg, achDispSeg, tResamp, 'linear', 'extrap');
        end
        recordedR = chForRegion > 0;
        if ~any(isfinite(achBoutAUC(recordedR))); continue; end

        fmeAUC_sil = [fmeAUC_sil; boutAUC]; %#ok
        fmeTr_sil  = [fmeTr_sil;  fmeResamp]; %#ok
        pupTr_sil  = [pupTr_sil;  pupResamp]; %#ok
        for r = 1:nR
            achAUC_sil{r} = [achAUC_sil{r}; achBoutAUC(r)];
            achTr_sil{r}  = [achTr_sil{r};  achBoutTr(r,:)];
        end
    end

    if isempty(fmeAUC_sil)
        fprintf('  Mouse %d: no bouts passed FME AUC threshold (%.2f).\n', m, fmeThreshold);
        continue;
    end
    fprintf('  Mouse %d: %d bouts used (AUC range %.2f-%.2f, threshold %.2f).\n', ...
        m, numel(fmeAUC_sil), min(fmeAUC_sil), max(fmeAUC_sil), fmeThreshold);

    % -----------------------------------------------------------------------
    %  Evoked condition: loop over stimulus types
    % -----------------------------------------------------------------------
    for s = 1:nS
        if ~isfield(fme_stim{m}, stimTypes{s}); continue; end

        fmeEv = fme_stim{m}.(stimTypes{s}).traces;
        if iscell(fmeEv); fmeEv = fmeEv{1}; end
        if isempty(fmeEv); continue; end

        % -----------------------------------------------------------------
        %  modOnly: select trials at moderate intensity only
        %  Uses nearest-neighbour match to target intensity (same logic as
        %  plot_ach_selectivity_magnitude.m)
        % -----------------------------------------------------------------
        if modOnly
            intVec = fme_stim{m}.(stimTypes{s}).(intFields{s});
            intVec = intVec(:);
            uniqueInts  = unique(intVec);
            tgt         = targetInt(stimTypes{s});
            [~, nearest] = min(abs(uniqueInts - tgt));
            modIntVal   = uniqueInts(nearest);
            modTrialIdx = intVec == modIntVal;   % logical [nTrials x 1]
            fprintf('    %s modOnly: target=%.2f, using %.2f (%d/%d trials)\n', ...
                stimTypes{s}, tgt, modIntVal, sum(modTrialIdx), numel(modTrialIdx));

            fmeEv = fmeEv(modTrialIdx, :);
        else
            modTrialIdx = true(size(fmeEv, 1), 1);   % all trials
        end

        if isempty(fmeEv); continue; end

        % Session-level z-score of evoked FME trials — same normalisation as
        % the continuous trace (fmeGlobalMu / fmeGlobalSd computed above).
        % Each trial is median-filtered row-wise before z-scoring.
        % IMPORTANT: evoked FME trials are at 30 Hz (nNative=241 over 8 s),
        % so the filter kernel must be computed at 30 Hz, not at fsAch_cont.
        fsFmeEv         = nNative / 8;                          % ~30.125 Hz
        medFiltSamp_ev  = max(3, round(boutMedFilt * fsFmeEv)); % e.g. ~15 samples
        fmeEvFilt = fmeEv;
        for ii = 1:size(fmeEv, 1)
            row = medfilt1(fmeEv(ii,:), medFiltSamp_ev);
            row(isnan(fmeEv(ii,:))) = NaN;
            fmeEvFilt(ii,:) = (row - fmeGlobalMu) / fmeGlobalSd;
        end

        % AUC and traces are onset-aligned: t=0 = first post-stimulus sample
        % where session-z-scored FME exceeds thresh (same threshold used for
        % spontaneous bout detection).  tNative covers -2 to 6 s re stimulus.
        postStimMask = tNative >= 0;   % only look for onset after stimulus
        fsFmeEv_eff  = nNative / 8;   % effective sample rate of trial traces (~30 Hz)
        aucWinSampEv = round(boutAUCWinSec * fsFmeEv_eff);

        % Per-trial onset index (into tNative), AUC, and display trace
        fmeAUC_ev    = nan(size(fmeEvFilt, 1), 1);
        fmeOnsetIdx  = nan(size(fmeEvFilt, 1), 1);   % column index into tNative
        for ii = 1:size(fmeEvFilt, 1)
            row = fmeEvFilt(ii,:);

            if risingEdgeOnly
                % Symmetric with spontaneous detection: onset must be a
                % RISING EDGE (previous sample below threshold), so the
                % pre-bout window is guaranteed sub-threshold in both
                % conditions. Removes the definitional asymmetry that
                % otherwise inflates pre-bout FME for evoked bouts.
                aboveRow = row > thresh;
                aboveRow(isnan(row)) = false;
                risingRow = [false, aboveRow(2:end) & ~aboveRow(1:end-1)];
                postIdx   = find(postStimMask & risingRow, 1, 'first');
            else
                % Default: first supra-threshold sample after stimulus onset
                postIdx = find(postStimMask & row > thresh, 1, 'first');
            end

            if isempty(postIdx); continue; end
            aucEnd = min(postIdx + aucWinSampEv, numel(tNative));
            fmeAUC_ev(ii)   = trapz(tNative(postIdx:aucEnd), row(postIdx:aucEnd));
            fmeOnsetIdx(ii) = postIdx;
        end

        % Only keep evoked trials where FME onset was detected and AUC > threshold
        evOk = ~isnan(fmeAUC_ev) & fmeAUC_ev > fmeThreshold;
        if sum(evOk) == 0; continue; end

        fmeAUC_ev_ok  = fmeAUC_ev(evOk);
        fmeOnsetOk    = fmeOnsetIdx(evOk);   % onset sample indices for kept trials
        nEv = numel(fmeAUC_ev_ok);

        % Diagnostic: collect intensity and FME AUC for all, above-threshold,
        % and matched trials (intensities stored after matching below)
        intVecFull = fme_stim{m}.(stimTypes{s}).(intFields{s})(:);
        if modOnly
            intVecMod = intVecFull(modTrialIdx);
        else
            intVecMod = intVecFull;
        end
        % All mod-intensity trials
        fmeAUC_allTrials = fmeAUC_ev;   % includes NaN (no onset) and below-thresh
        diagData{s}.intAll  = [diagData{s}.intAll;  intVecMod(isfinite(fmeAUC_ev))];
        diagData{s}.fmeAll  = [diagData{s}.fmeAll;  fmeAUC_ev(isfinite(fmeAUC_ev))];
        % Above-threshold trials
        diagData{s}.intAboveThresh = [diagData{s}.intAboveThresh; intVecMod(evOk)];
        diagData{s}.fmeAboveThresh = [diagData{s}.fmeAboveThresh; fmeAUC_ev_ok];

        if risingEdgeOnly
            fprintf('    %s: rising-edge rule — %d/%d trials had a valid onset\n', ...
                stimTypes{s}, sum(~isnan(fmeOnsetIdx)), numel(fmeOnsetIdx));
        end
        fprintf('    %s: %d evoked trials above threshold  evFME AUC range %.2f–%.2f  silFME AUC range %.2f–%.2f\n', ...
            stimTypes{s}, nEv, min(fmeAUC_ev_ok), max(fmeAUC_ev_ok), ...
            min(fmeAUC_sil), max(fmeAUC_sil));

        % ACh for evoked trials — skip missing regions, don't abort
        % Also collect z-scored traces resampled to nNative points
        achAUC_ev_ok  = cell(nR, 1);
        achTr_ev_ok   = cell(nR, 1);
        if ~isfield(ach_stim{m}, stimTypes{s}); continue; end

        % Full ACh trial matrix for this stimType before subsetting
        for r = 1:nR
            fiberIdx = find(strcmp(all_fiber_names_ach{m}, regions{r}));
            if isempty(fiberIdx)
                achAUC_ev_ok{r} = nan(nEv, 1);
                achTr_ev_ok{r}  = nan(nEv, nNative);
                continue;
            end
            evTracesFull = ach_stim{m}.(stimTypes{s}).traces{fiberIdx};

            % Apply same modOnly subsetting to ACh traces
            if modOnly
                % evTracesFull rows must align with original fme_stim trial count
                evTracesAll = evTracesFull;
                if size(evTracesAll, 1) ~= numel(modTrialIdx)
                    achAUC_ev_ok{r} = nan(nEv, 1);
                    achTr_ev_ok{r}  = nan(nEv, nNative);
                    continue;
                end
                evTraces = evTracesAll(modTrialIdx, :);
            else
                evTraces = evTracesFull;
            end

            if size(evTraces,1) ~= size(fmeEv,1)
                achAUC_ev_ok{r} = nan(nEv, 1);
                achTr_ev_ok{r}  = nan(nEv, nNative);
                continue;
            end
            achAll = compute_ach_trial_auc_onset(evTraces(evOk,:), fmeOnsetOk, tNative, boutAUCWinSec);
            achAUC_ev_ok{r} = achAll;
            % Resample each ACh trial trace onset-aligned to tResamp
            % tResamp = linspace(-preFmeSec, boutAUCWinSec, nNative), t=0=onset
            tAchTr    = linspace(-2, 6, size(evTraces,2));
            blMaskTr  = tAchTr >= -2 & tAchTr < 0;
            trMat     = evTraces(evOk, :);
            blMu      = mean(trMat(:, blMaskTr), 2, 'omitnan');
            trMat     = trMat - blMu;
            blSd      = std(trMat(:, blMaskTr), 0, 2, 'omitnan');
            blSd(blSd == 0) = 1;
            trMat     = trMat ./ blSd;
            achTr_ev_ok{r} = nan(nEv, nNative);
            for ii = 1:nEv
                tOnsetSec = tNative(fmeOnsetOk(ii));   % onset time re stimulus (s)
                tShifted  = tAchTr - tOnsetSec;        % t=0 at onset
                achTr_ev_ok{r}(ii,:) = interp1(tShifted, trMat(ii,:), tResamp, ...
                    'linear', 'extrap');
            end
        end

        % Evoked FME display traces — onset-aligned, same tResamp as spontaneous
        fmeTr_ev_ok = nan(nEv, nNative);
        fmeEvOkMat  = fmeEvFilt(evOk, :);
        for ii = 1:nEv
            tOnsetSec = tNative(fmeOnsetOk(ii));
            tShifted  = tNative - tOnsetSec;   % t=0 at onset
            fmeTr_ev_ok(ii,:) = interp1(tShifted, fmeEvOkMat(ii,:), tResamp, ...
                'linear', 'extrap');
        end

        % Evoked pupil display traces — onset-aligned, same tResamp
        pupTr_ev_ok = nan(nEv, nNative);
        if hasPupil && numel(pupil_stim) >= m && ~isempty(pupil_stim{m}) ...
                && isfield(pupil_stim{m}, stimTypes{s})
            pupEv = pupil_stim{m}.(stimTypes{s}).traces;
            if iscell(pupEv); pupEv = pupEv{1}; end
            if ~isempty(pupEv)
                if modOnly && exist('modTrialIdx','var') && ...
                        size(pupEv,1) == numel(modTrialIdx)
                    pupEv = pupEv(modTrialIdx, :);
                end
                if size(pupEv,1) == size(fmeEvFilt,1)
                    pupEvOkMat = pupEv(evOk, :);
                    for ii = 1:nEv
                        tOnsetSec = tNative(fmeOnsetOk(ii));
                        tShifted  = tNative - tOnsetSec;
                        pupTr_ev_ok(ii,:) = interp1(tShifted, pupEvOkMat(ii,:), ...
                            tResamp, 'linear', 'extrap');
                    end
                end
            end
        end

        % Magnitude-balanced matching: for each evoked trial, find candidates
        % above and below within ±matchTol. Pair each above candidate with
        % the below candidate of closest magnitude (|spont - evoked|). From
        % the resulting magnitude-matched pairs, randomly select one pair
        % then randomly choose the above or below member with 50/50 probability.
        % This guarantees both sign balance AND magnitude balance, so the
        % distribution of (spont - evoked) differences is centred at zero
        % with symmetric spread.
        matchIdx   = zeros(nEv, 1);
        matchValid = false(nEv, 1);
        for i = 1:nEv
            tolWindow = matchTol * abs(fmeAUC_ev_ok(i));
            allCands  = find(abs(fmeAUC_sil - fmeAUC_ev_ok(i)) <= tolWindow);
            if isempty(allCands); continue; end

            aboveIdx = allCands(fmeAUC_sil(allCands) > fmeAUC_ev_ok(i));
            belowIdx = allCands(fmeAUC_sil(allCands) < fmeAUC_ev_ok(i));

            % Require candidates on both sides
            if isempty(aboveIdx) || isempty(belowIdx); continue; end

            aboveDist = fmeAUC_sil(aboveIdx) - fmeAUC_ev_ok(i);  % positive
            belowDist = fmeAUC_ev_ok(i) - fmeAUC_sil(belowIdx);  % positive magnitude

            % For each above candidate, find the below candidate with
            % closest magnitude — build nAbove matched pairs
            nA = numel(aboveIdx);
            pairedAbove = zeros(nA,1);
            pairedBelow = zeros(nA,1);
            for ai = 1:nA
                [~, bi] = min(abs(belowDist - aboveDist(ai)));
                pairedAbove(ai) = aboveIdx(ai);
                pairedBelow(ai) = belowIdx(bi);
            end

            % Randomly select one magnitude-matched pair
            pairSel = randi(nA);

            % Randomly choose above or below member with 50/50 probability
            if rand() < 0.5
                matchIdx(i) = pairedAbove(pairSel);
            else
                matchIdx(i) = pairedBelow(pairSel);
            end
            matchValid(i) = true;
        end
        nMatched   = sum(matchValid);
        nUnmatched = nEv - nMatched;
        fprintf('    %s: %d / %d evoked trials matched within %.0f%% tolerance, balanced sub-pool (%d discarded)\n', ...
            stimTypes{s}, nMatched, nEv, matchTol*100, nUnmatched);
        if nMatched == 0; continue; end

        % Valid match: within tolerance AND at least one region has finite ACh
        validMatch = false(nEv, 1);
        for i = 1:nEv
            if ~matchValid(i); continue; end
            for r = 1:nR
                if isfinite(achAUC_ev_ok{r}(i)) && isfinite(achAUC_sil{r}(matchIdx(i)))
                    validMatch(i) = true; break;
                end
            end
        end
        fprintf('    %s: %d / %d matched trials with finite ACh\n', ...
            stimTypes{s}, sum(validMatch), nMatched);
        if sum(validMatch) == 0; continue; end

        % Diagnostic: store intensities and FME AUC of matched trials
        % validMatch indexes into the evOk subset, so expand back to full
        % mod-intensity trial space for intensity lookup
        evOkIdx = find(evOk);   % indices into mod-intensity trials that passed threshold
        matchedFullIdx = evOkIdx(validMatch);   % indices of matched trials in full mod space
        diagData{s}.intMatched  = [diagData{s}.intMatched;  intVecMod(matchedFullIdx)];
        diagData{s}.fmeMatched  = [diagData{s}.fmeMatched;  fmeAUC_ev_ok(validMatch)];

        % Store pooled across animals — one row per matched trial (for scatter)
        fmeAUC_evoked{s}    = [fmeAUC_evoked{s};    fmeAUC_ev_ok(validMatch)];
        fmeAUC_matched{s}   = [fmeAUC_matched{s};   fmeAUC_sil(matchIdx(validMatch))];
        fmeTrace_evoked{s}  = [fmeTrace_evoked{s};  fmeTr_ev_ok(validMatch,:)];
        fmeTrace_matched{s} = [fmeTrace_matched{s}; fmeTr_sil(matchIdx(validMatch),:)];

        % --- Pre-bout FME state (t = -preFmeSec to 0), per-animal mean ---
        % Tests whether evoked and matched spontaneous bouts arose from
        % comparable pre-bout arousal contexts (addresses R3 Major 2c).
        preMaskFme = tResamp >= -preFmeSec & tResamp < 0;
        evPre = mean(fmeTr_ev_ok(validMatch, preMaskFme), 2, 'omitnan');
        spPre = mean(fmeTr_sil(matchIdx(validMatch), preMaskFme), 2, 'omitnan');
        finPre = isfinite(evPre) & isfinite(spPre);
        if sum(finPre) > 0
            preFme_evoked_mn{s}  = [preFme_evoked_mn{s};  mean(evPre(finPre), 'omitnan')];
            preFme_matched_mn{s} = [preFme_matched_mn{s}; mean(spPre(finPre), 'omitnan')];
        end

        % Same for pupil (if available)
        if hasPupil && ~isempty(pupTr_sil)
            evPreP = mean(pupTr_ev_ok(validMatch, preMaskFme), 2, 'omitnan');
            spPreP = mean(pupTr_sil(matchIdx(validMatch), preMaskFme), 2, 'omitnan');
            finPreP = isfinite(evPreP) & isfinite(spPreP);
            if sum(finPreP) > 0
                prePup_evoked_mn{s}  = [prePup_evoked_mn{s};  mean(evPreP(finPreP), 'omitnan')];
                prePup_matched_mn{s} = [prePup_matched_mn{s}; mean(spPreP(finPreP), 'omitnan')];
            end
        end
        for r = 1:nR
            achAUC_evoked{s,r}    = [achAUC_evoked{s,r};    achAUC_ev_ok{r}(validMatch)];
            achAUC_matched{s,r}   = [achAUC_matched{s,r};   achAUC_sil{r}(matchIdx(validMatch))];
            achTrace_evoked{s,r}  = [achTrace_evoked{s,r};  achTr_ev_ok{r}(validMatch,:)];
            achTrace_matched{s,r} = [achTrace_matched{s,r}; achTr_sil{r}(matchIdx(validMatch),:)];

            % Per-trial difference traces: evoked minus matched spontaneous
            evTr = achTr_ev_ok{r}(validMatch,:);
            spTr = achTr_sil{r}(matchIdx(validMatch),:);
            finPairsTr = all(isfinite(evTr),2) & all(isfinite(spTr),2);
            achTrace_diff{s,r} = [achTrace_diff{s,r}; evTr(finPairsTr,:) - spTr(finPairsTr,:)];

            % Per-animal mean (one value per mouse) — used for boxplots/stats
            evTrials = achAUC_ev_ok{r}(validMatch);
            spTrials = achAUC_sil{r}(matchIdx(validMatch));
            finPairs = isfinite(evTrials) & isfinite(spTrials);
            if sum(finPairs) > 0
                achAUC_evoked_mn{s,r}  = [achAUC_evoked_mn{s,r};  mean(evTrials(finPairs), 'omitnan')];
                achAUC_matched_mn{s,r} = [achAUC_matched_mn{s,r}; mean(spTrials(finPairs), 'omitnan')];

                % Capture the raw vectors so the match can be re-drawn later
                bootIn{end+1} = struct( ...
                    'stim',      s, ...
                    'region',    r, ...
                    'mouse',     m, ...
                    'fmeEv',     fmeAUC_ev_ok(:), ...
                    'fmeSil',    fmeAUC_sil(:), ...
                    'achEv',     achAUC_ev_ok{r}(:), ...
                    'achSil',    achAUC_sil{r}(:)); %#ok

                % Per-animal mean traces (for peak comparison)
                evTrM = achTr_ev_ok{r}(validMatch,:);
                spTrM = achTr_sil{r}(matchIdx(validMatch),:);
                finTrM = all(isfinite(evTrM),2) & all(isfinite(spTrM),2);
                if sum(finTrM) > 0
                    achTrace_evoked_mn{s,r}  = [achTrace_evoked_mn{s,r}; ...
                        mean(evTrM(finTrM,:), 1, 'omitnan')];
                    achTrace_matched_mn{s,r} = [achTrace_matched_mn{s,r}; ...
                        mean(spTrM(finTrM,:), 1, 'omitnan')];
                    % Per-animal difference: mean of the paired per-trial
                    % differences, so pairing is preserved within animal
                    achTrace_diff_mn{s,r}    = [achTrace_diff_mn{s,r}; ...
                        mean(evTrM(finTrM,:) - spTrM(finTrM,:), 1, 'omitnan')];
                end
            end
        end
    end
end

%% =========================================================================
%  Pre-bout arousal state check (addresses R3 Major Comment 2c)
%  Compares mean FME and pupil in the pre-onset window (-preFmeSec to 0 s)
%  between stimulus-evoked bouts and their matched spontaneous counterparts.
%  If comparable, the two groups arose from equivalent arousal contexts and
%  explicit pre-stimulus state matching is unnecessary.
%  Unit of analysis: per-animal means. Paired Wilcoxon signed-rank.
%% =========================================================================
preSigLabels = {'Pre-bout FME (z)', 'Pre-bout pupil'};
preEvData    = {preFme_evoked_mn,  prePup_evoked_mn};
preSpData    = {preFme_matched_mn, prePup_matched_mn};
nPreSig      = 1 + double(hasPupil);

figure('Units','normalized','Position',[0.10 0.15 0.65 0.70]);
sgtitle(['Pre-bout Arousal State: Evoked vs Matched Spontaneous' titleSuffix], ...
    'FontSize', 12, 'FontWeight', 'bold');

preStatsRows = {};

for sg = 1:nPreSig
    for s = 1:nS
        ax = subplot(nPreSig, nS, (sg-1)*nS + s);
        hold(ax, 'on');

        evP = preEvData{sg}{s};
        spP = preSpData{sg}{s};
        nP  = min(numel(evP), numel(spP));
        if nP < 3
            title(ax, stimLabels{s}, 'FontSize', 10, 'FontWeight','bold');
            continue;
        end
        evP = evP(1:nP); spP = spP(1:nP);
        fin = isfinite(evP) & isfinite(spP);
        evP = evP(fin); spP = spP(fin);
        nP  = numel(evP);
        if nP < 3
            title(ax, stimLabels{s}, 'FontSize', 10, 'FontWeight','bold');
            continue;
        end

        cEv = stimColors(s,:);
        cSp = [0.55 0.55 0.55];

        for k = 1:nP
            plot(ax, [1 2], [evP(k) spP(k)], '-', ...
                'Color', [0.7 0.7 0.7], 'LineWidth', 0.9);
        end
        scatter(ax, ones(nP,1),   evP, 34, cEv, 'filled', ...
            'MarkerFaceAlpha', 0.8, 'MarkerEdgeColor','none');
        scatter(ax, 2*ones(nP,1), spP, 34, cSp, 'filled', ...
            'MarkerFaceAlpha', 0.8, 'MarkerEdgeColor','none');

        errorbar(ax, 0.78, mean(evP,'omitnan'), std(evP,0,'omitnan')/sqrt(nP), ...
            'o','Color',cEv,'MarkerFaceColor',cEv,'MarkerSize',7, ...
            'LineWidth',2.0,'CapSize',5,'LineStyle','none');
        errorbar(ax, 2.22, mean(spP,'omitnan'), std(spP,0,'omitnan')/sqrt(nP), ...
            'o','Color',cSp,'MarkerFaceColor',cSp,'MarkerSize',7, ...
            'LineWidth',2.0,'CapSize',5,'LineStyle','none');

        pPre = signrank(evP, spP);
        preStatsRows{end+1} = {preSigLabels{sg}, stimLabels{s}, nP, ...
            mean(evP,'omitnan'), std(evP,0,'omitnan')/sqrt(nP), ...
            mean(spP,'omitnan'), std(spP,0,'omitnan')/sqrt(nP), ...
            mean(evP-spP,'omitnan'), pPre}; %#ok

        set(ax,'XTick',[1 2],'XTickLabel',{'Evoked','Spont.'}, ...
            'TickDir','out','FontSize',9);
        xlim(ax,[0.6 2.4]);
        if s == 1; ylabel(ax, preSigLabels{sg}, 'FontSize', 10); end
        if sg == 1; title(ax, stimLabels{s}, 'FontSize', 10, 'FontWeight','bold'); end
        box off;
    end
end

annotation('textbox',[0.01 0 0.99 0.05], ...
    'String',['Mean FME and pupil in the pre-onset window (-0.5 to 0 s relative to bout onset).  ' ...
              'One point per mouse; grey lines connect paired observations.  ' ...
              'Large markers = mean \pm SEM.  Paired Wilcoxon signed-rank; ' ...
              'p-values printed to command window and saved to CSV.'], ...
    'EdgeColor','none','FontSize',8,'Color',[0.4 0.4 0.4], ...
    'HorizontalAlignment','center');

fprintf('\n============================================================\n');
fprintf('Pre-bout arousal state: evoked vs matched spontaneous\n');
fprintf('Unit: per-animal means.  Paired Wilcoxon signed-rank.\n');
fprintf('============================================================\n');
for k = 1:numel(preStatsRows)
    rw = preStatsRows{k};
    fprintf('\n%s  |  %s  (n = %d mice)\n', rw{1}, rw{2}, rw{3});
    fprintf('  Evoked:      %+.4f +/- %.4f (mean +/- SEM)\n', rw{4}, rw{5});
    fprintf('  Spontaneous: %+.4f +/- %.4f\n', rw{6}, rw{7});
    fprintf('  Difference:  %+.4f\n', rw{8});
    if rw{9} < 0.05
        fprintf('  Wilcoxon p = %.4f  <-- DIFFERENT: consider matching pre-bout state\n', rw{9});
    else
        fprintf('  Wilcoxon p = %.4f  (n.s. - pre-bout state comparable)\n', rw{9});
    end
end
fprintf('\n============================================================\n');

if ~isempty(preStatsRows)
    csvPre = fullfile(pwd, sprintf('prebout_state_comparison%s.csv', onsetTag));
    fidP = fopen(csvPre, 'w');
    if fidP == -1
        [pDir, pName, pExt] = fileparts(csvPre);
        csvPre = fullfile(pDir, sprintf('%s_%s%s', pName, ...
            datestr(now,'yyyymmdd_HHMMSS'), pExt));
        fidP = fopen(csvPre, 'w');
    end
    if fidP == -1
        warning(['Could not write pre-bout state CSV. ' ...
                 'Check that the file is not open in another program.']);
    else
    fprintf(fidP, 'signal,stimulus,n_mice,mean_evoked,sem_evoked,mean_spont,sem_spont,mean_diff,p_Wilcoxon_paired\n');
    for k = 1:numel(preStatsRows)
        rw = preStatsRows{k};
        fprintf(fidP, '%s,%s,%d,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\n', ...
            rw{1}, rw{2}, rw{3}, rw{4}, rw{5}, rw{6}, rw{7}, rw{8}, rw{9});
    end
    fclose(fidP);
    fprintf('Pre-bout state comparison saved to: %s\n', csvPre);
    end
end

%% =========================================================================
%  Diagnostic figure: FME AUC match quality
%  One panel per stimulus type — scatter of evoked vs matched spontaneous
%  FME AUC, with identity line. Tight pairing = points along diagonal.
%% =========================================================================
figure('Units','normalized','Position',[0.05 0.55 0.60 0.35]);
sgtitle(sprintf('FME AUC Match Quality (balanced \\pm%.0f%% tolerance)%s', matchTol*100, titleSuffix), ...
    'FontSize', 12, 'FontWeight', 'bold');
for s = 1:nS
    ax = subplot(1, nS, s);
    hold(ax, 'on');
    fEv = fmeAUC_evoked{s};
    fSp = fmeAUC_matched{s};
    n   = min(numel(fEv), numel(fSp));
    if n < 1
        title(ax, stimLabels{s}); continue;
    end
    scatter(ax, fEv(1:n), fSp(1:n), 20, stimColors(s,:), 'filled', ...
        'MarkerFaceAlpha', 0.5);
    lims = [min([fEv(1:n); fSp(1:n)]), max([fEv(1:n); fSp(1:n)])];
    plot(ax, lims, lims, '--k', 'LineWidth', 1.0);
    % ±10% tolerance lines
    plot(ax, lims, lims*(1+matchTol), ':', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.8);
    plot(ax, lims, lims*(1-matchTol), ':', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.8);
    xlabel(ax, 'Evoked FME AUC (z·s)', 'FontSize', 9);
    ylabel(ax, 'Matched spont. FME AUC (z·s)', 'FontSize', 9);
    title(ax, sprintf('%s  (n=%d)', stimLabels{s}, n), 'FontSize', 10, 'FontWeight', 'bold');
    axis(ax, 'square'); box off;
    set(ax, 'TickDir', 'out', 'FontSize', 9);
end

%% =========================================================================
%  FME AUC bias figure: distribution of (evoked − matched spontaneous)
%  per stimulus type. Half-violin + mean ± SEM + Wilcoxon vs zero.
%  Tests whether evoked FME is systematically larger or smaller than its
%  matched spontaneous counterpart within the tolerance window.
%% =========================================================================
figure('Units','normalized','Position',[0.05 0.55 0.65 0.38]);
sgtitle(sprintf('FME AUC Matching Bias (evoked \\minus spont)%s', titleSuffix), ...
    'FontSize', 12, 'FontWeight', 'bold');

nKde2      = 200;
cloudWidth2= 0.30;
dotSpread2 = 0.10;
cSpont2    = [0.6 0.6 0.6];

for s = 1:nS
    ax = subplot(1, nS, s);
    hold(ax, 'on');

    fEv = fmeAUC_evoked{s};
    fSp = fmeAUC_matched{s};
    n   = min(numel(fEv), numel(fSp));
    if n < 5
        title(ax, stimLabels{s}, 'FontSize', 10, 'FontWeight', 'bold');
        continue;
    end

    diffs = fEv(1:n) - fSp(1:n);
    diffs = diffs(isfinite(diffs));
    nd    = numel(diffs);
    if nd < 5; continue; end

    % Half-violin (kernel density) to the right
    [f, yi] = ksdensity(diffs, 'NumPoints', nKde2);
    f = f / max(f) * cloudWidth2;
    patch(ax, [1 + f, ones(1,nKde2)], [yi, fliplr(yi)], ...
        stimColors(s,:), 'FaceAlpha', 0.30, 'EdgeColor', 'none');
    plot(ax, 1 + f, yi, '-', 'Color', stimColors(s,:), 'LineWidth', 1.2);

    % Jittered dots to the left
    jit = dotSpread2 * (rand(nd,1) - 0.5);
    scatter(ax, 0.92 + jit, diffs, 10, stimColors(s,:), 'filled', ...
        'MarkerFaceAlpha', 0.25, 'MarkerEdgeColor', 'none');

    % Mean ± SEM
    muD  = mean(diffs, 'omitnan');
    semD = std(diffs, 0, 'omitnan') / sqrt(nd);
    errorbar(ax, 1, muD, semD, 'o', ...
        'Color', stimColors(s,:), 'MarkerFaceColor', stimColors(s,:), ...
        'MarkerSize', 7, 'LineWidth', 2.0, 'CapSize', 6, 'LineStyle', 'none');

    % Wilcoxon signed-rank vs zero
    pW = signrank(diffs);
    if pW < 0.001;    pStr = 'p<0.001';
    elseif pW < 0.01; pStr = sprintf('p=%.3f', pW);
    elseif pW < 0.05; pStr = sprintf('p=%.3f', pW);
    else;              pStr = sprintf('p=%.2f', pW);
    end

    yline(ax, 0, '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 1.0);
    text(ax, 0.60, 0.94, pStr, 'Units', 'normalized', ...
        'FontSize', 8.5, 'Color', 'k', 'HorizontalAlignment', 'center');
    text(ax, 0.60, 0.85, sprintf('n=%d pairs', nd), 'Units', 'normalized', ...
        'FontSize', 8, 'Color', [0.5 0.5 0.5], 'HorizontalAlignment', 'center');
    text(ax, 0.60, 0.76, sprintf('mean=%.3f', muD), 'Units', 'normalized', ...
        'FontSize', 8, 'Color', stimColors(s,:), 'HorizontalAlignment', 'center');

    set(ax, 'XTick', [], 'TickDir', 'out', 'FontSize', 9);
    xlim(ax, [0.55 1.45]);
    ylabel(ax, '\DeltaFME AUC (evoked \minus spont, z\cdots)', 'FontSize', 9);
    title(ax, stimLabels{s}, 'FontSize', 10, 'FontWeight', 'bold');
    box off;
end

annotation('textbox', [0.01 0.00 0.99 0.04], ...
    'String', ['\DeltaFME AUC = evoked \minus matched spontaneous, per pair.  ' ...
               'Half-violin = kernel density; dots = individual pairs (jittered).  ' ...
               'Circle + bar = mean \pm SEM.  Wilcoxon signed-rank vs zero.  ' ...
               sprintf('Matching tolerance: \\pm%.0f%%.', matchTol*100)], ...
    'EdgeColor', 'none', 'FontSize', 8, 'Color', [0.4 0.4 0.4], ...
    'HorizontalAlignment', 'center');

%% =========================================================================
%  Example FME trace pairs for AudRipple
%  Shows the match quality scatter with two highlighted example pairs:
%  one from the 60th–75th percentile of matched FME AUC (higher end of mean)
%  one from the 25th–40th percentile (lower end of mean).
%  Separate panels show onset-aligned FME traces for each example pair.
%% =========================================================================
sRipple = find(strcmp(stimTypes, 'AudRipple'));
if ~isempty(sRipple)
    fEvR = fmeAUC_evoked{sRipple};
    fSpR = fmeAUC_matched{sRipple};
    TrEv = fmeTrace_evoked{sRipple};    % [nPairs x nNative]
    TrSp = fmeTrace_matched{sRipple};   % [nPairs x nNative]
    n    = min([numel(fEvR), numel(fSpR), size(TrEv,1), size(TrSp,1)]);

    if n >= 10
        fEvR = fEvR(1:n); fSpR = fSpR(1:n);
        TrEv = TrEv(1:n,:); TrSp = TrSp(1:n,:);

        % Manually specified example pair indices
        hiPair = 101;
        loPair = 248;

        % Validate indices are within range and have finite traces
        validPairs = find(isfinite(fEvR) & isfinite(fSpR) & ...
                          all(isfinite(TrEv),2) & all(isfinite(TrSp),2));
        if ~ismember(hiPair, validPairs) || ~ismember(loPair, validPairs)
            warning('Specified pair indices not valid — falling back to percentile selection.');
            pct65 = prctile(fEvR(validPairs), 65);
            pct35 = prctile(fEvR(validPairs), 35);
            [~, hi] = min(abs(fEvR(validPairs) - pct65));
            hiPair  = validPairs(hi);
            [~, lo] = min(abs(fEvR(validPairs) - pct35));
            loPair  = validPairs(lo);
        end

        % Colors
        cEv   = stimColors(sRipple,:);   % evoked — ripple colour
        cSp   = [0.55 0.55 0.55];         % spontaneous — grey
        hiMark= [0.85 0.15 0.15];         % red marker for high example
        loMark= [0.15 0.35 0.85];         % blue marker for low example

        figure('Units','normalized','Position',[0.05 0.10 0.85 0.70]);
        sgtitle(['AudRipple: Example Matched FME Pairs' titleSuffix], ...
            'FontSize',12,'FontWeight','bold');

        % --- Panel 1: Match quality scatter with highlighted examples ---
        ax1 = subplot(2, 3, [1 4]);
        hold(ax1,'on');
        scatter(ax1, fEvR, fSpR, 18, cEv, 'filled','MarkerFaceAlpha',0.35,'MarkerEdgeColor','none');
        lims = [min([fEvR; fSpR]), max([fEvR; fSpR])];
        plot(ax1, lims, lims, '--k','LineWidth',1.0);
        plot(ax1, lims, lims*(1+matchTol), ':','Color',[0.6 0.6 0.6],'LineWidth',0.8);
        plot(ax1, lims, lims*(1-matchTol), ':','Color',[0.6 0.6 0.6],'LineWidth',0.8);
        % Highlight selected pairs
        scatter(ax1, fEvR(hiPair), fSpR(hiPair), 80, hiMark, 'filled', ...
            'MarkerEdgeColor','k','LineWidth',1.2);
        scatter(ax1, fEvR(loPair), fSpR(loPair), 80, loMark, 'filled', ...
            'MarkerEdgeColor','k','LineWidth',1.2);
        text(ax1, fEvR(hiPair), fSpR(hiPair)+0.05*range(lims), 'A', ...
            'Color',hiMark,'FontSize',10,'FontWeight','bold','HorizontalAlignment','center');
        text(ax1, fEvR(loPair), fSpR(loPair)+0.05*range(lims), 'B', ...
            'Color',loMark,'FontSize',10,'FontWeight','bold','HorizontalAlignment','center');
        xlabel(ax1,'Evoked FME AUC (z·s)','FontSize',9);
        ylabel(ax1,'Matched spont. FME AUC (z·s)','FontSize',9);
        title(ax1,'Match quality scatter','FontSize',10,'FontWeight','bold');
        axis(ax1,'square'); box off;
        set(ax1,'TickDir','out','FontSize',9);

        % --- Panels 2 & 5: High example (A) — FME traces ---
        for ei = 1:2
            if ei == 1
                pairIdx = hiPair; mc = hiMark; lbl = 'A (higher AUC)';
                axFme = subplot(2, 3, 2); axAuc = subplot(2, 3, 5);
            else
                pairIdx = loPair; mc = loMark; lbl = 'B (lower AUC)';
                axFme = subplot(2, 3, 3); axAuc = subplot(2, 3, 6);
            end

            % FME traces
            hold(axFme,'on');
            trEv = TrEv(pairIdx,:);
            trSp = TrSp(pairIdx,:);
            hSp = plot(axFme, tResamp, trSp, '-','Color',cSp,'LineWidth',1.8);
            hEv = plot(axFme, tResamp, trEv, '-','Color',cEv,'LineWidth',2.0);
            xline(axFme, 0,'--k','Alpha',0.4,'LineWidth',0.8);
            yline(axFme, 0,'-','Color',[0.85 0.85 0.85],'LineWidth',0.7);
            % Shade AUC window
            yL = ylim(axFme);
            patch(axFme, [0 boutAUCWinSec boutAUCWinSec 0], ...
                [yL(1) yL(1) yL(2) yL(2)], ...
                [0.9 0.9 0.9],'FaceAlpha',0.20,'EdgeColor','none');
            xlim(axFme,[tResamp(1) tResamp(end)]);
            ylabel(axFme,'FME (z)','FontSize',9);
            title(axFme, sprintf('Example %s', lbl),'FontSize',9,'FontWeight','bold','Color',mc);
            legend(axFme, [hEv hSp], {'Evoked','Spontaneous'}, ...
                'FontSize',7.5,'Box','off','Location','best');
            box off; set(axFme,'TickDir','out','FontSize',8,'XTickLabel',[]);

            % AUC value annotations below
            hold(axAuc,'on');
            aucEv = fEvR(pairIdx);
            aucSp = fSpR(pairIdx);
            bar(axAuc, 1, aucEv, 0.4,'FaceColor',cEv,'EdgeColor','none','FaceAlpha',0.8);
            bar(axAuc, 2, aucSp, 0.4,'FaceColor',cSp,'EdgeColor','none','FaceAlpha',0.8);
            text(axAuc, 1, aucEv+0.02, sprintf('%.2f',aucEv), ...
                'HorizontalAlignment','center','FontSize',8,'FontWeight','bold','Color',cEv);
            text(axAuc, 2, aucSp+0.02, sprintf('%.2f',aucSp), ...
                'HorizontalAlignment','center','FontSize',8,'FontWeight','bold','Color',cSp);
            set(axAuc,'XTick',[1 2],'XTickLabel',{'Evoked','Spont.'}, ...
                'TickDir','out','FontSize',8);
            ylabel(axAuc,'FME AUC (z·s)','FontSize',8);
            xlim(axAuc,[0.4 2.6]); ylim(axAuc,[0 inf]);
            xlabel(axAuc,'Time re FME onset (s)','FontSize',8);
            box off;
        end

        fprintf('  Example pair A (high): pair %d, evAUC=%.2f, spAUC=%.2f\n', ...
            hiPair, fEvR(hiPair), fSpR(hiPair));
        fprintf('  Example pair B (low):  pair %d, evAUC=%.2f, spAUC=%.2f\n', ...
            loPair, fEvR(loPair), fSpR(loPair));
    end
end
%  Layout: 1 row x nS columns, each with grouped boxplots per intensity.
%  Threshold marked as horizontal dashed line.
%% =========================================================================
figure('Units','normalized','Position',[0.02 0.25 0.92 0.45]);
sgtitle(['FME AUC Distributions by Stimulus & Intensity' titleSuffix], ...
    'FontSize', 13, 'FontWeight', 'bold');

for s = 1:nS
    ax = subplot(1, nS, s);
    hold(ax,'on');

    dd     = diagData{s};
    uInts2 = sort(unique(dd.intAll));
    nLev2  = numel(uInts2);

    yMaxAll = fmeThreshold;

    for ui = 1:nLev2
        mask = dd.intAboveThresh == uInts2(ui);
        vals = dd.fmeAboveThresh(mask);
        if isempty(vals); continue; end

        % Boxplot at x position = ui
        bp = boxplot(ax, vals, 'Positions', ui, 'Widths', 0.55, ...
            'Colors', stimColors(s,:), 'Symbol', '', 'Whisker', 1.5);
        set(bp,'LineWidth',1.4);
        hBox = findobj(bp,'Tag','Box');
        if ~isempty(hBox)
            patch(ax, get(hBox,'XData'), get(hBox,'YData'), ...
                stimColors(s,:), 'FaceAlpha', 0.35, ...
                'EdgeColor', stimColors(s,:), 'LineWidth', 1.4);
        end

        yMaxAll = max(yMaxAll, prctile(vals, 95));

        % n annotation below x axis
        nAllI = sum(dd.intAll == uInts2(ui));
        text(ax, ui, -0.15, sprintf('n=%d', numel(vals)), ...
            'HorizontalAlignment','center','FontSize',7,'Color',[0.4 0.4 0.4],...
            'Units','data','Clipping','off');
        text(ax, ui, -0.30, sprintf('(/%d)', nAllI), ...
            'HorizontalAlignment','center','FontSize',7,'Color',[0.65 0.65 0.65],...
            'Units','data','Clipping','off');
    end

    % Threshold line
    yline(ax, fmeThreshold, '--k', 'LineWidth', 1.0, 'Alpha', 0.6, ...
        'Label', sprintf('threshold=%.2f', fmeThreshold), 'FontSize', 7);

    % x-axis labels
    if strcmp(intFields{s}, 'int')
        xLabs2 = arrayfun(@(v) sprintf('%g dB', v), uInts2, 'UniformOutput', false);
    else
        xLabs2 = arrayfun(@(v) sprintf('%.2f', v), uInts2, 'UniformOutput', false);
    end
    set(ax,'XTick',1:nLev2,'XTickLabel',xLabs2,'TickDir','out','FontSize',8);
    xtickangle(ax,20);
    xlim(ax,[0.4 nLev2+0.6]);
    ylim(ax,[0 yMaxAll*1.15]);
    if s==1; ylabel(ax,'FME AUC (z·s)','FontSize',9); end
    xlabel(ax,'Intensity','FontSize',9);
    title(ax, stimLabels{s},'FontSize',10,'FontWeight','bold');
    box off;
end

annotation('textbox',[0.01 0 0.99 0.04],'String', ...
    ['Boxplots of FME AUC for trials with a detected FME event, pooled across mice.  ' ...
     'n = trials above threshold / total trials at that intensity.  ' ...
     'Whiskers = 1.5×IQR.  ' ...
     sprintf('Dashed line = FME threshold (%.2f z\xB7s).', fmeThreshold)], ...
    'EdgeColor','none','FontSize',8,'Color',[0.4 0.4 0.4],'HorizontalAlignment','center');

%% =========================================================================
%  Figure 1: Scatter — FME AUC vs ACh AUC
%  Layout: nS rows x nR columns
%  Evoked = filled circles, Spontaneous = open circles
%% =========================================================================
figure('Units','normalized','Position',[0.02 0.05 0.95 0.85]);
sgtitle(['FME AUC vs ACh AUC: Evoked vs Matched Spontaneous' titleSuffix], ...
    'FontSize', 13, 'FontWeight', 'bold');

for s = 1:nS
    for r = 1:nR
        ax = subplot(nS, nR, (s-1)*nR + r);
        hold(ax, 'on');

        ev  = achAUC_evoked{s,r};
        sp  = achAUC_matched{s,r};
        fEv = fmeAUC_evoked{s};
        fSp = fmeAUC_matched{s};

        n = min([numel(ev), numel(sp), numel(fEv), numel(fSp)]);
        ev = ev(1:n); sp = sp(1:n);
        fEv = fEv(1:n); fSp = fSp(1:n);

        % Keep only finite pairs
        finPairs = isfinite(ev) & isfinite(sp);
        ev = ev(finPairs); sp = sp(finPairs);
        fEv = fEv(finPairs); fSp = fSp(finPairs);
        n = numel(ev);

        if n < 3
            title(ax, sprintf('%s – %s', regions{r}, stimLabels{s}), ...
                'FontSize', 9, 'Color', regionColors(r,:));
            continue;
        end

        c = regionColors(r,:);

        % Evoked — each point is one animal's mean
        scatter(ax, fEv, ev, markerSz+20, stimColors(s,:), 'filled', ...
            'MarkerFaceAlpha', 0.75, 'MarkerEdgeColor', 'none');

        % Matched spontaneous
        scatter(ax, fSp, sp, markerSz+20, [0.6 0.6 0.6], ...
            'MarkerEdgeColor', [0.4 0.4 0.4], 'LineWidth', 0.8);

        % Regression lines
        if n >= 5
            xRange = linspace(min([fEv(1:n); fSp(1:n)]), ...
                               max([fEv(1:n); fSp(1:n)]), 100);
            pEv = polyfit(fEv(1:n), ev(1:n), 1);
            pSp = polyfit(fSp(1:n), sp(1:n), 1);
            plot(ax, xRange, polyval(pEv, xRange), '-', ...
                'Color', stimColors(s,:), 'LineWidth', 1.5);
            plot(ax, xRange, polyval(pSp, xRange), '-', ...
                'Color', [0.5 0.5 0.5], 'LineWidth', 1.5);

            % Pearson r annotations
            [rEv, pEvV] = corr(fEv(1:n), ev(1:n), 'Type', 'Pearson');
            [rSp, pSpV] = corr(fSp(1:n), sp(1:n), 'Type', 'Pearson');
            text(ax, 0.05, 0.95, sprintf('Evoked: r=%.2f (p=%.3f)', rEv, pEvV), ...
                'Units','normalized','FontSize',7.5,'Color',stimColors(s,:));
            text(ax, 0.05, 0.85, sprintf('Spont: r=%.2f (p=%.3f)', rSp, pSpV), ...
                'Units','normalized','FontSize',7.5,'Color',[0.4 0.4 0.4]);
        end

        xline(ax, 0, '-', 'Color', [0.85 0.85 0.85], 'LineWidth', 0.8);
        yline(ax, 0, '-', 'Color', [0.85 0.85 0.85], 'LineWidth', 0.8);

        if r == 1; ylabel(ax, {stimLabels{s}; 'ACh AUC (z·s)'}, 'FontSize', 9); end
        if s == nS; xlabel(ax, 'FME AUC (z·s)', 'FontSize', 9); end
        title(ax, regions{r}, 'FontSize', 10, 'FontWeight', 'bold', ...
            'Color', regionColors(r,:));
        box off;
        set(ax, 'TickDir', 'out', 'FontSize', 9);

        if s == 1 && r == 1
            legH = [scatter(nan,nan,markerSz,stimColors(1,:),'filled'), ...
                    scatter(nan,nan,markerSz,[0.6 0.6 0.6], ...
                        'MarkerEdgeColor',[0.4 0.4 0.4],'LineWidth',0.8)];
            legend(legH, {'Evoked','Spontaneous (matched)'}, ...
                'Location','best','FontSize',8,'Box','off');
        end
    end
end

%% =========================================================================
%  Figure 2: Boxplots — ACh AUC evoked vs spontaneous per region
%  poolTrials=false (default): one value per animal, Wilcoxon vs zero on
%                              per-animal (evoked − spont) differences.
%  poolTrials=true            : all matched trial pairs pooled, paired
%                              Wilcoxon signed-rank on evoked vs spont.
statsRowsEvSpont = {};   % collect for CSV export
%  Layout: nS rows x nR columns.
%% =========================================================================
figure('Units','normalized','Position',[0.02 0.05 0.95 0.85]);
sgtitle(['ACh AUC: Stimulus-Evoked vs Matched Spontaneous FME' titleSuffix], ...
    'FontSize', 13, 'FontWeight', 'bold');

for s = 1:nS
    for r = 1:nR
        ax = subplot(nS, nR, (s-1)*nR + r);
        hold(ax, 'on');

        if poolTrials
            ev = achAUC_evoked{s,r};
            sp = achAUC_matched{s,r};
        else
            ev = achAUC_evoked_mn{s,r};
            sp = achAUC_matched_mn{s,r};
        end

        n = min(numel(ev), numel(sp));
        if n < 3
            title(ax, sprintf('%s – %s', regions{r}, stimLabels{s}), ...
                'FontSize', 9, 'Color', regionColors(r,:));
            continue;
        end

        ev = ev(1:n); sp = sp(1:n);
        finPairs = isfinite(ev) & isfinite(sp);
        ev = ev(finPairs); sp = sp(finPairs);
        n  = numel(ev);

        if n < 3
            title(ax, sprintf('%s – %s', regions{r}, stimLabels{s}), ...
                'FontSize', 9, 'Color', regionColors(r,:));
            continue;
        end

        c = stimColors(s,:);

        % Evoked boxplot at x=1
        bpEv = boxplot(ax, ev, 'Positions', 1, 'Widths', 0.45, ...
            'Colors', c, 'Symbol', '', 'Whisker', 1.5);
        set(bpEv, 'LineWidth', 1.5);
        hBox = findobj(bpEv, 'Tag', 'Box');
        if ~isempty(hBox)
            patch(ax, get(hBox,'XData'), get(hBox,'YData'), ...
                c, 'FaceAlpha', 0.35, 'EdgeColor', c, 'LineWidth', 1.5);
        end

        % Spontaneous boxplot at x=2
        bpSp = boxplot(ax, sp, 'Positions', 2, 'Widths', 0.45, ...
            'Colors', [0.5 0.5 0.5], 'Symbol', '', 'Whisker', 1.5);
        set(bpSp, 'LineWidth', 1.5);
        hBox = findobj(bpSp, 'Tag', 'Box');
        if ~isempty(hBox)
            patch(ax, get(hBox,'XData'), get(hBox,'YData'), ...
                [0.7 0.7 0.7], 'FaceAlpha', 0.35, ...
                'EdgeColor', [0.5 0.5 0.5], 'LineWidth', 1.5);
        end

        % Per-animal data points and connecting lines are intentionally NOT
        % drawn here — boxes only, to keep the panels clean for use as insets.
        % Set showPoints = true to restore them.
        showPoints = false;
        if showPoints
            if poolTrials
                nShow = min(n, 30);
                for i = 1:nShow
                    plot(ax, [1 2], [ev(i) sp(i)], '-', ...
                        'Color', [0.7 0.7 0.7 0.3], 'LineWidth', 0.6);
                end
            else
                for i = 1:n
                    plot(ax, [1 2], [ev(i) sp(i)], '-', ...
                        'Color', [0.5 0.5 0.5 0.5], 'LineWidth', 0.8);
                end
                scatter(ax, ones(n,1),   ev, 28, c,             'filled', 'MarkerFaceAlpha', 0.7);
                scatter(ax, 2*ones(n,1), sp, 28, [0.5 0.5 0.5], 'filled', 'MarkerFaceAlpha', 0.7);
            end
        end

        % Stats — BOTH tails reported.
        % Two-sided is the default/primary test. One-sided (evoked > spont) is
        % also reported because the manuscript's hypothesis is directional, and
        % because the two-sided test has a hard floor at small n: with n=5 the
        % smallest attainable two-sided p is 0.0625, so it CANNOT reach 0.05
        % however consistent the effect. One-sided halves that floor to 0.031.
        % CAUTION: a one-sided test in the "evoked > spont" direction will miss
        % an effect running the other way. PFC differences are consistently
        % NEGATIVE, so read pW_one there as uninformative, not as a null.
        if poolTrials
            pW     = signrank(ev, sp);
            pW_one = signrank(ev, sp, 'tail', 'right');
            nLabel = sprintf('n=%d pairs', n);
        else
            pW     = signrank(ev - sp);
            pW_one = signrank(ev - sp, 0, 'tail', 'right');
            nLabel = sprintf('n=%d mice', n);
        end
        % smallest attainable p at this n, for transparency
        pFloor2 = 2 / 2^n;
        pFloor1 = 1 / 2^n;
        if pW < 0.001;    stars = '***';
        elseif pW < 0.01; stars = '**';
        elseif pW < 0.05; stars = '*';
        else;              stars = 'ns';
        end

        % Store for CSV export
        diffVec = ev - sp;
        statsRowsEvSpont{end+1} = {stimLabels{s}, regions{r}, n, ...
            mean(ev,'omitnan'),   std(ev,0,'omitnan')/sqrt(n), ...
            mean(sp,'omitnan'),   std(sp,0,'omitnan')/sqrt(n), ...
            mean(diffVec,'omitnan'), std(diffVec,0,'omitnan')/sqrt(n), ...
            pW, pW_one, pFloor2, pFloor1}; %#ok

        yline(ax, 0, '--k', 'LineWidth', 0.8, 'Alpha', 0.4);
        drawnow;
        if showPoints
            yTop = max([ev; sp]) + 0.05 * range(ax.YLim);
        else
            % Boxes only: place the bracket just above the upper whiskers
            wUp  = findobj(ax, 'Tag', 'Upper Whisker');
            wTop = -Inf;
            for wi = 1:numel(wUp)
                yd = get(wUp(wi), 'YData');
                wTop = max(wTop, max(yd));
            end
            if ~isfinite(wTop); wTop = max([ev; sp]); end
            yTop = wTop + 0.05 * range(ax.YLim);
        end
        plot(ax, [1 2], [yTop yTop], '-k', 'LineWidth', 1.0);
        text(ax, 1.5, yTop + 0.02*range(ax.YLim), stars, ...
            'HorizontalAlignment','center','FontSize',11,...
            'FontWeight','bold','Color','k');
        text(ax, 1.5, ax.YLim(1) - 0.08*range(ax.YLim), ...
            sprintf('W: p=%.3f  %s', pW, nLabel), ...
            'HorizontalAlignment','center','FontSize',7.5,...
            'Color','k','Clipping','off');

        set(ax, 'XTick', [1 2], 'XTickLabel', {'Evoked','Spont.'}, ...
            'TickDir', 'out', 'FontSize', 9);
        xlim(ax, [0.5 2.5]);
        if r == 1; ylabel(ax, {stimLabels{s}; 'ACh AUC (z·s)'}, 'FontSize', 9); end
        title(ax, regions{r}, 'FontSize', 10, 'FontWeight', 'bold', ...
            'Color', regionColors(r,:));
        box off;
    end
end

if poolTrials
    footStr = ['Pooled matched trial pairs across animals.  ' ...
               sprintf('Spontaneous bouts matched to evoked on FME AUC (\\pm%.0f%% tolerance).  ', matchTol*100) ...
               'W = Wilcoxon signed-rank test on matched pairs.  ' ...
               sprintf('FME threshold: AUC > %.2f z\xB7s', fmeThreshold)];
else
    footStr = ['Per-animal means of matched trials.  ' ...
               sprintf('Spontaneous bouts matched to evoked on FME AUC (\\pm%.0f%% tolerance).  ', matchTol*100) ...
               'W vs 0 = Wilcoxon signed-rank on per-animal (evoked \x2212 spont) differences.  ' ...
               sprintf('FME threshold: AUC > %.2f z\xB7s', fmeThreshold)];
end
annotation('textbox', [0.01 0.00 0.99 0.04], 'String', footStr, ...
    'EdgeColor', 'none', 'FontSize', 8, 'Color', [0.4 0.4 0.4], ...
    'HorizontalAlignment', 'center');

%% =========================================================================
%  CSV export: evoked vs spontaneous ACh AUC statistics
%  FDR correction (Benjamini-Hochberg) across all stim × region combinations
%% =========================================================================
if ~isempty(statsRowsEvSpont)
    nT        = numel(statsRowsEvSpont);
    pRaws     = cellfun(@(r) r{10}, statsRowsEvSpont);   % two-sided
    [ps, idx] = sort(pRaws(:));
    adj       = ps .* nT ./ (1:nT)';
    adj       = flipud(cummin(flipud(adj)));
    adj       = min(adj, 1);
    pAdj      = nan(nT,1); pAdj(idx) = adj;

    % Preserve for the combined AUC + peak printout further below
    % FDR for the one-sided family as well
    pRaws1    = cellfun(@(r) r{11}, statsRowsEvSpont);
    [ps1,id1] = sort(pRaws1(:));
    adj1      = ps1 .* nT ./ (1:nT)';
    adj1      = flipud(cummin(flipud(adj1)));
    adj1      = min(adj1, 1);
    pAdj1     = nan(nT,1); pAdj1(id1) = adj1;

    aucRowsSaved = statsRowsEvSpont;
    aucPAdjSaved = pAdj;

    if poolTrials
        unitStr = 'pairs';
    else
        unitStr = 'mice';
    end
    csvPath = fullfile(pwd, sprintf('ach_evoked_vs_spont_stats%s.csv', onsetTag));
    fid = fopen(csvPath, 'w');
    if fid == -1
        % File may be locked (e.g. open in Excel) — try a timestamped name
        [pDir, pName, pExt] = fileparts(csvPath);
        csvPath = fullfile(pDir, sprintf('%s_%s%s', pName, ...
            datestr(now,'yyyymmdd_HHMMSS'), pExt));
        fid = fopen(csvPath, 'w');
    end
    if fid == -1
        warning(['Could not write evoked-vs-spontaneous CSV. ' ...
                 'Check that the file is not open in another program ' ...
                 'and that %s is writable.'], pwd);
    else
    fprintf(fid, ['stimulus,region,n_%s,' ...
                  'mean_evoked,sem_evoked,' ...
                  'mean_spont,sem_spont,' ...
                  'mean_diff,sem_diff,' ...
                  'p_raw_twosided,p_FDR_twosided,sig_FDR05_twosided,' ...
                  'p_raw_onesided,p_FDR_onesided,sig_FDR05_onesided,' ...
                  'p_floor_twosided,p_floor_onesided\n'], unitStr);
    for k = 1:nT
        row = statsRowsEvSpont{k};
        fprintf(fid, ['%s,%s,%d,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,' ...
                      '%.4f,%.4f,%d,%.4f,%.4f,%d,%.4f,%.4f\n'], ...
            row{1}, row{2}, row{3}, ...
            row{4}, row{5}, row{6}, row{7}, row{8}, row{9}, ...
            row{10}, pAdj(k),  pAdj(k)<0.05, ...
            row{11}, pAdj1(k), pAdj1(k)<0.05, ...
            row{12}, row{13});
    end
    fclose(fid);
    fprintf('Evoked vs spontaneous ACh stats saved to: %s\n', csvPath);
    end
end


%% =========================================================================
%  Matching attrition: which mice reached the matching stage, and why
%  ------------------------------------------------------------------------
%  The matching analysis uses fewer mice than the trial-structured cohort.
%  This reports exactly who survived, per stimulus type and region, with the
%  bout counts behind it, so the exclusion can be described accurately in the
%  Methods rather than inferred.
%
%  Gates a mouse must pass, in order:
%    1. a *_traces.mat file exists in save_path
%    2. >=1 spontaneous FME bout survives detection
%       (>85th pct, duration < 1 s, FME AUC > fmeThreshold)
%    3. the region appears in all_fiber_names_ach{m}
%    4. >=1 evoked trial has a detectable FME onset above threshold
%       (with risingEdgeOnly, the onset must also be a RISING EDGE)
%    5. >=1 evoked trial finds spontaneous partners BOTH above and below it
%       within +/-matchTol  <-- usually the binding constraint
%    6. the paired ACh AUCs are finite
%
%  Reading the output: gates 4-5 depend only on FME, which is a WHOLE-ANIMAL
%  signal. So a mouse that reaches matching in one region but not another, for
%  the same stimulus, cannot have failed on FME — its failure is gate 6, a bad
%  ACh channel in the failing region. Only a mouse that drops out of ALL its
%  regions is a candidate for an FME-side cause.
%
%  Where gate 5 IS the cause, note that it is not random with respect to the
%  data: an animal whose evoked movements were consistently LARGER than any of
%  its spontaneous movements has no partner above, so every trial fails — and
%  that is precisely the animal with the strongest FME response to sound.
%  Treat that as a selection effect, not as neutral attrition, when writing up.
%% =========================================================================
if ~isempty(bootIn)
    fprintf('\n==================================================================\n');
    fprintf('Mice reaching the matching stage\n');
    fprintf('==================================================================\n');

    for s = 1:nS
        for r = 1:nR
            sel = find(cellfun(@(b) b.stim==s && b.region==r, bootIn));
            if isempty(sel)
                fprintf('\n%-13s %-4s : no mice reached matching\n', ...
                    stimLabels{s}, regions{r});
                continue;
            end
            fprintf('\n%-13s %-4s : %d mice\n', stimLabels{s}, regions{r}, numel(sel));
            fprintf('   %6s %10s %10s %12s %12s\n', ...
                'mouse','nEvoked','nSpont','med evFME','med spFME');
            for k = sel(:)'
                B = bootIn{k};
                fprintf('   %6d %10d %10d %12.3f %12.3f\n', ...
                    B.mouse, numel(B.fmeEv), numel(B.fmeSil), ...
                    median(B.fmeEv,'omitnan'), median(B.fmeSil,'omitnan'));
            end

            % Who had the region but did not reach matching?
            had = false(1, numel(ach_stim));
            for m = 1:numel(ach_stim)
                if isempty(ach_stim{m}); continue; end
                if numel(all_fiber_names_ach) < m; continue; end
                had(m) = any(strcmp(all_fiber_names_ach{m}, regions{r}));
            end
            reached = arrayfun(@(k) bootIn{k}.mouse, sel);
            missing = find(had & ~ismember(1:numel(ach_stim), reached));
            if ~isempty(missing)
                fprintf('   had %s but did NOT reach matching: %s\n', ...
                    regions{r}, mat2str(missing));
                % Distinguish the two failure modes. Matching (gates 4-5) uses
                % only FME, which is a WHOLE-ANIMAL signal — so if a mouse
                % reached matching for a different region with the same
                % stimulus, its FME must have matched fine, and the failure has
                % to be gate 6 (non-finite ACh AUCs) i.e. a bad channel in THIS
                % region. Only a mouse that failed in every one of its regions
                % is a candidate for an FME-side cause.
                for mm = missing
                    otherOK = false;
                    for rr = 1:nR
                        if rr == r; continue; end
                        if any(cellfun(@(b) b.stim==s && b.region==rr && ...
                                b.mouse==mm, bootIn))
                            otherOK = true; break;
                        end
                    end
                    if otherOK
                        fprintf('     mouse %d: reached matching in another region ' , mm);
                        fprintf('for this stimulus\n');
                        fprintf('       -> FME matched fine; failure is the ACh channel\n');
                        fprintf('          in %s (gate 6, non-finite ACh AUC)\n', regions{r});
                    else
                        fprintf('     mouse %d: failed in all its regions for this stimulus\n', mm);
                        fprintf('       -> consistent with an FME-side cause (gate 4 or 5):\n');
                        fprintf('          too few supra-threshold evoked bouts, or evoked\n');
                        fprintf('          bouts with no spontaneous partner above them\n');
                    end
                end
            end
        end
    end
    fprintf('\n==================================================================\n');
end

%% =========================================================================
%  Bootstrap of the matching procedure (R3 Major Comment 2c)
%  ------------------------------------------------------------------------
%  "If using a random selection, bootstrapping the matching procedure would
%   strengthen the result."
%
%  The magnitude-balanced match draws randomly: for each evoked bout it picks
%  one magnitude-matched pair at random, then the above- or below-member with
%  equal probability. Under a fixed seed the reported result is therefore ONE
%  realisation. This repeats the whole matching nBootMatch times with
%  different draws and reports the distribution of the per-animal
%  evoked-minus-spontaneous ACh AUC difference.
%
%  What this does and does not show: it establishes that the result does not
%  depend on the particular random draw. It is NOT a bootstrap over animals,
%  so it does not widen the confidence interval for sampling of mice — the
%  paired Wilcoxon on per-animal means remains the inferential test.
%
%  Nothing here feeds back into the reported analysis; it is a stability check.
%% =========================================================================
nBootMatch = 500;

if ~isempty(bootIn)
    fprintf('\nBootstrapping the matching procedure (%d re-draws)...\n', nBootMatch);

    % diffBoot{s,r} = [nBoot x nMice] per-animal mean difference per re-draw
    diffBoot = cell(nS, nR);
    for s = 1:nS
        for r = 1:nR
            diffBoot{s,r} = [];
        end
    end

    rngState = rng;          % preserve the caller's stream
    rng(12345);              % own reproducible stream for the bootstrap

    for b = 1:nBootMatch
        % accumulate one per-animal difference per (stim, region)
        acc = cell(nS, nR);
        for k = 1:numel(bootIn)
            B   = bootIn{k};
            nEvB = numel(B.fmeEv);
            dVals = nan(nEvB,1);

            for i = 1:nEvB
                tolWindow = matchTol * abs(B.fmeEv(i));
                allCands  = find(abs(B.fmeSil - B.fmeEv(i)) <= tolWindow);
                if isempty(allCands); continue; end
                aboveIdx = allCands(B.fmeSil(allCands) > B.fmeEv(i));
                belowIdx = allCands(B.fmeSil(allCands) < B.fmeEv(i));
                if isempty(aboveIdx) || isempty(belowIdx); continue; end

                aboveDist = B.fmeSil(aboveIdx) - B.fmeEv(i);
                belowDist = B.fmeEv(i) - B.fmeSil(belowIdx);

                nA = numel(aboveIdx);
                pairedAbove = zeros(nA,1); pairedBelow = zeros(nA,1);
                for ai = 1:nA
                    [~, bi] = min(abs(belowDist - aboveDist(ai)));
                    pairedAbove(ai) = aboveIdx(ai);
                    pairedBelow(ai) = belowIdx(bi);
                end

                pairSel = randi(nA);
                if rand() < 0.5
                    mIdx = pairedAbove(pairSel);
                else
                    mIdx = pairedBelow(pairSel);
                end

                if i <= numel(B.achEv) && mIdx <= numel(B.achSil)
                    dVals(i) = B.achEv(i) - B.achSil(mIdx);
                end
            end

            dVals = dVals(isfinite(dVals));
            if ~isempty(dVals)
                acc{B.stim, B.region}(end+1) = mean(dVals, 'omitnan');
            end
        end

        for s = 1:nS
            for r = 1:nR
                if ~isempty(acc{s,r})
                    % store as a row; pad with NaN so a re-draw that loses a
                    % mouse (no valid match) does not silently become a zero
                    row = nan(1, nMice);
                    row(1:numel(acc{s,r})) = acc{s,r};
                    diffBoot{s,r} = [diffBoot{s,r}; row];
                end
            end
        end
    end

    % trim the NaN padding to the widest actually used
    for s = 1:nS
        for r = 1:nR
            if isempty(diffBoot{s,r}); continue; end
            used = find(any(isfinite(diffBoot{s,r}), 1), 1, 'last');
            diffBoot{s,r} = diffBoot{s,r}(:, 1:used);
        end
    end

    rng(rngState);   % restore

    % ---- Report ----
    fprintf('\n==================================================================\n');
    fprintf('Matching bootstrap: per-animal ACh AUC difference (evoked - spont)\n');
    fprintf('%d independent re-draws of the magnitude-balanced match\n', nBootMatch);
    fprintf('==================================================================\n');
    fprintf('%-13s %-5s %6s %10s %10s %10s %8s %8s %9s %8s %9s\n', ...
        'Stimulus','Reg','nMice','mean diff','2.5 pct','97.5 pct', ...
        '%% same','SE ratio','%% 2-sided','%% 1-sided','min p');

    bootRows = {};
    for s = 1:nS
        for r = 1:nR
            D = diffBoot{s,r};
            if isempty(D); continue; end
            % mean across animals for each re-draw
            perDraw = mean(D, 2, 'omitnan');
            perDraw = perDraw(isfinite(perDraw));
            if numel(perDraw) < 10; continue; end

            % Wilcoxon on each re-draw's per-animal differences
            pDraw  = nan(size(D,1),1);
            pDraw1 = nan(size(D,1),1);
            nUse   = 0;
            for b = 1:size(D,1)
                v = D(b,:); v = v(isfinite(v));
                if numel(v) >= 3
                    pDraw(b)  = signrank(v);
                    pDraw1(b) = signrank(v, 0, 'tail', 'right');
                    nUse = numel(v);
                end
            end
            pDraw   = pDraw(isfinite(pDraw));
            pDraw1  = pDraw1(isfinite(pDraw1));
            fracSig  = 100 * mean(pDraw  < 0.05);
            fracSig1 = 100 * mean(pDraw1 < 0.05);

            % --- Two self-contained stability measures -------------------
            % (1) SIGN CONSISTENCY: proportion of re-draws in which the mean
            %     difference kept the same sign as the full-cohort mean. This
            %     needs no reference range — "the sign was identical in 500 of
            %     500 re-draws" is directly evaluable, and it is what the
            %     reviewer's question ("does the result depend on the random
            %     match?") actually asks.
            signRef  = sign(mean(perDraw, 'omitnan'));
            pctSign  = 100 * mean(sign(perDraw) == signRef);

            % (2) SE RATIO: both terms are uncertainties in the SAME
            %     quantity — the cohort mean difference — so the ratio is
            %     like-for-like.
            %
            %       sdRedraw : SD of the cohort mean across re-draws.
            %                  D is [nRedraws x nMice]; mean across animals
            %                  gives one cohort mean per draw; this is its SD.
            %                  = uncertainty from the RANDOM MATCH.
            %       seAnimals: SEM across animals within a draw, averaged over
            %                  draws (SD across animals / sqrt(n)).
            %                  = uncertainty from WHICH ANIMALS were sampled.
            %
            %     An earlier version divided by the raw SD across animals.
            %     That compared a standard error against a standard deviation,
            %     so the ratio carried a 1/sqrt(n) factor unrelated to matching
            %     stability (0.33 at n = 9 from that alone). The SEM version
            %     removes it: a ratio well below 1 means re-drawing the match
            %     perturbs the cohort mean far less than the choice of animals
            %     does.
            %
            %     sdAnimalsRaw is retained for the CSV, since the raw spread
            %     between animals is worth reporting in its own right.
            sdRedraw     = std(perDraw, 0, 'omitnan');
            sdAnimalsRaw = mean(std(D, 0, 2, 'omitnan'), 'omitnan');
            nAnimHere    = size(D, 2);
            seAnimals    = sdAnimalsRaw / sqrt(max(nAnimHere,1));
            if isfinite(seAnimals) && seAnimals > 0
                spreadRatio = sdRedraw / seAnimals;
            else
                spreadRatio = NaN;
            end
            floor2   = 2 / 2^max(nUse,1);
            floor1   = 1 / 2^max(nUse,1);

            fprintf('%-13s %-5s %6d %10.4f %10.4f %10.4f %7.1f%% %8.2f %8.1f%% %7.1f%% %9.4f%s\n', ...
                stimLabels{s}, regions{r}, size(D,2), ...
                mean(perDraw), prctile(perDraw,2.5), prctile(perDraw,97.5), ...
                pctSign, spreadRatio, fracSig, fracSig1, floor2, ...
                char(' ' + (floor2 > 0.05)*('!'-' ')));

            bootRows{end+1} = {stimLabels{s}, regions{r}, size(D,2), ...
                mean(perDraw), prctile(perDraw,2.5), prctile(perDraw,97.5), ...
                fracSig, median(pDraw), fracSig1, median(pDraw1), floor2, ...
                pctSign, spreadRatio, sdRedraw, seAnimals, sdAnimalsRaw}; %#ok
        end
    end
    fprintf('==================================================================\n');
    fprintf('"%% significant" = proportion of re-draws in which the paired\n');
    fprintf('  Wilcoxon across animals reached p<0.05. Near 100%% means the\n');
    fprintf('  conclusion does not depend on the particular random match.\n');
    fprintf('"%% same" = re-draws keeping the sign of the full-cohort mean.\n');
    fprintf('  Self-contained: needs no reference range, and is the direct\n');
    fprintf('  answer to "does the result depend on the random match?"\n');
    fprintf('"SE ratio" = (SD of the cohort mean across re-draws) / (SEM\n');
    fprintf('  across animals). Both are uncertainties in the SAME quantity,\n');
    fprintf('  so the comparison is like-for-like. Well below 1 means\n');
    fprintf('  re-drawing the match perturbs the cohort mean far less than\n');
    fprintf('  the choice of animals does. (Dividing by the raw SD instead\n');
    fprintf('  would carry a 1/sqrt(n) factor unrelated to stability.)\n');
    fprintf('"min p" = smallest TWO-SIDED p attainable at that n. A "!" marks\n');
    fprintf('  rows where it exceeds 0.05, so the two-sided test CANNOT reach\n');
    fprintf('  significance however consistent the effect (n=5 -> 0.0625).\n');
    fprintf('One-sided tests "evoked > spont". PFC differences are NEGATIVE, so\n');
    fprintf('  the one-sided column is uninformative there, not a null result.\n');
    fprintf('==================================================================\n');

    % ---- CSV ----
    if ~isempty(bootRows)
        csvB = fullfile(pwd, sprintf('ach_matching_bootstrap%s.csv', onsetTag));
        fidB = fopen(csvB, 'w');
        if fidB == -1
            warning('Could not write matching bootstrap CSV.');
        else
            fprintf(fidB, ['stimulus,region,n_mice,n_redraws,' ...
                           'mean_diff,ci_lo_2p5,ci_hi_97p5,' ...
                           'pct_redraws_same_sign,se_ratio_redraw_over_animals,' ...
                           'sd_across_redraws,sem_across_animals,sd_across_animals,' ...
                           'pct_redraws_sig_twosided,median_p_twosided,' ...
                           'pct_redraws_sig_onesided,median_p_onesided,' ...
                           'min_attainable_p_twosided\n']);
            for k = 1:numel(bootRows)
                rw = bootRows{k};
                fprintf(fidB, ['%s,%s,%d,%d,%.6f,%.6f,%.6f,' ...
                               '%.2f,%.4f,%.6f,%.6f,%.6f,' ...
                               '%.2f,%.6f,%.2f,%.6f,%.6f\n'], ...
                    rw{1}, rw{2}, rw{3}, nBootMatch, ...
                    rw{4}, rw{5}, rw{6}, ...
                    rw{12}, rw{13}, rw{14}, rw{15}, rw{16}, ...
                    rw{7}, rw{8}, rw{9}, rw{10}, rw{11});
            end
            fclose(fidB);
            fprintf('Matching bootstrap saved to: %s\n', csvB);
        end
    end
end

%% =========================================================================
%  Sanity check figure: mean ACh traces, grouped by stimulus modality
%  Layout: 2 rows x nR columns.
%    Row 1 = auditory stimuli (ripples solid, tones dashed)
%    Row 2 = visual stimuli   (gratings solid)
%  Grouping by modality keeps at most two evoked traces per panel, which
%  stays legible at publication size (the previous 1 x nR layout put three
%  same-hue traces in each panel, distinguishable only by line style).
%  Spontaneous: grey mean+/-SEM pooled across the stimulus types in that row.
%  Pre-onset baseline correction applied (display only).
%  UNITS DIFFER BETWEEN PANELS BY DESIGN:
%    traces here  -> ALL TRIALS pooled across mice; band = SEM across TRIALS
%    AUC boxplots -> PER ANIMAL; all statistics are per animal
%  The traces show the data; the boxplots carry the inference. State BOTH n's
%  in the figure legend, or a reader will assume the bands and the p-values
%  refer to the same unit.
%  Significance labels (added later) come from the peak-amplitude test:
%  per-animal mean within a fixed window around the grand-average peak,
%  paired Wilcoxon, FDR-corrected across stimulus type x region.
%% =========================================================================
tResampPlot = tResamp;   % t=0 = bout onset
lineStyles  = {'-', '--', ':'};
cSpont      = [0.6 0.6 0.6];

% Assign each stimulus type to a modality row by name
audIdx = find(contains(stimTypes, {'AudRipple','Tones'}));
visIdx = find(contains(stimTypes, {'VisGrating'}));
rowStimIdx   = {audIdx, visIdx};
rowStimNames = {'Auditory', 'Visual'};
rowStimIdx   = rowStimIdx(~cellfun(@isempty, rowStimIdx));
rowStimNames = rowStimNames(1:numel(rowStimIdx));
nRowsFig     = numel(rowStimIdx);

figure('Units','normalized','Position',[0.05 0.10 0.80 0.78]);
sgtitle(['ACh Traces: Evoked vs Matched Spontaneous' titleSuffix], ...
    'FontSize', 13, 'FontWeight', 'bold');

axTrace = gobjects(nRowsFig, nR);
preOnsetMask = tResampPlot >= -preFmeSec & tResampPlot < 0;

for rowI = 1:nRowsFig
    sList = rowStimIdx{rowI};

    for r = 1:nR
        ax = subplot(nRowsFig, nR, (rowI-1)*nR + r);
        axTrace(rowI,r) = ax;
        hold(ax, 'on');

        legHandles = gobjects(0);
        legLabels  = {};

        % --- Spontaneous: per-animal mean traces, AVERAGED over this row's
        %     stimulus types ---
        % Uses achTrace_matched_mn (one row per MOUSE). The row's stim types
        % are averaged WITHIN mouse, not concatenated: matching is within
        % mouse, so the same animals supply the matched spontaneous bouts for
        % every stimulus type. Concatenating counted each mouse once per
        % stimulus type (n=16 from 8 mice in the auditory row) and shrank the
        % SEM band by ~sqrt(2). Averaging keeps n = number of mice, matching
        % the evoked traces and the unit of analysis used for all statistics.
        % TRIAL-LEVEL: every matched spontaneous bout in this row's stimulus
        % types, pooled across mice. n = trials, so the band is SEM across
        % trials. The AUC boxplots and all statistics remain PER ANIMAL — the
        % two panels deliberately use different units and the legend must say so.
        % Note: a given spontaneous bout can be matched to trials of more than
        % one stimulus type, so rows pooled across stim types are not fully
        % independent.
        allSpont = [];
        for si = 1:numel(sList)
            Msp = achTrace_matched{sList(si), r};
            if isempty(Msp); continue; end
            Msp = Msp(all(isfinite(Msp),2), :);
            if isempty(Msp); continue; end
            blSp = mean(Msp(:, preOnsetMask), 2, 'omitnan');
            allSpont = [allSpont; Msp - blSp]; %#ok
        end
        if size(allSpont,1) >= 2
            muSp  = mean(allSpont, 1, 'omitnan');
            semSp = std(allSpont,  0, 1, 'omitnan') / sqrt(size(allSpont,1));
            fill(ax, [tResampPlot, fliplr(tResampPlot)], ...
                [muSp+semSp, fliplr(muSp-semSp)], ...
                cSpont, 'FaceAlpha', 0.20, 'EdgeColor', 'none');
            hSp = plot(ax, tResampPlot, muSp, '-', 'Color', cSpont, 'LineWidth', 1.8);
            legHandles(end+1) = hSp; %#ok
            legLabels{end+1}  = sprintf('Spontaneous (n=%d trials)', size(allSpont,1)); %#ok
        end

        % --- Evoked traces for the stim types in THIS row ---
        evC = regionColors(r,:);
        for k = 1:numel(sList)
            s = sList(k);
            Mev = achTrace_evoked{s,r};       % TRIAL-LEVEL
            Msp = achTrace_matched{s,r};
            if isempty(Mev) || isempty(Msp); continue; end
            nBoth = min(size(Mev,1), size(Msp,1));
            Mev = Mev(1:nBoth,:); Msp = Msp(1:nBoth,:);
            finBoth = all(isfinite(Mev),2) & all(isfinite(Msp),2);
            Mev = Mev(finBoth,:);
            n = size(Mev,1);
            if n < 2; continue; end
            blEv = mean(Mev(:, preOnsetMask), 2, 'omitnan');
            Mev  = Mev - blEv;

            muEv  = mean(Mev, 1, 'omitnan');
            semEv = std(Mev,  0, 1, 'omitnan') / sqrt(n);

            fill(ax, [tResampPlot, fliplr(tResampPlot)], ...
                [muEv+semEv, fliplr(muEv-semEv)], ...
                evC, 'FaceAlpha', 0.15, 'EdgeColor', 'none');
            hEv = plot(ax, tResampPlot, muEv, lineStyles{min(k,numel(lineStyles))}, ...
                'Color', evC, 'LineWidth', 2.0);

            legHandles(end+1) = hEv; %#ok
            legLabels{end+1}  = sprintf('%s (n=%d trials)', stimLabels{s}, n); %#ok
        end

        xline(ax, 0, '--k', 'Alpha', 0.4, 'LineWidth', 0.8);
        yline(ax, 0, '-',  'Color', [0.85 0.85 0.85], 'LineWidth', 0.7);
        xlim(ax, [tResampPlot(1) tResampPlot(end)]);

        if rowI == nRowsFig
            xlabel(ax, 'Time re FME onset (s)', 'FontSize', 9);
        end
        if r == 1
            ylabel(ax, {rowStimNames{rowI}; 'ACh (z)'}, ...
                'FontSize', 9, 'FontWeight', 'bold');
        end
        if rowI == 1
            title(ax, regions{r}, 'FontSize', 11, 'FontWeight', 'bold', ...
                'Color', regionColors(r,:));
        end
        if ~isempty(legHandles)
            legend(ax, legHandles, legLabels, 'FontSize', 7.5, 'Box', 'off', ...
                'Location', 'best');
        end
        box off;
        set(ax, 'TickDir', 'out', 'FontSize', 9);
    end
end

% Shared y-axis within each region column so rows are directly comparable
drawnow;
for r = 1:nR
    yMinC = Inf; yMaxC = -Inf;
    for rowI = 1:nRowsFig
        if ~isgraphics(axTrace(rowI,r)); continue; end
        yl = ylim(axTrace(rowI,r));
        yMinC = min(yMinC, yl(1)); yMaxC = max(yMaxC, yl(2));
    end
    for rowI = 1:nRowsFig
        if ~isgraphics(axTrace(rowI,r)); continue; end
        ylim(axTrace(rowI,r), [yMinC yMaxC]);
    end
end

%% =========================================================================
%  Peak amplitude comparison: evoked vs matched spontaneous
%  ---------------------------------------------------------------------
%  The peak window is defined from the GRAND AVERAGE of both conditions
%  combined (per stimulus type x region), then each animal's mean trace is
%  averaged within that same fixed window. Using a common, condition-blind
%  window avoids the upward bias that arises from taking each animal's own
%  maximum: the max of a noisy trace grows with noise, and the spontaneous
%  condition averages over many more bouts than the evoked condition, so
%  per-animal maxima would not be comparable between conditions.
%  Unit of analysis: per-animal means. Paired Wilcoxon signed-rank.
%  FDR (Benjamini-Hochberg) across all stimulus type x region comparisons.
%% =========================================================================
peakHalfWin  = 0.10;   % +/- s around grand-average peak time
earlyAUCWin  = 0.25;   % early AUC window: 0 to this many seconds from onset
                       % (captures the early divergence visible in the traces,
                       %  which the full 0-1 s window dilutes)

peakRows   = {};
peakSigStr = cell(nS, nR);   % for annotating the trace figure
for s = 1:nS
    for r = 1:nR
        peakSigStr{s,r} = '';
    end
end

for s = 1:nS
    for r = 1:nR
        Ev = achTrace_evoked_mn{s,r};
        Sp = achTrace_matched_mn{s,r};
        if isempty(Ev) || isempty(Sp); continue; end
        nMouse = min(size(Ev,1), size(Sp,1));
        if nMouse < 3; continue; end
        Ev = Ev(1:nMouse,:); Sp = Sp(1:nMouse,:);

        % Pre-onset baseline correction (matches the display convention)
        preM = tResamp >= -preFmeSec & tResamp < 0;
        Ev = Ev - mean(Ev(:,preM), 2, 'omitnan');
        Sp = Sp - mean(Sp(:,preM), 2, 'omitnan');

        % Grand average across BOTH conditions -> condition-blind peak time
        grand = mean([Ev; Sp], 1, 'omitnan');
        postM = tResamp >= 0 & tResamp <= boutAUCWinSec;
        grandPost = grand;
        grandPost(~postM) = -Inf;
        [~, pkIdx] = max(grandPost);
        tPeak = tResamp(pkIdx);

        % Fixed window around that peak, applied identically to both
        pkWin = tResamp >= (tPeak - peakHalfWin) & tResamp <= (tPeak + peakHalfWin);
        if sum(pkWin) < 2; continue; end

        evPk = mean(Ev(:, pkWin), 2, 'omitnan');
        spPk = mean(Sp(:, pkWin), 2, 'omitnan');
        fin  = isfinite(evPk) & isfinite(spPk);
        evPk = evPk(fin); spPk = spPk(fin);
        nMouse = numel(evPk);
        if nMouse < 3; continue; end

        pPk = signrank(evPk, spPk);

        % --- Early-window AUC (0 to earlyAUCWin) on the same animals ---
        earlyM = tResamp >= 0 & tResamp <= earlyAUCWin;
        if sum(earlyM) >= 2
            evEa = trapz(tResamp(earlyM), Ev(:, earlyM), 2);
            spEa = trapz(tResamp(earlyM), Sp(:, earlyM), 2);
            evEa = evEa(fin); spEa = spEa(fin);
            pEa  = signrank(evEa, spEa);
            muEvEa = mean(evEa,'omitnan');  seEvEa = std(evEa,0,'omitnan')/sqrt(nMouse);
            muSpEa = mean(spEa,'omitnan');  seSpEa = std(spEa,0,'omitnan')/sqrt(nMouse);
            muDfEa = mean(evEa-spEa,'omitnan');
        else
            pEa = NaN; muEvEa = NaN; seEvEa = NaN;
            muSpEa = NaN; seSpEa = NaN; muDfEa = NaN;
        end

        peakRows{end+1} = {stimLabels{s}, regions{r}, nMouse, tPeak, ...
            mean(evPk,'omitnan'), std(evPk,0,'omitnan')/sqrt(nMouse), ...
            mean(spPk,'omitnan'), std(spPk,0,'omitnan')/sqrt(nMouse), ...
            mean(evPk-spPk,'omitnan'), pPk, s, r, ...
            muEvEa, seEvEa, muSpEa, seSpEa, muDfEa, pEa}; %#ok
    end
end

% FDR across all stim x region comparisons
if ~isempty(peakRows)
    nPk       = numel(peakRows);
    pRawPk    = cellfun(@(x) x{10}, peakRows);
    [ps, idx] = sort(pRawPk(:));
    adj       = ps .* nPk ./ (1:nPk)';
    adj       = flipud(cummin(flipud(adj)));
    adj       = min(adj, 1);
    pAdjPk    = nan(nPk,1); pAdjPk(idx) = adj;

    % FDR for the early-AUC family (its own family of stim x region tests)
    pRawEa = cellfun(@(x) x{18}, peakRows);
    finEa  = isfinite(pRawEa);
    pAdjEa = nan(nPk,1);
    if any(finEa)
        pe        = pRawEa(finEa);
        nEa       = numel(pe);
        [psE,idE] = sort(pe(:));
        adjE      = psE .* nEa ./ (1:nEa)';
        adjE      = flipud(cummin(flipud(adjE)));
        adjE      = min(adjE, 1);
        tmpE      = nan(nEa,1); tmpE(idE) = adjE;
        pAdjEa(finEa) = tmpE;
    end

    for k = 1:nPk
        sIdx = peakRows{k}{11};
        rIdx = peakRows{k}{12};
        if     pAdjPk(k) < 0.001; st = '***';
        elseif pAdjPk(k) < 0.01;  st = '**';
        elseif pAdjPk(k) < 0.05;  st = '*';
        else;                      st = 'n.s.';
        end
        peakSigStr{sIdx,rIdx} = st;
    end

    % ------------------------------------------------------------------
    %  Combined printout: AUC (0-1 s) and PEAK amplitude side by side
    %  Both are per-animal means, paired Wilcoxon, FDR-corrected within
    %  their own family of stimulus type x region comparisons.
    % ------------------------------------------------------------------
    fprintf('\n==========================================================================\n');
    fprintf('ACh: EVOKED vs MATCHED SPONTANEOUS  (unit = per-animal means)\n');
    fprintf('  AUC  : integrated 0-%.1f s from bout onset\n', boutAUCWinSec);
    fprintf('  AUC-e: integrated 0-%.2f s from bout onset (early window)\n', earlyAUCWin);
    fprintf('  PEAK : mean within +/-%.2f s of grand-average peak (condition-blind)\n', peakHalfWin);
    fprintf('  Paired Wilcoxon signed-rank; raw and FDR-corrected p both shown.\n');
    fprintf('==========================================================================\n');
    fprintf('%-13s %-4s %3s | %-34s | %-34s | %-34s\n', ...
        'Stimulus','Reg','n', ...
        sprintf('AUC 0-%.2gs ev/sp/diff (praw/pFDR)', boutAUCWinSec), ...
        sprintf('AUC 0-%.2gs ev/sp/diff (praw/pFDR)', earlyAUCWin), ...
        'PEAK       ev/sp/diff (praw/pFDR)');
    fprintf('%s\n', repmat('-', 1, 130));

    for k = 1:nPk
        rw   = peakRows{k};
        sIdx = rw{11};  rIdx = rw{12};

        % Find the matching AUC row for this stimulus type x region
        aucStr = '               --                 ';
        if exist('aucRowsSaved','var')
            for j = 1:numel(aucRowsSaved)
                if strcmp(aucRowsSaved{j}{1}, rw{1}) && strcmp(aucRowsSaved{j}{2}, rw{2})
                    aw = aucRowsSaved{j};
                    aucStr = sprintf('%+5.2f/%+5.2f/%+5.2f (2s %.3f | 1s %.3f)%s', ...
                        aw{4}, aw{6}, aw{8}, aucPAdjSaved(j), aw{11}, ...
                        char(' ' + (aucPAdjSaved(j)<0.05)*('*'-' ')));
                    break;
                end
            end
        end

        pkStr = sprintf('%+5.2f/%+5.2f/%+5.2f (%.3f/%.3f)%s', ...
            rw{5}, rw{7}, rw{9}, rw{10}, pAdjPk(k), ...
            char(' ' + (pAdjPk(k)<0.05)*('*'-' ')));

        if isfinite(rw{18})
            eaStr = sprintf('%+5.2f/%+5.2f/%+5.2f (%.3f/%.3f)%s', ...
                rw{13}, rw{15}, rw{17}, rw{18}, pAdjEa(k), ...
                char(' ' + (pAdjEa(k)<0.05)*('*'-' ')));
        else
            eaStr = '               --                 ';
        end

        fprintf('%-13s %-4s %3d | %-34s | %-34s | %-34s\n', ...
            rw{1}, rw{2}, rw{3}, aucStr, eaStr, pkStr);
    end

    fprintf('%s\n', repmat('-', 1, 130));
    fprintf('Cells show: evoked / spontaneous / difference  (p_raw / p_FDR)\n');
    fprintf('* = p_FDR < 0.05 (FDR corrected within each metric separately).\n');
    fprintf('Peak times (s from onset):\n');
    for k = 1:nPk
        rw = peakRows{k};
        fprintf('    %-14s %-5s  %+.3f s\n', rw{1}, rw{2}, rw{4});
    end
    fprintf('==========================================================================\n');

    % CSV
    csvPk = fullfile(pwd, sprintf('ach_peak_and_earlyAUC_evoked_vs_spont%s.csv', onsetTag));
    fidK  = fopen(csvPk, 'w');
    if fidK == -1
        [dK,nK,eK] = fileparts(csvPk);
        csvPk = fullfile(dK, sprintf('%s_%s%s', nK, datestr(now,'yyyymmdd_HHMMSS'), eK));
        fidK  = fopen(csvPk, 'w');
    end
    if fidK == -1
        warning('Could not write peak comparison CSV.');
    else
        fprintf(fidK, ['stimulus,region,n_mice,peak_time_s,' ...
                       'peak_mean_evoked,peak_sem_evoked,peak_mean_spont,peak_sem_spont,' ...
                       'peak_mean_diff,peak_p_raw,peak_p_FDR,peak_sig_FDR05,' ...
                       'earlyAUC_win_s,' ...
                       'earlyAUC_mean_evoked,earlyAUC_sem_evoked,' ...
                       'earlyAUC_mean_spont,earlyAUC_sem_spont,' ...
                       'earlyAUC_mean_diff,earlyAUC_p_raw,earlyAUC_p_FDR,earlyAUC_sig_FDR05\n']);
        for k = 1:nPk
            rw = peakRows{k};
            fprintf(fidK, ['%s,%s,%d,%.4f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%d,' ...
                           '%.4f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%d\n'], ...
                rw{1}, rw{2}, rw{3}, rw{4}, rw{5}, rw{6}, rw{7}, rw{8}, ...
                rw{9}, rw{10}, pAdjPk(k), pAdjPk(k)<0.05, ...
                earlyAUCWin, rw{13}, rw{14}, rw{15}, rw{16}, rw{17}, rw{18}, ...
                pAdjEa(k), pAdjEa(k)<0.05);
        end
        fclose(fidK);
        fprintf('Peak comparison saved to: %s\n', csvPk);
    end

    % --- Annotate peak-comparison significance on the trace figure ---
    % One short line per stimulus type in each panel, e.g. "Ripple **".
    if exist('axTrace','var') && exist('rowStimIdx','var')
        for rowI = 1:size(axTrace,1)
            sList = rowStimIdx{rowI};
            for r = 1:nR
                ax = axTrace(rowI,r);
                if ~isgraphics(ax); continue; end
                lbl = {};
                for k = 1:numel(sList)
                    s = sList(k);
                    if isempty(peakSigStr{s,r}); continue; end
                    lbl{end+1} = sprintf('%s %s', stimLabels{s}, peakSigStr{s,r}); %#ok
                end
                if isempty(lbl); continue; end
                text(ax, 0.03, 0.97, lbl, 'Units','normalized', ...
                    'FontSize', 7.5, 'Color', [0.25 0.25 0.25], ...
                    'VerticalAlignment','top', 'HorizontalAlignment','left');
            end
        end
    end
end

%% =========================================================================
%  Sanity check figure — collapsed across stimulus types
%  Layout: 1 row x nR columns.
%  Evoked: all stimulus types pooled into one mean±SEM trace per region.
%  Spontaneous: same single grey trace as above.
%  Pre-onset baseline correction applied (display only).
%% =========================================================================
figure('Units','normalized','Position',[0.05 0.20 0.80 0.45]);
sgtitle(['ACh Traces: Evoked (all stim pooled) vs Matched Spontaneous' titleSuffix], ...
    'FontSize', 13, 'FontWeight', 'bold');

for r = 1:nR
    ax = subplot(1, nR, r);
    hold(ax, 'on');

    legHandles = gobjects(0);
    legLabels  = {};

    % --- Spontaneous: per-animal means, averaged across ALL stim types ---
    % Per-mouse (achTrace_matched_mn), averaged within mouse rather than
    % concatenated, so n = mice and each animal is counted once.
    allSpont = [];
    for s = 1:nS
        Msp = achTrace_matched{s,r};        % TRIAL-LEVEL
        if isempty(Msp); continue; end
        Msp = Msp(all(isfinite(Msp),2), :);
        if isempty(Msp); continue; end
        blSp = mean(Msp(:, preOnsetMask), 2, 'omitnan');
        allSpont = [allSpont; Msp - blSp]; %#ok
    end
    if size(allSpont,1) >= 2
        muSp  = mean(allSpont, 1, 'omitnan');
        semSp = std(allSpont,  0, 1, 'omitnan') / sqrt(size(allSpont,1));
        fill(ax, [tResampPlot, fliplr(tResampPlot)], ...
            [muSp+semSp, fliplr(muSp-semSp)], ...
            cSpont, 'FaceAlpha', 0.20, 'EdgeColor', 'none');
        hSp = plot(ax, tResampPlot, muSp, '-', 'Color', cSpont, 'LineWidth', 1.8);
        legHandles(end+1) = hSp; %#ok
        legLabels{end+1}  = sprintf('Spontaneous (n=%d trials)', size(allSpont,1)); %#ok
    end

    % --- Evoked: per-animal means, averaged across ALL stim types ---
    allEvoked = [];
    for s = 1:nS
        Mev = achTrace_evoked{s,r};         % TRIAL-LEVEL
        if isempty(Mev); continue; end
        Mev = Mev(all(isfinite(Mev),2), :);
        if isempty(Mev); continue; end
        blEv = mean(Mev(:, preOnsetMask), 2, 'omitnan');
        allEvoked = [allEvoked; Mev - blEv]; %#ok
    end
    if size(allEvoked,1) >= 2
        muEv  = mean(allEvoked, 1, 'omitnan');
        semEv = std(allEvoked,  0, 1, 'omitnan') / sqrt(size(allEvoked,1));
        evC   = regionColors(r,:);
        fill(ax, [tResampPlot, fliplr(tResampPlot)], ...
            [muEv+semEv, fliplr(muEv-semEv)], ...
            evC, 'FaceAlpha', 0.20, 'EdgeColor', 'none');
        hEv = plot(ax, tResampPlot, muEv, '-', 'Color', evC, 'LineWidth', 2.2);
        legHandles(end+1) = hEv; %#ok
        legLabels{end+1}  = sprintf('Evoked (n=%d trials)', size(allEvoked,1)); %#ok
    end

    xline(ax, 0, '--k', 'Alpha', 0.4, 'LineWidth', 0.8);
    yline(ax, 0, '-', 'Color', [0.85 0.85 0.85], 'LineWidth', 0.7);
    xlim(ax, [tResampPlot(1) tResampPlot(end)]);
    xlabel(ax, 'Time re FME onset (s)', 'FontSize', 9);
    if r == 1; ylabel(ax, 'ACh (z)', 'FontSize', 9); end
    title(ax, regions{r}, 'FontSize', 11, 'FontWeight', 'bold', ...
        'Color', regionColors(r,:));
    if ~isempty(legHandles)
        legend(ax, legHandles, legLabels, 'FontSize', 8, 'Box', 'off', ...
            'Location', 'best');
    end
    box off;
    set(ax, 'TickDir', 'out', 'FontSize', 9);
end

%% =========================================================================
%  Difference figure: mean(evoked − spontaneous) grouped by stimulus type
%  Layout: 1 row x nS columns.
%  Each panel shows all brain regions overlaid in region colours.
%  Line style distinguishes regions (solid / dashed / dotted).
%  Pre-onset baseline correction applied (display only).
%  Traces are PER-ANIMAL differences (each animal's mean of its paired
%  evoked-minus-spontaneous differences), so the band is SEM across animals
%  and reflects exactly the quantity the paired statistics evaluate.
%% =========================================================================
figure('Units','normalized','Position',[0.05 0.30 0.80 0.45]);
sgtitle(['ACh Difference (Evoked \minus Spontaneous) by Stimulus Type' titleSuffix], ...
    'FontSize', 13, 'FontWeight', 'bold');

lineStylesR = {'-', '--', ':'};   % one per region

for s = 1:nS
    ax = subplot(1, nS, s);
    hold(ax, 'on');

    legHandles = gobjects(0);
    legLabels  = {};

    for r = 1:nR
        Md = achTrace_diff_mn{s,r};   % per-animal difference traces
        if isempty(Md); continue; end
        Md = Md(all(isfinite(Md),2),:);
        n  = size(Md,1);
        if n < 2; continue; end

        % Pre-onset baseline correction (display only)
        blMd = mean(Md(:, preOnsetMask), 2, 'omitnan');
        Md   = Md - blMd;

        muD  = mean(Md, 1, 'omitnan');
        semD = std(Md,  0, 1, 'omitnan') / sqrt(n);

        fill(ax, [tResampPlot, fliplr(tResampPlot)], ...
            [muD+semD, fliplr(muD-semD)], ...
            regionColors(r,:), 'FaceAlpha', 0.15, 'EdgeColor', 'none');
        h = plot(ax, tResampPlot, muD, lineStylesR{r}, ...
            'Color', regionColors(r,:), 'LineWidth', 2.0);

        legHandles(end+1) = h; %#ok
        legLabels{end+1}  = sprintf('%s (n=%d mice)', regions{r}, n); %#ok
    end

    % AUC window shading
    yL = ylim(ax);
    patch(ax, [0 boutAUCWinSec boutAUCWinSec 0], ...
        [yL(1) yL(1) yL(2) yL(2)], ...
        [0.9 0.9 0.9], 'FaceAlpha', 0.18, 'EdgeColor', 'none');

    xline(ax, 0, '--k', 'Alpha', 0.4, 'LineWidth', 0.8);
    yline(ax, 0, '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 1.0);

    xlim(ax, [tResampPlot(1) tResampPlot(end)]);
    xlabel(ax, 'Time re FME onset (s)', 'FontSize', 9);
    if s == 1; ylabel(ax, '\DeltaACh (z)', 'FontSize', 9); end
    title(ax, stimLabels{s}, 'FontSize', 11, 'FontWeight', 'bold');
    if ~isempty(legHandles)
        legend(ax, legHandles, legLabels, 'FontSize', 8, 'Box', 'off', ...
            'Location', 'best');
    end
    box off;
    set(ax, 'TickDir', 'out', 'FontSize', 9);
end

annotation('textbox', [0.01 0.00 0.99 0.04], ...
    'String', ['\DeltaACh = evoked \minus matched spontaneous, onset-aligned per pair, ' ...
               'averaged within animal (n = mice; band = SEM across animals).  ' ...
               'Pre-onset baseline correction applied for display.  ' ...
               'Grey shading = AUC window (0\–1 s).  ' ...
               'Zero = no stimulus-specific ACh beyond FME-matched baseline.' titleSuffix], ...
    'EdgeColor', 'none', 'FontSize', 8, 'Color', [0.4 0.4 0.4], ...
    'HorizontalAlignment', 'center');

%% =========================================================================
%  Old region-grouped difference figure (suppressed — code retained)
%% =========================================================================
if false %#ok<UNRCH>

figure('Units','normalized','Position',[0.05 0.30 0.80 0.38]);
sgtitle(['ACh Difference (Evoked \minus Spontaneous)' titleSuffix], ...
    'FontSize', 13, 'FontWeight', 'bold');

for r = 1:nR
    ax = subplot(1, nR, r);
    hold(ax, 'on');

    legHandles = gobjects(0);
    legLabels  = {};

    for s = 1:nS
        Md = achTrace_diff{s,r};
        if isempty(Md); continue; end
        Md = Md(all(isfinite(Md),2),:);
        n  = size(Md,1);
        if n < 2; continue; end

        muD  = mean(Md, 1, 'omitnan');
        semD = std(Md,  0, 1, 'omitnan') / sqrt(n);

        fill(ax, [tResampPlot, fliplr(tResampPlot)], ...
            [muD+semD, fliplr(muD-semD)], ...
            regionColors(r,:), 'FaceAlpha', 0.12, 'EdgeColor', 'none');
        h = plot(ax, tResampPlot, muD, lineStyles{s}, ...
            'Color', regionColors(r,:), 'LineWidth', 2.0);

        legHandles(end+1) = h; %#ok
        legLabels{end+1}  = sprintf('%s (n=%d)', stimLabels{s}, n); %#ok
    end

    yL = ylim(ax);
    patch(ax, [0 boutAUCWinSec boutAUCWinSec 0], ...
        [yL(1) yL(1) yL(2) yL(2)], ...
        [0.9 0.9 0.9], 'FaceAlpha', 0.18, 'EdgeColor', 'none');
    xline(ax, 0, '--k', 'Alpha', 0.4, 'LineWidth', 0.8);
    yline(ax, 0, '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 1.0);
    xlim(ax, [tResampPlot(1) tResampPlot(end)]);
    xlabel(ax, 'Time re FME onset (s)', 'FontSize', 9);
    if r == 1; ylabel(ax, '\DeltaACh (z)', 'FontSize', 9); end
    title(ax, regions{r}, 'FontSize', 11, 'FontWeight', 'bold', 'Color', regionColors(r,:));
    if ~isempty(legHandles)
        legend(ax, legHandles, legLabels, 'FontSize', 8, 'Box', 'off', 'Location', 'best');
    end
    box off; set(ax, 'TickDir', 'out', 'FontSize', 9);
end

end   % if false — old region-grouped difference figure

%% =========================================================================
%  Schematic figure: example continuous FME + AC ACh traces
%  (suppressed — code retained for reference)
%% =========================================================================
if false %#ok<UNRCH>
%  For each qualifying mouse (AC recorded, >=minBoutsForSchematic matched
%  spontaneous bouts), show:
%    Left  — ~30 s of spontaneous FME + AC ACh with shaded included bouts
%    Right — single representative AudRipple trial whose evoked FME AUC
%            matches a shown spontaneous bout (FME on top, ACh below)
%  Multiple mice are plotted; user picks the best one.
%% =========================================================================
minBoutsForSchematic = 5;    % minimum matched spont bouts to qualify
winSec_schema        = 30;   % continuous trace window length (s)
acRegionIdx          = find(strcmp(regions, 'AC'));
stimTypeEvoked       = 'AudRipple';
stimIdxEvoked        = find(strcmp(stimTypes, stimTypeEvoked));
shadeFaceAlpha       = 0.25;
shadeEdgeColor       = 'none';
shadeColor           = [0.85 0.85 0.85];  % grey for included bouts
acColor              = regionColors(acRegionIdx, :);
fmeColor             = [0.3 0.3 0.3];

nSchema = 0;   % count of qualifying mice plotted

for m = 1:nMice
    if isempty(ach_stim{m}) || isempty(fme_stim{m}); continue; end
    if m > nSessions; continue; end

    % Check AC is recorded for this mouse
    tmp2    = load(fullfile(sigFiles(m).folder, sigFiles(m).name), 'traces');
    tr2     = tmp2.traces;
    chForR2 = tr2.chForRegion;
    if chForR2(acRegionIdx) == 0; continue; end

    % Check enough matched spont bouts exist for AC in AudRipple
    if numel(achAUC_evoked_mn{stimIdxEvoked, acRegionIdx}) < 1; continue; end

    % -----------------------------------------------------------------------
    %  Re-detect spontaneous bouts for this mouse (same pipeline as above)
    % -----------------------------------------------------------------------
    fmeCont2    = tr2.fme(:)';
    achCh2      = {tr2.achCh1, tr2.achCh2};
    fsAch2      = tr2.fsAch;
    fsFme2      = fsAch2;
    nContFme2   = numel(fmeCont2);
    tContFme2   = (0:nContFme2-1) / fsFme2;
    tContAch2   = (0:numel(achCh2{1})-1) / fsAch2;

    medFiltSamp2 = round(boutMedFilt * fsFme2);
    fmeFilt2     = medfilt1(fmeCont2, medFiltSamp2);
    fmeFilt2(isnan(fmeCont2)) = NaN;
    fmeGMu2 = mean(fmeFilt2(~isnan(fmeFilt2)), 'omitnan');
    fmeGSd2 = std( fmeFilt2(~isnan(fmeFilt2)), 'omitnan');
    if fmeGSd2 == 0; fmeGSd2 = 1; end
    fmeFiltZ2 = (fmeFilt2 - fmeGMu2) ./ fmeGSd2;

    validS2   = fmeFiltZ2(~isnan(fmeFiltZ2));
    thresh2   = prctile(validS2, boutThreshPct) * boutThreshMult;

    above2 = fmeFiltZ2 > thresh2;
    above2(isnan(fmeFiltZ2)) = false;
    above2 = above2(:)';
    minGapSamp2 = round(minGapSec * fsFme2);
    d2 = diff([1, above2, 1]);
    falls2 = find(d2 == -1); rises2 = find(d2 == 1);
    for gi = 1:numel(falls2)
        if gi <= numel(rises2) && (rises2(gi)-falls2(gi)) < minGapSamp2
            above2(falls2(gi):rises2(gi)-1) = true;
        end
    end
    d2    = diff([0, above2, 0]);
    bOns2 = find(d2 ==  1);
    bOffs2= find(d2 == -1) - 1;

    % Filter: duration <=1s, AUC > fmeThreshold, AC ACh finite
    aucWinSamp2 = round(boutAUCWinSec * fsFme2);
    validBouts  = false(size(bOns2));
    boutFmeAUC2 = nan(size(bOns2));
    for ib = 1:numel(bOns2)
        dur = (bOffs2(ib) - bOns2(ib)) / fsFme2;
        if dur > boutAUCWinSec; continue; end
        iOn2  = bOns2(ib);
        iOff2 = min(iOn2 + aucWinSamp2, nContFme2);
        auc2  = trapz(tContFme2(iOn2:iOff2), fmeFiltZ2(iOn2:iOff2));
        if auc2 <= fmeThreshold; continue; end
        % Check AC ACh baseline available
        tOns2   = tContFme2(iOn2);
        blOn2   = find(tContAch2 >= tOns2 - 2, 1, 'first');
        blOff2  = find(tContAch2 <= tOns2,     1, 'last');
        if isempty(blOn2) || isempty(blOff2) || blOff2 <= blOn2; continue; end
        ch2 = chForR2(acRegionIdx);
        if ch2 == 0; continue; end
        validBouts(ib)  = true;
        boutFmeAUC2(ib) = auc2;
    end

    nValid = sum(validBouts);
    if nValid < minBoutsForSchematic; continue; end

    validIdx = find(validBouts);

    % -----------------------------------------------------------------------
    %  Find a ~30 s window containing >=2 valid bouts
    % -----------------------------------------------------------------------
    winSampFme = round(winSec_schema * fsFme2);
    foundWin   = false;
    winStart   = 1;
    for ib = 1:numel(validIdx)-1
        i1 = validIdx(ib); i2 = validIdx(ib+1);
        % window centred loosely around the pair
        wS = max(1, bOns2(i1) - round(5 * fsFme2));
        wE = wS + winSampFme;
        if wE > nContFme2; continue; end
        % count valid bouts inside window
        boutsInWin = validIdx(tContFme2(bOns2(validIdx)) >= tContFme2(wS) & ...
                               tContFme2(bOns2(validIdx)) <= tContFme2(wE));
        if numel(boutsInWin) >= 2
            winStart  = wS;
            foundWin  = true;
            boutsShown = boutsInWin;
            break;
        end
    end
    if ~foundWin; continue; end

    winEnd   = winStart + winSampFme;
    tWin     = tContFme2(winStart:winEnd) - tContFme2(winStart);
    fmeWin   = fmeFiltZ2(winStart:winEnd);
    ch2      = chForR2(acRegionIdx);
    achRaw2  = achCh2{ch2};

    % Z-score AC ACh for display using first bout's baseline
    tOns_first = tContFme2(bOns2(boutsShown(1)));
    blOn_d  = find(tContAch2 >= tOns_first - 2, 1, 'first');
    blOff_d = find(tContAch2 <= tOns_first,     1, 'last');
    blMu_d  = mean(achRaw2(blOn_d:blOff_d), 'omitnan');
    blSd_d  = std( achRaw2(blOn_d:blOff_d), 'omitnan');
    if blSd_d == 0; blSd_d = 1; end
    achContZ = (achRaw2 - blMu_d) / blSd_d;

    achWinIdx = find(tContAch2 >= tContFme2(winStart) & ...
                     tContAch2 <= tContFme2(winEnd));
    tAchWin   = tContAch2(achWinIdx) - tContFme2(winStart);
    achWin    = achContZ(achWinIdx);

    % -----------------------------------------------------------------------
    %  Select evoked AudRipple trial: moderate intensity, FME AUC closest
    %  to the median of shown spontaneous bout AUC values
    % -----------------------------------------------------------------------
    medSpAUC = median(boutFmeAUC2(boutsShown), 'omitnan');

    % Get evoked data for this mouse
    if ~isfield(fme_stim{m}, stimTypeEvoked) || ~isfield(ach_stim{m}, stimTypeEvoked)
        continue;
    end
    fmeEvRaw = fme_stim{m}.(stimTypeEvoked).traces;
    if iscell(fmeEvRaw); fmeEvRaw = fmeEvRaw{1}; end
    intVecEv = fme_stim{m}.(stimTypeEvoked).(intFields{stimIdxEvoked});
    intVecEv = intVecEv(:);
    uInts    = unique(intVecEv);
    [~,ni]   = min(abs(uInts - targetInt(stimTypeEvoked)));
    modMask  = intVecEv == uInts(ni);
    fmeEvMod = fmeEvRaw(modMask, :);

    % Session-z-score and median-filter evoked FME (same as main pipeline)
    fsFmeEv2     = nNative / 8;
    mfSamp_ev2   = max(3, round(boutMedFilt * fsFmeEv2));
    fmeEvFilt2   = fmeEvMod;
    for ii = 1:size(fmeEvMod,1)
        row = medfilt1(fmeEvMod(ii,:), mfSamp_ev2);
        row(isnan(fmeEvMod(ii,:))) = NaN;
        fmeEvFilt2(ii,:) = (row - fmeGMu2) / fmeGSd2;
    end

    % Detect onset and compute AUC per evoked trial
    postMask2    = tNative >= 0;
    aucWinSampEv2= round(boutAUCWinSec * fsFmeEv2);
    evFmeAUC2    = nan(size(fmeEvFilt2,1),1);
    evOnsetIdx2  = nan(size(fmeEvFilt2,1),1);
    for ii = 1:size(fmeEvFilt2,1)
        row   = fmeEvFilt2(ii,:);
        pIdx  = find(postMask2 & row > thresh2, 1, 'first');
        if isempty(pIdx); continue; end
        aEnd  = min(pIdx + aucWinSampEv2, numel(tNative));
        evFmeAUC2(ii)   = trapz(tNative(pIdx:aEnd), row(pIdx:aEnd));
        evOnsetIdx2(ii) = pIdx;
    end

    % Pick trial with FME AUC closest to median spontaneous
    evOk2 = ~isnan(evFmeAUC2) & evFmeAUC2 > fmeThreshold;
    if sum(evOk2) == 0; continue; end
    okIdx    = find(evOk2);
    [~, bi]  = min(abs(evFmeAUC2(okIdx) - medSpAUC));
    trialIdx = okIdx(bi);

    % Evoked FME trace (session z-scored, onset-aligned)
    tOnsetEvSec = tNative(evOnsetIdx2(trialIdx));
    tFmeEv2     = tNative - tOnsetEvSec;
    fmeEvTr     = fmeEvFilt2(trialIdx,:);

    % Evoked ACh trace (per-trial z-score, onset-aligned)
    fiberIdx2 = find(strcmp(all_fiber_names_ach{m}, 'AC'));
    if isempty(fiberIdx2); continue; end
    achEvFull = ach_stim{m}.(stimTypeEvoked).traces{fiberIdx2};
    achEvMod  = achEvFull(modMask, :);
    tAchEv2   = linspace(-2, 6, size(achEvMod,2));
    achEvTr   = achEvMod(trialIdx,:);
    blMaskEv2 = tAchEv2 >= -2 & tAchEv2 < 0;
    blMuEv2   = mean(achEvTr(blMaskEv2), 'omitnan');
    blSdEv2   = std( achEvTr(blMaskEv2), 'omitnan');
    if blSdEv2 == 0; blSdEv2 = 1; end
    achEvTrZ  = (achEvTr - blMuEv2) / blSdEv2;
    tAchEv2_sh = tAchEv2 - tOnsetEvSec;   % onset-aligned

    % -----------------------------------------------------------------------
    %  Plot schematic figure for this mouse
    % -----------------------------------------------------------------------
    nSchema = nSchema + 1;
    figure('Units','normalized','Position',[0.03 0.05 0.90 0.55], ...
        'Name', sprintf('Schematic Mouse %d', m));
    sgtitle(sprintf('Example Traces — Mouse %d%s', m, titleSuffix), ...
        'FontSize', 12, 'FontWeight', 'bold');

    % Column 1: spontaneous | Column 2: evoked
    % Row 1: FME | Row 2: ACh
    axFmeSp  = subplot(2, 2, 1);
    axFmeEv  = subplot(2, 2, 2);
    axAchSp  = subplot(2, 2, 3);
    axAchEv  = subplot(2, 2, 4);

    % --- Spontaneous FME ---
    hold(axFmeSp, 'on');
    plot(axFmeSp, tWin, fmeWin, '-', 'Color', fmeColor, 'LineWidth', 0.8);
    for ib = boutsShown
        tOn_b  = tContFme2(bOns2(ib))  - tContFme2(winStart);
        tOff_b = tContFme2(bOffs2(ib)) - tContFme2(winStart);
        % shade full detected bout duration
        patch(axFmeSp, [tOn_b tOff_b tOff_b tOn_b], ...
            [min(fmeWin) min(fmeWin) max(fmeWin) max(fmeWin)], ...
            shadeColor, 'FaceAlpha', shadeFaceAlpha, 'EdgeColor', shadeEdgeColor);
    end
    yline(axFmeSp, thresh2, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.8);
    ylabel(axFmeSp, 'FME (z)', 'FontSize', 9);
    title(axFmeSp, 'Spontaneous', 'FontSize', 10, 'FontWeight', 'bold');
    xlim(axFmeSp, [tWin(1) tWin(end)]);
    box off; set(axFmeSp, 'TickDir','out','FontSize',9,'XTickLabel',[]);

    % --- Spontaneous ACh ---
    hold(axAchSp, 'on');
    plot(axAchSp, tAchWin, achWin, '-', 'Color', acColor, 'LineWidth', 0.9);
    for ib = boutsShown
        tOn_b  = tContFme2(bOns2(ib))  - tContFme2(winStart);
        tOff_b = tContFme2(bOffs2(ib)) - tContFme2(winStart);
        yL = ylim(axAchSp);
        patch(axAchSp, [tOn_b tOff_b tOff_b tOn_b], ...
            [yL(1) yL(1) yL(2) yL(2)], ...
            shadeColor, 'FaceAlpha', shadeFaceAlpha, 'EdgeColor', shadeEdgeColor);
    end
    yline(axAchSp, 0, '-', 'Color', [0.85 0.85 0.85], 'LineWidth', 0.7);
    xlabel(axAchSp, 'Time (s)', 'FontSize', 9);
    ylabel(axAchSp, 'AC ACh (z)', 'FontSize', 9);
    xlim(axAchSp, [tWin(1) tWin(end)]);
    box off; set(axAchSp, 'TickDir','out','FontSize',9);

    % --- Evoked FME (onset-aligned) ---
    hold(axFmeEv, 'on');
    % shade the detected evoked bout duration
    evBoutDurSec = min((bOffs2(boutsShown(1)) - bOns2(boutsShown(1))) / fsFme2, boutAUCWinSec);
    patch(axFmeEv, [0 evBoutDurSec evBoutDurSec 0], ...
        [min(fmeEvTr) min(fmeEvTr) max(fmeEvTr) max(fmeEvTr)], ...
        shadeColor, 'FaceAlpha', shadeFaceAlpha, 'EdgeColor', shadeEdgeColor);
    plot(axFmeEv, tFmeEv2, fmeEvTr, '-', 'Color', fmeColor, 'LineWidth', 0.8);
    xline(axFmeEv, 0, '--k', 'Alpha', 0.4, 'LineWidth', 0.8);
    yline(axFmeEv, thresh2, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.8);
    title(axFmeEv, sprintf('Evoked (AudRipple, %d dB)', uInts(ni)), ...
        'FontSize', 10, 'FontWeight', 'bold');
    xlim(axFmeEv, [-2 4]);
    box off; set(axFmeEv, 'TickDir','out','FontSize',9,'XTickLabel',[]);

    % --- Evoked ACh (onset-aligned) ---
    hold(axAchEv, 'on');
    patch(axAchEv, [0 evBoutDurSec evBoutDurSec 0], ...
        [-10 -10 10 10], ...   % will be clipped by ylim
        shadeColor, 'FaceAlpha', shadeFaceAlpha, 'EdgeColor', shadeEdgeColor);
    plot(axAchEv, tAchEv2_sh, achEvTrZ, '-', 'Color', acColor, 'LineWidth', 1.2);
    xline(axAchEv, 0, '--k', 'Alpha', 0.4, 'LineWidth', 0.8);
    yline(axAchEv, 0, '-', 'Color', [0.85 0.85 0.85], 'LineWidth', 0.7);
    xlabel(axAchEv, 'Time re FME onset (s)', 'FontSize', 9);
    xlim(axAchEv, [-2 4]);
    box off; set(axAchEv, 'TickDir','out','FontSize',9);

    % Link y-axes within columns so FME panels share scale, ACh panels share scale
    linkaxes([axFmeSp, axFmeEv], 'y');
    linkaxes([axAchSp, axAchEv], 'y');
    % Link x-axes within rows so spontaneous and evoked FME have same time scale
    % (spontaneous x = absolute time within window, evoked x = onset-relative)
    % — not linked since scales differ intentionally

    fprintf('  Schematic mouse %d: %d valid bouts, window %.1f–%.1f s, trial %d (evFME AUC=%.2f, medSpAUC=%.2f)\n', ...
        m, nValid, tContFme2(winStart), tContFme2(winEnd), trialIdx, evFmeAUC2(trialIdx), medSpAUC);
end

if nSchema == 0
    fprintf('  No qualifying mice found for schematic figure (need AC recording + >=%d matched bouts).\n', ...
        minBoutsForSchematic);
end

end   % if false — schematic figure

%% =========================================================================
end   % main function

%  Local subfunctions
%% =========================================================================

function aucVec = compute_trial_auc(traces, tOrig, blMask, aucMask)
mat   = traces;
bl    = mean(mat(:, blMask), 2, 'omitnan');
mat   = mat - bl;
blSD  = std(mat(:, blMask), 0, 2, 'omitnan');
blSD(blSD == 0) = 1;
mat   = mat ./ blSD;
tSeg  = tOrig(aucMask);
aucVec = trapz(tSeg, mat(:, aucMask), 2);

end

function aucVec = compute_ach_trial_auc_onset(traces, onsetIdxVec, ~, winSec)
% Compute per-trial ACh AUC from FME-onset to onset+winSec.
% onsetIdxVec: [nTrials x 1] index into tNative (241-pt) for each trial's
%              FME onset — converted here to the ACh time vector (8139-pt).
% Baseline: -2 to 0 s re stimulus (pre-stimulus, independent of onset).
tOrig   = linspace(-2, 6, size(traces, 2));   % ACh time vector at actual rate
blMask  = tOrig >= -2 & tOrig < 0;
tNative = linspace(-2, 6, 241);               % FME time vector for onset lookup
nTrials = size(traces, 1);
aucVec  = nan(nTrials, 1);
for ii = 1:nTrials
    if isnan(onsetIdxVec(ii)); continue; end
    % Convert FME onset time to ACh sample index
    tOnsetSec = tNative(onsetIdxVec(ii));
    iOn  = find(tOrig >= tOnsetSec, 1, 'first');
    iOff = find(tOrig <= tOnsetSec + winSec, 1, 'last');
    if isempty(iOn) || isempty(iOff) || iOff <= iOn; continue; end
    row  = traces(ii,:);
    blMu = mean(row(blMask), 'omitnan');
    blSd = std( row(blMask), 'omitnan');
    if blSd == 0; blSd = 1; end
    rowZ = (row - blMu) / blSd;
    aucVec(ii) = trapz(tOrig(iOn:iOff), rowZ(iOn:iOff));
end
end


function aucVec = compute_ach_trial_auc(traces)
tOrig   = linspace(-2, 6, size(traces,2));
blMask  = tOrig >= -2 & tOrig < 0;
aucMask = tOrig >= 0  & tOrig <= 1;
mat     = traces;
bl      = mean(mat(:, blMask), 2, 'omitnan');
mat     = mat - bl;
blSD    = std(mat(:, blMask), 0, 2, 'omitnan');
blSD(blSD == 0) = 1;
mat     = mat ./ blSD;
tSeg    = tOrig(aucMask);
aucVec  = trapz(tSeg, mat(:, aucMask), 2);
end
