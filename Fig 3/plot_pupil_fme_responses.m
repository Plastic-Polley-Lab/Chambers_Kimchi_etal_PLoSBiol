function plot_pupil_fme_responses(fme_stim, pupil_stim, fme_silent, pupil_silent, stimTypes_in)
% PLOT_PUPIL_FME_RESPONSES  Plot FME and pupil responses by stimulus type.
%
% Usage:
%   plot_pupil_fme_responses(fme_stim, pupil_stim)
%   plot_pupil_fme_responses(fme_stim, pupil_stim, {'Tones','VisGrating'})
%   plot_pupil_fme_responses(fme_stim, pupil_stim, {'Tones'})
%
% fme_silent, pupil_silent: separate cell arrays with .Silent.traces for noise floor.
% stimTypes_in: optional cell array of stimulus type names.
%   Default: {'AudRipple','VisGrating'}
%
% Noise floor: computed from Silent trials (fme_stim{m}.Silent / pupil_stim{m}.Silent)
%   using the same 0-1 s AUC window and pre-stimulus baseline subtraction.
%   Shown as mean ± SEM horizontal shaded band across the AUC panel.
%
% Figure 1: FME  — 1 row x (nS+1) columns: one trace per stim + combined AUC
% Figure 2: Pupil — same layout, independent y-axis per trace panel
% Figure 3: FME  diagnostic (per-mouse means, highest intensity, no processing)
% Figure 4: Pupil diagnostic (same)

nMice       = numel(fme_stim);
nNative     = 241;
tNative     = linspace(-2, 6, nNative);
baselineIdx = tNative >= -2 & tNative < 0;
aucIdx      = tNative >= 0  & tNative <= 1;

% Stimulus type configuration
if nargin < 5 || isempty(stimTypes_in)
    % All three stimulus types by default. Previously {'AudRipple','VisGrating'}
    % only, so tones were silently absent from the figure and the stats table
    % while every other function in the pipeline included them. Order matches
    % the acquisition table so rows line up across spreadsheets.
    stimTypes = {'AudRipple', 'Tones', 'VisGrating'};
else
    stimTypes = stimTypes_in;
end
nS = numel(stimTypes);

intFieldMap = containers.Map({'AudRipple','Tones','VisGrating'}, ...
                              {'int',      'int',  'contrast'});
labelMap    = containers.Map({'AudRipple','Tones','VisGrating'}, ...
                              {'Aud. Ripple','Tones','Vis. Grating'});
xLabelMap   = containers.Map({'AudRipple','Tones','VisGrating'}, ...
                              {'Intensity (dB)','Intensity (dB)','Contrast'});

intFields  = cellfun(@(s) intFieldMap(s), stimTypes, 'UniformOutput', false);
stimLabels = cellfun(@(s) labelMap(s),    stimTypes, 'UniformOutput', false);

signalColors = [
    0.15, 0.15, 0.15;   % near-black -> FME
    0.15, 0.15, 0.15;   % near-black -> Pupil
];
signalLabels = {'FME', 'Pupil'};
nSig = 2;

grayLo = 0.72;   % lightest shade (lowest intensity)
grayHi = 0.20;   % darkest shade  (highest intensity)

%% =========================================================================
%  Noise floor from Silent trials
%  For each signal (FME, pupil), collect per-mouse AUC values from Silent
%  trials using the same baseline subtraction and 0-1 s AUC window.
%  Summary: mean ± SEM across mice.
%% =========================================================================
noiseFloorAUC     = cell(nSig, 1);   % {sigIdx} = [nMice x 1] per-mouse silent AUC mean
noiseFloorMouseIdx= cell(nSig, 1);   % {sigIdx} = [nMice x 1] mouse index
rawSilent         = {fme_silent, pupil_silent};

for sigIdx = 1:nSig
    src = rawSilent{sigIdx};
    mouseNF    = [];
    mouseNFIdx = [];
    for m = 1:nMice
        if isempty(src{m}); continue; end
        if ~isfield(src{m}, 'Silent'); continue; end
        silTraces = src{m}.Silent.traces;
        if iscell(silTraces); silTraces = silTraces{1}; end
        if isempty(silTraces); continue; end

        % Baseline subtraction consistent with stimulus trials
        bl       = mean(silTraces(:, baselineIdx), 2, 'omitnan');
        silProc  = silTraces - bl;

        % AUC over 0-1 s window, mean across silent trials for this mouse
        aucPerTrial = trapz(tNative(aucIdx), silProc(:, aucIdx), 2);
        mouseNF(end+1)    = mean(aucPerTrial, 'omitnan'); %#ok
        mouseNFIdx(end+1) = m; %#ok
    end
    noiseFloorAUC{sigIdx}      = mouseNF(:);
    noiseFloorMouseIdx{sigIdx} = mouseNFIdx(:);
end

% Summarise noise floor: mean ± SEM across mice
nfMu  = cellfun(@(v) mean(v, 'omitnan'),                         noiseFloorAUC);
nfSem = cellfun(@(v) std(v, 0, 'omitnan') / max(sqrt(numel(v)),1), noiseFloorAUC);

fprintf('Noise floor (FME):   mean=%.4f  SEM=%.4f  (n=%d mice)\n', ...
    nfMu(1), nfSem(1), numel(noiseFloorAUC{1}));
fprintf('Noise floor (Pupil): mean=%.4f  SEM=%.4f  (n=%d mice)\n', ...
    nfMu(2), nfSem(2), numel(noiseFloorAUC{2}));

%% =========================================================================
%  Data collection
%% =========================================================================
allInts      = cell(nS, 1);
shadeFactors = cell(nS, 1);
mouseMeans   = cell(nS, nSig);
aucVals      = cell(nS, nSig);
aucMouseIdx  = cell(nS, nSig);   % mouse index corresponding to each aucVals entry

for si = 1:nS
    ivals = [];
    for m = 1:nMice
        if isempty(fme_stim{m}); continue; end
        ivals = union(ivals, unique(fme_stim{m}.(stimTypes{si}).(intFields{si})));
    end
    allInts{si}      = sort(ivals);
    nLev             = numel(allInts{si});
    shadeFactors{si} = linspace(0.65, 0, nLev);

    for sigIdx = 1:nSig
        mouseMeans{si,sigIdx} = cell(nLev, 1);
        aucVals{si,sigIdx}    = cell(nLev, 1);
        aucMouseIdx{si,sigIdx}= cell(nLev, 1);
        for iInt = 1:nLev
            mouseMeans{si,sigIdx}{iInt}  = [];
            aucVals{si,sigIdx}{iInt}     = [];
            aucMouseIdx{si,sigIdx}{iInt} = [];
        end
    end

    for m = 1:nMice
        if isempty(fme_stim{m}) || isempty(pupil_stim{m}); continue; end

        fmeTraces = fme_stim{m}.(stimTypes{si}).traces;
        pupTraces = pupil_stim{m}.(stimTypes{si}).traces;
        ints      = fme_stim{m}.(stimTypes{si}).(intFields{si});

        if iscell(fmeTraces); fmeTraces = fmeTraces{1}; end
        if iscell(pupTraces); pupTraces = pupTraces{1}; end

        rawData = {fmeTraces, pupTraces};

        for iInt = 1:nLev
            trialIdx = ints == allInts{si}(iInt);
            if sum(trialIdx) == 0; continue; end

            for sigIdx = 1:nSig
                mat      = rawData{sigIdx}(trialIdx, :);
                baseline = mean(mat(:, baselineIdx), 2, 'omitnan');
                mat      = mat - baseline;
                muTrial  = mean(mat, 1, 'omitnan');
                mouseMeans{si,sigIdx}{iInt} = [mouseMeans{si,sigIdx}{iInt}; muTrial];
                auc = trapz(tNative(aucIdx), muTrial(aucIdx));
                aucVals{si,sigIdx}{iInt}     = [aucVals{si,sigIdx}{iInt};     auc];
                aucMouseIdx{si,sigIdx}{iInt} = [aucMouseIdx{si,sigIdx}{iInt}; m];
            end
        end
    end
end

%% =========================================================================
%  Ripple reference (only when called with {'Tones'} alone)
%% =========================================================================
showRippleRef = isequal(stimTypes, {'Tones'});
rippleMeans   = cell(nSig, 1);
rippleAUC     = cell(nSig, 1);
rippleInts    = [];

if showRippleRef
    for m = 1:nMice
        if isempty(fme_stim{m}); continue; end
        if ~isfield(fme_stim{m}, 'AudRipple'); continue; end
        rippleInts = union(rippleInts, unique(fme_stim{m}.AudRipple.int));
    end
    rippleInts = sort(rippleInts);
    nLevRip    = numel(rippleInts);

    for sigIdx = 1:nSig
        rippleMeans{sigIdx} = cell(nLevRip, 1);
        rippleAUC{sigIdx}   = cell(nLevRip, 1);
        for iInt = 1:nLevRip
            rippleMeans{sigIdx}{iInt} = [];
            rippleAUC{sigIdx}{iInt}   = [];
        end
    end

    for m = 1:nMice
        if isempty(fme_stim{m}) || isempty(pupil_stim{m}); continue; end
        if ~isfield(fme_stim{m}, 'AudRipple'); continue; end

        fmeRip = fme_stim{m}.AudRipple.traces;
        pupRip = pupil_stim{m}.AudRipple.traces;
        if iscell(fmeRip); fmeRip = fmeRip{1}; end
        if iscell(pupRip); pupRip = pupRip{1}; end
        intsRip = fme_stim{m}.AudRipple.int;
        rawRip  = {fmeRip, pupRip};

        for iInt = 1:nLevRip
            trialIdx = intsRip == rippleInts(iInt);
            if sum(trialIdx) == 0; continue; end
            for sigIdx = 1:nSig
                mat     = rawRip{sigIdx}(trialIdx, :);
                bl      = mean(mat(:, baselineIdx), 2, 'omitnan');
                mat     = mat - bl;
                muTrial = mean(mat, 1, 'omitnan');
                rippleMeans{sigIdx}{iInt} = [rippleMeans{sigIdx}{iInt}; muTrial];
                auc = trapz(tNative(aucIdx), muTrial(aucIdx));
                rippleAUC{sigIdx}{iInt} = [rippleAUC{sigIdx}{iInt}; auc];
            end
        end
    end
end

%% =========================================================================
%  Y-limits
%% =========================================================================
fme_upper = []; fme_lower = [];
if showRippleRef && ~isempty(rippleInts)
    for iInt = 1:numel(rippleInts)
        if isempty(rippleMeans{1}{iInt}); continue; end
        mu  = mean(rippleMeans{1}{iInt}, 1, 'omitnan');
        sem = std( rippleMeans{1}{iInt}, 0, 1, 'omitnan') ./ ...
              sqrt(size(rippleMeans{1}{iInt}, 1));
        fme_upper = [fme_upper; mu + sem];
        fme_lower = [fme_lower; mu - sem];
    end
else
    nLev_aud = numel(allInts{1});
    for iInt = 1:nLev_aud
        if isempty(mouseMeans{1,1}{iInt}); continue; end
        mu  = mean(mouseMeans{1,1}{iInt}, 1, 'omitnan');
        sem = std( mouseMeans{1,1}{iInt}, 0, 1, 'omitnan') ./ ...
              sqrt(size(mouseMeans{1,1}{iInt}, 1));
        fme_upper = [fme_upper; mu + sem];
        fme_lower = [fme_lower; mu - sem];
    end
end
fme_pad   = 0.15 * (max(fme_upper(:)) - min(fme_lower(:)));
yLims_fme = [min(fme_lower(:)) - fme_pad, max(fme_upper(:)) + fme_pad];

pupil_yLims = cell(nS, 1);
for si = 1:nS
    p_upper = []; p_lower = [];
    if showRippleRef && ~isempty(rippleInts)
        for iInt = 1:numel(rippleInts)
            if isempty(rippleMeans{2}{iInt}); continue; end
            mu  = mean(rippleMeans{2}{iInt}, 1, 'omitnan');
            sem = std( rippleMeans{2}{iInt}, 0, 1, 'omitnan') ./ ...
                  sqrt(size(rippleMeans{2}{iInt}, 1));
            p_upper = [p_upper; mu + sem];
            p_lower = [p_lower; mu - sem];
        end
    end
    nLev = numel(allInts{si});
    for iInt = 1:nLev
        if isempty(mouseMeans{si,2}{iInt}); continue; end
        mu  = mean(mouseMeans{si,2}{iInt}, 1, 'omitnan');
        sem = std( mouseMeans{si,2}{iInt}, 0, 1, 'omitnan') ./ ...
              sqrt(size(mouseMeans{si,2}{iInt}, 1));
        p_upper = [p_upper; mu + sem];
        p_lower = [p_lower; mu - sem];
    end
    if isempty(p_upper)
        pupil_yLims{si} = [-1 1];
    else
        p_pad = 0.15 * (max(p_upper(:)) - min(p_lower(:)));
        pupil_yLims{si} = [min(p_lower(:)) - p_pad, max(p_upper(:)) + p_pad];
    end
end

%% =========================================================================
%  Helper: draw noise floor shading on AUC panel
%% =========================================================================
    function draw_noise_floor(ax, xLo, xHi, mu, sem, col)
        patch(ax, [xLo xHi xHi xLo], ...
            [mu-sem, mu-sem, mu+sem, mu+sem], ...
            col, 'FaceAlpha', 0.18, 'EdgeColor', 'none');
        h = plot(ax, [xLo xHi], [mu mu], '--', ...
            'Color', col, 'LineWidth', 1.0);
        h.Color(4) = 0.55;
    end

%% =========================================================================
%  Helper: paired Wilcoxon test of stimulus AUC vs noise floor AUC
%  Matches mice present in both vectors by mouse index, then tests the
%  paired differences. Returns p-value and star string.
%% =========================================================================
    function [pVal, stars, pTwo, mDiff] = noise_floor_test(stimAUC, stimMIdx, nfAUC, nfMIdx)
        % Find mice present in both
        commonMice = intersect(stimMIdx, nfMIdx);
        if numel(commonMice) < 3
            pVal  = NaN;
            stars = '';
            pTwo  = NaN;
            mDiff = NaN;
            return;
        end
        stimVec = arrayfun(@(m) stimAUC(stimMIdx == m), commonMice);
        nfVec   = arrayfun(@(m) nfAUC(nfMIdx   == m), commonMice);

        % One-sided (stimulus > noise floor) drives the asterisks, since the
        % hypothesis for a sensory response is directional.
        % TWO-SIDED is also returned: the one-sided test cannot detect a
        % response BELOW the noise floor, and at least one condition is
        % expected to go that way — visual gratings constrict the pupil via
        % the light reflex. Both are written to the CSV.
        pVal  = signrank(stimVec(:), nfVec(:), 'tail', 'right');
        pTwo  = signrank(stimVec(:), nfVec(:));
        mDiff = mean(stimVec(:) - nfVec(:), 'omitnan');
        if pVal < 0.001;    stars = '***';
        elseif pVal < 0.01; stars = '**';
        elseif pVal < 0.05; stars = '*';
        else;                stars = '';
        end
    end

nfColor    = [0.55 0.55 0.55];   % neutral grey for noise floor
statsTable = {};   % collect p-values for CSV export

%% =========================================================================
%  Figure 1: FME — 1 x (nS+1)
%% =========================================================================
sig = 1;
c   = signalColors(sig,:);
figure('Units','normalized','Position',[0.05 0.20 0.88 0.42]);

for si = 1:nS
    nLev  = numel(allInts{si});
    ivals = allInts{si};
    ax    = subplot(1, nS+1, si);
    hold(ax, 'on');

    if showRippleRef && ~isempty(rippleInts)
        ripAllMat = [];
        for iRip = 1:numel(rippleInts)
            if ~isempty(rippleMeans{sig}{iRip})
                ripAllMat = [ripAllMat; rippleMeans{sig}{iRip}]; %#ok
            end
        end
        if ~isempty(ripAllMat)
            muR  = mean(ripAllMat, 1, 'omitnan');
            semR = std( ripAllMat, 0, 1, 'omitnan') ./ sqrt(size(ripAllMat, 1));
            shadedErrorBar(tNative, muR, semR, ...
                'lineProps', {'Color', [0.85 0.78 0.65], 'LineWidth', 1.5}, ...
                'transparent', 1, 'patchSaturation', 0.25);
        end
    end

    for iInt = 1:nLev
        if isempty(mouseMeans{si,sig}{iInt}); continue; end
        mu  = mean(mouseMeans{si,sig}{iInt}, 1, 'omitnan');
        sem = std( mouseMeans{si,sig}{iInt}, 0, 1, 'omitnan') ./ ...
              sqrt(size(mouseMeans{si,sig}{iInt}, 1));
        gv     = linspace(grayLo, grayHi, nLev);
        cShade = repmat(gv(iInt), 1, 3);
        shadedErrorBar(tNative, mu, sem, ...
            'lineProps', {'Color', cShade, 'LineWidth', 1.8}, ...
            'transparent', 1, 'patchSaturation', 0.15);
    end

    xline(0, '--k', 'LineWidth', 1.0, 'Alpha', 0.5);
    xlabel('Time from onset (s)', 'FontSize', 11);
    if si == 1; ylabel('baseline sub. (a.u.)', 'FontSize', 11); end
    title(stimLabels{si}, 'FontSize', 12, 'FontWeight', 'bold');
    xlim([-2 6]); ylim(yLims_fme);
    box off; set(ax, 'TickDir', 'out', 'FontSize', 11);

    if si == 1
        legHandles = gobjects(nLev, 1);
        legLabels  = cell(nLev, 1);
        for iInt = 1:nLev
            gv     = linspace(grayLo, grayHi, nLev);
            cShade = repmat(gv(iInt), 1, 3);
            legHandles(iInt) = plot(nan, nan, '-', 'Color', cShade, 'LineWidth', 2.2);
            legLabels{iInt}  = sprintf('%g dB', ivals(iInt));
        end
        legend(legHandles, legLabels, 'Location', 'best', 'FontSize', 9, 'Box', 'off');
    end
end

% AUC panel — FME
ax = subplot(1, nS+1, nS+1);
hold(ax, 'on');
nLevMax = max(cellfun(@numel, allInts));
xLo = 0.5; xHi = nLevMax + 0.5;

% Noise floor band (drawn first, behind data)
if ~isnan(nfMu(sig)) && ~isnan(nfSem(sig))
    draw_noise_floor(ax, xLo, xHi, nfMu(sig), nfSem(sig), nfColor);
end

if showRippleRef && ~isempty(rippleInts)
    nLevRip   = numel(rippleInts);
    muRipAll  = nan(nLevRip, 1);
    semRipAll = nan(nLevRip, 1);
    validRip  = false(nLevRip, 1);
    for iInt = 1:nLevRip
        vals = rippleAUC{sig}{iInt};
        vals = vals(isfinite(vals));
        if numel(vals) < 1; continue; end
        validRip(iInt)  = true;
        muRipAll(iInt)  = mean(vals, 'omitnan');
        semRipAll(iInt) = std(vals, 0, 'omitnan') ./ sqrt(max(numel(vals),1));
    end
    xRip = find(validRip);
    lRip = plot(ax, xRip, muRipAll(validRip), '-', ...
        'Color', [0.85 0.78 0.65], 'LineWidth', 1.8);
    lRip.Color(4) = 0.55;
    for iInt = xRip'
        ebRip = errorbar(ax, iInt, muRipAll(iInt), semRipAll(iInt), ...
            'o', 'Color', [0.85 0.78 0.65], ...
            'MarkerFaceColor', [0.85 0.78 0.65], ...
            'MarkerSize', 8, 'LineWidth', 1.5, 'CapSize', 5, ...
            'LineStyle', 'none');
        ebRip.Color(4) = 0.55;
    end
end

% One style per stimulus type. Was hardcoded to two entries, which broke as
% soon as a third stimulus type was passed. Cycled so any number works.
baseStyles = {'-', '--', ':', '-.'};
lineStyles = baseStyles(mod(0:numel(stimTypes)-1, numel(baseStyles)) + 1);
for si = 1:nS
    nLev     = numel(allInts{si});
    muAUC    = zeros(nLev, 1);
    semAUC   = zeros(nLev, 1);
    validInt = false(nLev, 1);
    for iInt = 1:nLev
        vals = aucVals{si,sig}{iInt};
        vals = vals(isfinite(vals));
        if numel(vals) < 2; continue; end
        validInt(iInt) = true;
        muAUC(iInt)    = mean(vals, 'omitnan');
        semAUC(iInt)   = std(vals, 0, 'omitnan') ./ sqrt(numel(vals));
    end
    xVals = find(validInt);
    plot(xVals, muAUC(validInt), lineStyles{si}, 'Color', c, 'LineWidth', 1.5);
    errorbar(xVals, muAUC(validInt), semAUC(validInt), ...
        'o', 'Color', c, 'MarkerFaceColor', c, ...
        'MarkerSize', 7, 'LineWidth', 1.5, 'CapSize', 5, 'LineStyle', 'none');

    % Significance labels: one-sided Wilcoxon vs noise floor, per intensity level
    drawnow;
    yLimNow = ylim(ax);
    yRange  = range(yLimNow);
    for iInt = xVals'
        stimV = aucVals{si,sig}{iInt};
        stimV = stimV(isfinite(stimV));
        mIdx  = aucMouseIdx{si,sig}{iInt};
        mIdx  = mIdx(isfinite(aucVals{si,sig}{iInt}));
        [pV, stars, pV2, mDiff] = noise_floor_test(stimV, mIdx, ...
            noiseFloorAUC{sig}, noiseFloorMouseIdx{sig});
        if ~isnan(pV)
            statsTable{end+1} = {signalLabels{sig}, stimLabels{si}, iInt, ...
                numel(stimV), mean(stimV,'omitnan'), ...
                std(stimV,0,'omitnan')/sqrt(numel(stimV)), pV, pV2, mDiff}; %#ok
        end
        if ~isempty(stars)
            yTop = muAUC(iInt) + semAUC(iInt) + 0.04 * yRange;
            text(ax, iInt, yTop, stars, ...
                'HorizontalAlignment', 'center', 'FontSize', 11, ...
                'FontWeight', 'bold', 'Color', 'k');
        end
    end
end

% Legend — include noise floor
legH = [plot(nan,nan,'-', 'Color',c,'LineWidth',1.8), ...
        plot(nan,nan,'--','Color',c,'LineWidth',1.8)];
legLab = stimLabels;
if showRippleRef
    legH(end+1) = plot(nan,nan,'o','Color',[0.85 0.78 0.65], ...
        'MarkerFaceColor',[0.85 0.78 0.65],'MarkerSize',8,'LineWidth',1.5);
    legLab{end+1} = 'Ripple (ref)';
end
legH(end+1) = patch(nan, nan, nfColor, 'FaceAlpha', 0.18, 'EdgeColor', 'none');
legLab{end+1} = 'Noise floor (silent)';
legend(legH, legLab, 'Location', 'best', 'FontSize', 9, 'Box', 'off');

yline(0, '-', 'Color', [0.75 0.75 0.75], 'LineWidth', 0.8);
xlim(ax, [xLo, xHi]); set(ax, 'XTick', 1:nLevMax);
xlabel('Intensity level', 'FontSize', 11);
ylabel('AUC (0-1 s, a.u.)', 'FontSize', 11);
title('AUC by Intensity', 'FontSize', 12, 'FontWeight', 'bold');
axis(ax, 'square'); box off; set(ax, 'TickDir', 'out', 'FontSize', 11);
sgtitle(sprintf('FME Responses – %s', strjoin(stimLabels, ' & ')), ...
    'FontSize', 14, 'FontWeight', 'bold');

%% =========================================================================
%  Figure 2: Pupil — 1 x (nS+1)
%% =========================================================================
sig = 2;
c   = signalColors(sig,:);
figure('Units','normalized','Position',[0.05 0.20 0.88 0.42]);

for si = 1:nS
    nLev  = numel(allInts{si});
    ivals = allInts{si};
    ax    = subplot(1, nS+1, si);
    hold(ax, 'on');

    if showRippleRef && ~isempty(rippleInts)
        ripAllMat = [];
        for iRip = 1:numel(rippleInts)
            if ~isempty(rippleMeans{sig}{iRip})
                ripAllMat = [ripAllMat; rippleMeans{sig}{iRip}]; %#ok
            end
        end
        if ~isempty(ripAllMat)
            muR  = mean(ripAllMat, 1, 'omitnan');
            semR = std( ripAllMat, 0, 1, 'omitnan') ./ sqrt(size(ripAllMat, 1));
            shadedErrorBar(tNative, muR, semR, ...
                'lineProps', {'Color', [0.85 0.78 0.65], 'LineWidth', 1.5}, ...
                'transparent', 1, 'patchSaturation', 0.25);
        end
    end

    for iInt = 1:nLev
        if isempty(mouseMeans{si,sig}{iInt}); continue; end
        mu  = mean(mouseMeans{si,sig}{iInt}, 1, 'omitnan');
        sem = std( mouseMeans{si,sig}{iInt}, 0, 1, 'omitnan') ./ ...
              sqrt(size(mouseMeans{si,sig}{iInt}, 1));
        gv     = linspace(grayLo, grayHi, nLev);
        cShade = repmat(gv(iInt), 1, 3);
        shadedErrorBar(tNative, mu, sem, ...
            'lineProps', {'Color', cShade, 'LineWidth', 1.8}, ...
            'transparent', 1, 'patchSaturation', 0.15);
    end

    xline(0, '--k', 'LineWidth', 1.0, 'Alpha', 0.5);
    xlabel('Time from onset (s)', 'FontSize', 11);
    if si == 1; ylabel('baseline sub. (a.u.)', 'FontSize', 11); end
    title(stimLabels{si}, 'FontSize', 12, 'FontWeight', 'bold');
    xlim([-2 6]); ylim(pupil_yLims{si});
    box off; set(ax, 'TickDir', 'out', 'FontSize', 11);

    if si == 1
        legHandles = gobjects(nLev, 1);
        legLabels  = cell(nLev, 1);
        for iInt = 1:nLev
            gv     = linspace(grayLo, grayHi, nLev);
            cShade = repmat(gv(iInt), 1, 3);
            legHandles(iInt) = plot(nan, nan, '-', 'Color', cShade, 'LineWidth', 2.2);
            if strcmp(intFields{si}, 'int')
                legLabels{iInt} = sprintf('%g dB', ivals(iInt));
            else
                legLabels{iInt} = sprintf('%g', ivals(iInt));
            end
        end
        legend(legHandles, legLabels, 'Location', 'best', 'FontSize', 9, 'Box', 'off');
    end
end

% AUC panel — Pupil
ax = subplot(1, nS+1, nS+1);
hold(ax, 'on');
nLevMax = max(cellfun(@numel, allInts));
xLo = 0.5; xHi = nLevMax + 0.5;

% Noise floor band (drawn first, behind data)
if ~isnan(nfMu(sig)) && ~isnan(nfSem(sig))
    draw_noise_floor(ax, xLo, xHi, nfMu(sig), nfSem(sig), nfColor);
end

if showRippleRef && ~isempty(rippleInts)
    nLevRip   = numel(rippleInts);
    muRipAll  = nan(nLevRip, 1);
    semRipAll = nan(nLevRip, 1);
    validRip  = false(nLevRip, 1);
    for iInt = 1:nLevRip
        vals = rippleAUC{sig}{iInt};
        vals = vals(isfinite(vals));
        if numel(vals) < 1; continue; end
        validRip(iInt)  = true;
        muRipAll(iInt)  = mean(vals, 'omitnan');
        semRipAll(iInt) = std(vals, 0, 'omitnan') ./ sqrt(max(numel(vals),1));
    end
    xRip = find(validRip);
    lRip = plot(ax, xRip, muRipAll(validRip), '-', ...
        'Color', [0.85 0.78 0.65], 'LineWidth', 1.8);
    lRip.Color(4) = 0.55;
    for iInt = xRip'
        ebRip = errorbar(ax, iInt, muRipAll(iInt), semRipAll(iInt), ...
            'o', 'Color', [0.85 0.78 0.65], ...
            'MarkerFaceColor', [0.85 0.78 0.65], ...
            'MarkerSize', 8, 'LineWidth', 1.5, 'CapSize', 5, ...
            'LineStyle', 'none');
        ebRip.Color(4) = 0.55;
    end
end

for si = 1:nS
    nLev     = numel(allInts{si});
    muAUC    = zeros(nLev, 1);
    semAUC   = zeros(nLev, 1);
    validInt = false(nLev, 1);
    for iInt = 1:nLev
        vals = aucVals{si,sig}{iInt};
        vals = vals(isfinite(vals));
        if numel(vals) < 2; continue; end
        validInt(iInt) = true;
        muAUC(iInt)    = mean(vals, 'omitnan');
        semAUC(iInt)   = std(vals, 0, 'omitnan') ./ sqrt(numel(vals));
    end
    xVals = find(validInt);
    plot(xVals, muAUC(validInt), lineStyles{si}, 'Color', c, 'LineWidth', 1.5);
    errorbar(xVals, muAUC(validInt), semAUC(validInt), ...
        'o', 'Color', c, 'MarkerFaceColor', c, ...
        'MarkerSize', 7, 'LineWidth', 1.5, 'CapSize', 5, 'LineStyle', 'none');

    % Significance labels: one-sided Wilcoxon vs noise floor, per intensity level
    drawnow;
    yLimNow = ylim(ax);
    yRange  = range(yLimNow);
    for iInt = xVals'
        stimV = aucVals{si,sig}{iInt};
        stimV = stimV(isfinite(stimV));
        mIdx  = aucMouseIdx{si,sig}{iInt};
        mIdx  = mIdx(isfinite(aucVals{si,sig}{iInt}));
        [pV, stars, pV2, mDiff] = noise_floor_test(stimV, mIdx, ...
            noiseFloorAUC{sig}, noiseFloorMouseIdx{sig});
        if ~isnan(pV)
            statsTable{end+1} = {signalLabels{sig}, stimLabels{si}, iInt, ...
                numel(stimV), mean(stimV,'omitnan'), ...
                std(stimV,0,'omitnan')/sqrt(numel(stimV)), pV, pV2, mDiff}; %#ok
        end
        if ~isempty(stars)
            yTop = muAUC(iInt) + semAUC(iInt) + 0.04 * yRange;
            text(ax, iInt, yTop, stars, ...
                'HorizontalAlignment', 'center', 'FontSize', 11, ...
                'FontWeight', 'bold', 'Color', 'k');
        end
    end
end

% Legend — include noise floor
legH = [plot(nan,nan,'-', 'Color',c,'LineWidth',1.8), ...
        plot(nan,nan,'--','Color',c,'LineWidth',1.8)];
legLab = stimLabels;
if showRippleRef
    legH(end+1) = plot(nan,nan,'o','Color',[0.85 0.78 0.65], ...
        'MarkerFaceColor',[0.85 0.78 0.65],'MarkerSize',8,'LineWidth',1.5);
    legLab{end+1} = 'Ripple (ref)';
end
legH(end+1) = patch(nan, nan, nfColor, 'FaceAlpha', 0.18, 'EdgeColor', 'none');
legLab{end+1} = 'Noise floor (silent)';
legend(legH, legLab, 'Location', 'best', 'FontSize', 9, 'Box', 'off');

yline(0, '-', 'Color', [0.75 0.75 0.75], 'LineWidth', 0.8);
xlim(ax, [xLo, xHi]); set(ax, 'XTick', 1:nLevMax);
xlabel('Intensity level', 'FontSize', 11);
ylabel('AUC (0-1 s, a.u.)', 'FontSize', 11);
title('AUC by Intensity', 'FontSize', 12, 'FontWeight', 'bold');
axis(ax, 'square'); box off; set(ax, 'TickDir', 'out', 'FontSize', 11);
sgtitle(sprintf('Pupil Responses – %s', strjoin(stimLabels, ' & ')), ...
    'FontSize', 14, 'FontWeight', 'bold');

%% =========================================================================
%  Figures 3 & 4: Diagnostic — per-mouse means, highest intensity, no processing
%% =========================================================================
rawSources  = {'FME – per-mouse means, highest intensity (no processing)', ...
               'Pupil – per-mouse means, highest intensity (no processing)'};
rawData_all = {fme_stim, pupil_stim};
mouseColorMap = lines(nMice);

for sig = 1:nSig
    src = rawData_all{sig};
    figure('Units','normalized','Position',[0.05 0.08 0.70 0.42]);

    for si = 1:nS
        ax         = subplot(1, nS, si);
        hold(ax, 'on');
        ivals      = allInts{si};
        highestInt = ivals(end);
        nDrawn     = 0;

        for m = 1:nMice
            if isempty(src{m}) || ~isfield(src{m}, stimTypes{si}); continue; end
            traces = src{m}.(stimTypes{si}).traces;
            ints   = src{m}.(stimTypes{si}).(intFields{si});
            if iscell(traces); traces = traces{1}; end
            trialIdx = ints == highestInt;
            if sum(trialIdx) < 1; continue; end
            muMouse = mean(traces(trialIdx, :), 1, 'omitnan');
            plot(ax, tNative, muMouse, '-', ...
                'Color', mouseColorMap(m,:), 'LineWidth', 1.8);
            nDrawn = nDrawn + 1;
        end

        xline(ax, 0, '--k', 'LineWidth', 1.0, 'Alpha', 0.5);
        xlabel('Time from onset (s)', 'FontSize', 11);
        if si == 1; ylabel('Raw signal (a.u., unprocessed)', 'FontSize', 11); end
        if strcmp(intFields{si}, 'int')
            intStr = sprintf('%g dB', highestInt);
        else
            intStr = sprintf('contrast %g', highestInt);
        end
        title(sprintf('%s  -  %s  (n=%d mice)', stimLabels{si}, intStr, nDrawn), ...
            'FontSize', 11, 'FontWeight', 'bold');
        xlim([-2 6]); box off; set(ax, 'TickDir', 'out', 'FontSize', 11);

        legHandles = []; legLabels = {};
        for m = 1:nMice
            if isempty(src{m}) || ~isfield(src{m}, stimTypes{si}); continue; end
            ints = src{m}.(stimTypes{si}).(intFields{si});
            if sum(ints == highestInt) < 1; continue; end
            legHandles(end+1) = plot(nan, nan, '-', ...
                'Color', mouseColorMap(m,:), 'LineWidth', 2.0); %#ok
            legLabels{end+1}  = sprintf('Mouse %d', m); %#ok
        end
        legend(legHandles, legLabels, 'Location', 'best', 'FontSize', 9, 'Box', 'off');
    end
    sgtitle(rawSources{sig}, 'FontSize', 13, 'FontWeight', 'bold');
end

%% =========================================================================
%  Save noise-floor comparison stats to CSV
%  FDR correction (Benjamini-Hochberg) across all signal × stim × intensity
%  Unit: per-mouse mean AUC matched to noise-floor (paired Wilcoxon, one-sided)
%% =========================================================================
if ~isempty(statsTable)
    pRaws     = cellfun(@(r) r{7}, statsTable);
    nT        = numel(pRaws);
    [ps, idx] = sort(pRaws(:));
    adj       = ps .* nT ./ (1:nT)';
    adj       = flipud(cummin(flipud(adj)));
    adj       = min(adj, 1);
    pAdj      = nan(nT,1); pAdj(idx) = adj;

    csvPath = fullfile(pwd, 'fme_pupil_noise_floor_stats.csv');
    fid = fopen(csvPath, 'w');
    % Test identity written on every row so the stats spreadsheet is
    % self-describing. NOTE the FDR family here is broader than in Fig 2:
    % correction is across ALL signal x stimulus x intensity combinations in
    % one family, not within stimulus type. That is a deliberate choice but
    % differs from the convention elsewhere, so it is recorded explicitly.
    % FDR for the two-sided family, computed separately — each tail is its
    % own family of tests and correcting them together would be wrong.
    pRaw2 = cellfun(@(r) r{8}, statsTable);
    pAdj2 = nan(numel(pRaw2),1);
    fin2  = isfinite(pRaw2);
    if any(fin2)
        pf = pRaw2(fin2); nf = numel(pf);
        [ps,ix] = sort(pf(:));
        adj = ps .* nf ./ (1:nf)';
        adj = flipud(cummin(flipud(adj)));
        adj = min(adj,1);
        tmp = nan(nf,1); tmp(ix) = adj;
        pAdj2(fin2) = tmp;
    end

    fprintf(fid, ['cohort,test,comparison,MC_correction,' ...
                  'correction_family,unit_of_analysis,' ...
                  'signal,stimulus,intensity_level,n_mice,mean_AUC,sem_AUC,' ...
                  'mean_diff_vs_noise,' ...
                  'p_raw_onesided,p_FDR_onesided,sig_FDR05_onesided,' ...
                  'p_raw_twosided,p_FDR_twosided,sig_FDR05_twosided\n']);
    for k = 1:numel(statsTable)
        row = statsTable{k};
        % row{1} is the signal name (FME or pupil) and doubles as the cohort.
        % Both tails reported: the one-sided test drives the figure asterisks
        % but cannot detect a response BELOW the noise floor (e.g. pupil
        % constriction to visual gratings). mean_diff_vs_noise gives the
        % direction.
        fprintf(fid, ['%s,Wilcoxon signed-rank (paired),' ...
                      'vs silent-trial noise floor,Benjamini-Hochberg FDR,' ...
                      'all signals x stimuli x intensities (each tail its own family),' ...
                      'per mouse,' ...
                      '%s,%s,%d,%d,%.4f,%.4f,%.4f,%.4f,%.4f,%d,%.4f,%.4f,%d\n'], ...
            row{1}, row{1}, row{2}, row{3}, row{4}, row{5}, row{6}, ...
            row{9}, ...
            row{7}, pAdj(k),  pAdj(k)<0.05, ...
            row{8}, pAdj2(k), pAdj2(k)<0.05);
    end
    fclose(fid);
    fprintf('Noise-floor stats saved to: %s\n', csvPath);
end

end
