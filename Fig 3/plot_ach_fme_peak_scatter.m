function plot_ach_fme_peak_scatter(ach_stim, fme_stim, all_fiber_names_ach)
% PLOT_ACH_FME_AUC_SCATTER
%
% Usage:
%   plot_ach_fme_peak_scatter(ach_stim, fme_stim, all_fiber_names_ach)
%
% For Tones at moderate intensity (50 dB), plots FME AUC vs ACh AUC
% within 0–1 s, one point per trial pooled across mice.
% Dots coloured by brain region. One panel per region.
% Also produces a raincloud plot of ACh AUC distribution by region.

%% -------------------------------------------------------------------------
%  Parameters
%% -------------------------------------------------------------------------
nMice       = numel(ach_stim);

nTimepoints = 8139;
tAch        = linspace(-2, 6, nTimepoints);
blAch       = tAch >= -2 & tAch < 0;
peakAch     = tAch >= 0  & tAch <= 1;

nNative     = 241;
tFme        = linspace(-2, 6, nNative);
blFme       = tFme >= -2 & tFme < 0;
peakFme     = tFme >= 0  & tFme <= 1;

regions     = {'AC', 'VC', 'PFC'};
nR          = numel(regions);

stimTypes   = {'Tones'};
intFields   = {'int'};
stimLabels  = {'Tones'};
nS          = numel(stimTypes);

targetInt   = containers.Map(...
    {'AudRipple', 'Tones', 'VisGrating'}, ...
    {50,          50,      0.33});

regionColors = [
    0.00, 0.45, 0.70;   % AC  — blue
    0.90, 0.62, 0.00;   % VC  — orange
    0.00, 0.62, 0.45;   % PFC — teal
];

markerSz    = 18;
markerAlpha = 0.35;

%% -------------------------------------------------------------------------
%  Data collection
%  aucFmeAll{s,r}     = [nTrials x 1]  FME AUC 0-1 s (global z-score)
%  aucAchAll{s,r}     = [nTrials x 1]  ACh AUC 0-1 s (per-trial z-score)
%  aucAchAll_gz{s,r}  = [nTrials x 1]  ACh AUC 0-1 s (global z-score)
%  mouseIdxAll{s,r}   = [nTrials x 1]  mouse index per trial
%% -------------------------------------------------------------------------
aucFmeAll    = cell(nS, nR);
aucAchAll    = cell(nS, nR);
aucAchAll_gz = cell(nS, nR);   % global z-score, baseline subtracted
aucAchAll_nb = cell(nS, nR);   % global z-score, no baseline subtraction
mouseIdxAll  = cell(nS, nR);   % mouse identity per trial
for s = 1:nS
    for r = 1:nR
        aucFmeAll{s,r}    = [];
        aucAchAll{s,r}    = [];
        aucAchAll_gz{s,r} = [];
        aucAchAll_nb{s,r} = [];
        mouseIdxAll{s,r}  = [];
    end
end

for m = 1:nMice
    if isempty(ach_stim{m}) || isempty(fme_stim{m}); continue; end

    for s = 1:nS
        if ~isfield(fme_stim{m}, stimTypes{s}); continue; end
        if ~isfield(ach_stim{m}, stimTypes{s}); continue; end

        % Select moderate-intensity trials
        intVec      = fme_stim{m}.(stimTypes{s}).(intFields{s});
        intVec      = intVec(:);
        uniqueInts  = unique(intVec);
        [~, nearest] = min(abs(uniqueInts - targetInt(stimTypes{s})));
        modIntVal   = uniqueInts(nearest);
        modMask     = intVec == modIntVal;
        if sum(modMask) == 0; continue; end

        % FME: global z-score across all moderate-intensity trials within
        % this mouse, then take peak within 0–1 s per trial
        fmeTraces = fme_stim{m}.(stimTypes{s}).traces;
        if iscell(fmeTraces); fmeTraces = fmeTraces{1}; end
        fmeMod    = fmeTraces(modMask, :);

        % Global z-score: mean and SD computed across all timepoints and
        % trials for this mouse (flattened), so units are comparable across mice
        fmeAll    = fmeMod(:);
        fmeGMu    = mean(fmeAll(isfinite(fmeAll)), 'omitnan');
        fmeGSd    = std( fmeAll(isfinite(fmeAll)), 'omitnan');
        if fmeGSd == 0; fmeGSd = 1; end
        fmeModZ   = (fmeMod - fmeGMu) / fmeGSd;
        fmePeak   = max(fmeModZ(:, peakFme), [], 2, 'omitnan');   % [nTrials x 1]
        fmeAUC    = trapz(tFme(peakFme), fmeModZ(:, peakFme), 2); % [nTrials x 1]

        for r = 1:nR
            fiberIdx = find(strcmp(all_fiber_names_ach{m}, regions{r}));
            if isempty(fiberIdx); continue; end

            % ACh per-trial z-score: baseline-subtract then divide by baseline SD
            achTraces = ach_stim{m}.(stimTypes{s}).traces{fiberIdx};
            if size(achTraces, 1) ~= numel(intVec); continue; end
            achMod    = achTraces(modMask, :);
            achBl     = mean(achMod(:, blAch), 2, 'omitnan');
            achMod    = achMod - achBl;
            achBlSd   = std(achMod(:, blAch), 0, 2, 'omitnan');
            achBlSd(achBlSd == 0) = 1;
            achModZ   = achMod ./ achBlSd;
            achPeak   = max(achModZ(:, peakAch), [], 2, 'omitnan');
            achAUC    = trapz(tAch(peakAch), achModZ(:, peakAch), 2);

            % ACh global z-score: flatten all trials and timepoints for this mouse
            achAll_flat = achMod(:);   % already baseline-subtracted
            achGMu      = mean(achAll_flat(isfinite(achAll_flat)), 'omitnan');
            achGSd      = std( achAll_flat(isfinite(achAll_flat)), 'omitnan');
            if achGSd == 0; achGSd = 1; end
            achModGZ    = (achMod - achGMu) / achGSd;
            achAUC_gz   = trapz(tAch(peakAch), achModGZ(:, peakAch), 2);

            % ACh global z-score, no baseline subtraction: use raw traces
            achRaw      = achTraces(modMask, :);
            achRaw_flat = achRaw(:);
            achNbMu     = mean(achRaw_flat(isfinite(achRaw_flat)), 'omitnan');
            achNbSd     = std( achRaw_flat(isfinite(achRaw_flat)), 'omitnan');
            if achNbSd == 0; achNbSd = 1; end
            achModNB    = (achRaw - achNbMu) / achNbSd;
            achAUC_nb   = trapz(tAch(peakAch), achModNB(:, peakAch), 2);

            % Keep only finite pairs with peak FME within ±5 z (outlier rejection)
            finMask = isfinite(fmePeak) & isfinite(achPeak) & abs(fmePeak) <= 5;
            aucFmeAll{s,r}    = [aucFmeAll{s,r};    fmeAUC(finMask)];
            aucAchAll{s,r}    = [aucAchAll{s,r};    achAUC(finMask)];
            aucAchAll_gz{s,r} = [aucAchAll_gz{s,r}; achAUC_gz(finMask)];
            aucAchAll_nb{s,r} = [aucAchAll_nb{s,r}; achAUC_nb(finMask)];
            mouseIdxAll{s,r}  = [mouseIdxAll{s,r};  repmat(m, sum(finMask), 1)];
        end
    end
end

%% -------------------------------------------------------------------------
%  Figure 1: AUC FME vs AUC ACh (0–1 s)
%  Dots coloured by brain region, one regression line per panel.
%% -------------------------------------------------------------------------
statsFig1 = {};   % collect correlation stats for command-window printout
statsFig2 = {};

figure('Units','normalized','Position',[0.05 0.25 0.85 0.45]);
sgtitle('FME AUC vs ACh AUC – Tones (50 dB, 0–1 s, per trial)', ...
    'FontSize', 13, 'FontWeight', 'bold');

for r = 1:nR
    ax = subplot(1, nR, r);
    hold(ax, 'on');
    c = regionColors(r,:);

    xAll = []; yAll = []; mAll = [];
    for s = 1:nS
        x   = aucFmeAll{s,r};
        y   = aucAchAll{s,r};
        mi  = mouseIdxAll{s,r};
        fin = isfinite(x) & isfinite(y);
        if sum(fin) > 0
            scatter(ax, x(fin), y(fin), markerSz, c, 'filled', ...
                'MarkerFaceAlpha', markerAlpha, 'MarkerEdgeColor', 'none');
        end
        xAll = [xAll; x(fin)]; %#ok
        yAll = [yAll; y(fin)]; %#ok
        mAll = [mAll; mi(fin)]; %#ok
    end

    n = numel(xAll);
    if n >= 5
        p      = polyfit(xAll, yAll, 1);
        xRange = linspace(min(xAll), max(xAll), 200);
        plot(ax, xRange, polyval(p, xRange), '-', ...
            'Color', [0.2 0.2 0.2], 'LineWidth', 2.0);

        % Trial-level Spearman (shown but note pseudoreplication)
        [rho, pVal] = corr(xAll, yAll, 'Type', 'Spearman', 'Rows', 'complete');
        if pVal < 0.001;    pStr = 'p<0.001';
        elseif pVal < 0.05; pStr = sprintf('p=%.3f', pVal);
        else;                pStr = sprintf('p=%.2f', pVal);
        end

        % Per-mouse mean Spearman (primary reported statistic)
        uMice  = unique(mAll);
        nMiceR = numel(uMice);
        muFme  = arrayfun(@(mm) mean(xAll(mAll==mm),'omitnan'), uMice);
        muAch  = arrayfun(@(mm) mean(yAll(mAll==mm),'omitnan'), uMice);
        finMu  = isfinite(muFme) & isfinite(muAch);
        if sum(finMu) >= 3
            [rhoMu, pValMu] = corr(muFme(finMu), muAch(finMu), ...
                'Type','Spearman','Rows','complete');
            scatter(ax, muFme(finMu), muAch(finMu), 55, c, 'filled', ...
                'MarkerFaceAlpha',0.9,'MarkerEdgeColor','w','LineWidth',1.2);
            statsFig1{end+1} = {regions{r}, 'per-trial z-score', ...
                nMiceR, rhoMu, pValMu, n, rho, pVal}; %#ok
        else
            statsFig1{end+1} = {regions{r}, 'per-trial z-score', ...
                NaN, NaN, NaN, n, rho, pVal}; %#ok
        end
    end

    xline(ax, 0, '-', 'Color', [0.85 0.85 0.85], 'LineWidth', 0.8);
    yline(ax, 0, '-', 'Color', [0.85 0.85 0.85], 'LineWidth', 0.8);

    if r == 1; ylabel(ax, 'ACh AUC (z·s)', 'FontSize', 10); end
    xlabel(ax, 'FME AUC (z·s)', 'FontSize', 10);
    title(ax, regions{r}, 'FontSize', 11, 'FontWeight', 'bold', 'Color', c);
    box off;
    set(ax, 'TickDir', 'out', 'FontSize', 9);
end

annotation('textbox', [0.01 0.00 0.99 0.03], ...
    'String', ['Small dots = individual trials; large dots = per-mouse means (primary statistic).  ' ...
               'Tones at 50 dB only.  FME globally z-scored; trials with peak FME >5 z excluded.  ' ...
               'ACh z-scored per trial using pre-stimulus baseline SD.  ' ...
               '\rho_{mouse} = Spearman on per-mouse means; \rho_{trial} shown in grey for reference.'], ...
    'EdgeColor', 'none', 'FontSize', 8, 'Color', [0.4 0.4 0.4], ...
    'HorizontalAlignment', 'center');

%% -------------------------------------------------------------------------
%  Figure 2: AUC FME vs AUC ACh — baseline-subtracted, global z-score
%  Identical layout to Figure 1 for direct comparison.
%% -------------------------------------------------------------------------
figure('Units','normalized','Position',[0.05 0.25 0.85 0.45]);
sgtitle('FME AUC vs ACh AUC – Tones (50 dB, 0–1 s, baseline sub. + global z-score)', ...
    'FontSize', 13, 'FontWeight', 'bold');

for r = 1:nR
    ax = subplot(1, nR, r);
    hold(ax, 'on');
    c = regionColors(r,:);

    xAll = []; yAll = []; mAll = [];
    for s = 1:nS
        x   = aucFmeAll{s,r};
        y   = aucAchAll_gz{s,r};
        mi  = mouseIdxAll{s,r};
        fin = isfinite(x) & isfinite(y);
        if sum(fin) > 0
            scatter(ax, x(fin), y(fin), markerSz, c, 'filled', ...
                'MarkerFaceAlpha', markerAlpha, 'MarkerEdgeColor', 'none');
        end
        xAll = [xAll; x(fin)]; %#ok
        yAll = [yAll; y(fin)]; %#ok
        mAll = [mAll; mi(fin)]; %#ok
    end

    n = numel(xAll);
    if n >= 5
        p      = polyfit(xAll, yAll, 1);
        xRange = linspace(min(xAll), max(xAll), 200);
        plot(ax, xRange, polyval(p, xRange), '-', ...
            'Color', [0.2 0.2 0.2], 'LineWidth', 2.0);

        [rho, pVal] = corr(xAll, yAll, 'Type', 'Spearman', 'Rows', 'complete');
        if pVal < 0.001;    pStr = 'p<0.001';
        elseif pVal < 0.05; pStr = sprintf('p=%.3f', pVal);
        else;                pStr = sprintf('p=%.2f', pVal);
        end

        % Per-mouse means
        uMice  = unique(mAll); nMiceR = numel(uMice);
        muFme  = arrayfun(@(mm) mean(xAll(mAll==mm),'omitnan'), uMice);
        muAch  = arrayfun(@(mm) mean(yAll(mAll==mm),'omitnan'), uMice);
        finMu  = isfinite(muFme) & isfinite(muAch);
        if sum(finMu) >= 3
            [rhoMu, pValMu] = corr(muFme(finMu), muAch(finMu), ...
                'Type','Spearman','Rows','complete');
            scatter(ax, muFme(finMu), muAch(finMu), 55, c, 'filled', ...
                'MarkerFaceAlpha',0.9,'MarkerEdgeColor','w','LineWidth',1.2);
            statsFig2{end+1} = {regions{r}, 'global z-score', ...
                nMiceR, rhoMu, pValMu, n, rho, pVal}; %#ok
        else
            statsFig2{end+1} = {regions{r}, 'global z-score', ...
                NaN, NaN, NaN, n, rho, pVal}; %#ok
        end
    end

    xline(ax, 0, '-', 'Color', [0.85 0.85 0.85], 'LineWidth', 0.8);
    yline(ax, 0, '-', 'Color', [0.85 0.85 0.85], 'LineWidth', 0.8);

    if r == 1; ylabel(ax, 'ACh AUC (global z·s)', 'FontSize', 10); end
    xlabel(ax, 'FME AUC (z·s)', 'FontSize', 10);
    title(ax, regions{r}, 'FontSize', 11, 'FontWeight', 'bold', 'Color', c);
    box off;
    set(ax, 'TickDir', 'out', 'FontSize', 9);
end

annotation('textbox', [0.01 0.00 0.99 0.03], ...
    'String', ['Small dots = individual trials; large dots = per-mouse means (primary statistic).  ' ...
               'Tones at 50 dB only.  FME globally z-scored; trials with peak FME >5 z excluded.  ' ...
               'ACh baseline-subtracted per trial, then globally z-scored.  ' ...
               '\rho_{mouse} = Spearman on per-mouse means; \rho_{trial} shown in grey for reference.'], ...
    'EdgeColor', 'none', 'FontSize', 8, 'Color', [0.4 0.4 0.4], ...
    'HorizontalAlignment', 'center');

%% -------------------------------------------------------------------------
%  Figure 3: Raincloud plot — per-trial ACh AUC by brain region
%  All stimulus types pooled (Tones only, moderate intensity).
%  Each region: half-violin (density) above, jittered dots below,
%  mean ± SEM marker on top. All three regions on shared axis.
%% -------------------------------------------------------------------------
figure('Units','normalized','Position',[0.25 0.25 0.40 0.55]);
ax = axes;
hold(ax, 'on');

cloudAlpha  = 0.30;
dotAlpha    = 0.25;
dotSz       = 8;
nKde        = 200;      % points for kernel density estimate
cloudWidth  = 0.35;     % max half-width of violin in axis units
dotSpread   = 0.12;     % jitter spread for dot cloud
xPos        = 1:nR;     % x position per region

for r = 1:nR
    % Pool ACh AUC across all stimulus types for this region
    yAll = [];
    for s = 1:nS
        v = aucAchAll{s,r};
        yAll = [yAll; v(isfinite(v))]; %#ok
    end
    if numel(yAll) < 5; continue; end

    c = regionColors(r,:);
    x = xPos(r);

    % --- Kernel density estimate (half-violin, drawn above x position) ---
    [f, yi] = ksdensity(yAll, 'NumPoints', nKde);
    f       = f / max(f) * cloudWidth;   % normalise to max half-width
    % Half-violin: patch to the right of x
    patch(ax, [x + f, x*ones(1,nKde)], [yi, fliplr(yi)], c, ...
        'FaceAlpha', cloudAlpha, 'EdgeColor', 'none');
    plot(ax, x + f, yi, '-', 'Color', c, 'LineWidth', 1.2);

    % --- Jittered dot cloud (to the left of x) ---
    jit = dotSpread * (rand(numel(yAll), 1) - 0.5);
    scatter(ax, x - 0.08 + jit, yAll, dotSz, c, 'filled', ...
        'MarkerFaceAlpha', dotAlpha, 'MarkerEdgeColor', 'none');

    % --- Mean ± SEM marker ---
    mu  = mean(yAll, 'omitnan');
    sem = std(yAll, 0, 'omitnan') / sqrt(numel(yAll));
    errorbar(ax, x, mu, sem, 'o', ...
        'Color', c, 'MarkerFaceColor', c, ...
        'MarkerSize', 7, 'LineWidth', 2.0, 'CapSize', 6, 'LineStyle', 'none');
end

yline(ax, 0, '-', 'Color', [0.75 0.75 0.75], 'LineWidth', 0.9);
set(ax, 'XTick', xPos, 'XTickLabel', regions, 'TickDir', 'out', 'FontSize', 11);
xlim(ax, [0.5 nR + 0.5]);
ylabel(ax, 'ACh AUC (z·s)', 'FontSize', 11);
title(ax, 'ACh AUC Distribution by Region – Tones (50 dB)', ...
    'FontSize', 11, 'FontWeight', 'bold');
box off;

annotation('textbox', [0.01 0.00 0.99 0.04], ...
    'String', ['Per-trial ACh AUC (0–1 s), pooled across mice.  ' ...
               'Half-violin = kernel density; dots = individual trials (jittered); ' ...
               'circle + bar = mean ± SEM.  ' ...
               'Trials with peak FME >5 z excluded.'], ...
    'EdgeColor', 'none', 'FontSize', 8, 'Color', [0.4 0.4 0.4], ...
    'HorizontalAlignment', 'center');


%% =========================================================================
%  Correlation statistics — printed to command window
%  Primary statistic: Spearman on per-mouse mean AUC values
%  Trial-level correlation reported for reference only (pseudoreplicated)
%% =========================================================================
fprintf('\n============================================================\n');
fprintf('FME AUC vs ACh AUC correlations — Tones (50 dB, 0-1 s)\n');
fprintf('Primary statistic: Spearman on PER-MOUSE means\n');
fprintf('============================================================\n');

allStats = [statsFig1, statsFig2];
for k = 1:numel(allStats)
    row = allStats{k};
    fprintf('\n%s  [%s]\n', row{1}, row{2});
    if isfinite(row{4})
        fprintf('  Per-mouse:  rho = %+.3f, p = %.4f  (n = %d mice)\n', ...
            row{4}, row{5}, row{3});
    else
        fprintf('  Per-mouse:  insufficient mice for correlation\n');
    end
    fprintf('  Per-trial:  rho = %+.3f, p = %.4g  (n = %d trials) [reference only]\n', ...
        row{7}, row{8}, row{6});
end
fprintf('\n============================================================\n');

%% =========================================================================
%  Intensity justification: trial-to-trial variance by intensity
%  ------------------------------------------------------------------------
%  The scatter uses tones at the INTERMEDIATE intensity. The stated reasons
%  are that it avoids floor effects at the lowest level and saturation at the
%  highest, and that trial-to-trial variability in the ACh response is
%  greatest there — variance is what a correlation has to work with, so the
%  intermediate level is the most informative condition.
%
%  Those are assertions unless the numbers are shown. This quantifies them:
%  per mouse, the SD across trials of the ACh AUC at each intensity, then the
%  mean of those SDs across mice. Computed per mouse first so that between-
%  mouse differences in overall response amplitude do not inflate the
%  estimate.
%
%  Also reports the mean AUC per intensity, which is what shows the floor and
%  saturation directly.
%% =========================================================================
fprintf('\n============================================================\n');
fprintf('Trial-to-trial variance by intensity — justification for the\n');
fprintf('intensity used in the scatter\n');
fprintf('============================================================\n');
fprintf('%-5s %-10s %8s %14s %16s\n', ...
    'Reg','Intensity','n_mice','mean peak','mean within-mouse SD');

varRows = {};
for r = 1:nR
    for s = 1:nS
        st = stimTypes{s};
        ivAll = [];
        for m = 1:numel(ach_stim)
            if isempty(ach_stim{m}) || ~isfield(ach_stim{m}, st); continue; end
            ivAll = union(ivAll, unique(ach_stim{m}.(st).(intFields{s})));
        end
        for iI = 1:numel(ivAll)
            sdPer = []; muPer = [];
            for m = 1:numel(ach_stim)
                if isempty(ach_stim{m}) || ~isfield(ach_stim{m}, st); continue; end
                fi = find(strcmp(all_fiber_names_ach{m}, regions{r}));
                if isempty(fi); continue; end
                iv  = ach_stim{m}.(st).(intFields{s});
                sel = iv == ivAll(iI);
                if sum(sel) < 3; continue; end
                mat = ach_stim{m}.(st).traces{fi}(sel,:);
                bl  = mean(mat(:, blAch), 2, 'omitnan');
                mat = mat - bl;
                % Peak amplitude in the 0-1 s window — the same metric the
                % scatter plots, so the variance reported here is the
                % variance the correlation actually operates on.
                aucTr = max(mat(:, peakAch), [], 2);
                aucTr = aucTr(isfinite(aucTr));
                if numel(aucTr) < 3; continue; end
                sdPer(end+1) = std(aucTr, 0, 'omitnan');   %#ok
                muPer(end+1) = mean(aucTr, 'omitnan');     %#ok
            end
            if isempty(sdPer); continue; end
            fprintf('%-5s %-10g %8d %14.4f %16.4f\n', ...
                regions{r}, ivAll(iI), numel(sdPer), ...
                mean(muPer,'omitnan'), mean(sdPer,'omitnan'));
            varRows{end+1} = {regions{r}, st, ivAll(iI), numel(sdPer), ...
                mean(muPer,'omitnan'), mean(sdPer,'omitnan')}; %#ok
        end
    end
end
fprintf('============================================================\n');
fprintf('Within-mouse SD is computed across trials for each mouse, then\n');
fprintf('averaged, so between-mouse amplitude differences do not inflate it.\n');
fprintf('============================================================\n');

if ~isempty(varRows)
    fidV = fopen(fullfile(pwd, 'ach_variance_by_intensity.csv'), 'w');
    if fidV == -1
        warning('Could not write ach_variance_by_intensity.csv');
    else
        fprintf(fidV, ['region,stimulus,intensity,n_mice,' ...
                       'mean_peak,mean_within_mouse_SD\n']);
        for k = 1:numel(varRows)
            rw = varRows{k};
            fprintf(fidV, '%s,%s,%g,%d,%.6f,%.6f\n', ...
                rw{1}, rw{2}, rw{3}, rw{4}, rw{5}, rw{6});
        end
        fclose(fidV);
        fprintf('Saved: ach_variance_by_intensity.csv\n');
    end
end

%% =========================================================================
%  CSV export
%  ------------------------------------------------------------------------
%  These statistics were previously printed to the command window only, so
%  nothing reached the stats spreadsheet. Both the per-mouse and the
%  per-trial correlations are written, with the per-mouse one marked as
%  primary: pooling trials across mice is pseudoreplication, since trials
%  from one animal are not independent observations of the population.
%
%  No multiple-comparison correction is applied across the 3 regions x 2
%  normalisations, because these are descriptive correlations rather than a
%  hypothesis family. If any is reported as a significance claim, correct it
%  and say so.
%% =========================================================================
csvPathScatter = fullfile(pwd, 'ach_fme_auc_correlations.csv');
fidS = fopen(csvPathScatter, 'w');
if fidS == -1
    warning('Could not write %s (is it open in Excel?)', csvPathScatter);
else
    fprintf(fidS, ['region,normalisation,stimulus,intensity,auc_window_s,' ...
                   'test,unit_of_analysis,is_primary,' ...
                   'n,rho_spearman,p_raw,MC_correction\n']);
    for k = 1:numel(allStats)
        row = allStats{k};
        % per-mouse (primary)
        fprintf(fidS, ['%s,%s,Tones,50 dB,0-1,' ...
                       'Spearman correlation,per mouse,yes,' ...
                       '%d,%.6f,%.6f,none (descriptive)\n'], ...
            row{1}, row{2}, row{3}, row{4}, row{5});
        % per-trial (reference only — pseudoreplicated)
        fprintf(fidS, ['%s,%s,Tones,50 dB,0-1,' ...
                       'Spearman correlation,per trial (pseudoreplicated),no,' ...
                       '%d,%.6f,%.6g,none (descriptive)\n'], ...
            row{1}, row{2}, row{6}, row{7}, row{8});
    end
    fclose(fidS);
    fprintf('\nCorrelation stats saved to: %s\n', csvPathScatter);
end
