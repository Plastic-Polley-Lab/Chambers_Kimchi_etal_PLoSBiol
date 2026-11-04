function plot_ach_responses(ach_stim, all_fiber_names_ach, ...
                            ach_stim_ref, gfp_stim, all_fiber_names_gfp, ...
                            stimTypes_in, aucWinSec)
% PLOT_ACH_RESPONSES  — Manuscript Figure 2
%
% Usage:
%   plot_ach_responses(ach_stim, all_fiber_names_ach, ...
%       ach_stim_ref, gfp_stim, all_fiber_names_gfp)
%   plot_ach_responses(..., stimTypes_in)
%
% Figure 1: trace panels (nS rows x nR columns), -0.5 to +1 s, shared y-axis
% Figure 2: AUC errorbar panels (nS rows x 2 cols: ACh+405nm | GFP)
%           with significance markers below x-axis
% CSV:      ach_responses_auc_stats[_<window>].csv
%
% aucWinSec (optional, default 1.0): end of the AUC window in seconds from
%   stimulus onset. Window is always 0 to aucWinSec.
%
%   Set to 0.4 for the duration-control check requested by R2 (tones were
%   0.4 s, auditory and visual gratings 2 s). IMPORTANT CAVEAT: a 0-0.4 s
%   window does NOT equate the stimuli for duration. Given GRAB-ACh3.0
%   binding kinetics and the latency of ACh release, the first 400 ms is
%   dominated by the ONSET response, which all three stimuli share; what
%   differs with duration is the sustained and offset portion, which falls
%   largely outside this window for every stimulus. So the correct claim is
%   "the regional pattern is present in the onset response, before duration
%   can plausibly differentiate the stimuli" — NOT "we controlled for
%   duration". Duration is confounded with complexity in this design and
%   cannot be fully separated post hoc; say so in the text.

%% -------------------------------------------------------------------------
%  Parameters
%% -------------------------------------------------------------------------
nMiceAch = numel(ach_stim);
nMiceRef = numel(ach_stim_ref);
nMiceGfp = numel(gfp_stim);

% ---- Normalise fibre-name nesting -------------------------------------
% all_fiber_names_gfp arrives with an extra cell layer in some datasets:
% each entry is a cell of 1x1 cells rather than a cell of char vectors, so
% strcmp(names, 'AC') compares cell against char and returns false for every
% region. Every GFP mouse is then skipped SILENTLY — no error, just empty
% panels and no rows in the stats table. Flattened here so the function works
% with either nesting.
for q = 1:numel(all_fiber_names_gfp)
    if iscell(all_fiber_names_gfp{q})
        inner = all_fiber_names_gfp{q};
        for w = 1:numel(inner)
            while iscell(inner{w}) && ~isempty(inner{w})
                inner{w} = inner{w}{1};
            end
        end
        all_fiber_names_gfp{q} = inner;
    end
end
for q = 1:numel(all_fiber_names_ach)
    if iscell(all_fiber_names_ach{q})
        inner = all_fiber_names_ach{q};
        for w = 1:numel(inner)
            while iscell(inner{w}) && ~isempty(inner{w})
                inner{w} = inner{w}{1};
            end
        end
        all_fiber_names_ach{q} = inner;
    end
end

if nargin < 7 || isempty(aucWinSec)
    aucWinSec = 1.0;                 % published window
end

nTimepoints = 8139;
tAch        = linspace(-2, 6, nTimepoints);
baselineIdx = tAch >= -2  & tAch < 0;
plotIdx     = tAch >= -0.5 & tAch <= max(1.0, aucWinSec);
aucIdx      = tAch >= 0    & tAch <= aucWinSec;

% Tag outputs by window so runs do not overwrite each other
if abs(aucWinSec - 1.0) < 1e-9
    winTag = '';
else
    winTag = sprintf('_%gms', round(aucWinSec*1000));
end
fprintf('ACh AUC window: 0-%.3f s\n', aucWinSec);
tPlot       = tAch(plotIdx);

if nargin < 6 || isempty(stimTypes_in)
    stimTypes = {'AudRipple', 'VisGrating', 'Tones'};
else
    stimTypes = stimTypes_in;
end
nS = numel(stimTypes);
gfpUsed = cell(nS, 3);   % {stim,region} = GFP mouse indices contributing
refUsed = cell(nS, 3);   % {stim,region} = 405 nm mouse indices contributing

intFieldMap = containers.Map({'AudRipple','Tones','VisGrating'}, ...
                              {'int',      'int',  'contrast'});
labelMap    = containers.Map({'AudRipple','Tones','VisGrating'}, ...
                              {'Aud. Ripple','Tones','Vis. Grating'});
intFields  = cellfun(@(s) intFieldMap(s), stimTypes, 'UniformOutput', false);
stimLabels = cellfun(@(s) labelMap(s),    stimTypes, 'UniformOutput', false);

regions     = {'AC', 'VC', 'PFC'};
nR          = numel(regions);

regionColors = [
    0.00, 0.45, 0.70;
    0.90, 0.62, 0.00;
    0.00, 0.62, 0.45;
];

%% -------------------------------------------------------------------------
%  Data collection
%  mmTracesAch{s,r}{iInt} = [nMice x nPlotPts]  per stim type AND region
%  aucAch{s,r}{iInt}     = [nMice x 1]          per region
%  aucRef{s,r}{iInt}     = [nMice x 1]
%  aucGfp{s,r}{iInt}     = [nMice x 1]
%  allIntsAll{s}         = intensity vector
%% -------------------------------------------------------------------------
mmTracesAch = cell(nS,nR);   % {s,r}{iInt} = [nMice x nPlotPts]
allIntsAll  = cell(nS,1);
aucAch      = cell(nS,nR);
aucRef      = cell(nS,nR);
aucGfp      = cell(nS,nR);

% Utility: process one trace matrix, return per-mouse mean trace and AUC
processTrials = @(tr, tidx) deal( ...
    mean(bsxfun(@rdivide, ...
         bsxfun(@minus, tr(tidx,:), mean(tr(tidx,baselineIdx),2,'omitnan')), ...
         max(std(bsxfun(@minus, tr(tidx,:), mean(tr(tidx,baselineIdx),2,'omitnan')), ...
                 0, 2, 'omitnan'), eps)), ...
    1, 'omitnan'));   % returns mean z-scored trace [1 x nT]

for s = 1:nS
    stimType = stimTypes{s};
    intField = intFields{s};

    % Find all intensity levels
    ivals = [];
    for m = 1:nMiceAch
        if isempty(ach_stim{m}); continue; end
        if ~isfield(ach_stim{m}, stimType); continue; end
        ivals = union(ivals, unique(ach_stim{m}.(stimType).(intField)));
    end
    allInts = sort(ivals);
    nLev    = numel(allInts);
    allIntsAll{s} = allInts;

    for r = 1:nR
        mmTracesAch{s,r} = cell(nLev,1);
        for iInt = 1:nLev; mmTracesAch{s,r}{iInt} = []; end
    end
    for r = 1:nR
        aucAch{s,r} = cell(nLev,1);
        aucRef{s,r} = cell(nLev,1);
        aucGfp{s,r} = cell(nLev,1);
        for iInt = 1:nLev
            aucAch{s,r}{iInt} = [];
            aucRef{s,r}{iInt} = [];
            aucGfp{s,r}{iInt} = [];
        end
    end

    % 405 nm skip counters — must be initialised BEFORE the ACh/405 loop
    % below, not alongside the GFP counters further down.
    refSkip = struct('noEntry',0,'noStim',0,'noFibre',0,'noIntVec',0,'sizeMismatch',0);
    for rr = 1:nR; refUsed{s,rr} = []; end

    % ACh and 405nm (same mice)
    for m = 1:nMiceAch
        if isempty(ach_stim{m}); continue; end
        if ~isfield(ach_stim{m}, stimType); continue; end
        intVec = ach_stim{m}.(stimType).(intField)(:);
        for r = 1:nR
            fIdx = find(strcmp(all_fiber_names_ach{m}, regions{r}));
            if isempty(fIdx); continue; end

            % ACh
            tr = ach_stim{m}.(stimType).traces{fIdx};
            if size(tr,1) ~= numel(intVec); continue; end
            for iInt = 1:nLev
                tidx = intVec == allInts(iInt);
                if sum(tidx)==0; continue; end
                mat  = tr(tidx,:);
                bl   = mean(mat(:,baselineIdx),2,'omitnan');
                mat  = mat - bl;
                blSd = std(mat(:,baselineIdx),0,2,'omitnan');
                blSd(blSd==0) = 1;
                mat  = mat ./ blSd;
                muTr = mean(mat(:,plotIdx),1,'omitnan');
                mmTracesAch{s,r}{iInt} = [mmTracesAch{s,r}{iInt}; muTr];
                aucAch{s,r}{iInt}    = [aucAch{s,r}{iInt}; ...
                    mean(trapz(tAch(aucIdx), mat(:,aucIdx), 2),'omitnan')];
            end

            % 405nm
            % The reference has its OWN failure points, independent of ACh:
            % a mouse can have ACh for this stimulus but no reference entry,
            % no reference field for this stimulus, or a reference trace with
            % a different trial count. Each previously skipped silently, which
            % is why n differs between signal and reference and between
            % stimulus types. Each is now counted and reported.
            if numel(ach_stim_ref) < m || isempty(ach_stim_ref{m})
                refSkip.noEntry = refSkip.noEntry + 1;
            elseif ~isfield(ach_stim_ref{m}, stimType)
                refSkip.noStim = refSkip.noStim + 1;
            elseif numel(ach_stim_ref{m}.(stimType).traces) < fIdx
                refSkip.noFibre = refSkip.noFibre + 1;
            elseif ~isfield(ach_stim_ref{m}.(stimType), intField)
                refSkip.noIntVec = refSkip.noIntVec + 1;
            elseif size(ach_stim_ref{m}.(stimType).traces{fIdx},1) ~= ...
                   numel(ach_stim_ref{m}.(stimType).(intField))
                % Compared against the REFERENCE's own intensity vector.
                % Previously this was compared against the ACh intensity
                % vector, which is wrong: ach_stim_ref carries its own
                % trial list, and a session where the signal lost a trial
                % during extraction but the reference did not (119 vs 120)
                % was discarded even though the reference data were valid.
                refSkip.sizeMismatch = refSkip.sizeMismatch + 1;
                if refSkip.sizeMismatch <= 5
                    fprintf(['  405 %s mouse %d %s: ref has %d trials but ' ...
                             'its own intensity vector has %d — skipped\n'], ...
                        stimType, m, regions{r}, ...
                        size(ach_stim_ref{m}.(stimType).traces{fIdx},1), ...
                        numel(ach_stim_ref{m}.(stimType).(intField)));
                end
            end
            if numel(ach_stim_ref) >= m && ~isempty(ach_stim_ref{m}) && ...
               isfield(ach_stim_ref{m}, stimType) && ...
               numel(ach_stim_ref{m}.(stimType).traces) >= fIdx
                tr2 = ach_stim_ref{m}.(stimType).traces{fIdx};
                % The reference has its own trial list — use it, not the ACh one.
                if isfield(ach_stim_ref{m}.(stimType), intField)
                    intVecRef = ach_stim_ref{m}.(stimType).(intField)(:);
                else
                    intVecRef = [];
                end
                if ~isempty(intVecRef) && size(tr2,1) == numel(intVecRef)
                    refUsed{s,r}(end+1) = m;   % contributed
                    for iInt = 1:nLev
                        tidx = intVecRef == allInts(iInt);
                        if sum(tidx)==0; continue; end
                        mat2  = tr2(tidx,:);
                        bl2   = mean(mat2(:,baselineIdx),2,'omitnan');
                        mat2  = mat2 - bl2;
                        blSd2 = std(mat2(:,baselineIdx),0,2,'omitnan');
                        blSd2(blSd2==0) = 1;
                        mat2  = mat2 ./ blSd2;
                        aucRef{s,r}{iInt} = [aucRef{s,r}{iInt}; ...
                            mean(trapz(tAch(aucIdx), mat2(:,aucIdx), 2),'omitnan')];
                    end
                end
            end
        end
    end

    % GFP (separate mice)
    % Diagnostic: every guard below silently skips a mouse, so when the GFP
    % panels come out empty there is no indication of which gate rejected it.
    % Each rejection is now reported once per stimulus type.
    gfpSkip = struct('empty',0,'noStim',0,'noFibre',0,'sizeMismatch',0,'ok',0);
    for rr = 1:nR; gfpUsed{s,rr} = []; end
    for m = 1:nMiceGfp
        if isempty(gfp_stim{m}); gfpSkip.empty = gfpSkip.empty + 1; continue; end
        if ~isfield(gfp_stim{m}, stimType)
            gfpSkip.noStim = gfpSkip.noStim + 1; continue;
        end
        intVecG = gfp_stim{m}.(stimType).(intField)(:);
        for r = 1:nR
            fIdx = find(strcmp(all_fiber_names_gfp{m}, regions{r}));
            if isempty(fIdx); gfpSkip.noFibre = gfpSkip.noFibre + 1; continue; end
            trG = gfp_stim{m}.(stimType).traces{fIdx};
            if size(trG,1) ~= numel(intVecG)
                gfpSkip.sizeMismatch = gfpSkip.sizeMismatch + 1;
                if gfpSkip.sizeMismatch <= 3
                    fprintf(['  GFP %s mouse %d %s: traces have %d rows but ' ...
                             '%s has %d entries — skipped\n'], ...
                        stimType, m, regions{r}, size(trG,1), intField, numel(intVecG));
                end
                continue;
            end
            gfpSkip.ok = gfpSkip.ok + 1;
            gfpUsed{s,r}(end+1) = m;   % which mice contributed, per stim x region
            for iInt = 1:nLev
                tidx = intVecG == allInts(iInt);
                if sum(tidx)==0; continue; end
                matG  = trG(tidx,:);
                blG   = mean(matG(:,baselineIdx),2,'omitnan');
                matG  = matG - blG;
                blSdG = std(matG(:,baselineIdx),0,2,'omitnan');
                blSdG(blSdG==0) = 1;
                matG  = matG ./ blSdG;
                aucGfp{s,r}{iInt} = [aucGfp{s,r}{iInt}; ...
                    mean(trapz(tAch(aucIdx), matG(:,aucIdx), 2),'omitnan')];
            end
        end
    end
    fprintf(['  405 %s: skipped %d no ref entry, %d no stim field, ' ...
             '%d no fibre, %d no intensity vector, %d size mismatch\n'], ...
        stimType, refSkip.noEntry, refSkip.noStim, ...
        refSkip.noFibre, refSkip.noIntVec, refSkip.sizeMismatch);
    for rr = 1:nR
        fprintf('      405 %-4s n=%d  mice: %s\n', regions{rr}, ...
            numel(refUsed{s,rr}), mat2str(refUsed{s,rr}));
    end

    fprintf(['  GFP %s: %d mouse-regions used | skipped: %d empty, ' ...
             '%d no stim field, %d no fibre in region, %d size mismatch ' ...
             '(of %d GFP mice)\n'], ...
        stimType, gfpSkip.ok, gfpSkip.empty, gfpSkip.noStim, ...
        gfpSkip.noFibre, gfpSkip.sizeMismatch, nMiceGfp);
    for rr = 1:nR
        fprintf('      %-4s n=%d  mice: %s\n', regions{rr}, ...
            numel(gfpUsed{s,rr}), mat2str(gfpUsed{s,rr}));
    end
end

%% =========================================================================
%  Figure 1: Trace panels (nS rows x nR columns)
%% =========================================================================
figure('Units','normalized','Position',[0.02 0.05 0.92 0.88]);
sgtitle('ACh Responses by Stimulus Type, Intensity & Brain Region', ...
    'FontSize',13,'FontWeight','bold');
axT = gobjects(nS,nR);

for s = 1:nS
    nLev = numel(allIntsAll{s});
    for r = 1:nR
        ax = subplot(nS, nR, (s-1)*nR + r);
        axT(s,r) = ax;
        hold(ax,'on');
        for iInt = 1:nLev
            mat = mmTracesAch{s,r}{iInt};
            if isempty(mat); continue; end
            mat = mat(all(isfinite(mat),2),:);
            if size(mat,1)<2; continue; end
            mu  = mean(mat,1,'omitnan');
            sem = std(mat,0,1,'omitnan')/sqrt(size(mat,1));
            t   = (iInt-1)/max(nLev-1,1);
            col = regionColors(r,:)*t + (regionColors(r,:)*0.25+[1 1 1]*0.75)*(1-t);
            shadedErrorBar(tPlot, mu, sem, ...
                'lineProps',{'Color',col,'LineWidth',1.6}, ...
                'transparent',1,'patchSaturation',0.12);
        end
        xline(ax,0,'--k','LineWidth',0.8,'Alpha',0.4);
        yline(ax,0,'-','Color',[0.82 0.82 0.82],'LineWidth',0.7);
        xlim(ax,[-0.5 1.0]);
        if r==1; ylabel(ax,{stimLabels{s};'ACh (z)'},'FontSize',9,'FontWeight','bold'); end
        if s==1; title(ax,regions{r},'FontSize',11,'FontWeight','bold','Color',regionColors(r,:)); end
        if s==nS; xlabel(ax,'Time from onset (s)','FontSize',8); end
        if s==1 && r==nR
            legH=gobjects(nLev,1); legLab=cell(nLev,1);
            for iInt=1:nLev
                t=(iInt-1)/max(nLev-1,1);
                col=regionColors(r,:)*t+(regionColors(r,:)*0.25+[1 1 1]*0.75)*(1-t);
                legH(iInt)=plot(nan,nan,'-','Color',col,'LineWidth',2);
                if strcmp(intFields{s},'int'); legLab{iInt}=sprintf('%g dB',allIntsAll{s}(iInt));
                else; legLab{iInt}=sprintf('%.2f',allIntsAll{s}(iInt)); end
            end   % for iInt legend
            legend(ax,legH,legLab,'Location','best','FontSize',7.5,'Box','off');
        end
        box off; set(ax,'TickDir','out','FontSize',8);
    end   % for r
end   % for s (Figure 1)
drawnow;
yMin=Inf; yMax=-Inf;
for s=1:nS; for r=1:nR; yl=ylim(axT(s,r)); yMin=min(yMin,yl(1)); yMax=max(yMax,yl(2)); end; end
for s=1:nS; for r=1:nR; ylim(axT(s,r),[yMin yMax]); end; end
annotation('textbox',[0.01 0 0.99 0.03], ...
    'String',['Mean \pm SEM.  Baseline subtracted, z-scored per trial.  ' ...
              'Light\x2013dark = low\x2013high intensity.  X: \minus0.5 to +1 s.'], ...
    'EdgeColor','none','FontSize',8,'Color',[0.4 0.4 0.4],'HorizontalAlignment','center');

%% =========================================================================
%  Significance (one-sided Wilcoxon vs zero, FDR within stim type)
%% =========================================================================
sigAch  = cell(nS,nR);
pAdjAll = cell(nS,nR);
for s = 1:nS
    nLev  = numel(allIntsAll{s});
    % Initialise
    for r = 1:nR
        sigAch{s,r}  = false(nLev,1);
        pAdjAll{s,r} = nan(nLev,1);
    end
    % Collect raw p-values in order: r1/iInt1, r1/iInt2, ... r2/iInt1, ...
    pRaws = nan(nR*nLev,1);
    k = 0;
    for r = 1:nR
        for iInt = 1:nLev
            k = k+1;
            v = aucAch{s,r}{iInt}; v = v(isfinite(v));
            if numel(v)<3; pRaws(k)=1; else; pRaws(k)=signrank(v,0,'tail','right'); end
        end
    end
    % FDR within stim type
    nT        = numel(pRaws);
    [ps, idx] = sort(pRaws(:));
    adj       = ps.*nT./(1:nT)';
    adj       = flipud(cummin(flipud(adj)));
    adj       = min(adj,1);
    pAdj      = nan(nT,1); pAdj(idx) = adj;
    % Assign back
    k = 0;
    for r = 1:nR
        for iInt = 1:nLev
            k = k+1;
            sigAch{s,r}(iInt)  = pAdj(k) < 0.05;
            pAdjAll{s,r}(iInt) = pAdj(k);
        end
    end
end

%% =========================================================================
%  Figure 2: AUC errorbar panels (nS rows x 2 cols: ACh+405nm | GFP)
%% =========================================================================
% ---- Sizing for direct export (no manual adjustment needed) -------------
% Figure is sized in centimetres so exported dimensions are predictable.
% Target: ~2-column supplementary figure, legible at print size.
FS_TITLE  = 15;   % sgtitle
FS_PANEL  = 13;   % panel titles
FS_LABEL  = 13;   % axis labels
FS_TICK   = 11;   % tick labels
FS_LEG    = 10;   % legend
FS_STAR   = 13;   % significance asterisks
LW_LINE   = 2.0;  % ACh errorbar lines
LW_REF    = 1.4;  % 405nm / GFP dashed lines
MS_MAIN   = 7;    % ACh markers
MS_REF    = 5;    % reference markers
AX_LW     = 1.0;  % axis line width

figure('Units','centimeters','Position',[2 2 17 22], ...
       'Color','w','PaperPositionMode','auto');
sgtitle(sprintf('ACh AUC (0-%g s) per Intensity', aucWinSec), ...
    'FontSize',FS_TITLE,'FontWeight','bold');
axAUC = gobjects(nS,2);

for s = 1:nS
    nLev  = numel(allIntsAll{s});
    xVals = 1:nLev;
    if strcmp(intFields{s},'int')
        xLabs = arrayfun(@(v) sprintf('%g',v), allIntsAll{s},'UniformOutput',false);
        xUnit = 'Intensity (dB)';
    else
        xLabs = arrayfun(@(v) sprintf('%.2f',v), allIntsAll{s},'UniformOutput',false);
        xUnit = 'Contrast';
    end

    % --- ACh + 405nm ---
    ax = subplot(nS,2,(s-1)*2+1);
    axAUC(s,1) = ax; hold(ax,'on');
    yTopAch = -Inf;
    for r = 1:nR
        muV=nan(nLev,1); seV=nan(nLev,1);
        muR=nan(nLev,1); seR=nan(nLev,1);
        for iInt=1:nLev
            v=aucAch{s,r}{iInt}; v=v(isfinite(v));
            if numel(v)>=2
                muV(iInt)=mean(v,'omitnan');
                seV(iInt)=std(v,0,'omitnan')/sqrt(numel(v));
                yTopAch = max(yTopAch, muV(iInt)+seV(iInt));
            end
            v=aucRef{s,r}{iInt}; v=v(isfinite(v));
            if numel(v)>=2; muR(iInt)=mean(v,'omitnan'); seR(iInt)=std(v,0,'omitnan')/sqrt(numel(v)); end
        end
        c=regionColors(r,:);
        errorbar(ax,xVals,muV,seV,'-o','Color',c,'MarkerFaceColor',c, ...
            'MarkerSize',MS_MAIN,'LineWidth',LW_LINE,'CapSize',5);
        plot(ax,xVals,muR,'--o','Color',c,'MarkerFaceColor','w', ...
            'MarkerSize',MS_REF,'LineWidth',LW_REF);
    end
    yline(ax,0,'-','Color',[0.8 0.8 0.8],'LineWidth',AX_LW);
    set(ax,'XTick',xVals,'XTickLabel',xLabs,'TickDir','out', ...
        'FontSize',FS_TICK,'LineWidth',AX_LW);
    xlim(ax,[0.5 nLev+0.5]);
    if s==1; title(ax,'ACh (solid) + 405nm (dashed)','FontSize',FS_PANEL,'FontWeight','bold'); end
    ylabel(ax,{stimLabels{s}; sprintf('AUC 0-%gs (z\\cdots)', aucWinSec)}, ...
        'FontSize',FS_LABEL,'FontWeight','bold');
    if s==nS; xlabel(ax,xUnit,'FontSize',FS_LABEL); end

    % Coloured asterisks above errorbars — one row per region, stacked upward
    drawnow; yL=ylim(ax);
    yStep = 0.07*range(yL);
    % Place legend in bottom-right to avoid asterisk conflict in row 1
    if s==1
        legH=gobjects(nR*2,1); legLab=cell(nR*2,1);
        for r=1:nR
            legH(r)=plot(nan,nan,'-o','Color',regionColors(r,:),'MarkerFaceColor',regionColors(r,:),'MarkerSize',MS_REF,'LineWidth',LW_LINE);
            legLab{r}=regions{r};
            legH(nR+r)=plot(nan,nan,'--o','Color',regionColors(r,:),'MarkerFaceColor','w','MarkerSize',MS_REF,'LineWidth',LW_REF);
            legLab{nR+r}=[regions{r} ' 405nm'];
        end
        legend(ax,legH,legLab,'Location','southeast','FontSize',FS_LEG,'Box','off');
    end
    % Asterisks above the top of the errorbars
    yAsterBase = yTopAch + 0.08*range(yL);
    for r=1:nR
        for iInt=1:nLev
            p = pAdjAll{s,r}(iInt);
            if p < 0.001;     stars = '***';
            elseif p < 0.01;  stars = '**';
            elseif p < 0.05;  stars = '*';
            else;              stars = '';
            end
            if ~isempty(stars)
                text(ax, iInt, yAsterBase+(r-1)*yStep, stars, ...
                    'Color',regionColors(r,:),'FontSize',FS_STAR,'FontWeight','bold', ...
                    'HorizontalAlignment','center','VerticalAlignment','bottom');
            end
        end
    end
    % Expand ylim to show asterisks
    ylim(ax,[yL(1), max(yL(2), yAsterBase+nR*yStep)]);
    box off;

    % --- GFP: significance ---
    nLev = numel(allIntsAll{s});
    pRawsG = nan(nR*nLev,1); k=0;
    for r=1:nR
        for iInt=1:nLev
            k=k+1;
            v=aucGfp{s,r}{iInt}; v=v(isfinite(v));
            if numel(v)<3; pRawsG(k)=1; else; pRawsG(k)=signrank(v,0,'tail','right'); end
        end
    end
    nT=[]; nT=numel(pRawsG);
    [psG,idxG]=sort(pRawsG(:));
    adjG=psG.*nT./(1:nT)'; adjG=flipud(cummin(flipud(adjG))); adjG=min(adjG,1);
    pAdjG=nan(nT,1); pAdjG(idxG)=adjG;
    sigGfp=cell(nR,1); k=0;
    for r=1:nR
        sigGfp{r}=false(nLev,1);
        for iInt=1:nLev
            k=k+1; sigGfp{r}(iInt)=pAdjG(k)<0.05;
        end
    end

    % --- GFP panel ---
    ax=subplot(nS,2,(s-1)*2+2);
    axAUC(s,2)=ax; hold(ax,'on');
    yTopGfp=-Inf;
    for r=1:nR
        muG=nan(nLev,1); seG=nan(nLev,1);
        for iInt=1:nLev
            v=aucGfp{s,r}{iInt}; v=v(isfinite(v));
            if numel(v)>=2
                muG(iInt)=mean(v,'omitnan');
                seG(iInt)=std(v,0,'omitnan')/sqrt(numel(v));
                yTopGfp=max(yTopGfp,muG(iInt)+seG(iInt));
            end
        end
        errorbar(ax,xVals,muG,seG,'--o','Color',regionColors(r,:),'MarkerFaceColor',regionColors(r,:), ...
            'MarkerSize',MS_MAIN,'LineWidth',LW_LINE,'CapSize',5);
    end
    yline(ax,0,'-','Color',[0.8 0.8 0.8],'LineWidth',AX_LW);
    set(ax,'XTick',xVals,'XTickLabel',xLabs,'TickDir','out', ...
        'FontSize',FS_TICK,'LineWidth',AX_LW);
    xlim(ax,[0.5 nLev+0.5]);
    if s==1; title(ax,'GFP control','FontSize',FS_PANEL,'FontWeight','bold'); end
    if s==nS; xlabel(ax,xUnit,'FontSize',FS_LABEL); end
    % Asterisks above errorbars
    drawnow; yLG=ylim(ax);
    yStepG=0.07*range(yLG);
    yAsterBaseG=yTopGfp+0.08*range(yLG);
    for r=1:nR
        for iInt=1:nLev
            pVal = pAdjG((r-1)*nLev+iInt);
            if pVal < 0.001;     starsG = '***';
            elseif pVal < 0.01;  starsG = '**';
            elseif pVal < 0.05;  starsG = '*';
            else;                 starsG = '';
            end
            if ~isempty(starsG)
                text(ax,iInt,yAsterBaseG+(r-1)*yStepG,starsG, ...
                    'Color',regionColors(r,:),'FontSize',FS_STAR,'FontWeight','bold', ...
                    'HorizontalAlignment','center','VerticalAlignment','bottom');
            end
        end
    end
    ylim(ax,[yLG(1), max(yLG(2), yAsterBaseG+nR*yStepG)]);
    box off;
end
for s=1:nS; linkaxes([axAUC(s,1) axAUC(s,2)],'y'); end
% Apply single shared y-axis across all errorbar panels
drawnow;
yMinAll=Inf; yMaxAll=-Inf;
for s=1:nS; for c=1:2
    yl=ylim(axAUC(s,c));
    yMinAll=min(yMinAll,yl(1)); yMaxAll=max(yMaxAll,yl(2));
end; end
for s=1:nS; for c=1:2
    ylim(axAUC(s,c),[yMinAll yMaxAll]);
end; end
annotation('textbox',[0.01 0 0.99 0.05], ...
    'String',['Mean \pm SEM across mice.  ' ...
              'Asterisks: * p<0.05, ** p<0.01, *** p<0.001 ' ...
              '(one-sided Wilcoxon vs zero, FDR-corrected within stimulus type); ' ...
              'colour indicates brain region, absent = n.s.'], ...
    'EdgeColor','none','FontSize',FS_LEG,'Color',[0.35 0.35 0.35], ...
    'HorizontalAlignment','center');

%% =========================================================================
%  CSV export
%% =========================================================================
tableRows={};
ctrlRowsFig2={};   % 405 nm and GFP rows, FDR-corrected below
for s=1:nS
    nLev=numel(allIntsAll{s}); stimRows={};
    for r=1:nR
        for iInt=1:nLev
            v=aucAch{s,r}{iInt}; v=v(isfinite(v));
            if numel(v)<3; continue; end
            if strcmp(intFields{s},'int'); intStr=sprintf('%g dB',allIntsAll{s}(iInt));
            else; intStr=sprintf('%.2f contrast',allIntsAll{s}(iInt)); end
            stimRows{end+1}={'ACh',stimLabels{s},regions{r},intStr,numel(v), ...
                mean(v,'omitnan'),std(v,0,'omitnan')/sqrt(numel(v)), ...
                signrank(v,0,'tail','right'),pAdjAll{s,r}(iInt),sigAch{s,r}(iInt)}; %#ok

            % --- Controls: same test, same tail, same FDR family structure --
            % 405 nm and GFP were previously tested only to place asterisks on
            % the figure (GFP) or not at all (405 nm); neither reached the CSV.
            % The figure and the stats table must agree, so all three cohorts
            % are now written with an identical one-sided Wilcoxon signed-rank
            % against zero. FDR is applied within cohort x stimulus type below.
            vR = aucRef{s,r}{iInt}; vR = vR(isfinite(vR));
            if numel(vR) >= 3
                ctrlRowsFig2{end+1} = {'405nm',stimLabels{s},regions{r},intStr, ...
                    numel(vR), mean(vR,'omitnan'), ...
                    std(vR,0,'omitnan')/sqrt(numel(vR)), ...
                    signrank(vR,0,'tail','right')}; %#ok
            end
            vG = aucGfp{s,r}{iInt}; vG = vG(isfinite(vG));
            if numel(vG) >= 3
                ctrlRowsFig2{end+1} = {'GFP',stimLabels{s},regions{r},intStr, ...
                    numel(vG), mean(vG,'omitnan'), ...
                    std(vG,0,'omitnan')/sqrt(numel(vG)), ...
                    signrank(vG,0,'tail','right')}; %#ok
            end
        end   % for iInt
    end   % for r
    for k=1:numel(stimRows); tableRows{end+1}=stimRows{k}; end %#ok
end   % for s (CSV)
% ---- FDR for the control cohorts -------------------------------------
% Same family structure as ACh: corrected within cohort x stimulus type,
% so the 9 region x intensity tests of one cohort/stimulus form one family.
if ~isempty(ctrlRowsFig2)
    cohorts   = unique(cellfun(@(x) x{1}, ctrlRowsFig2, 'UniformOutput', false));
    stimsSeen = unique(cellfun(@(x) x{2}, ctrlRowsFig2, 'UniformOutput', false));
    pAdjCtrlF2 = nan(numel(ctrlRowsFig2),1);
    famSize    = nan(numel(ctrlRowsFig2),1);
    for ci = 1:numel(cohorts)
        for si = 1:numel(stimsSeen)
            sel = find(cellfun(@(x) strcmp(x{1},cohorts{ci}) && ...
                                     strcmp(x{2},stimsSeen{si}), ctrlRowsFig2));
            if isempty(sel); continue; end
            pf = cellfun(@(k) ctrlRowsFig2{k}{8}, num2cell(sel));
            nf = numel(pf);
            [ps,ix] = sort(pf(:));
            adj = ps .* nf ./ (1:nf)';
            adj = flipud(cummin(flipud(adj)));
            adj = min(adj,1);
            tmp = nan(nf,1); tmp(ix) = adj;
            pAdjCtrlF2(sel) = tmp;
            famSize(sel)    = nf;
        end
    end
    for k = 1:numel(ctrlRowsFig2)
        rw = ctrlRowsFig2{k};
        tableRows{end+1} = {rw{1}, rw{2}, rw{3}, rw{4}, rw{5}, rw{6}, rw{7}, ...
            rw{8}, pAdjCtrlF2(k), pAdjCtrlF2(k) < 0.05}; %#ok
    end
end

headers={'Cohort','Test','Test_tail','Comparison','MC_correction','Correction_family', ...
         'Stimulus','Region','Intensity','auc_win_end_s','n_mice','Mean_AUC','SEM_AUC', ...
         'p_raw','p_FDR','Significant_FDR05'};
csvName = sprintf('ach_responses_auc_stats%s.csv', winTag);
fid=fopen(fullfile(pwd,csvName),'w');
fprintf(fid,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n',headers{:});
for k=1:numel(tableRows)
    row=tableRows{k};
    % Test identity written on every row so the stats spreadsheet is
    % self-describing and does not depend on remembering which function
    % produced which file.
    % Cohort now comes from the row (ACh, 405nm or GFP) — all three are
    % tested identically, so the table matches the asterisks on the figure.
    fprintf(fid,['%s,Wilcoxon signed-rank,one-sided (greater),vs zero,' ...
                 'Benjamini-Hochberg FDR,within cohort x stimulus type,' ...
                 '%s,%s,%s,%.3f,%d,%.4f,%.4f,%.4f,%.4f,%d\n'], ...
        row{1},row{2},row{3},row{4},aucWinSec,row{5},row{6},row{7}, ...
        row{8},row{9},row{10});
end
fclose(fid);
fprintf('Saved: %s  (AUC window 0-%.3f s)\n', csvName, aucWinSec);
