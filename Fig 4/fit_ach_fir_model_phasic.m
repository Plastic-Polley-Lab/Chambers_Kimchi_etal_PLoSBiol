function fit_ach_fir_model_phasic(ach_stim, fme_stim, pupil_stim, all_fiber_names_ach, logScale)
% FIT_ACH_FIR_MODEL_PHASIC
%
% Usage:
%   fit_ach_fir_model_phasic(ach_stim, fme_stim, pupil_stim, all_fiber_names_ach)
%   fit_ach_fir_model_phasic(ach_stim, fme_stim, pupil_stim, all_fiber_names_ach, logScale)
%
% Inputs:
%   ach_stim            — trial-structured ACh data
%   fme_stim            — trial-structured FME data
%   pupil_stim          — trial-structured pupil data
%   all_fiber_names_ach — fiber region names per mouse, cell{m}
%   logScale            — (optional, default false) if true, Plot 2 variance
%                         explained bars are shown on a log10 y-axis
%
% Fits a FIR GLM predicting the phasic ACh response (0 to +1 s post-stimulus)
% using stimulus, FME, and pupil predictors with lags spanning the full
% trial (-2 to +6 s), allowing detection of long-latency arousal-ACh
% relationships that gate the phasic response.
%
% Uses ridge regression with k-fold cross-validation.
% Plots: (1) FIR kernels, (2) variance explained per stimulus type,
%        (3) variance explained collapsed across stimuli with stats,
%        (4) observed vs predicted phasic traces.

%% -------------------------------------------------------------------------
%  Parse optional argument
% -------------------------------------------------------------------------
if nargin < 5 || isempty(logScale)
    logScale = false;
end

%% -------------------------------------------------------------------------
%  Parameters
% -------------------------------------------------------------------------
nMice      = numel(ach_stim);
regions    = {'AC', 'VC', 'PFC'};
stimTypes  = {'AudRipple', 'Tones', 'VisGrating'};
stimLabels = {'Aud. Ripple', 'Tones', 'Vis. Grating'};

fsTarget      = 25;
trialDur      = 8;
nTarget       = round(fsTarget * trialDur);   % full trial length at target fs
tDS           = linspace(-2, 6, nTarget);

% Phasic response window
phasicStart   = 0;
phasicEnd     = 1;
phasicIdx     = tDS >= phasicStart & tDS <= phasicEnd;
nPhasic       = sum(phasicIdx);
tPhasic       = tDS(phasicIdx);
onsetInPhasic = find(tPhasic >= phasicStart, 1, 'first');

% FIR lag parameters
firWinCausal  = 6;
firWinAcausal = 2;
nTapsCausal   = round(fsTarget * firWinCausal);
nTapsAcausal  = round(fsTarget * firWinAcausal);
% Stimulus taps: a delta at onset can only be observed at lags inside the
% phasic window, so taps beyond nPhasic are structurally zero. They add
% nothing to the fit but DO count in size(X,2), which sets the ridge penalty
% (lambda = 1e-4*trace(X'X)/size(X,2)) — 125 dead columns made lambda ~1.3x
% smaller than intended. Tested against the previous nTapsCausal setting:
% R2 values shift only slightly and the regional orderings are unchanged.
nTapsStim     = nPhasic;
nTapsPred     = nTapsCausal + nTapsAcausal;

tKernelStim   = (0:nTapsStim-1) / fsTarget;
tKernelPred   = (-nTapsAcausal:nTapsCausal-1) / fsTarget;

% Native FME/pupil parameters
nNative  = 241;
tNative  = linspace(-2, 6, nNative);
blIdx    = tNative >= -2 & tNative < 0;

kFolds   = 5;

% Seed the random stream. Fold assignment below uses randperm, so without a
% seed the cross-validated R^2 and the predicted traces in Plot 4 differ on
% every run. (Unique R^2 is unaffected — it comes from a full-data fit — but
% the CV numbers should be reproducible.)
rng(42);

% Colors
predColors = [
    0.00, 0.45, 0.70;   % blue        → stimulus
    0.90, 0.62, 0.00;   % orange      → FME
    0.80, 0.47, 0.65;   % purple/pink → pupil
];
predLabels   = {'Stimulus', 'FME', 'Pupil'};
regionColors = [
    0.00, 0.45, 0.70;
    0.90, 0.62, 0.00;
    0.00, 0.62, 0.45;
];

%% -------------------------------------------------------------------------
%  Helper functions
% -------------------------------------------------------------------------
    function X = build_fir_causal(pred, nTaps, nTime)
        X = zeros(nTime, nTaps);
        for lag = 1:nTaps
            X(lag:end, lag) = pred(1:end-lag+1);
        end
    end

    function X = build_fir_acausal(pred, nTaps, nTime)
        X = zeros(nTime, nTaps);
        for lag = 1:nTaps
            X(1:end-lag, lag) = pred(lag+1:end);
        end
    end

    function X = build_fir_full(pred, nTapsCaus, nTapsAcaus, nTime)
        % build_fir_acausal returns columns in order lag = +1, +2, ... samples
        % AHEAD, i.e. from -0.04 s to -2 s on the kernel axis. tKernelPred runs
        % -2 s -> 0, so the acausal block must be FLIPPED to match, otherwise
        % the pre-onset half of every FME/pupil kernel is plotted time-reversed.
        % (Column order does not affect the fit — ridge with lambda*I is
        % invariant to permutation — so R^2, predictions and CV are unchanged.
        % This only corrects the lag labelling of the plotted kernels.)
        Xa = build_fir_acausal(pred, nTapsAcaus, nTime);
        X  = [fliplr(Xa), build_fir_causal(pred, nTapsCaus, nTime)];
    end

    function X = build_stim_matrix(nTime, onsetSample, nTaps)
        stim              = zeros(nTime, 1);
        stim(onsetSample) = 1;
        X = build_fir_causal(stim, nTaps, nTime);
    end

    function b = ridge_fit(X, y)
        lambda = 1e-4 * trace(X'*X) / size(X,2);
        b      = (X'*X + lambda * eye(size(X,2))) \ (X'*y);
    end

    function r2 = compute_r2(y, yhat)
        ss_res = sum((y - yhat).^2);
        ss_tot = sum((y - mean(y)).^2);
        r2     = max(0, 1 - ss_res/ss_tot);
    end

    function mat = detrend_baseline(traces, tLocal, blMask)
        tBl = tLocal(blMask);
        mat = zeros(size(traces));
        for tr = 1:size(traces,1)
            p         = polyfit(tBl, traces(tr,blMask), 1);
            mat(tr,:) = traces(tr,:) - polyval(p, tLocal);
        end
    end

    function vec = interpolate_nans(vec, tLocal)
        nanMask = isnan(vec);
        if any(nanMask) && sum(~nanMask) > 1
            vec(nanMask) = interp1(tLocal(~nanMask), vec(~nanMask), ...
                tLocal(nanMask), 'linear', 'extrap');
        end
    end

    function mat = clean_traces(mat, tLocal)
        for tr = 1:size(mat,1)
            mat(tr,:) = interpolate_nans(mat(tr,:), tLocal);
        end
        mat(any(~isfinite(mat),2),:) = NaN;
    end

    function mat = zscore_rows(mat)
        mu        = mean(mat, 2, 'omitnan');
        sd        = std(mat,  0, 2, 'omitnan');
        sd(sd==0) = 1;
        mat       = (mat - mu) ./ sd;
    end

    function [pAdjOut, pRawOut] = region_pairs_stats(vals)
        pRawOut = nan(3,1);
        pr = [1 2; 1 3; 2 3];
        for q = 1:3
            a = vals{pr(q,1)}; b = vals{pr(q,2)};
            if numel(a) >= 3 && numel(b) >= 3
                pRawOut(q) = ranksum(a, b);
            end
        end
        pAdjOut = nan(3,1);
        fin = isfinite(pRawOut);
        if any(fin)
            pf = pRawOut(fin); nT = numel(pf);
            [ps,ix] = sort(pf(:));
            adj = ps .* nT ./ (1:nT)';
            adj = flipud(cummin(flipud(adj)));
            tmp = nan(nT,1); tmp(ix) = min(adj,1);
            pAdjOut(fin) = tmp;
        end
    end

    function draw_brackets(ax, xPos, vals, pAdjB, yTop)
        pr = [1 2; 1 3; 2 3];
        step = 0.085 * max(yTop, eps);
        lvl  = 0;
        for q = 1:3
            if ~isfinite(pAdjB(q)) || pAdjB(q) >= 0.05; continue; end
            lvl = lvl + 1;
            yB  = yTop + lvl*step;
            x1 = xPos(pr(q,1)); x2 = xPos(pr(q,2));
            plot(ax, [x1 x1 x2 x2], [yB-step*0.18 yB yB yB-step*0.18], ...
                '-', 'Color', [0.3 0.3 0.3], 'LineWidth', 0.9);
            if     pAdjB(q) < 0.001; st = '***';
            elseif pAdjB(q) < 0.01;  st = '**';
            else;                     st = '*';
            end
            text(ax, (x1+x2)/2, yB + step*0.08, st, 'HorizontalAlignment','center', ...
                'FontSize', 11, 'FontWeight','bold', 'Color',[0.2 0.2 0.2]);
        end
        yl = ylim(ax);
        ylim(ax, [min(0,yl(1)), max(yTop + (lvl+1.2)*step, yl(2))]);
    end

    function pAdj = fdr_bh(pVals)
        n    = numel(pVals);
        pAdj = ones(size(pVals));
        if n == 0; return; end
        [pSorted, sortIdx] = sort(pVals(:));
        ranks = (1:n)';
        bh    = pSorted .* n ./ ranks;
        for i = n-1:-1:1
            bh(i) = min(bh(i), bh(i+1));
        end
        bh = min(bh, 1);
        pAdj(sortIdx) = bh;
        pAdj = reshape(pAdj, size(pVals));
    end

    function s = format_pval(p)
        if p < 0.001;    s = 'p<0.001';
        elseif p < 0.05; s = sprintf('p=%.3f', p);
        else;             s = sprintf('p=%.2f',  p);
        end
    end

    function sig = safe_smooth(sig, ~)
        nSig = numel(sig);
        if nSig < 3; return; end
        kernLen = min(12, max(3, floor(nSig / 5)));
        if mod(kernLen, 2) == 0; kernLen = kernLen - 1; end
        if kernLen < 3; return; end
        kG  = gausswin(kernLen, 3);
        kG  = kG ./ sum(kG);
        pad = floor(kernLen / 2);
        sigPad = [sig(pad:-1:1), sig, sig(end:-1:end-pad+1)];
        sigConv = conv(sigPad, kG, 'valid');
        nOut = numel(sigConv);
        if nOut == nSig
            sig = sigConv;
        elseif nOut > nSig
            trim = floor((nOut - nSig) / 2);
            sig  = sigConv(trim+1 : trim+nSig);
        end
    end

%% -------------------------------------------------------------------------
%  Storage
% -------------------------------------------------------------------------
nR = numel(regions);
nS = numel(stimTypes);
nP = 3;

% Combined-stimulus fit: one model per mouse x region using ALL trials from
% all stimulus types pooled, so each animal contributes ONE value.
% (Pooling r2unique across stimulus types instead would enter each animal
% three times, which inflates n and understates the SEM.)
% Auditory-only variant of the combined fit (sensitivity check).
% Visual gratings change retinal illumination, so on those trials part of the
% pupil signal is a pupillary light reflex rather than an arousal measure.
% The light reflex is stimulus-locked, so most of that variance is shared with
% the stimulus regressor and is already excluded from pupil's UNIQUE R^2 — but
% this variant tests that directly by refitting on auditory trials only.
% Reported as a sensitivity check; the main panels use all stimulus types.
r2full_aud   = cell(nR, 1);
r2unique_aud = cell(nR, nP);

r2full_all   = cell(nR, 1);
r2unique_all = cell(nR, nP);
nTrialsAll   = cell(nR, 1);
for r = 1:nR
    r2full_all{r} = [];
    nTrialsAll{r} = [];
    r2full_aud{r} = [];
    for p = 1:nP
        r2unique_all{r,p} = [];
        r2unique_aud{r,p} = [];
    end
end

kernels       = cell(nR, nS, nP);   % PARTIAL: from the full model
% MARGINAL kernels: each predictor fitted ALONE (predictor + intercept).
% The full-model kernel is a partial coefficient — what that predictor adds
% beyond the others — and its sign can invert when predictors are collinear
% (a suppressor effect). FME and pupil are strongly coupled here, so the
% marginal fit is the check on whether a negative partial kernel reflects a
% real negative relationship or is an artefact of the other regressors.
kernelsSingle = cell(nR, nS, nP);
adjCorrAll    = cell(nR, nS);   % [nMice x nP] adjacent-tap correlation
r2unique      = cell(nR, nS, nP);
r2full        = cell(nR, nS);
cvR2          = cell(nR, nS);
obsTraces     = cell(nR, nS);
predTraces    = cell(nR, nS);
partialTraces = cell(nR, nS, nP);

% Single-trial observed and CV-predicted traces, pooled across animals.
% Needed because trial-AVERAGED plots destroy exactly the variance that FME
% and pupil explain: their contributions differ trial to trial and cancel on
% averaging, while the stimulus contribution is identical on every trial and
% survives intact. The single-trial view is what makes the R^2 ordering legible.
stObs  = cell(nR, nS);            % [nTrials x nPhasic] pooled
stObsM  = cell(nR, nS);           % {mouse} = [nTrials x nPhasic] per mouse
stFullM = cell(nR, nS);
stPartM = cell(nR, nS, nP);
stFull = cell(nR, nS);
stPart = cell(nR, nS, nP);
% Single-predictor R² for commonality analysis
% r2single{r,s,p} = [nMice x 1] R² from model with predictor p alone
r2single      = cell(nR, nS, nP);
% Pairwise-predictor R² for commonality analysis
% r2pair{r,s,1}=S+FME, {r,s,2}=S+Pupil, {r,s,3}=FME+Pupil
r2pair        = cell(nR, nS, 3);

% Per-mouse means averaged across stimulus types — used for collapsed stats
% so each mouse contributes exactly one observation regardless of nS.
% r2unique_mn{r,p}  = [nMice x 1]
% r2full_mn{r}      = [nMice x 1]
% r2single_mn{r,p}  = [nMice x 1]
% r2pair_mn{r,pp}   = [nMice x 1]
r2unique_mn   = cell(nR, nP);
r2full_mn     = cell(nR, 1);
r2single_mn   = cell(nR, nP);
r2pair_mn     = cell(nR, 3);
for r = 1:nR
    r2full_mn{r} = [];
    for p = 1:nP
        r2unique_mn{r,p}  = [];
        r2single_mn{r,p}  = [];
    end
    for pp = 1:3
        r2pair_mn{r,pp} = [];
    end
end

for r = 1:nR; for s = 1:nS
    r2full{r,s}     = [];
    cvR2{r,s}       = [];
    adjCorrAll{r,s} = [];
    obsTraces{r,s}  = [];
    predTraces{r,s} = [];
    stObs{r,s}      = [];
    stFull{r,s}     = [];
    stObsM{r,s}     = {};
    stFullM{r,s}    = {};
    for p = 1:nP
        stPart{r,s,p}  = [];
        stPartM{r,s,p} = {};
    end
    for p = 1:nP
        kernels{r,s,p}       = [];
        kernelsSingle{r,s,p} = [];
        r2unique{r,s,p}      = [];
        partialTraces{r,s,p} = [];
        r2single{r,s,p}      = [];
    end
    for pp = 1:3
        r2pair{r,s,pp} = [];
    end
end; end
% -------------------------------------------------------------------------
fprintf('Fitting phasic FIR models with %d-fold CV...\n', kFolds);
fprintf('Phasic window: %.1f to %.1f s  (%d samples at %d Hz)\n', ...
    phasicStart, phasicEnd, nPhasic, fsTarget);
fprintf('FIR lags: -%d to +%d s (%d taps for FME/pupil)\n', ...
    firWinAcausal, firWinCausal, nTapsPred);

for m = 1:nMice
    if isempty(ach_stim{m}) || isempty(fme_stim{m}) || isempty(pupil_stim{m})
        continue;
    end
    fprintf('  Mouse %d/%d\n', m, nMice);

    for s = 1:nS
        stimType = stimTypes{s};

        fmeRaw  = clean_traces(fme_stim{m}.(stimType).traces, tNative);
        fmeProc = detrend_baseline(fmeRaw, tNative, blIdx);
        fmeProc = zscore_rows(fmeProc);
        fmeDS   = zeros(size(fmeProc,1), nTarget);
        for tr = 1:size(fmeProc,1)
            if any(isnan(fmeProc(tr,:))); fmeDS(tr,:) = NaN;
            else; fmeDS(tr,:) = interp1(tNative, fmeProc(tr,:), tDS, 'linear','extrap');
            end
        end

        pupRaw  = clean_traces(pupil_stim{m}.(stimType).traces, tNative);
        pupProc = detrend_baseline(pupRaw, tNative, blIdx);
        pupProc = zscore_rows(pupProc);
        pupDS   = zeros(size(pupProc,1), nTarget);
        for tr = 1:size(pupProc,1)
            if any(isnan(pupProc(tr,:))); pupDS(tr,:) = NaN;
            else; pupDS(tr,:) = interp1(tNative, pupProc(tr,:), tDS, 'linear','extrap');
            end
        end

        nTrials = size(fmeDS,1);

        for r = 1:nR
            region   = regions{r};
            fiberIdx = find(strcmp(all_fiber_names_ach{m}, region));
            if isempty(fiberIdx); continue; end

            achRaw  = ach_stim{m}.(stimType).traces{fiberIdx};
            tOrig   = linspace(-2, 6, size(achRaw,2));
            achRaw  = clean_traces(achRaw, tOrig);
            achProc = detrend_baseline(achRaw, tOrig, tOrig >= -2 & tOrig < 0);
            achProc = zscore_rows(achProc);
            achDS   = zeros(nTrials, nTarget);
            for tr = 1:nTrials
                if any(isnan(achProc(tr,:))); achDS(tr,:) = NaN;
                else; achDS(tr,:) = interp1(tOrig, achProc(tr,:), tDS, 'linear','extrap');
                end
            end

            validTrials = ~any(isnan(achDS),2) & ...
                          ~any(isnan(fmeDS),2) & ...
                          ~any(isnan(pupDS),2);
            if sum(validTrials) < kFolds + 2
                fprintf('    Skipping %s %s mouse %d: only %d valid trials\n', ...
                    region, stimType, m, sum(validTrials));
                continue;
            end

            achV   = achDS(validTrials,:);
            fmeV   = fmeDS(validTrials,:);
            pupV   = pupDS(validTrials,:);
            nValid = sum(validTrials);

            % --- Provenance diagnostic -----------------------------------
            % Unique R^2 is a deterministic function of the trials that enter
            % the fit, so if these counts differ between runs the R^2 values
            % will differ too. Printing them makes the input to each fit
            % auditable and is the first thing to compare when values move
            % unexpectedly between versions of the pipeline.
            nAchBad = sum(any(isnan(achDS),2));
            nFmeBad = sum(any(isnan(fmeDS),2));
            nPupBad = sum(any(isnan(pupDS),2));
            fprintf(['    %-4s %-12s mouse %2d: %3d/%3d trials used ' ...
                     '(dropped: ACh %d, FME %d, pupil %d)\n'], ...
                region, stimType, m, nValid, size(achDS,1), ...
                nAchBad, nFmeBad, nPupBad);

            nTotal  = nValid * nPhasic;
            Y       = zeros(nTotal,1);
            X_s     = zeros(nTotal, nTapsStim);
            X_f     = zeros(nTotal, nTapsPred);
            X_p     = zeros(nTotal, nTapsPred);

            for tr = 1:nValid
                idx        = (tr-1)*nPhasic + (1:nPhasic);
                Y(idx)     = achV(tr, phasicIdx)';
                X_s(idx,:) = build_stim_matrix(nPhasic, onsetInPhasic, nTapsStim);
                X_f_full   = build_fir_full(fmeV(tr,:)', nTapsCausal, nTapsAcausal, nTarget);
                X_p_full   = build_fir_full(pupV(tr,:)', nTapsCausal, nTapsAcausal, nTarget);
                X_f(idx,:) = X_f_full(phasicIdx,:);
                X_p(idx,:) = X_p_full(phasicIdx,:);
            end

            % --- Identifiability diagnostic ---------------------------------
            % How correlated are neighbouring lag columns? Adjacent-tap
            % correlation governs whether individual tap weights can be
            % estimated: if lag k and lag k+1 are near-identical columns, the
            % fit cannot attribute weight between them and splits it
            % arbitrarily. This measures it from the actual design matrix
            % rather than assuming it from how smooth a signal looks.
            % Collected here and summarised after the loop.
            adjCorr = nan(1,nP);
            for pp = 1:nP
                switch pp
                    case 1; Xd = X_s;
                    case 2; Xd = X_f;
                    case 3; Xd = X_p;
                end
                cc = nan(1, size(Xd,2)-1);
                for kk = 1:size(Xd,2)-1
                    a = Xd(:,kk); b = Xd(:,kk+1);
                    if std(a) > 0 && std(b) > 0
                        cc(kk) = abs(corr(a, b));
                    end
                end
                adjCorr(pp) = mean(cc, 'omitnan');
            end
            adjCorrAll{r,s} = [adjCorrAll{r,s}; adjCorr];

            X_full     = [X_s, X_f, X_p, ones(nTotal,1)];
            X_no_stim  = [     X_f, X_p, ones(nTotal,1)];
            X_no_fme   = [X_s,      X_p, ones(nTotal,1)];
            X_no_pupil = [X_s, X_f,      ones(nTotal,1)];

            if any(~isfinite(X_full(:))); continue; end

            b_full = ridge_fit(X_full, Y);
            if any(~isfinite(b_full)); continue; end

            k_stim  = b_full(1:nTapsStim);
            k_fme   = b_full(nTapsStim+1           : nTapsStim+nTapsPred);
            k_pupil = b_full(nTapsStim+nTapsPred+1 : nTapsStim+2*nTapsPred);

            kernels{r,s,1} = [kernels{r,s,1}; k_stim'];
            kernels{r,s,2} = [kernels{r,s,2}; k_fme'];
            kernels{r,s,3} = [kernels{r,s,3}; k_pupil'];

            % --- Marginal kernels: each predictor fitted alone ---
            Xs_only = [X_s, ones(nTotal,1)];
            Xf_only = [X_f, ones(nTotal,1)];
            Xp_only = [X_p, ones(nTotal,1)];
            bS1 = ridge_fit(Xs_only, Y);
            bF1 = ridge_fit(Xf_only, Y);
            bP1 = ridge_fit(Xp_only, Y);
            if all(isfinite(bS1)) && all(isfinite(bF1)) && all(isfinite(bP1))
                kernelsSingle{r,s,1} = [kernelsSingle{r,s,1}; bS1(1:nTapsStim)'];
                kernelsSingle{r,s,2} = [kernelsSingle{r,s,2}; bF1(1:nTapsPred)'];
                kernelsSingle{r,s,3} = [kernelsSingle{r,s,3}; bP1(1:nTapsPred)'];
            end

            r2_full    = compute_r2(Y, X_full    * b_full);
            r2_no_stim = compute_r2(Y, X_no_stim  * ridge_fit(X_no_stim,  Y));
            r2_no_fme  = compute_r2(Y, X_no_fme   * ridge_fit(X_no_fme,   Y));
            r2_no_pup  = compute_r2(Y, X_no_pupil * ridge_fit(X_no_pupil, Y));

            r2full{r,s}     = [r2full{r,s};     r2_full];
            r2unique{r,s,1} = [r2unique{r,s,1}; r2_full - r2_no_stim];
            r2unique{r,s,2} = [r2unique{r,s,2}; r2_full - r2_no_fme];
            r2unique{r,s,3} = [r2unique{r,s,3}; r2_full - r2_no_pup];

            % ---------------------------------------------------------------
            %  Commonality analysis: single-predictor and pairwise models
            % ---------------------------------------------------------------
            X_only_s   = [X_s,             ones(nTotal,1)];
            X_only_f   = [     X_f,         ones(nTotal,1)];
            X_only_p   = [          X_p,    ones(nTotal,1)];
            X_sf       = [X_s, X_f,         ones(nTotal,1)];
            X_sp       = [X_s,       X_p,   ones(nTotal,1)];
            X_fp       = [     X_f,  X_p,   ones(nTotal,1)];

            r2s  = compute_r2(Y, X_only_s * ridge_fit(X_only_s, Y));
            r2f  = compute_r2(Y, X_only_f * ridge_fit(X_only_f, Y));
            r2p  = compute_r2(Y, X_only_p * ridge_fit(X_only_p, Y));
            r2sf = compute_r2(Y, X_sf     * ridge_fit(X_sf,     Y));
            r2sp = compute_r2(Y, X_sp     * ridge_fit(X_sp,     Y));
            r2fp = compute_r2(Y, X_fp     * ridge_fit(X_fp,     Y));

            r2single{r,s,1} = [r2single{r,s,1}; r2s];
            r2single{r,s,2} = [r2single{r,s,2}; r2f];
            r2single{r,s,3} = [r2single{r,s,3}; r2p];
            r2pair{r,s,1}   = [r2pair{r,s,1};   r2sf];
            r2pair{r,s,2}   = [r2pair{r,s,2};   r2sp];
            r2pair{r,s,3}   = [r2pair{r,s,3};   r2fp];

            foldAssign    = mod(randperm(nValid), kFolds) + 1;
            cvPred        = zeros(nValid, nPhasic);
            cvPartial     = zeros(nValid, nPhasic, nP);
            cvR2_folds    = zeros(kFolds,1);

            for fold = 1:kFolds
                testIdx  = find(foldAssign == fold);
                trainIdx = find(foldAssign ~= fold);

                nTr = numel(trainIdx) * nPhasic;
                nTe = numel(testIdx)  * nPhasic;

                Ytr   = zeros(nTr,1);
                Xtr_s = zeros(nTr, nTapsStim);
                Xtr_f = zeros(nTr, nTapsPred);
                Xtr_p = zeros(nTr, nTapsPred);

                Yte   = zeros(nTe,1);
                Xte_s = zeros(nTe, nTapsStim);
                Xte_f = zeros(nTe, nTapsPred);
                Xte_p = zeros(nTe, nTapsPred);

                for ii = 1:numel(trainIdx)
                    tr  = trainIdx(ii);
                    idx = (ii-1)*nPhasic + (1:nPhasic);
                    Ytr(idx)     = achV(tr,phasicIdx)';
                    Xtr_s(idx,:) = build_stim_matrix(nPhasic, onsetInPhasic, nTapsStim);
                    Xf_full      = build_fir_full(fmeV(tr,:)', nTapsCausal, nTapsAcausal, nTarget);
                    Xp_full      = build_fir_full(pupV(tr,:)', nTapsCausal, nTapsAcausal, nTarget);
                    Xtr_f(idx,:) = Xf_full(phasicIdx,:);
                    Xtr_p(idx,:) = Xp_full(phasicIdx,:);
                end

                for ii = 1:numel(testIdx)
                    tr  = testIdx(ii);
                    idx = (ii-1)*nPhasic + (1:nPhasic);
                    Yte(idx)     = achV(tr,phasicIdx)';
                    Xte_s(idx,:) = build_stim_matrix(nPhasic, onsetInPhasic, nTapsStim);
                    Xf_full      = build_fir_full(fmeV(tr,:)', nTapsCausal, nTapsAcausal, nTarget);
                    Xp_full      = build_fir_full(pupV(tr,:)', nTapsCausal, nTapsAcausal, nTarget);
                    Xte_f(idx,:) = Xf_full(phasicIdx,:);
                    Xte_p(idx,:) = Xp_full(phasicIdx,:);
                end

                X_train = [Xtr_s, Xtr_f, Xtr_p, ones(nTr,1)];
                X_test  = [Xte_s, Xte_f, Xte_p, ones(nTe,1)];

                if any(~isfinite(X_train(:))); continue; end
                b_tr = ridge_fit(X_train, Ytr);
                if any(~isfinite(b_tr)); continue; end

                yhat_te = X_test * b_tr;
                cvR2_folds(fold) = compute_r2(Yte, yhat_te);

                b_s = b_tr(1:nTapsStim);
                b_f = b_tr(nTapsStim+1           : nTapsStim+nTapsPred);
                b_p = b_tr(nTapsStim+nTapsPred+1 : nTapsStim+2*nTapsPred);
                b_i = b_tr(end);

                for ii = 1:numel(testIdx)
                    tr  = testIdx(ii);
                    idx = (ii-1)*nPhasic + (1:nPhasic);
                    cvPred(tr,:)      = yhat_te(idx);
                    cvPartial(tr,:,1) = Xte_s(idx,:)*b_s + b_i;
                    cvPartial(tr,:,2) = Xte_f(idx,:)*b_f + b_i;
                    cvPartial(tr,:,3) = Xte_p(idx,:)*b_p + b_i;
                end
            end

            mousePred  = mean(cvPred, 1, 'omitnan');
            mouseObs   = mean(achV(:,phasicIdx), 1, 'omitnan');
            mouseCV_R2 = mean(cvR2_folds(cvR2_folds > 0), 'omitnan');

            if any(~isfinite(mousePred)); continue; end

            predTraces{r,s}  = [predTraces{r,s};  mousePred];
            obsTraces{r,s}   = [obsTraces{r,s};   mouseObs];

            % single trials, pooled across animals
            okTr = all(isfinite(cvPred), 2);
            % per-mouse blocks: correlations must be computed WITHIN animal,
            % otherwise between-mouse differences in mean ACh and mean FME
            % create correlation that is not trial-to-trial.
            stObsM{r,s}{end+1}  = achV(okTr, phasicIdx);
            stFullM{r,s}{end+1} = cvPred(okTr,:);
            for pp = 1:nP
                stPartM{r,s,pp}{end+1} = cvPartial(okTr,:,pp);
            end
            stObs{r,s}  = [stObs{r,s};  achV(okTr, phasicIdx)];
            stFull{r,s} = [stFull{r,s}; cvPred(okTr,:)];
            for p = 1:nP
                stPart{r,s,p} = [stPart{r,s,p}; cvPartial(okTr,:,p)];
            end
            cvR2{r,s}        = [cvR2{r,s};        mouseCV_R2];
            for p = 1:nP
                partialTraces{r,s,p} = [partialTraces{r,s,p}; ...
                    mean(cvPartial(:,:,p),1,'omitnan')];
            end
        end
    end

    % -------------------------------------------------------------------
    %  After all s and r loops for this mouse: compute per-mouse means
    %  averaged across stimulus types for each region × predictor.
    %  This ensures each mouse contributes exactly one observation to
    %  cross-region and cross-predictor collapsed statistics.
    % -------------------------------------------------------------------
    for r = 1:nR
        % Collect this mouse's values across stim types for this region
        uMat = nan(nS, nP);   % [nS x nP] unique R²
        fMat = nan(nS, 1);    % [nS x 1]  full R²
        sMat = nan(nS, nP);   % [nS x nP] single-predictor R²
        pMat = nan(nS, 3);    % [nS x 3]  pairwise R²

        for s = 1:nS
            % Check if this mouse contributed to this r,s combination
            % by checking if the last entry was just added (same length)
            nU = numel(r2unique{r,s,1});
            nF = numel(r2full{r,s});
            if nF == 0; continue; end

            % The most recently appended value is the last element
            fMat(s) = r2full{r,s}(end);
            for p = 1:nP
                if numel(r2unique{r,s,p}) >= nU && nU > 0
                    uMat(s,p) = r2unique{r,s,p}(end);
                end
                if numel(r2single{r,s,p}) > 0
                    sMat(s,p) = r2single{r,s,p}(end);
                end
            end
            for pp = 1:3
                if numel(r2pair{r,s,pp}) > 0
                    pMat(s,pp) = r2pair{r,s,pp}(end);
                end
            end
        end

        % Only append if at least one stim type contributed
        if all(isnan(fMat)); continue; end

        r2full_mn{r} = [r2full_mn{r}; mean(fMat, 'omitnan')];
        for p = 1:nP
            r2unique_mn{r,p} = [r2unique_mn{r,p}; mean(uMat(:,p), 'omitnan')];
            r2single_mn{r,p} = [r2single_mn{r,p}; mean(sMat(:,p), 'omitnan')];
        end
        for pp = 1:3
            r2pair_mn{r,pp} = [r2pair_mn{r,pp}; mean(pMat(:,pp), 'omitnan')];
        end
    end
end

%% -------------------------------------------------------------------------
%  Combined-stimulus fit — one model per mouse x region, ALL trials pooled
%  ------------------------------------------------------------------------
%  Each animal contributes ONE value per region, rather than one per stimulus
%  type. This is a single fit on all the animal's trials, not an average of
%  three separate fits.
%
%  The stimulus regressor is a delta at onset, identical for every trial
%  regardless of type, so the combined stimulus term means "a stimulus
%  occurred" and does not distinguish modality. Stimulus-specific effects are
%  therefore still read from the per-stimulus fits above; the combined fit is
%  used for the full model, FME and pupil terms.
%
%  Note the trial gate applies to the POOLED count, so a mouse with too few
%  trials in any single stimulus type may still qualify here. n for the
%  combined panels can therefore exceed n for the per-stimulus panels.
%
%  Preprocessing is repeated here rather than cached from the loop above, so
%  this block is self-contained; it calls the same nested helpers, so the two
%  cannot silently diverge.
%% -------------------------------------------------------------------------
% Run the combined fit twice: once over all stimulus types (the reported
% panels), once over auditory types only (sensitivity check for the pupillary
% light reflex on visual trials).
audOnlyIdx = find(~strcmpi(stimTypes, 'VisGrating'));
fitPasses  = {1:nS, audOnlyIdx};
passLabel  = {'all stimulus types', 'auditory only'};

for pass = 1:2
fprintf('\nFitting combined-stimulus models (%s)...\n', passLabel{pass});

for m = 1:nMice
    if isempty(ach_stim{m}) || isempty(fme_stim{m}) || isempty(pupil_stim{m})
        continue;
    end

    for r = 1:nR
        region   = regions{r};
        fiberIdx = find(strcmp(all_fiber_names_ach{m}, region));
        if isempty(fiberIdx); continue; end

        Yall = []; Xs_all = []; Xf_all = []; Xp_all = [];

        for s = fitPasses{pass}
            stimType = stimTypes{s};
            if ~isfield(fme_stim{m}, stimType);   continue; end
            if ~isfield(pupil_stim{m}, stimType); continue; end
            if ~isfield(ach_stim{m}, stimType);   continue; end

            fmeProc = zscore_rows(detrend_baseline( ...
                clean_traces(fme_stim{m}.(stimType).traces, tNative), tNative, blIdx));
            pupProc = zscore_rows(detrend_baseline( ...
                clean_traces(pupil_stim{m}.(stimType).traces, tNative), tNative, blIdx));

            fmeDS = nan(size(fmeProc,1), nTarget);
            pupDS = nan(size(pupProc,1), nTarget);
            for tr = 1:size(fmeProc,1)
                if ~any(isnan(fmeProc(tr,:)))
                    fmeDS(tr,:) = interp1(tNative, fmeProc(tr,:), tDS, 'linear','extrap');
                end
            end
            for tr = 1:size(pupProc,1)
                if ~any(isnan(pupProc(tr,:)))
                    pupDS(tr,:) = interp1(tNative, pupProc(tr,:), tDS, 'linear','extrap');
                end
            end

            achRaw  = ach_stim{m}.(stimType).traces{fiberIdx};
            tOrig   = linspace(-2, 6, size(achRaw,2));
            achProc = zscore_rows(detrend_baseline( ...
                clean_traces(achRaw, tOrig), tOrig, tOrig >= -2 & tOrig < 0));
            achDS   = nan(size(achProc,1), nTarget);
            for tr = 1:size(achProc,1)
                if ~any(isnan(achProc(tr,:)))
                    achDS(tr,:) = interp1(tOrig, achProc(tr,:), tDS, 'linear','extrap');
                end
            end

            nTr = min([size(achDS,1), size(fmeDS,1), size(pupDS,1)]);
            if nTr < 1; continue; end
            achDS = achDS(1:nTr,:); fmeDS = fmeDS(1:nTr,:); pupDS = pupDS(1:nTr,:);

            vt = ~any(isnan(achDS),2) & ~any(isnan(fmeDS),2) & ~any(isnan(pupDS),2);
            if ~any(vt); continue; end

            achV = achDS(vt,:); fmeV = fmeDS(vt,:); pupV = pupDS(vt,:);
            for tr = 1:size(achV,1)
                Yall   = [Yall;   achV(tr, phasicIdx)'];                     %#ok
                Xs_all = [Xs_all; build_stim_matrix(nPhasic, onsetInPhasic, nTapsStim)]; %#ok
                Xff    = build_fir_full(fmeV(tr,:)', nTapsCausal, nTapsAcausal, nTarget);
                Xpf    = build_fir_full(pupV(tr,:)', nTapsCausal, nTapsAcausal, nTarget);
                Xf_all = [Xf_all; Xff(phasicIdx,:)];                          %#ok
                Xp_all = [Xp_all; Xpf(phasicIdx,:)];                          %#ok
            end
        end

        nTrPooled = numel(Yall) / nPhasic;
        if nTrPooled < kFolds + 2
            continue;
        end

        nTot   = numel(Yall);
        Xfull  = [Xs_all, Xf_all, Xp_all, ones(nTot,1)];
        XnoS   = [        Xf_all, Xp_all, ones(nTot,1)];
        XnoF   = [Xs_all,         Xp_all, ones(nTot,1)];
        XnoP   = [Xs_all, Xf_all,         ones(nTot,1)];
        if any(~isfinite(Xfull(:))); continue; end

        bF = ridge_fit(Xfull, Yall);
        if any(~isfinite(bF)); continue; end

        r2F   = compute_r2(Yall, Xfull * bF);
        r2nS  = compute_r2(Yall, XnoS * ridge_fit(XnoS, Yall));
        r2nF  = compute_r2(Yall, XnoF * ridge_fit(XnoF, Yall));
        r2nP  = compute_r2(Yall, XnoP * ridge_fit(XnoP, Yall));

        if pass == 1
            r2full_all{r}     = [r2full_all{r};     r2F];
            r2unique_all{r,1} = [r2unique_all{r,1}; r2F - r2nS];
            r2unique_all{r,2} = [r2unique_all{r,2}; r2F - r2nF];
            r2unique_all{r,3} = [r2unique_all{r,3}; r2F - r2nP];
            nTrialsAll{r}     = [nTrialsAll{r};     nTrPooled];
        else
            r2full_aud{r}     = [r2full_aud{r};     r2F];
            r2unique_aud{r,1} = [r2unique_aud{r,1}; r2F - r2nS];
            r2unique_aud{r,2} = [r2unique_aud{r,2}; r2F - r2nF];
            r2unique_aud{r,3} = [r2unique_aud{r,3}; r2F - r2nP];
        end
    end
end
end   % for pass

for r = 1:nR
    fprintf('  %-4s combined fit: %d mice (all stim), %d mice (auditory only)\n', ...
        regions{r}, numel(r2full_all{r}), numel(r2full_aud{r}));
end

%% -------------------------------------------------------------------------
%  Identifiability of individual tap weights
%  ------------------------------------------------------------------------
%  Measured, not assumed. The question is how correlated NEIGHBOURING lag
%  columns of the design matrix are. Where that correlation is high, lag k and
%  lag k+1 are near-identical predictors, the fit cannot attribute weight
%  between them, and individual tap weights are arbitrary even though their
%  sum is well determined. Where it is low, tap weights mean something.
%
%  This does NOT follow from how smooth a signal looks overall: what matters
%  is similarity between copies shifted by one 40 ms sample.
%% -------------------------------------------------------------------------
fprintf('\n==================================================================\n');
fprintf('Adjacent-tap correlation in the design matrix (mean |r|)\n');
fprintf('Higher = individual tap weights less identifiable\n');
fprintf('==================================================================\n');
fprintf('%-13s %-5s', 'Stimulus','Reg');
for p = 1:nP; fprintf('%12s', predLabels{p}); end
fprintf('\n');

for s = 1:nS
    for r = 1:nR
        A = adjCorrAll{r,s};
        if isempty(A); continue; end
        fprintf('%-13s %-5s', stimLabels{s}, regions{r});
        for p = 1:nP
            fprintf('%12.3f', mean(A(:,p), 'omitnan'));
        end
        fprintf('\n');
    end
end
fprintf('------------------------------------------------------------------\n');
allA = [];
for s = 1:nS; for r = 1:nR
    if ~isempty(adjCorrAll{r,s}); allA = [allA; adjCorrAll{r,s}]; end
end; end
if ~isempty(allA)
    fprintf('%-19s', 'POOLED');
    for p = 1:nP; fprintf('%12.3f', mean(allA(:,p), 'omitnan')); end
    fprintf('\n');
end
fprintf('==================================================================\n');
fprintf('Guide: |r| below ~0.5 -> tap weights carry information.\n');
fprintf('       |r| above ~0.9 -> tap weights are essentially arbitrary;\n');
fprintf('       interpret only sums (cumulative kernels) and unique R^2.\n');
fprintf('==================================================================\n');

fprintf('Done. Plotting...\n');

%% -------------------------------------------------------------------------
%  Smoothing kernel
% -------------------------------------------------------------------------
kSmoothSig = 2;
kGauss     = gausswin(6*kSmoothSig, 3);
kGauss     = kGauss ./ sum(kGauss);
zeroIdxPred = nTapsAcausal + 1;

%% =========================================================================
%  Plot 1: FIR kernels — RE-ENABLED (R1 comment 2c: "it would be
%  interesting to see the filters obtained for the different regressors").
%  Note the lag axes differ by predictor: stimulus is causal only (0 to +6 s),
%  FME and pupil span -2 to +6 s relative to the predicted ACh sample.
% =========================================================================
if false %#ok<UNRCH>   % raw tap weights — superseded by Plot 1b (cumulative)
figure('Units','normalized','Position',[0.02 0.05 0.96 0.85]);
sgtitle('FIR Kernels – Phasic ACh Response', ...
    'FontSize', 14, 'FontWeight', 'bold');

for s = 1:nS
    for r = 1:nR
        ax = subplot(nS, nR, (s-1)*nR + r);
        hold(ax, 'on');

        for p = 1:nP
            K = kernels{r,s,p};
            if isempty(K); continue; end
            K = K(all(isfinite(K),2),:);
            if size(K,1) < 2; continue; end

            mu  = mean(K,1,'omitnan');
            sem = std(K,0,1,'omitnan') ./ sqrt(size(K,1));
            if any(~isfinite(mu)) || any(~isfinite(sem)); continue; end

            if p == 1
                tK   = tKernelStim;
                mu   = safe_smooth(mu,  kGauss);
                sem  = safe_smooth(sem, kGauss);
            else
                tK   = tKernelPred;
                zIdx = zeroIdxPred;
                mu  = [safe_smooth(mu(1:zIdx-1),  kGauss), safe_smooth(mu(zIdx:end),  kGauss)];
                sem = [safe_smooth(sem(1:zIdx-1), kGauss), safe_smooth(sem(zIdx:end), kGauss)];
            end

            shadedErrorBar(tK, mu, sem, ...
                'lineProps',{'Color',predColors(p,:),'LineWidth',1.8}, ...
                'transparent',1,'patchSaturation',0.15);
        end

        xline(0,'--k','LineWidth',1.0,'Alpha',0.5);
        yline(0,'-','Color',[0.7 0.7 0.7],'LineWidth',0.8);
        xlabel('Lag (s)','FontSize',11);
        ylabel('Coefficient (z-scored)','FontSize',11);
        title(sprintf('%s – %s', regions{r}, stimLabels{s}), ...
            'FontSize',11,'FontWeight','bold','Color',regionColors(r,:));
        xlim([-firWinAcausal firWinCausal]);
        box off;
        set(ax,'TickDir','out','FontSize',11);

        if s == 1 && r == 1
            legH = gobjects(nP,1);
            for p = 1:nP
                legH(p) = plot(nan,nan,'-','Color',predColors(p,:),'LineWidth',2);
            end
            legend(legH, predLabels,'Location','best','FontSize',9,'Box','off');
        end
    end
end

end   % if false — Plot 1 (raw tap weights)

%% =========================================================================
%  Plot 2: Variance explained — 3 figures (one per stim type)
%  Each figure has 3 subplots (one per predictor), each with nR bars
%  (one per brain region) and an independent y-axis.
% =========================================================================
if false %#ok<UNRCH>
scaleStr = 'Linear';
if logScale; scaleStr = 'Log_{10}'; end

barWidth    = 0.55;
regionPairs = [1 2; 1 3; 2 3];
nRP         = size(regionPairs, 1);

for s = 1:nS
    figure('Units','normalized','Position',[0.05 0.15 0.75 0.55]);
    sgtitle(sprintf('Unique R² – %s  (%s scale)', stimLabels{s}, scaleStr), ...
        'FontSize', 13, 'FontWeight', 'bold');

    % Pre-compute unpaired Wilcoxon between regions per predictor
    pRaw2 = ones(nRP, nP);
    for p = 1:nP
        for rp = 1:nRP
            a = r2unique{regionPairs(rp,1), s, p};
            b = r2unique{regionPairs(rp,2), s, p};
            a = a(isfinite(a)); b = b(isfinite(b));
            if numel(a) >= 3 && numel(b) >= 3
                pRaw2(rp,p) = ranksum(a, b);
            end
        end
    end
    pAdj2 = reshape(fdr_bh(pRaw2(:)), nRP, nP);

    for p = 1:nP
        ax = subplot(1, nP, p);
        hold(ax, 'on');

        barTops  = nan(nR,1);
        fullVals = nan(nR,1);

        for r = 1:nR
            vals = r2unique{r,s,p};
            if isempty(vals); continue; end
            vals = vals(isfinite(vals));
            if isempty(vals); continue; end

            muR2  = mean(vals, 'omitnan');
            semR2 = std(vals, 0, 'omitnan') / sqrt(numel(vals));

            if logScale
                muPlot = log10(max(muR2, 1e-6));
                loPl   = muPlot - log10(max(muR2 - semR2, 1e-6));
                hiPl   = log10(max(muR2 + semR2, 1e-6)) - muPlot;
                bar(ax, r, muPlot, barWidth, ...
                    'FaceColor', regionColors(r,:), 'EdgeColor', 'none', 'FaceAlpha', 0.75);
                errorbar(ax, r, muPlot, loPl, hiPl, 'k', ...
                    'LineWidth', 1.5, 'CapSize', 5, 'LineStyle', 'none');
                barTops(r) = muPlot + hiPl;
            else
                bar(ax, r, muR2, barWidth, ...
                    'FaceColor', regionColors(r,:), 'EdgeColor', 'none', 'FaceAlpha', 0.75);
                errorbar(ax, r, muR2, semR2, 'k', ...
                    'LineWidth', 1.5, 'CapSize', 5, 'LineStyle', 'none');
                barTops(r) = muR2 + semR2;
            end

            % Store full model value for this region
            fv = r2full{r,s}; fv = fv(isfinite(fv));
            if ~isempty(fv)
                muFull = mean(fv,'omitnan');
                if logScale; muFull = log10(max(muFull,1e-6)); end
                fullVals(r) = muFull;
            end
        end

        % P-value annotations between regions (raw, uncorrected)
        drawnow;
        yRng   = range(ax.YLim);
        bracketY = max(barTops(isfinite(barTops))) + 0.08*yRng;
        for rp = 1:nRP
            pRaw  = pRaw2(rp,p);
            if pRaw < 0.001
                pStr = 'p<0.001';
            else
                pStr = sprintf('p=%.3f', pRaw);
            end
            bY = bracketY + (rp-1)*0.14*yRng;
            r1 = regionPairs(rp,1); r2idx = regionPairs(rp,2);
            plot(ax,[r1 r2idx],[bY bY],'-k','LineWidth',0.8);
            text(ax,(r1+r2idx)/2, bY+0.02*yRng, pStr, ...
                'HorizontalAlignment','center','FontSize',7.5,'Color','k');
        end

        % Set y-axis: for AudRipple, cut at bar level + bracket room;
        % full model lines are drawn with Clipping off so they remain visible
        % even if they exceed ylim, and are labelled.
        if s == 1   % AudRipple — cut y-axis
            yBarMax = max(barTops(isfinite(barTops)));
            yFullMax = max(fullVals(isfinite(fullVals)));
            % Set ylim to cover bars + brackets, not full model
            yTop = yBarMax + 0.45*yBarMax;   % 45% headroom above bars
            ylim(ax, [0, yTop]);
            drawnow;
            % Draw full model lines; if above ylim add text annotation
            for r = 1:nR
                if isnan(fullVals(r)); continue; end
                fv = fullVals(r);
                if fv <= yTop
                    yline(ax, fv, ':', 'Color', regionColors(r,:), ...
                        'LineWidth', 1.2, 'Alpha', 0.7, ...
                        'Label', sprintf('%s full', regions{r}), 'FontSize', 7);
                else
                    % Clip at top and annotate with arrow + text
                    plot(ax, [0.5 nR+0.5], [yTop*0.97 yTop*0.97], ':', ...
                        'Color', regionColors(r,:), 'LineWidth', 1.2);
                    text(ax, nR+0.55, yTop*0.97, sprintf('%s full R²=%.3f', regions{r}, fv), ...
                        'FontSize', 7, 'Color', regionColors(r,:), ...
                        'VerticalAlignment','middle','Clipping','off');
                end
            end
        else
            ylim(ax, [0 inf]);
            for r = 1:nR
                if isnan(fullVals(r)); continue; end
                yline(ax, fullVals(r), '--', 'Color', regionColors(r,:), ...
                    'LineWidth', 1.0, 'Alpha', 0.5, ...
                    'Label', sprintf('%s full R²', regions{r}), 'FontSize', 8);
            end
        end

        set(ax, 'XTick', 1:nR, 'XTickLabel', regions, ...
            'TickDir', 'out', 'FontSize', 11);
        xlim(ax, [0.4 nR + 0.6]);

        if logScale
            if p == 1; ylabel(ax, 'Unique R² (log_{10})', 'FontSize', 11); end
        else
            if p == 1; ylabel(ax, 'Unique R²', 'FontSize', 11); end
        end

        title(ax, predLabels{p}, 'FontSize', 11, 'FontWeight', 'bold', ...
            'Color', predColors(p,:));
        box off;
    end

    if s == 1
        annotation('textbox',[0.01 0 0.99 0.04],'String', ...
            ['Unpaired Wilcoxon between regions (FDR-corrected).  ' ...
             'Y-axis cut to bar level; full model R² shown as dotted lines (annotated if above axis).'], ...
            'EdgeColor','none','FontSize',8,'Color',[0.4 0.4 0.4],'HorizontalAlignment','center');
    else
        annotation('textbox',[0.01 0 0.99 0.04],'String', ...
            'Unpaired Wilcoxon between regions (FDR-corrected).  Dotted lines = full model R² per region.', ...
            'EdgeColor','none','FontSize',8,'Color',[0.4 0.4 0.4],'HorizontalAlignment','center');
    end
end

end   % if false — Plot 2
%% =========================================================================
%  Plot 1b: Cumulative FIR kernels — one predictor per row, region per column
%  ------------------------------------------------------------------------
%  Layout: nP rows (stimulus, FME, pupil) x nR columns (AC, VC, PFC), with the
%  three stimulus types overlaid within each panel.
%
%  Why this layout: normalisation now happens WITHIN a predictor, so curves
%  that share units are compared with each other. Grouping by stimulus type
%  instead forced three predictors with different units, different reliability
%  and different estimable lag ranges onto one axis, where pupil's noise
%  obscured the rest.
%
%  Normalisation is PER PANEL, by the largest absolute value across the three
%  stimulus curves — not per curve. This preserves their relative magnitudes,
%  so a panel shows which stimulus type drives the larger contribution rather
%  than forcing all three to equal height.
%
%  The STIMULUS row shows RAW kernel weights; the FME and pupil rows show
%  CUMULATIVE weights. Cumulating is a remedy for collinearity, which the
%  stimulus regressor does not suffer from (adjacent-tap |r| = 0.04 — its
%  design matrix is the identity), so its weights are already identifiable and
%  cumulating them only obscured the shape. The two rows answer different
%  questions: the stimulus row is a response profile, the behavioural rows are
%  accumulated influence across lags.
%
%  Interpreting the lag axis (after the acausal-block fix): NEGATIVE lag =
%  predictor AFTER the ACh sample (acausal taps); POSITIVE lag = predictor
%  BEFORE it (causal taps). The value at lag 0 is therefore the total acausal
%  weight, and the rise from 0 to the right edge is the total causal weight.
%  Those two block sums are the readable quantities; individual lag weights are
%  not identifiable (adjacent-tap |r| = 0.85 for FME, 0.98 for pupil).
%
%  Weights are referenced to the mean of the fitted window, not to baseline, so
%  early samples of a rising response carry negative weight. The stimulus
%  kernel's negative excursion is a property of that parameterisation and does
%  NOT indicate a decrease in ACh (cf. Fig 2).
%% =========================================================================
binWidthSec = [0 0 0.2];           % [stimulus FME pupil], seconds; 0 = no binning
binTapsVec  = max(1, round(binWidthSec * fsTarget));

lagMaxPlot  = 3.0;                 % FME/pupil: causal taps are zero beyond ~3 s
lagMaxStim  = (nPhasic - 1) / fsTarget;   % stimulus: only nPhasic taps are live

stimLines   = {'-','--',':'};      % one line style per stimulus type

figure('Units','centimeters','Position',[2 2 22 16], ...
       'Color','w','PaperPositionMode','auto');
sgtitle('Cumulative FIR kernels', 'FontSize', 13, 'FontWeight','bold');

for p = 1:nP
    for r = 1:nR
        ax = subplot(nP, nR, (p-1)*nR + r);
        hold(ax,'on');

        % ---- gather all stimulus types for this predictor x region ----
        muSet = cell(nS,1); semSet = cell(nS,1); tSet = cell(nS,1);
        for s = 1:nS
            K = kernels{r,s,p};
            if isempty(K); continue; end
            K = K(all(isfinite(K),2),:);
            if size(K,1) < 2; continue; end

            if p == 1
                tK = tKernelStim;  keep = tK <= lagMaxStim;
            else
                tK = tKernelPred;  keep = tK <= lagMaxPlot;
            end

            Kk = K(:, keep);
            tt = tK(keep);

            binTaps = binTapsVec(p);
            if binTaps > 1 && size(Kk,2) >= binTaps
                nB = floor(size(Kk,2) / binTaps);
                Kb = zeros(size(Kk,1), nB); tb = zeros(1, nB);
                for bIdx = 1:nB
                    cols       = (bIdx-1)*binTaps + (1:binTaps);
                    Kb(:,bIdx) = sum(Kk(:,cols), 2);
                    tb(bIdx)   = tt(cols(end));
                end
                Kk = Kb; tt = tb;
            end

            if p == 1
                Kc = Kk;              % stimulus: raw weights, already identifiable
            else
                Kc = cumsum(Kk, 2);   % FME/pupil: cumulate to escape collinearity
            end
            muSet{s}    = mean(Kc, 1, 'omitnan');
            semSet{s}   = std(Kc, 0, 1, 'omitnan') ./ sqrt(size(Kc,1));
            tSet{s}     = tt;
        end

        % ---- normalise by the largest |value| ACROSS stimulus types ----
        denom = 0;
        for s = 1:nS
            if ~isempty(muSet{s}); denom = max(denom, max(abs(muSet{s}))); end
        end
        if denom == 0 || ~isfinite(denom); denom = 1; end

        for s = 1:nS
            if isempty(muSet{s}); continue; end
            mu  = muSet{s}  ./ denom;
            sem = semSet{s} ./ denom;
            tt  = tSet{s};
            fill(ax, [tt, fliplr(tt)], [mu+sem, fliplr(mu-sem)], ...
                predColors(p,:), 'FaceAlpha', 0.10, 'EdgeColor','none');
            plot(ax, tt, mu, stimLines{s}, 'Color', predColors(p,:), ...
                'LineWidth', 1.8);
        end

        xline(ax, 0, '--k', 'Alpha', 0.45, 'LineWidth', 1.0);
        yline(ax, 0, '-', 'Color', [0.85 0.85 0.85], 'LineWidth', 0.8);
        if p == 1
            xlim(ax, [0 lagMaxStim]);
        else
            xlim(ax, [-firWinAcausal lagMaxPlot]);
        end
        ylim(ax, [-1.15 1.15]);

        if r == 1
            if p == 1; ylab = 'weight (norm.)'; else; ylab = 'cum. weight (norm.)'; end
            ylabel(ax, {predLabels{p}; ylab}, ...
                'FontSize', 10, 'FontWeight','bold', 'Color', predColors(p,:));
        end
        if p == nP; xlabel(ax, 'Lag (s)', 'FontSize', 10); end
        if p == 1
            title(ax, regions{r}, 'FontSize', 11, 'FontWeight','bold', ...
                'Color', regionColors(r,:));
        end
        box off; set(ax,'TickDir','out','FontSize',9);

        if p == 1 && r == nR
            lh = gobjects(nS,1);
            for s = 1:nS
                lh(s) = plot(ax, nan, nan, stimLines{s}, ...
                    'Color', [0.35 0.35 0.35], 'LineWidth', 1.8);
            end
            legend(ax, lh, stimLabels, 'FontSize', 8, 'Box','off', ...
                'Location','southeast');
        end
    end
end

annotation('textbox',[0.01 0 0.99 0.055], ...
    'String', sprintf([ ...
        'STIMULUS row: raw kernel weights (identifiable; adjacent-tap |r| = 0.04).  ' ...
        'FME and PUPIL rows: cumulative weights, because their adjacent taps are ' ...
        'collinear (|r| = 0.85 and 0.98).  Mean across mice (mean \\pm SEM), normalised PER PANEL ' ...
        'to the largest absolute value across stimulus types, so relative ' ...
        'magnitudes within a panel are preserved.  NEGATIVE lag = predictor ' ...
        'AFTER the ACh sample; POSITIVE lag = predictor BEFORE it.  The value at ' ...
        'lag 0 is the total acausal weight and the rise beyond 0 the total causal ' ...
        'weight; individual lag weights are not identifiable (adjacent-tap ' ...
        '|r| = 0.85 FME, 0.98 pupil), so pupil is binned at %g s.  Weights are ' ...
        'referenced to the fitted-window mean, not baseline, so early samples of ' ...
        'a rising response are negative — this does not indicate reduced ACh.'], ...
        binWidthSec(3)), ...
    'EdgeColor','none','FontSize',8.5,'Color',[0.4 0.4 0.4], ...
    'HorizontalAlignment','center');

%% =========================================================================
%  Plot 2c: FIGURE 4 PANEL — the four plots to report
%  ------------------------------------------------------------------------
%    1. Full model R^2      — all stimuli pooled, one point per animal
%    2. Pupil unique R^2    — all stimuli pooled, one point per animal
%    3. FME unique R^2      — all stimuli pooled, one point per animal
%    4. Stimulus unique R^2 — split by stimulus type (grouped bars)
%
%  Panels 1-3 use the combined fit, so each animal appears once. Panel 4 uses
%  the per-stimulus fits, because the combined stimulus regressor cannot
%  distinguish modality. n may therefore differ between panel 4 and panels
%  1-3, and is printed on each panel.
%
%  Stats: unpaired Wilcoxon rank-sum between regions (no mouse contributes to
%  all three), FDR-corrected within each panel. Brackets shown for surviving
%  comparisons; the full set goes to the CSV.
%% =========================================================================
figure('Units','centimeters','Position',[2 2 23 7.0], ...
       'Color','w','PaperPositionMode','auto');

panelData  = {r2full_all, r2unique_all(:,3), r2unique_all(:,2)};
panelNames = {'Full model', 'Pupil unique', 'FME unique'};
fig4Rows   = {};

barAlphaF = 0.45; dotSzF = 26; jitF = 0.11;

% ---- Panels 1-3: combined fit, one point per animal ----
for pnl = 1:3
    ax = subplot(1, 4, pnl);
    hold(ax,'on');

    D = panelData{pnl};
    vals = cell(nR,1);
    for r = 1:nR
        v = D{r}; vals{r} = v(isfinite(v));
    end

    yTop = -Inf;
    for r = 1:nR
        v = vals{r};
        if isempty(v); continue; end
        mu = mean(v,'omitnan'); se = std(v,0,'omitnan')/sqrt(numel(v));
        bar(ax, r, mu, 0.62, 'FaceColor', regionColors(r,:), ...
            'FaceAlpha', barAlphaF, 'EdgeColor', regionColors(r,:), 'LineWidth', 0.9);
        errorbar(ax, r, mu, se, 'Color',[0.25 0.25 0.25], ...
            'LineWidth',1.1,'CapSize',6,'LineStyle','none');
        jj = (rand(numel(v),1)-0.5)*2*jitF;
        scatter(ax, r+jj, v, dotSzF, regionColors(r,:), 'filled', ...
            'MarkerFaceAlpha',0.75,'MarkerEdgeColor','w','LineWidth',0.3);
        yTop = max([yTop; v(:); mu+se]);
    end

    [pAdjB, pRawB] = region_pairs_stats(vals);
    draw_brackets(ax, 1:nR, vals, pAdjB, yTop);

    pr = [1 2; 1 3; 2 3];
    for q = 1:3
        fig4Rows{end+1} = {panelNames{pnl}, 'all stimuli', ...
            regions{pr(q,1)}, regions{pr(q,2)}, ...
            numel(vals{pr(q,1)}), numel(vals{pr(q,2)}), pRawB(q), pAdjB(q)}; %#ok
    end

    set(ax,'XTick',1:nR,'XTickLabel',regions,'TickDir','out','FontSize',10,'LineWidth',0.9);
    xlim(ax,[0.4 nR+0.6]);
    ylabel(ax, 'R^2', 'FontSize', 11);
    title(ax, sprintf('%s  (n=%d/%d/%d)', panelNames{pnl}, ...
        numel(vals{1}), numel(vals{2}), numel(vals{3})), ...
        'FontSize', 10.5, 'FontWeight','bold');
    box off;
end

% ---- Panel 4: stimulus unique R^2, split by stimulus type ----
ax = subplot(1, 4, 4);
hold(ax,'on');
gw   = 0.26;                       % offset between region bars within a group
yTop4 = -Inf;
for s = 1:nS
    vals = cell(nR,1);
    for r = 1:nR
        v = r2unique{r,s,1}; vals{r} = v(isfinite(v));
    end
    xPos = s + (-1:1)*gw;
    for r = 1:nR
        v = vals{r};
        if isempty(v); continue; end
        mu = mean(v,'omitnan'); se = std(v,0,'omitnan')/sqrt(numel(v));
        bar(ax, xPos(r), mu, gw*0.85, 'FaceColor', regionColors(r,:), ...
            'FaceAlpha', barAlphaF, 'EdgeColor', regionColors(r,:), 'LineWidth', 0.8);
        errorbar(ax, xPos(r), mu, se, 'Color',[0.25 0.25 0.25], ...
            'LineWidth',1.0,'CapSize',4,'LineStyle','none');
        jj = (rand(numel(v),1)-0.5)*2*(gw*0.28);
        scatter(ax, xPos(r)+jj, v, dotSzF*0.8, regionColors(r,:), 'filled', ...
            'MarkerFaceAlpha',0.75,'MarkerEdgeColor','w','LineWidth',0.25);
        yTop4 = max([yTop4; v(:); mu+se]);
    end
    [pAdjB, pRawB] = region_pairs_stats(vals);
    draw_brackets(ax, xPos, vals, pAdjB, yTop4);
    pr = [1 2; 1 3; 2 3];
    for q = 1:3
        fig4Rows{end+1} = {'Stimulus unique', stimLabels{s}, ...
            regions{pr(q,1)}, regions{pr(q,2)}, ...
            numel(vals{pr(q,1)}), numel(vals{pr(q,2)}), pRawB(q), pAdjB(q)}; %#ok
    end
end
set(ax,'XTick',1:nS,'XTickLabel',stimLabels,'TickDir','out','FontSize',10,'LineWidth',0.9);
xlim(ax,[0.4 nS+0.6]);
ylabel(ax,'Unique R^2','FontSize',11);
title(ax,'Stimulus unique, by type','FontSize',10.5,'FontWeight','bold');
legH4 = gobjects(nR,1);
for r = 1:nR
    legH4(r) = bar(ax, nan, nan, 'FaceColor', regionColors(r,:), ...
        'FaceAlpha', barAlphaF, 'EdgeColor', regionColors(r,:));
end
legend(ax, legH4, regions, 'FontSize', 8.5, 'Box','off','Location','best');
box off;

annotation('textbox',[0.01 0 0.99 0.055], ...
    'String',['Bars: mean \pm SEM across mice; dots: individual mice.  ' ...
              'Panels 1-3 pool all stimulus types into a single fit per animal ' ...
              '(one point per animal); panel 4 uses the per-stimulus fits.  ' ...
              'Brackets: unpaired Wilcoxon rank-sum between regions, FDR-corrected ' ...
              'within panel (* p<0.05, ** p<0.01, *** p<0.001).'], ...
    'EdgeColor','none','FontSize',8.5,'Color',[0.4 0.4 0.4], ...
    'HorizontalAlignment','center');

% ---- CSV ----
if ~isempty(fig4Rows)
    csvF4 = fullfile(pwd, 'fig4_region_stats.csv');
    fidF4 = fopen(csvF4,'w');
    if fidF4 == -1
        warning('Could not write fig4_region_stats.csv');
    else
        fprintf(fidF4, 'measure,stimulus,region_1,region_2,n_1,n_2,p_raw_ranksum,p_FDR,sig_FDR05\n');
        for k = 1:numel(fig4Rows)
            rw = fig4Rows{k};
            fprintf(fidF4, '%s,%s,%s,%s,%d,%d,%.6f,%.6f,%d\n', ...
                rw{1}, rw{2}, rw{3}, rw{4}, rw{5}, rw{6}, rw{7}, rw{8}, rw{8} < 0.05);
        end
        fclose(fidF4);
        fprintf('Figure 4 region stats saved to: %s\n', csvF4);
    end
end

%% =========================================================================
%  Plot 2b: Region comparison bar plots — one figure per measure
%  ------------------------------------------------------------------------
%  Layout per figure: 1 row x nS panels (one per stimulus type). Each panel
%  has nR semi-transparent bars (one per brain region), mean +/- SEM, with
%  individual animals overlaid as jittered dots.
%
%  Significance: unpaired Wilcoxon rank-sum between regions (no mouse
%  contributes to all three), FDR-corrected across the 3 pairs within each
%  panel. Brackets are drawn only for comparisons that survive correction, to
%  keep the panels readable; the full set including non-significant
%  comparisons goes to the CSV.
%
%  The FULL MODEL R^2 gets its own figure here rather than being overlaid as
%  dotted reference lines on the unique-R^2 plots: full model R^2 is an order
%  of magnitude larger, so sharing an axis compressed the unique values into
%  the bottom of the plot.
%% =========================================================================
if false %#ok<UNRCH>
measureData  = cell(nP+1, 1);
measureNames = [predLabels, {'Full model'}];
for p = 1:nP
    % Build the nR x nS slice explicitly rather than with squeeze, which can
    % collapse a singleton dimension and break the {r,s} indexing below.
    slice = cell(nR, nS);
    for r = 1:nR
        for s = 1:nS
            slice{r,s} = r2unique{r,s,p};
        end
    end
    measureData{p} = slice;
end
measureData{nP+1} = r2full;                       % {r,s}

barAlpha   = 0.45;
dotSz      = 26;
dotAlpha   = 0.75;
jitterW    = 0.11;
regionPairsB = [1 2; 1 3; 2 3];
nRPB       = size(regionPairsB,1);

barStatRows = {};

for mi = 1:numel(measureData)
    M = measureData{mi};

    figure('Units','centimeters','Position',[2 2 19 7.5], ...
           'Color','w','PaperPositionMode','auto');
    sgtitle(sprintf('%s R^2', measureNames{mi}), ...
        'FontSize', 13, 'FontWeight','bold');

    axB = gobjects(nS,1);
    for s = 1:nS
        ax = subplot(1, nS, s);
        axB(s) = ax;
        hold(ax,'on');

        mu = nan(nR,1); se = nan(nR,1); vals = cell(nR,1);
        for r = 1:nR
            v = M{r,s};
            v = v(isfinite(v));
            vals{r} = v;
            if numel(v) >= 1
                mu(r) = mean(v,'omitnan');
                se(r) = std(v,0,'omitnan')/sqrt(max(numel(v),1));
            end
        end

        % bars + errorbars
        for r = 1:nR
            if ~isfinite(mu(r)); continue; end
            bar(ax, r, mu(r), 0.62, 'FaceColor', regionColors(r,:), ...
                'FaceAlpha', barAlpha, 'EdgeColor', regionColors(r,:), ...
                'LineWidth', 0.9);
            if isfinite(se(r))
                errorbar(ax, r, mu(r), se(r), 'Color', [0.25 0.25 0.25], ...
                    'LineWidth', 1.1, 'CapSize', 6, 'LineStyle','none');
            end
        end

        % individual animals
        for r = 1:nR
            v = vals{r};
            if isempty(v); continue; end
            jit = (rand(numel(v),1)-0.5) * 2 * jitterW;
            scatter(ax, r + jit, v, dotSz, regionColors(r,:), 'filled', ...
                'MarkerFaceAlpha', dotAlpha, 'MarkerEdgeColor','w', ...
                'LineWidth', 0.3);
        end

        % pairwise stats
        pRawB = nan(nRPB,1);
        for rp = 1:nRPB
            a = vals{regionPairsB(rp,1)};
            b = vals{regionPairsB(rp,2)};
            if numel(a) >= 3 && numel(b) >= 3
                pRawB(rp) = ranksum(a, b);
            end
        end
        fin  = isfinite(pRawB);
        pAdjB = nan(nRPB,1);
        if any(fin)
            pf = pRawB(fin); nT = numel(pf);
            [ps,ix] = sort(pf(:));
            adj = ps .* nT ./ (1:nT)';
            adj = flipud(cummin(flipud(adj)));
            tmp = nan(nT,1); tmp(ix) = min(adj,1);
            pAdjB(fin) = tmp;
        end

        for rp = 1:nRPB
            barStatRows{end+1} = {measureNames{mi}, stimLabels{s}, ...
                regions{regionPairsB(rp,1)}, regions{regionPairsB(rp,2)}, ...
                numel(vals{regionPairsB(rp,1)}), numel(vals{regionPairsB(rp,2)}), ...
                pRawB(rp), pAdjB(rp)}; %#ok
        end

        % significance brackets (only those surviving FDR)
        drawnow;
        yl = ylim(ax);
        allV = vertcat(vals{:});
        yTop = max([allV(:); mu+se], [], 'omitnan');
        if ~isfinite(yTop); yTop = yl(2); end
        step = 0.085 * max(yTop, eps);
        lvl  = 0;
        for rp = 1:nRPB
            if ~isfinite(pAdjB(rp)) || pAdjB(rp) >= 0.05; continue; end
            lvl = lvl + 1;
            yB  = yTop + lvl*step;
            r1 = regionPairsB(rp,1); r2b = regionPairsB(rp,2);
            plot(ax, [r1 r1 r2b r2b], [yB-step*0.18 yB yB yB-step*0.18], ...
                '-', 'Color', [0.3 0.3 0.3], 'LineWidth', 0.9);
            if     pAdjB(rp) < 0.001; st = '***';
            elseif pAdjB(rp) < 0.01;  st = '**';
            else;                      st = '*';
            end
            text(ax, (r1+r2b)/2, yB + step*0.10, st, ...
                'HorizontalAlignment','center', 'FontSize', 12, ...
                'FontWeight','bold', 'Color', [0.2 0.2 0.2]);
        end
        ylim(ax, [min(0, yl(1)), max(yTop + (lvl+1.2)*step, yl(2))]);

        set(ax,'XTick',1:nR,'XTickLabel',regions,'TickDir','out', ...
            'FontSize',10,'LineWidth',0.9);
        xlim(ax,[0.4 nR+0.6]);
        if s == 1
            ylabel(ax, sprintf('%s R^2', measureNames{mi}), 'FontSize', 11);
        end
        title(ax, stimLabels{s}, 'FontSize', 11, 'FontWeight','bold');
        box off;
    end

    annotation('textbox',[0.01 0 0.99 0.055], ...
        'String',['Bars: mean \pm SEM across mice; dots: individual mice.  ' ...
                  'Brackets: unpaired Wilcoxon rank-sum between regions, ' ...
                  'FDR-corrected across the 3 pairs within each panel ' ...
                  '(* p<0.05, ** p<0.01, *** p<0.001).  ' ...
                  'Non-significant comparisons are in the CSV.'], ...
        'EdgeColor','none','FontSize',8.5,'Color',[0.4 0.4 0.4], ...
        'HorizontalAlignment','center');
end

% ---- CSV for all region comparisons in these panels ----
if ~isempty(barStatRows)
    csvPathB = fullfile(pwd, 'fir_region_bars_stats.csv');
    fidB = fopen(csvPathB, 'w');
    if fidB == -1
        warning('Could not write fir_region_bars_stats.csv');
    else
        fprintf(fidB, ['measure,stimulus,region_1,region_2,n_1,n_2,' ...
                       'p_raw_ranksum,p_FDR,significant_FDR05\n']);
        for k = 1:numel(barStatRows)
            rw = barStatRows{k};
            fprintf(fidB, '%s,%s,%s,%s,%d,%d,%.6f,%.6f,%d\n', ...
                rw{1}, rw{2}, rw{3}, rw{4}, rw{5}, rw{6}, ...
                rw{7}, rw{8}, rw{8} < 0.05);
        end
        fclose(fidB);
        fprintf('Region comparison bar stats saved to: %s\n', csvPathB);
    end
end


end   % if false — Plot 2b
%% =========================================================================
%  Plot 3: Variance explained collapsed across stimuli with stats
% =========================================================================
if false %#ok<UNRCH>
r2collapsed = cell(nR, nP);
for r = 1:nR; for p = 1:nP
    pooled = [];
    for s = 1:nS
        vals = r2unique{r,s,p};
        pooled = [pooled; vals(isfinite(vals))];
    end
    r2collapsed{r,p} = pooled;
end; end

pRaw_t = ones(nR,nP);
pRaw_w = ones(nR,nP);
for r = 1:nR; for p = 1:nP
    vals = r2collapsed{r,p};
    if numel(vals) < 3 || all(vals==0); continue; end
    [~, pRaw_t(r,p)] = ttest(vals);
    pRaw_w(r,p)      = signrank(vals);
end; end

pAdj_t = reshape(fdr_bh(pRaw_t(:)), nR, nP);
pAdj_w = reshape(fdr_bh(pRaw_w(:)), nR, nP);

figure('Units','normalized','Position',[0.10 0.15 0.78 0.55]);
sgtitle(sprintf('Variance Explained Collapsed Across Stimuli – Phasic Window  (FDR-corrected)  [%s scale]', scaleStr), ...
    'FontSize', 13, 'FontWeight', 'bold');

for r = 1:nR
    ax = subplot(1, nR, r);
    hold(ax, 'on');

    yMaxAll = -Inf;

    for p = 1:nP
        vals = r2collapsed{r,p};
        if isempty(vals); continue; end

        muR2  = mean(vals,'omitnan');
        semR2 = std(vals,0,'omitnan') ./ sqrt(numel(vals));

        if logScale
            muPlot = log10(max(muR2, 1e-6));
            loPl   = muPlot - log10(max(muR2 - semR2, 1e-6));
            hiPl   = log10(max(muR2 + semR2, 1e-6)) - muPlot;
            bar(ax, p, muPlot, 0.6, 'FaceColor', predColors(p,:), 'EdgeColor', 'none', 'FaceAlpha', 0.75);
            errorbar(ax, p, muPlot, loPl, hiPl, 'k', 'LineWidth', 1.5, 'CapSize', 6, 'LineStyle', 'none');
            yMaxAll = max(yMaxAll, muPlot + hiPl);
        else
            bar(ax, p, muR2, 0.6, 'FaceColor', predColors(p,:), 'EdgeColor', 'none', 'FaceAlpha', 0.75);
            errorbar(ax, p, muR2, semR2, 'k', 'LineWidth', 1.5, 'CapSize', 6, 'LineStyle', 'none');
            yMaxAll = max(yMaxAll, muR2 + semR2);
        end
    end

    fullPooled = [];
    for s = 1:nS
        fv = r2full{r,s};
        fullPooled = [fullPooled; fv(isfinite(fv))];
    end
    if ~isempty(fullPooled)
        muFP = mean(fullPooled,'omitnan');
        if logScale; muFP = log10(max(muFP,1e-6)); end
        yline(ax, muFP,'--k','LineWidth',1.2,'Alpha',0.6,'Label','Full R²','FontSize',9);
    end

    hasAnySig = any(pAdj_t(r,:) < 0.05);
    ylim(ax, [0, yMaxAll + (0.06 * hasAnySig) + 0.02]);
    drawnow;

    for p = 1:nP
        vals = r2collapsed{r,p};
        if isempty(vals); continue; end

        muR2  = mean(vals,'omitnan');
        semR2 = std(vals,0,'omitnan') ./ sqrt(numel(vals));
        if logScale
            muPlot = log10(max(muR2,1e-6));
            hiPl   = log10(max(muR2+semR2,1e-6)) - muPlot;
            yTop   = muPlot + hiPl;
        else
            yTop = muR2 + semR2;
        end
        yOff = 0.025 * range(ax.YLim);

        pT = pAdj_t(r,p);
        if pT < 0.001;    stars = '***';
        elseif pT < 0.01; stars = '**';
        elseif pT < 0.05; stars = '*';
        else;              stars = '';
        end
        if ~isempty(stars)
            text(ax, p, yTop+yOff, stars, 'HorizontalAlignment','center',...
                'FontSize',13,'FontWeight','bold','Color','k');
        end

        statStr = sprintf('t: %s\nW: %s', format_pval(pAdj_t(r,p)), format_pval(pAdj_w(r,p)));
        text(ax, p, ax.YLim(1) - 0.005, statStr,...
            'HorizontalAlignment','center','VerticalAlignment','top',...
            'FontSize',7.5,'Color','k','Units','data','Clipping','off');
    end

    set(ax,'XTick',1:nP,'XTickLabel',predLabels,'TickDir','out','FontSize',11);
    if logScale
        ylabel('Unique R² (log_{10})','FontSize',12);
    else
        ylabel('Unique R²','FontSize',12);
    end
    title(regions{r},'FontSize',13,'FontWeight','bold','Color',regionColors(r,:));
    xlim([0.4 nP+0.6]);
    box off;
end

annotation('textbox',[0.01 0.00 0.99 0.04],...
    'String',['* p<0.05   ** p<0.01   *** p<0.001 (t-test, FDR-corrected, above bars)     '...
              't: one-sample t-test     W: Wilcoxon signed-rank (both FDR-corrected, below bars)'],...
    'EdgeColor','none','FontSize',8.5,'Color',[0.4 0.4 0.4],...
    'HorizontalAlignment','center');

end   % if false — Plot 3
%% =========================================================================
%  Plot 4: Observed vs predicted phasic traces
% =========================================================================
figure('Units','normalized','Position',[0.02 0.05 0.96 0.85]);
sgtitle(sprintf('Observed vs CV-Predicted Phasic ACh (%d-fold CV, window: %.1f–%.1f s)', ...
    kFolds, phasicStart, phasicEnd),'FontSize',14,'FontWeight','bold');

obsColor  = [0.2 0.2 0.2];
predColor = [0.5 0.5 0.5];

for s = 1:nS
    for r = 1:nR
        ax = subplot(nS, nR, (s-1)*nR + r);
        hold(ax, 'on');

        OBS = obsTraces{r,s};
        if ~isempty(OBS)
            OBS = OBS(all(isfinite(OBS),2),:);
            if size(OBS,1) >= 2
                mu  = mean(OBS,1,'omitnan');
                sem = std(OBS,0,1,'omitnan') ./ sqrt(size(OBS,1));
                mu  = safe_smooth(mu,  kGauss);
                sem = safe_smooth(sem, kGauss);
                shadedErrorBar(tPhasic, mu, sem,...
                    'lineProps',{'Color',obsColor,'LineWidth',2.0},...
                    'transparent',1,'patchSaturation',0.10);
            end
        end

        PRED = predTraces{r,s};
        if ~isempty(PRED)
            PRED = PRED(all(isfinite(PRED),2),:);
            if size(PRED,1) >= 2
                mu  = mean(PRED,1,'omitnan');
                sem = std(PRED,0,1,'omitnan') ./ sqrt(size(PRED,1));
                mu  = safe_smooth(mu,  kGauss);
                sem = safe_smooth(sem, kGauss);
                shadedErrorBar(tPhasic, mu, sem,...
                    'lineProps',{'Color',predColor,'LineWidth',1.8,'LineStyle','--'},...
                    'transparent',1,'patchSaturation',0.08);
            end
        end

        for p = 1:nP
            PP = partialTraces{r,s,p};
            if isempty(PP); continue; end
            PP = PP(all(isfinite(PP),2),:);
            if size(PP,1) < 2; continue; end
            mu  = mean(PP,1,'omitnan');
            sem = std(PP,0,1,'omitnan') ./ sqrt(size(PP,1));
            mu  = safe_smooth(mu,  kGauss);
            sem = safe_smooth(sem, kGauss);
            shadedErrorBar(tPhasic, mu, sem,...
                'lineProps',{'Color',predColors(p,:),'LineWidth',1.4,'LineStyle',':'},...
                'transparent',1,'patchSaturation',0.12);
        end

        xline(0,'--k','LineWidth',1.0,'Alpha',0.5);
        yline(0,'-','Color',[0.8 0.8 0.8],'LineWidth',0.8);
        xlabel('Time from stimulus onset (s)','FontSize',11);
        ylabel('ACh (z-scored)','FontSize',11);

        cvStr = '–';
        if ~isempty(cvR2{r,s})
            cv = cvR2{r,s}(isfinite(cvR2{r,s}));
            if ~isempty(cv); cvStr = sprintf('%.2f',mean(cv,'omitnan')); end
        end
        title(sprintf('%s – %s  (CV R²=%s)',regions{r},stimLabels{s},cvStr),...
            'FontSize',10,'FontWeight','bold','Color',regionColors(r,:));

        xlim([phasicStart phasicEnd]);
        box off;
        set(ax,'TickDir','out','FontSize',11);

        if s == 1 && r == 1
            legH = [
                plot(nan,nan,'-', 'Color',obsColor,  'LineWidth',2.0),...
                plot(nan,nan,'--','Color',predColor,  'LineWidth',1.8)
            ];
            legLab = {'Observed','Full model (CV)'};
            for p = 1:nP
                legH(end+1)   = plot(nan,nan,':','Color',predColors(p,:),'LineWidth',1.8); %#ok
                legLab{end+1} = [predLabels{p} ' only'];
            end
            legend(legH, legLab,'Location','best','FontSize',9,'Box','off');
        end
    end
end

%% =========================================================================
%  Sensitivity check: all stimulus types vs auditory only
%  ------------------------------------------------------------------------
%  Visual gratings change retinal illumination, so on those trials part of the
%  pupil signal is a pupillary light reflex rather than an arousal signal.
%  Because the light reflex is stimulus-locked, most of that variance is SHARED
%  with the stimulus regressor and is therefore already excluded from pupil's
%  unique R^2 — the partitioning handles it. This prints the comparison so the
%  claim can be checked rather than asserted.
%
%  Reported as a sensitivity check; the figure panels use all stimulus types,
%  so that every panel describes the same trials.
%% =========================================================================
fprintf('\n==================================================================\n');
fprintf('Sensitivity check: all stimulus types vs auditory only\n');
fprintf('(visual gratings drive a pupillary light reflex)\n');
fprintf('==================================================================\n');
fprintf('%-14s %-5s %8s %10s %10s %9s\n', ...
    'Measure','Reg','n all/aud','mean all','mean aud','%% change');

sensRows = {};
measAll = [{r2full_all}, {r2unique_all(:,1)}, {r2unique_all(:,2)}, {r2unique_all(:,3)}];
measAud = [{r2full_aud}, {r2unique_aud(:,1)}, {r2unique_aud(:,2)}, {r2unique_aud(:,3)}];
measNm  = {'Full model','Stimulus unique','FME unique','Pupil unique'};

for mi = 1:numel(measNm)
    for r = 1:nR
        a = measAll{mi}{r}; a = a(isfinite(a));
        b = measAud{mi}{r}; b = b(isfinite(b));
        if isempty(a) || isempty(b); continue; end
        mA = mean(a,'omitnan'); mB = mean(b,'omitnan');
        pctChange = 100 * (mB - mA) / max(abs(mA), eps);
        fprintf('%-14s %-5s %4d/%-4d %10.4f %10.4f %8.1f%%\n', ...
            measNm{mi}, regions{r}, numel(a), numel(b), mA, mB, pctChange);
        sensRows{end+1} = {measNm{mi}, regions{r}, numel(a), numel(b), ...
            mA, mB, pctChange}; %#ok
    end
end
fprintf('==================================================================\n');
fprintf('If pupil unique R^2 is similar with and without visual trials, the\n');
fprintf('light reflex is not driving it — state that in one sentence. A large\n');
fprintf('drop would instead belong in the Discussion.\n');
fprintf('==================================================================\n');

if ~isempty(sensRows)
    csvS = fullfile(pwd, 'fir_auditory_only_sensitivity.csv');
    fidS = fopen(csvS,'w');
    if fidS == -1
        warning('Could not write fir_auditory_only_sensitivity.csv');
    else
        fprintf(fidS, 'measure,region,n_all_stim,n_auditory_only,mean_all_stim,mean_auditory_only,pct_change\n');
        for k = 1:numel(sensRows)
            rw = sensRows{k};
            fprintf(fidS, '%s,%s,%d,%d,%.6f,%.6f,%.3f\n', ...
                rw{1}, rw{2}, rw{3}, rw{4}, rw{5}, rw{6}, rw{7});
        end
        fclose(fidS);
        fprintf('Sensitivity check saved to: %s\n', csvS);
    end
end

%% =========================================================================
%  Plot 11: Trial-by-trial correlation between prediction and data, by time bin
%  ------------------------------------------------------------------------
%  The trial-AVERAGED figure (Plot 4) shows that the stimulus-only prediction
%  tracks the observed mean closely while FME and pupil look flat. That is an
%  artefact of averaging: the stimulus contribution is identical on every trial
%  and survives averaging intact, whereas the behavioural contributions differ
%  trial to trial and cancel.
%
%  This figure removes that artefact. At each time bin, the correlation is
%  computed ACROSS TRIALS between observed and predicted ACh, so trials are
%  never averaged together and nothing washes out.
%
%  The stimulus-only prediction is CONSTANT across trials, so its correlation
%  is undefined (zero variance) and is plotted as zero — which is the point:
%  it cannot explain trial-to-trial structure, however well it reproduces the
%  mean. FME and pupil can.
%
%  WITHIN-MOUSE, then averaged. Pooling trials across animals would let
%  between-mouse differences in mean ACh and mean FME generate correlation that
%  is not trial-to-trial at all. Averaging uses the Fisher z transform, since
%  correlations are not additive.
%
%  NOTE: correlation is scale-free, so this shows WHETHER a predictor tracks
%  trial-to-trial structure, not HOW MUCH variance it accounts for. Magnitudes
%  are the unique R^2 values in Fig 4. Use this alongside Plot 4, not instead
%  of it — R1 2c asked for the timecourse of data and predictions, which the
%  averaged figure provides.
%% =========================================================================
figure('Units','centimeters','Position',[1 1 24 16], ...
       'Color','w','PaperPositionMode','auto');
sgtitle('Trial-by-trial correlation of prediction with observed ACh', ...
    'FontSize', 13, 'FontWeight','bold');

corrColors = [0.2 0.2 0.2; predColors];        % full model + 3 predictors
corrNames  = [{'Full model'}, predLabels];

for s = 1:nS
    for r = 1:nR
        ax = subplot(nS, nR, (s-1)*nR + r);
        hold(ax,'on');

        nMouseHere = numel(stObsM{r,s});
        if nMouseHere < 2
            axis(ax,'off'); continue;
        end

        for k = 1:4          % 1 = full model, 2:4 = each predictor alone
            Z = nan(nMouseHere, nPhasic);     % Fisher z per mouse per bin
            for mm = 1:nMouseHere
                O = stObsM{r,s}{mm};
                if k == 1
                    P = stFullM{r,s}{mm};
                else
                    P = stPartM{r,s,k-1}{mm};
                end
                if isempty(O) || isempty(P) || size(O,1) < 5; continue; end
                for bIdx = 1:nPhasic
                    a = O(:,bIdx); b = P(:,bIdx);
                    ok = isfinite(a) & isfinite(b);
                    if sum(ok) < 5; continue; end
                    if std(b(ok)) == 0 || std(a(ok)) == 0
                        Z(mm,bIdx) = 0;      % constant predictor: no trial variance
                        continue;
                    end
                    rr = corr(a(ok), b(ok));
                    rr = max(min(rr, 0.9999), -0.9999);
                    Z(mm,bIdx) = atanh(rr);              % Fisher z
                end
            end

            muZ  = mean(Z, 1, 'omitnan');
            semZ = std(Z, 0, 1, 'omitnan') ./ sqrt(sum(isfinite(Z),1));
            if all(~isfinite(muZ)); continue; end
            muR  = tanh(muZ);                            % back to r
            hiR  = tanh(muZ + semZ);
            loR  = tanh(muZ - semZ);

            ok = isfinite(muR);
            fill(ax, [tPhasic(ok), fliplr(tPhasic(ok))], [hiR(ok), fliplr(loR(ok))], ...
                corrColors(k,:), 'FaceAlpha', 0.12, 'EdgeColor','none');
            plot(ax, tPhasic(ok), muR(ok), '-', 'Color', corrColors(k,:), ...
                'LineWidth', 1.8);
        end

        yline(ax, 0, '-', 'Color', [0.8 0.8 0.8], 'LineWidth', 0.8);
        xlim(ax, [phasicStart phasicEnd]);
        ylim(ax, [-0.25 1]);
        if r == 1
            ylabel(ax, {stimLabels{s}; 'correlation (r)'}, ...
                'FontSize', 10, 'FontWeight','bold');
        end
        if s == nS; xlabel(ax, 'Time from onset (s)', 'FontSize', 10); end
        if s == 1
            title(ax, sprintf('%s (n=%d mice)', regions{r}, nMouseHere), ...
                'FontSize', 11, 'FontWeight','bold', 'Color', regionColors(r,:));
        end
        box off; set(ax,'TickDir','out','FontSize',9);

        if s == 1 && r == nR
            lh = gobjects(4,1);
            for k = 1:4
                lh(k) = plot(ax, nan, nan, '-', 'Color', corrColors(k,:), 'LineWidth', 1.8);
            end
            legend(ax, lh, corrNames, 'FontSize', 8, 'Box','off', 'Location','northeast');
        end
    end
end

annotation('textbox',[0.01 0 0.99 0.05], ...
    'String',['Correlation across trials between observed and predicted ACh at ' ...
              'each time bin, computed WITHIN mouse and averaged across mice via ' ...
              'the Fisher z transform (mean \pm SEM).  The stimulus-only ' ...
              'prediction is identical on every trial, so it carries no ' ...
              'trial-to-trial variance and its correlation is zero by ' ...
              'construction — it reproduces the trial mean (Fig SX-A) but cannot ' ...
              'explain variation between trials.  Correlation indicates whether ' ...
              'a predictor tracks trial-to-trial structure, not how much variance ' ...
              'it explains (see Fig 4).'], ...
    'EdgeColor','none','FontSize',8.5,'Color',[0.4 0.4 0.4], ...
    'HorizontalAlignment','center');

%% =========================================================================
%  Plot 9: Single-trial heatmaps — why stimulus reproduces the MEAN while
%          FME/pupil explain the VARIANCE
%  ------------------------------------------------------------------------
%  Trial-averaged observed-vs-predicted plots are the convention for models
%  whose regressors are event-locked, but they are actively misleading for
%  continuous behavioural regressors: averaging across trials preserves the
%  stimulus contribution intact (it is identical on every trial) while
%  cancelling the FME and pupil contributions (they differ trial to trial).
%  That is why the stimulus-only trace tracks the observed mean almost
%  perfectly despite contributing the LEAST unique R^2.
%
%  Read the columns: the stimulus-only panel is the SAME STRIPE REPEATED on
%  every row — it cannot express trial-to-trial structure. The observed and
%  FME-only panels do vary down the rows. Trials are sorted by observed
%  response amplitude, so any column that tracks the observed gradient is
%  capturing across-trial variance.
%% =========================================================================
heatStim = 1;                      % which stimulus type to display
nShowMax = 400;                    % row cap, for legibility

colsHM   = {'Observed', 'Full model', 'Stimulus only', 'FME only', 'Pupil only'};
nCol     = numel(colsHM);

figure('Units','centimeters','Position',[1 1 26 15], ...
       'Color','w','PaperPositionMode','auto');
sgtitle(sprintf('Single trials, pooled across mice — %s', stimLabels{heatStim}), ...
    'FontSize', 13, 'FontWeight','bold');

for r = 1:nR
    O = stObs{r,heatStim};
    if isempty(O); continue; end

    % sort trials by observed amplitude so across-trial structure is visible
    [~, ord] = sort(mean(O, 2, 'omitnan'), 'descend');
    if numel(ord) > nShowMax
        ord = ord(round(linspace(1, numel(ord), nShowMax)));
    end

    mats = {O(ord,:), stFull{r,heatStim}(ord,:), ...
            stPart{r,heatStim,1}(ord,:), stPart{r,heatStim,2}(ord,:), ...
            stPart{r,heatStim,3}(ord,:)};

    % common colour scale across columns within a region
    allV = cell2mat(cellfun(@(M) M(:), mats, 'UniformOutput', false)');
    cl   = prctile(abs(allV(isfinite(allV))), 99);
    if ~isfinite(cl) || cl == 0; cl = 1; end

    for cIdx = 1:nCol
        ax = subplot(nR, nCol, (r-1)*nCol + cIdx);
        imagesc(ax, tPhasic, 1:size(mats{cIdx},1), mats{cIdx}, [-cl cl]);
        set(ax, 'YDir','normal', 'TickDir','out', 'FontSize', 8);
        colormap(ax, parula);
        if r == 1;   title(ax, colsHM{cIdx}, 'FontSize', 10, 'FontWeight','bold'); end
        if cIdx == 1
            ylabel(ax, {regions{r}; 'trials (sorted)'}, 'FontSize', 9, ...
                'FontWeight','bold', 'Color', regionColors(r,:));
        else
            set(ax, 'YTickLabel', []);
        end
        if r == nR;  xlabel(ax, 'Time (s)', 'FontSize', 9); end
        if cIdx == nCol
            cb = colorbar(ax); cb.Label.String = 'ACh (z)'; cb.FontSize = 8;
        end
    end
end

annotation('textbox',[0.01 0 0.99 0.045], ...
    'String',['Trials sorted by observed amplitude, pooled across mice, common ' ...
              'colour scale within each region.  The stimulus-only column is ' ...
              'identical on every trial by construction, so it can reproduce the ' ...
              'trial MEAN but none of the across-trial variance — which is why it ' ...
              'contributes the least unique R^2 despite matching the averaged trace.'], ...
    'EdgeColor','none','FontSize',8.5,'Color',[0.4 0.4 0.4], ...
    'HorizontalAlignment','center');

%% =========================================================================
%  Plot 10: Observed vs predicted scatter — ties directly to the R^2 bars
%  R^2 is literally the fit of this scatter, so it is the most direct visual
%  counterpart to the bar plot. One point per trial-sample, pooled.
%% =========================================================================
figure('Units','centimeters','Position',[1 1 24 9], ...
       'Color','w','PaperPositionMode','auto');
sgtitle(sprintf('Observed vs CV-predicted, single trial-samples — %s', ...
    stimLabels{heatStim}), 'FontSize', 13, 'FontWeight','bold');

maxPts = 6000;
for r = 1:nR
    ax = subplot(1, nR, r);
    hold(ax,'on');

    O = stObs{r,heatStim};
    P = stFull{r,heatStim};
    if isempty(O); continue; end
    x = O(:); y = P(:);
    fin = isfinite(x) & isfinite(y);
    x = x(fin); y = y(fin);
    if numel(x) > maxPts
        sel = randperm(numel(x), maxPts);
        xs = x(sel); ys = y(sel);
    else
        xs = x; ys = y;
    end

    scatter(ax, xs, ys, 5, regionColors(r,:), 'filled', ...
        'MarkerFaceAlpha', 0.10, 'MarkerEdgeColor','none');

    lims = [min([x;y]) max([x;y])];
    plot(ax, lims, lims, '-', 'Color', [0.45 0.45 0.45], 'LineWidth', 1.2);
    bfit = polyfit(x, y, 1);
    plot(ax, lims, polyval(bfit, lims), '-', 'Color', regionColors(r,:), ...
        'LineWidth', 1.8);

    r2s = 1 - sum((x-y).^2) / sum((x-mean(x)).^2);
    text(ax, 0.05, 0.94, sprintf('R^2 = %.3f', max(r2s,0)), ...
        'Units','normalized', 'FontSize', 10, 'FontWeight','bold', ...
        'Color', regionColors(r,:));
    text(ax, 0.05, 0.86, sprintf('n = %d samples', numel(x)), ...
        'Units','normalized', 'FontSize', 8.5, 'Color', [0.45 0.45 0.45]);

    xlim(ax, lims); ylim(ax, lims); axis(ax,'square');
    xlabel(ax, 'Observed ACh (z)', 'FontSize', 10);
    if r == 1; ylabel(ax, 'Predicted ACh (z)', 'FontSize', 10); end
    title(ax, regions{r}, 'FontSize', 11, 'FontWeight','bold', ...
        'Color', regionColors(r,:));
    box off; set(ax,'TickDir','out','FontSize', 9);
end

annotation('textbox',[0.01 0 0.99 0.05], ...
    'String',['One point per trial-sample (sub-sampled for display; R^2 computed ' ...
              'on all points).  Grey line = unity, coloured line = least-squares fit.  ' ...
              'This is the scatter whose fit the R^2 bars summarise.'], ...
    'EdgeColor','none','FontSize',8.5,'Color',[0.4 0.4 0.4], ...
    'HorizontalAlignment','center');

%% =========================================================================
%  Helper: stars from p-value
% =========================================================================
    function stars = pstars(p)
        if p < 0.001;    stars = '***';
        elseif p < 0.01; stars = '**';
        elseif p < 0.05; stars = '*';
        else;             stars = 'ns';
        end
    end

%% =========================================================================
%  Collapse unique R² across stimulus types using per-mouse means
%  Each mouse contributes exactly one value per region × predictor.
% =========================================================================
r2coll = cell(nR, nP);
for r = 1:nR; for p = 1:nP
    r2coll{r,p} = r2unique_mn{r,p};
end; end

%% =========================================================================
%  Plot 5: Paired predictor comparisons within each region
%  Paired Wilcoxon across mice: S vs FME, S vs Pupil, FME vs Pupil
%  Uses collapsed-across-stim r2unique values.
%  Layout: 1 row x nR columns, each showing 3 pairwise comparisons.
% =========================================================================
if false %#ok<UNRCH>
pairNames  = {'S vs FME', 'S vs Pupil', 'FME vs Pupil'};
pairIdx    = [1 2; 1 3; 2 3];   % predictor index pairs
nPairs     = size(pairIdx, 1);

% Collect raw p-values for FDR correction (nPairs x nR)
pRaw5 = ones(nPairs, nR);
for r = 1:nR
    for pp = 1:nPairs
        a = r2coll{r, pairIdx(pp,1)};
        b = r2coll{r, pairIdx(pp,2)};
        [~, ia, ib] = intersect(... % approximate pairing by position
            (1:numel(a))', (1:numel(b))');
        if numel(ia) >= 3
            pRaw5(pp,r) = signrank(a(ia), b(ib));
        end
    end
end
pAdj5 = reshape(fdr_bh(pRaw5(:)), nPairs, nR);

figure('Units','normalized','Position',[0.05 0.15 0.85 0.50]);
sgtitle('Predictor Comparisons (paired Wilcoxon, FDR-corrected, collapsed across stim)', ...
    'FontSize', 12, 'FontWeight', 'bold');

for r = 1:nR
    ax = subplot(1, nR, r);
    hold(ax, 'on');

    yMaxAll = -Inf;
    for p = 1:nP
        vals = r2coll{r,p};
        if isempty(vals); continue; end
        muR2  = mean(vals,'omitnan');
        semR2 = std(vals,0,'omitnan') / sqrt(numel(vals));
        bar(ax, p, muR2, 0.6, 'FaceColor', predColors(p,:), ...
            'EdgeColor','none','FaceAlpha',0.75);
        errorbar(ax, p, muR2, semR2, 'k', 'LineWidth',1.5, ...
            'CapSize',5, 'LineStyle','none');
        yMaxAll = max(yMaxAll, muR2 + semR2);
    end

    drawnow; yRange = range(ax.YLim);
    % Significance brackets between pairs
    bracketY = yMaxAll + 0.08*yRange;
    for pp = 1:nPairs
        p1 = pairIdx(pp,1); p2 = pairIdx(pp,2);
        st = pstars(pAdj5(pp,r));
        bY = bracketY + (pp-1)*0.12*yRange;
        plot(ax,[p1 p2],[bY bY],'-k','LineWidth',0.8);
        text(ax,(p1+p2)/2, bY+0.02*yRange, st, ...
            'HorizontalAlignment','center','FontSize',10,'FontWeight','bold');
    end

    set(ax,'XTick',1:nP,'XTickLabel',predLabels,'TickDir','out','FontSize',10);
    xlim(ax,[0.4 nP+0.6]); ylim(ax,[0 inf]);
    if r==1; ylabel(ax,'Unique R²','FontSize',11); end
    title(ax, regions{r},'FontSize',12,'FontWeight','bold','Color',regionColors(r,:));
    box off;
end
annotation('textbox',[0.01 0 0.99 0.04],'String', ...
    'Paired Wilcoxon on pooled-across-stim unique R² values. FDR-corrected across all pairs × regions.', ...
    'EdgeColor','none','FontSize',8,'Color',[0.4 0.4 0.4],'HorizontalAlignment','center');

end   % if false — Plot 5 figures
%% =========================================================================
%  Plot 6: Between-region comparisons for each predictor
%  Unpaired Wilcoxon between region pairs (AC/VC, AC/PFC, VC/PFC)
%  collapsed across stimulus types.
%  Layout: 1 row x nP columns.
% =========================================================================
regionPairs  = [1 2; 1 3; 2 3];
regionPairLb = {'AC vs VC','AC vs PFC','VC vs PFC'};
nRP          = size(regionPairs,1);

% Compute p-values for all 3 predictors
pRaw6 = ones(nRP, nP);
for p = 1:nP
    for rp = 1:nRP
        a = r2coll{regionPairs(rp,1), p};
        b = r2coll{regionPairs(rp,2), p};
        if numel(a) >= 3 && numel(b) >= 3
            pRaw6(rp,p) = ranksum(a, b);
        end
    end
end
% FDR correct each predictor separately (3 pairs each)
pAdj6 = ones(nRP, nP);
for p = 1:nP
    pAdj6(:,p) = fdr_bh(pRaw6(:,p));
end

% One figure per predictor
predTitles6 = {'Stimulus Unique R²', 'FME Unique R²', 'Pupil Unique R²'};
for p = 1:nP
    if false %#ok<UNRCH>   % Plot 6 figures suppressed; stats + CSV retained
    figure('Units','normalized','Position',[0.28+0.02*(p-1) 0.25 0.38 0.50]);
    sgtitle(sprintf('%s Across Regions (unpaired Wilcoxon, FDR-corrected)', ...
        predTitles6{p}), 'FontSize', 12, 'FontWeight', 'bold');

    ax = axes;
    hold(ax, 'on');

    yMaxAll = -Inf;
    for r = 1:nR
        vals = r2coll{r,p};
        if isempty(vals); continue; end
        muR2  = mean(vals,'omitnan');
        semR2 = std(vals,0,'omitnan') / sqrt(numel(vals));
        bar(ax, r, muR2, 0.6, 'FaceColor', regionColors(r,:), ...
            'EdgeColor','none','FaceAlpha',0.75);
        errorbar(ax, r, muR2, semR2, 'k', 'LineWidth',1.5, ...
            'CapSize',5, 'LineStyle','none');
        yMaxAll = max(yMaxAll, muR2 + semR2);
    end

    drawnow; yRange = range(ax.YLim);
    bracketY = yMaxAll + 0.08*yRange;
    for rp = 1:nRP
        r1   = regionPairs(rp,1); r2idx = regionPairs(rp,2);
        st   = pstars(pAdj6(rp,p));
        bY   = bracketY + (rp-1)*0.12*yRange;
        plot(ax,[r1 r2idx],[bY bY],'-k','LineWidth',0.8);
        text(ax,(r1+r2idx)/2, bY+0.02*yRange, st, ...
            'HorizontalAlignment','center','FontSize',10,'FontWeight','bold');
    end

    set(ax,'XTick',1:nR,'XTickLabel',regions,'TickDir','out','FontSize',11);
    xlim(ax,[0.4 nR+0.6]); ylim(ax,[0 inf]);
    ylabel(ax, sprintf('Unique R² (%s)', predLabels{p}), 'FontSize',12);
    title(ax, predLabels{p},'FontSize',13,'FontWeight','bold','Color',predColors(p,:));
    box off;
    annotation('textbox',[0.01 0 0.99 0.04],'String', ...
        'Unpaired Wilcoxon rank-sum between regions. FDR-corrected across 3 pairs.', ...
        'EdgeColor','none','FontSize',8,'Color',[0.4 0.4 0.4],'HorizontalAlignment','center');
    end   % if false — Plot 6 figures
end

% CSV export: region comparisons for all predictors (Plot 6)
csvPath6 = fullfile(pwd, 'fir_region_comparison_stats.csv');
fid6 = fopen(csvPath6, 'w');
fprintf(fid6, 'predictor,region1,region2,n_region1,mean_region1,sem_region1,n_region2,mean_region2,sem_region2,p_raw_ranksum,p_FDR,significant_FDR05\n');
for p = 1:nP
    for rp = 1:nRP
        r1   = regionPairs(rp,1);
        r2   = regionPairs(rp,2);
        a    = r2coll{r1,p}; a = a(isfinite(a));
        b    = r2coll{r2,p}; b = b(isfinite(b));
        fprintf(fid6, '%s,%s,%s,%d,%.4f,%.4f,%d,%.4f,%.4f,%.4f,%.4f,%d\n', ...
            predLabels{p}, regions{r1}, regions{r2}, ...
            numel(a), mean(a,'omitnan'), std(a,0,'omitnan')/sqrt(numel(a)), ...
            numel(b), mean(b,'omitnan'), std(b,0,'omitnan')/sqrt(numel(b)), ...
            pRaw6(rp,p), pAdj6(rp,p), pAdj6(rp,p)<0.05);
    end
end
fclose(fid6);
fprintf('FIR region comparison stats saved to: %s\n', csvPath6);

%% =========================================================================
%  Plot 7: Between-stimulus-type comparisons
%  Paired Wilcoxon across mice (same mouse always has all 3 stim types)
%  within each predictor and region.
%  Layout: nP rows x nR columns.
% =========================================================================
if false %#ok<UNRCH>
stimPairs  = [1 2; 1 3; 2 3];
stimPairLb = {'Ripple vs Tones','Ripple vs Vis','Tones vs Vis'};
nSP        = size(stimPairs,1);

% Compute per-mouse means collapsed across brain regions for each
% predictor × stim type. Each mouse contributes the mean of whichever
% regions it has recorded.
% r2unique_byStim_mn{p,s} = [nMice x 1] mean across regions per mouse
r2unique_byStim_mn = cell(nP, nS);
for p = 1:nP
    for s = 1:nS
        % Collect all (mouse, region) observations, average within mouse
        % Use r2unique{r,s,p} entries aligned by position within region
        % Since we don't store mouse identity per entry, we use
        % r2unique_mn approach: average r2unique across regions per mouse
        % using the per-stim per-region per-mouse values stored in
        % r2unique{r,s,p}. We approximate by collecting min-length vectors
        % across regions and averaging row-wise.
        regionVecs = {};
        for r = 1:nR
            v = r2unique{r,s,p};
            v = v(isfinite(v));
            if ~isempty(v); regionVecs{end+1} = v; end %#ok
        end
        if isempty(regionVecs)
            r2unique_byStim_mn{p,s} = [];
            continue;
        end
        % Each region may have different numbers of mice; take the minimum
        % length and average row-wise (approximate — mice are not
        % individually tracked across regions here)
        minN = min(cellfun(@numel, regionVecs));
        mat  = cell2mat(cellfun(@(v) v(1:minN), regionVecs, 'UniformOutput', false)');
        r2unique_byStim_mn{p,s} = mean(mat, 2, 'omitnan');
    end
end

% Paired Wilcoxon across stim types, collapsed across regions, per predictor
pRaw7c = ones(nSP, nP);
for p = 1:nP
    for sp = 1:nSP
        s1 = stimPairs(sp,1); s2 = stimPairs(sp,2);
        a  = r2unique_byStim_mn{p,s1}; a = a(isfinite(a));
        b  = r2unique_byStim_mn{p,s2}; b = b(isfinite(b));
        n  = min(numel(a), numel(b));
        if n >= 3
            pRaw7c(sp,p) = signrank(a(1:n), b(1:n));
        end
    end
end
pAdj7c = reshape(fdr_bh(pRaw7c(:)), nSP, nP);

% One figure per predictor
for p = 1:nP
    figure('Units','normalized','Position',[0.28+0.02*(p-1) 0.15 0.38 0.55]);
    sgtitle(sprintf('%s – Stim Type Comparison (region-averaged, paired Wilcoxon)', ...
        predLabels{p}), 'FontSize',11,'FontWeight','bold');

    ax = axes;
    hold(ax,'on');

    yMaxAll = -Inf;
    gv = linspace(0.65, 0.15, nS);
    for s = 1:nS
        vals = r2unique_byStim_mn{p,s};
        vals = vals(isfinite(vals));
        if isempty(vals); continue; end
        col   = predColors(p,:) * gv(s) + [1 1 1]*(1-gv(s));
        muR2  = mean(vals,'omitnan');
        semR2 = std(vals,0,'omitnan') / sqrt(numel(vals));
        bar(ax, s, muR2, 0.6, 'FaceColor', col, ...
            'EdgeColor','none','FaceAlpha',0.85);
        errorbar(ax, s, muR2, semR2, 'k','LineWidth',1.5,'CapSize',5,'LineStyle','none');
        yMaxAll = max(yMaxAll, muR2 + semR2);
    end

    drawnow; yRange = range(ax.YLim);
    bracketY = yMaxAll + 0.08*yRange;
    for sp = 1:nSP
        s1 = stimPairs(sp,1); s2 = stimPairs(sp,2);
        st = pstars(pAdj7c(sp,p));
        bY = bracketY + (sp-1)*0.12*yRange;
        plot(ax,[s1 s2],[bY bY],'-k','LineWidth',0.8);
        text(ax,(s1+s2)/2, bY+0.02*yRange, st, ...
            'HorizontalAlignment','center','FontSize',10,'FontWeight','bold');
    end

    set(ax,'XTick',1:nS,'XTickLabel',stimLabels,'TickDir','out','FontSize',11);
    xtickangle(ax,20);
    xlim(ax,[0.4 nS+0.6]); ylim(ax,[0 inf]);
    ylabel(ax, sprintf('Unique R² (%s)', predLabels{p}),'FontSize',12);
    title(ax, predLabels{p},'FontSize',13,'FontWeight','bold','Color',predColors(p,:));
    box off;
    annotation('textbox',[0.01 0 0.99 0.04],'String', ...
        ['Unique R² averaged across recorded brain regions per mouse.  ' ...
         'Paired Wilcoxon across mice.  FDR-corrected across 3 stim pairs.'], ...
        'EdgeColor','none','FontSize',8,'Color',[0.4 0.4 0.4],'HorizontalAlignment','center');
end

end   % if false — Plot 7 figures
%% =========================================================================
%  Plot 8: Commonality analysis (suppressed — code retained for reference)
% =========================================================================
if false %#ok<UNRCH>

commColors = [
    predColors(1,:);                    % Unique S     — blue
    predColors(2,:);                    % Unique FME   — orange
    predColors(3,:);                    % Unique Pupil — purple
    mean(predColors([1 2],:),1)*0.85;  % Shared S+FME
    mean(predColors([1 3],:),1)*0.85;  % Shared S+Pupil
    mean(predColors([2 3],:),1)*0.85;  % Shared FME+Pupil
    [0.35 0.35 0.35];                  % Shared All
];
commLabels = {'Unique S','Unique FME','Unique Pupil', ...
              'Shared S+FME','Shared S+Pupil','Shared FME+Pupil','Shared All'};
nComm = 7;

% Pre-compute compMat for each region (one row per mouse, NOT floored)
compMats = cell(nR,1);
for r = 1:nR
    nMiceHere = numel(r2full_mn{r});
    cMat = nan(nMiceHere, nComm);
    for mi = 1:nMiceHere
        rs   = r2single_mn{r,1}(mi);
        rf   = r2single_mn{r,2}(mi);
        rp_v = r2single_mn{r,3}(mi);
        rsf  = r2pair_mn{r,1}(mi);
        rsp  = r2pair_mn{r,2}(mi);
        rfp  = r2pair_mn{r,3}(mi);
        rall = r2full_mn{r}(mi);
        uS   = rall - rfp;
        uF   = rall - rsp;
        uP   = rall - rsf;
        sAll = rs + rf + rp_v - rsf - rsp - rfp + rall;
        sSF  = rsf - rs - rf + sAll;
        sSP  = rsp - rs - rp_v + sAll;
        sFP  = rfp - rf - rp_v + sAll;
        cMat(mi,:) = [uS uF uP sSF sSP sFP sAll];
    end
    compMats{r} = cMat;
end

%% -------------------------------------------------------------------------
%  Plot 8a: Stacked bar (mean components floored at 0, SEM on total)
% -------------------------------------------------------------------------
figure('Units','normalized','Position',[0.05 0.55 0.85 0.38]);
sgtitle('Commonality Analysis – Stacked Bar (mean across mice)', ...
    'FontSize',13,'FontWeight','bold');

for r = 1:nR
    ax = subplot(1, nR, r);
    hold(ax,'on');
    cMat = compMats{r};
    if isempty(cMat)
        title(ax,regions{r},'FontSize',11,'FontWeight','bold','Color',regionColors(r,:));
        continue;
    end

    muComp = mean(max(cMat,0), 1, 'omitnan');
    bottom = 0; bw = 0.5;
    for ci = 1:nComm
        h = muComp(ci);
        if h <= 0; continue; end
        patch(ax, [1-bw 1+bw 1+bw 1-bw], ...
            [bottom bottom bottom+h bottom+h], ...
            commColors(ci,:), 'EdgeColor','none','FaceAlpha',0.85);
        bottom = bottom + h;
    end

    totalPerMouse = sum(cMat, 2, 'omitnan');
    totalMu  = mean(totalPerMouse,'omitnan');
    totalSem = std(totalPerMouse,'omitnan') / sqrt(size(cMat,1));
    errorbar(ax, 1, totalMu, totalSem, 'k','LineWidth',1.5,'CapSize',6,'LineStyle','none');

    if ~isempty(r2full_mn{r})
        yline(ax, mean(r2full_mn{r},'omitnan'), '--k','LineWidth',1.0,'Alpha',0.5,...
            'Label','Full R²','FontSize',8);
    end
    set(ax,'XTick',[],'TickDir','out','FontSize',10);
    if r==1; ylabel(ax,'R² (decomposed)','FontSize',11); end
    title(ax,regions{r},'FontSize',12,'FontWeight','bold','Color',regionColors(r,:));
    xlim(ax,[0.4 1.6]); ylim(ax,[0 inf]); box off;

    if r == nR
        legH = gobjects(nComm,1);
        for ci = 1:nComm
            legH(ci) = patch(nan,nan,commColors(ci,:),'EdgeColor','none','FaceAlpha',0.85);
        end
        legend(ax, legH, commLabels,'Location','eastoutside','FontSize',9,'Box','off');
    end
end
annotation('textbox',[0.01 0 0.99 0.04],'String', ...
    ['8a: Stacked bar = mean of each component (floored at 0 for display).  ' ...
     'Error bar = SEM of total R².  One value per mouse (averaged across stim types).'], ...
    'EdgeColor','none','FontSize',8,'Color',[0.4 0.4 0.4],'HorizontalAlignment','center');

%% -------------------------------------------------------------------------
%  Plots 8b and 8c: suppressed — code retained for reference
% -------------------------------------------------------------------------
if false %#ok<UNRCH>

%% -------------------------------------------------------------------------
%  Plot 8b: Grouped dot-and-errorbar per component, coloured by region
% -------------------------------------------------------------------------
figure('Units','normalized','Position',[0.02 0.05 0.96 0.50]);
sgtitle('Commonality Analysis – Component-wise Mean ± SEM by Region', ...
    'FontSize',13,'FontWeight','bold');

groupGap = 1.2; bwG = 0.55;
xPosG = zeros(nR, nComm);
for ci = 1:nComm
    for r = 1:nR
        xPosG(r,ci) = (ci-1)*(nR + groupGap) + r;
    end
end
xTickG = mean(xPosG,1);
xMaxG  = max(xPosG(:)) + 1;

ax = axes('Units','normalized','Position',[0.06 0.18 0.90 0.68]);
hold(ax,'on');

for ci = 1:nComm
    for r = 1:nR
        cMat = compMats{r};
        if isempty(cMat); continue; end
        vals = cMat(:,ci); vals = vals(isfinite(vals));
        if isempty(vals); continue; end
        muC  = mean(vals,'omitnan');
        semC = std(vals,0,'omitnan') / sqrt(numel(vals));
        xg   = xPosG(r,ci);
        bar(ax, xg, muC, bwG, 'FaceColor', commColors(ci,:), ...
            'EdgeColor', regionColors(r,:), 'LineWidth', 1.2, 'FaceAlpha', 0.75);
        errorbar(ax, xg, muC, semC, 'k','LineWidth',1.3,'CapSize',4,'LineStyle','none');
        jit = 0.07*randn(numel(vals),1);
        scatter(ax, xg+jit, vals, 18, regionColors(r,:), 'filled', ...
            'MarkerFaceAlpha',0.5,'MarkerEdgeColor','none');
    end
    if ci < nComm
        xDiv = (xPosG(nR,ci) + xPosG(1,ci+1)) / 2;
        xline(ax, xDiv, '-','Color',[0.88 0.88 0.88],'LineWidth',0.8);
    end
end

yline(ax, 0, '-','Color',[0.7 0.7 0.7],'LineWidth',0.8);
set(ax,'XTick',xTickG,'XTickLabel',commLabels,'TickDir','out','FontSize',9);
xtickangle(ax,20);
xlim(ax,[0.3 xMaxG]);
ylabel(ax,'R² component','FontSize',11); box off;

legH = gobjects(nR,1);
for r = 1:nR
    legH(r) = patch(nan,nan,'w','EdgeColor',regionColors(r,:),'LineWidth',2);
end
legend(ax, legH, regions,'Location','northeast','FontSize',10,'Box','off');

annotation('textbox',[0.01 0 0.99 0.04],'String', ...
    ['8b: Bars filled by component colour, edged by region colour.  ' ...
     'Dots = individual mice.  Mean ± SEM per component per region.'], ...
    'EdgeColor','none','FontSize',8,'Color',[0.4 0.4 0.4],'HorizontalAlignment','center');

%% -------------------------------------------------------------------------
%  Plot 8c: Connected dot profile — one line per mouse, mean ± SEM overlaid
% -------------------------------------------------------------------------
figure('Units','normalized','Position',[0.02 0.05 0.96 0.50]);
sgtitle('Commonality Analysis – Per-Mouse Profiles', ...
    'FontSize',13,'FontWeight','bold');

for r = 1:nR
    ax = subplot(1, nR, r);
    hold(ax,'on');
    cMat = compMats{r};
    if isempty(cMat)
        title(ax,regions{r},'FontSize',11,'FontWeight','bold','Color',regionColors(r,:));
        continue;
    end
    c = regionColors(r,:);
    nMiceR = size(cMat,1);

    for mi = 1:nMiceR
        plot(ax, 1:nComm, cMat(mi,:), '-o', ...
            'Color', [c 0.25], 'MarkerSize', 4, ...
            'MarkerFaceColor', c, 'MarkerEdgeColor','none','LineWidth',0.8);
    end

    muC  = mean(cMat, 1, 'omitnan');
    semC = std(cMat,  0, 1, 'omitnan') / sqrt(nMiceR);
    fill(ax, [1:nComm, fliplr(1:nComm)], [muC+semC, fliplr(muC-semC)], ...
        c, 'FaceAlpha',0.20,'EdgeColor','none');
    plot(ax, 1:nComm, muC, '-o', 'Color', c, 'LineWidth', 2.2, ...
        'MarkerSize', 6, 'MarkerFaceColor', c, 'MarkerEdgeColor','w');

    yline(ax, 0, '-','Color',[0.75 0.75 0.75],'LineWidth',0.9);
    set(ax,'XTick',1:nComm,'XTickLabel',commLabels,'TickDir','out','FontSize',8);
    xtickangle(ax,25);
    xlim(ax,[0.5 nComm+0.5]);
    if r==1; ylabel(ax,'R² component','FontSize',11); end
    title(ax,regions{r},'FontSize',11,'FontWeight','bold','Color',c);
    box off;
end
annotation('textbox',[0.01 0 0.99 0.04],'String', ...
    ['8c: Thin lines = individual mice; bold line = mean ± SEM (shaded).  ' ...
     'Negative values shown (not floored).  One value per mouse per component.'], ...
    'EdgeColor','none','FontSize',8,'Color',[0.4 0.4 0.4],'HorizontalAlignment','center');

end   % if false — Plots 8b and 8c

end   % if false — Plot 8 commonality

end
