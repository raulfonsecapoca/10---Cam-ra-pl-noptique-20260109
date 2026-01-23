%% Process a decoded/rectified microlens light field from .mat
% Robust for LF stored as 5-D array (e.g., 9x9x381x382x4).
%
% Output:
% - Central view (per-channel + combined)
% - Integrated image over (u,v) (per-channel + combined)
% - Mosaic of subaperture views (combined or selected channel)
% - Refocus (shift-and-sum) on combined LF
% - Channel diagnostics (correlation, energy)

clear; close all; clc;

%% ---------------- USER SETTINGS ----------------
matFile = "LF_01.mat";

% Display normalization (percentiles)
pLow  = 1;
pHigh = 99.9;

% How to combine the 4 channels into a single intensity per raxel:
% "mean"     : simple average across C
% "median"   : robust to outlier channels
% "weighted" : automatic weights from per-channel SNR proxy (recommended)
combineMode = "weighted";   % "mean" | "median" | "weighted"
epsW = 1e-12;

% Visualization
showPerChannelFigures = true;
showMosaic = true;

% Mosaic options
mosaicUse = "combined";     % "combined" or "channel"
mosaicChannel = 1;          % used only if mosaicUse="channel"

% Refocus options (shift-and-sum)
doRefocus = true;
alphaList = [-1.5 -0.8 0 0.8 1.5];   % try a wider range if needed
refocusUse = "combined";             % "combined" or "channel"
refocusChannel = 1;                  % used only if refocusUse="channel"

% OPTIONAL: if your LF is known to be already [U V Y X C], you can set this true
forceNoPermute = false;

%% ---------------- 1) LOAD + FIND 5D LF VARIABLE ----------------
vars = whos('-file', matFile);
disp("Variables in MAT file:");
for k = 1:numel(vars)
    szStr = sprintf('%dx', vars(k).size); szStr(end) = [];
    fprintf("  %s : [%s]  %s\n", vars(k).name, szStr, vars(k).class);
end

S = load(matFile);

lfName = "";
for k = 1:numel(vars)
    name = vars(k).name;
    val = S.(name);
    if isnumeric(val) && ndims(val) == 5
        lfName = name;
        break;
    end
end
if lfName == ""
    error("No 5-D numeric variable found in the MAT file.");
end

LFraw = S.(lfName);
fprintf("\nUsing LF variable: %s\n", lfName);
fprintf("LF size: %s, class: %s\n", mat2str(size(LFraw)), class(LFraw));

LF = double(LFraw);

%% ---------------- 2) STANDARDIZE DIM ORDER TO [U V Y X C] ----------------
% We try to detect U,V from size==9, channels from size==4, and spatial from 381/382.
sz = size(LF);

if ~forceNoPermute
    dims9  = find(sz == 9);
    dims4  = find(sz == 4);
    dimsYX = find((sz == 381) | (sz == 382));

    if numel(dims9) >= 2 && ~isempty(dims4) && numel(dimsYX) >= 2
        du = dims9(1);
        dv = dims9(2);

        dY = dimsYX(find(sz(dimsYX) == 381, 1, "first"));
        dX = dimsYX(find(sz(dimsYX) == 382, 1, "first"));
        if isempty(dY), dY = dimsYX(1); end
        if isempty(dX), dX = dimsYX(2); end

        dC = dims4(1);

        perm = [du dv dY dX dC];
        LF = permute(LF, perm);
        fprintf("Permute used: %s\n", mat2str(perm));
    else
        fprintf("Permute not applied (auto-detect failed). Using original order.\n");
    end
else
    fprintf("forceNoPermute=true -> using LF as-is.\n");
end

[U,V,Y,X,C] = size(LF);
fprintf("Standardized LF: U=%d V=%d Y=%d X=%d C=%d\n\n", U,V,Y,X,C);

if C ~= 4
    warning("C is %d (not 4). The code still works, but combineMode assumes multi-channel samples.", C);
end

%% ---------------- 3) NORMALIZE FOR DISPLAY ONLY ----------------
LFdisp = LF;  % keep LF as-is for math; LFdisp only for visualization

for c = 1:C
    tmp = LFdisp(:,:,:,:,c);
    lo = prctile(tmp(:), pLow);
    hi = prctile(tmp(:), pHigh);
    LFdisp(:,:,:,:,c) = (tmp - lo) / max(eps, (hi - lo));
end
LFdisp = max(0, min(1, LFdisp));

%% ---------------- 4) CHANNEL DIAGNOSTICS ----------------
% We compute simple stats to understand if channels are similar (samples) or color-like.
fprintf("Channel diagnostics:\n");
energy = zeros(C,1);
mu = zeros(C,1);
sig = zeros(C,1);

for c = 1:C
    tmp = LFdisp(:,:,:,:,c);
    mu(c) = mean(tmp(:));
    sig(c)= std(tmp(:));
    energy(c) = mean(tmp(:).^2);
    fprintf("  c=%d: mean=%.4f, std=%.4f, energy=%.4f\n", c, mu(c), sig(c), energy(c));
end

% Correlation between channels (using a subsample for speed)
nSample = min(200000, numel(LFdisp(:,:,:,:,1)));
idx = randperm(numel(LFdisp(:,:,:,:,1)), nSample);
chanVec = zeros(nSample, C);
for c = 1:C
    tmp = LFdisp(:,:,:,:,c);
    chanVec(:,c) = tmp(idx);
end
R = corrcoef(chanVec);
disp("Approx channel correlation matrix:");
disp(R);

figure('Name','Channel correlation matrix');
imagesc(R); axis image; colorbar;
title("Channel correlation (display-normalized)");
set(gca,'XTick',1:C,'YTick',1:C);

%% ---------------- 5) COMBINE CHANNELS THE RIGHT WAY ----------------
% We create LF1 = LFdisp combined over channel dimension -> [U V Y X]
% combineMode:
% - mean/median: straightforward
% - weighted: uses a simple SNR proxy: w_c ~ std(tmp) / (mean(tmp)+eps)
%   Then normalized weights. This tends to downweight very flat/noisy channels.
switch lower(combineMode)
    case "mean"
        w = ones(C,1) / C;

    case "median"
        w = []; % not used

    case "weighted"
        snrProxy = sig ./ max(epsW, mu);          % higher -> more contrast relative to mean
        snrProxy = max(epsW, snrProxy);
        w = snrProxy / sum(snrProxy);

    otherwise
        error("Unknown combineMode: %s", combineMode);
end

if lower(combineMode) == "median"
    LF1 = median(LFdisp, 5);
else
    LF1 = zeros(U,V,Y,X);
    for c = 1:C
        LF1 = LF1 + w(c) * LFdisp(:,:,:,:,c);
    end
end

% Print weights if used
if exist('w','var') && ~isempty(w)
    fprintf("\nChannel combine mode: %s\n", combineMode);
    fprintf("Weights: "); fprintf("%.3f ", w); fprintf("\n\n");
end

%% ---------------- 6) CENTRAL VIEW + INTEGRATED (PER-CHANNEL + COMBINED) ----------------
uc = ceil(U/2);
vc = ceil(V/2);

% Central combined view
centerCombined = squeeze(LF1(uc, vc, :, :)); % [Y X]
figure('Name','Central view (combined)');
imagesc(centerCombined); axis image off; colormap gray;
title(sprintf("Central view (combined) u=%d v=%d", uc, vc));

if showPerChannelFigures
    for c = 1:C
        centerC = squeeze(LFdisp(uc, vc, :, :, c));
        figure('Name',sprintf('Central view (channel %d)', c));
        imagesc(centerC); axis image off; colormap gray;
        title(sprintf("Central view (channel %d) u=%d v=%d", c, uc, vc));
    end
end

% Integrated over u,v (combined)
integratedCombined = squeeze(mean(mean(LF1, 1), 2)); % [Y X]
figure('Name','Integrated over (u,v) (combined)');
imagesc(integratedCombined); axis image off; colormap gray;
title("Integrated image: mean over (u,v) (combined)");

if showPerChannelFigures
    for c = 1:C
        img2D_c = squeeze(mean(mean(LFdisp(:,:,:,:,c), 1), 2)); % [Y X]
        figure('Name',sprintf('Integrated over (u,v) (channel %d)', c));
        imagesc(img2D_c); axis image off; colormap gray;
        title(sprintf("Integrated image: mean over (u,v) (channel %d)", c));
    end
end

%% ---------------- 7) MOSAIC OF SUBAPERTURE VIEWS ----------------
if showMosaic
    if lower(mosaicUse) == "combined"
        getTile = @(u,v) squeeze(LF1(u,v,:,:)); % [Y X]
        figName = "Mosaic 9x9 (combined)";
    else
        ch = mosaicChannel;
        getTile = @(u,v) squeeze(LFdisp(u,v,:,:,ch));
        figName = sprintf("Mosaic 9x9 (channel %d)", ch);
    end

    mosaic = zeros(U*Y, V*X);
    for u = 1:U
        for v = 1:V
            tile = getTile(u,v);
            r0 = (u-1)*Y + 1;
            c0 = (v-1)*X + 1;
            mosaic(r0:r0+Y-1, c0:c0+X-1) = tile;
        end
    end

    figure('Name',figName);
    imagesc(mosaic); axis image off; colormap gray;
    title(figName);
end

%% ---------------- 8) REFOCUS (SHIFT-AND-SUM) ----------------
if doRefocus
    uCoord = (1:U) - uc;
    vCoord = (1:V) - vc;

    if lower(refocusUse) == "combined"
        getView = @(u,v) squeeze(LF1(u,v,:,:));
        tag = "combined";
    else
        ch = refocusChannel;
        getView = @(u,v) squeeze(LFdisp(u,v,:,:,ch));
        tag = sprintf("channel %d", ch);
    end

    for a = alphaList
        acc = zeros(Y, X);

        for u = 1:U
            for v = 1:V
                img = getView(u,v);

                % Shift in pixels proportional to angular offset
                dy = a * uCoord(u);
                dx = a * vCoord(v);

                acc = acc + imtranslate(img, [dx dy], 'linear', 'FillValues', 0);
            end
        end

        refocused = acc / (U*V);

        figure('Name', sprintf('Refocus alpha=%.2f (%s)', a, tag));
        imagesc(refocused); axis image off; colormap gray;
        title(sprintf("Refocus (shift-and-sum) alpha=%.2f (%s)", a, tag));
    end
end

disp("Done.");
