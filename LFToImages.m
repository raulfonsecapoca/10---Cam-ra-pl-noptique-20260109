%% LF visualization + refocus sweep + all-in-focus fusion (RGB + Weight)
% Assumes LF is 5-D: [j, i, l, k, c], with c=1..3 RGB, c=4 Weight

clear; close all; clc;

%% ---------------- USER SETTINGS ----------------
ImgName = 'IMG_0005';
matFile = ['Output\LF_from_' ImgName '.mat'];
lfVarName = '';                     % leave empty to auto-detect 5D numeric var

% Visualization normalization (percentiles)
pLow  = 1;
pHigh = 99.9;

% Refocus alpha sweep (dense / quasi-continuous)
doRefocus = true;
alphaMin = -1.5;
alphaMax =  1.5;
nAlpha   = 40;
alphaList = linspace(alphaMin, alphaMax, nAlpha);

% Use weights in refocus
useWeights = true;
weightGamma = 1.0;      % >1 emphasizes high-confidence regions

% Fusion settings (all-in-focus)
doAllInFocus = true;
focusWin = 9;           % odd window size for local focus measure smoothing
laplacianEps = 1e-12;
fusionSoftness = 1.0;   % higher -> softer blending, lower -> more winner-take-all

% Mosaic options
showMosaicChannels = true;
showMosaicColor    = true;

%% ---------------- 1) LOAD + FIND LF ----------------
S = load(matFile);
vars = fieldnames(S);

LF = [];
if ~isempty(lfVarName) && isfield(S, lfVarName)
    LF = S.(lfVarName);
else
    for k = 1:numel(vars)
        v = S.(vars{k});
        if isnumeric(v) && ndims(v) == 5
            LF = v;
            fprintf("Using LF variable: %s\n", vars{k});
            break;
        end
    end
end
if isempty(LF)
    error("No 5-D numeric variable found in MAT file.");
end

LF = double(LF);
sz = size(LF);
fprintf("LF size: %s\n", mat2str(sz));

% Expect [U V Y X C] but Dansereau (toolbox) naming is [j i l k c]
U = sz(1); V = sz(2); Y = sz(3); X = sz(4); C = sz(5);

if C < 4
    error("Expected at least 4 channels (RGB + Weight). Got C=%d.", C);
end

% Split channels
LF_RGB = LF(:,:,:,:,1:3);
LF_W   = LF(:,:,:,:,4);

% Normalize for display only
[LFdisp_RGB, LFdisp_W] = normalizeForDisplay(LF_RGB, LF_W, pLow, pHigh);


%% ---------------- 2) MOSAIC 9x9 FOR EACH CHANNEL ----------------
%if showMosaicChannels
%    % Channel 1..3 + Weight channel 4
%    figure('Name','Mosaic 9x9 - R'); imagesc(buildMosaic(squeeze(LFdisp_RGB(:,:,:,:,1)))); axis image off; colormap gray; title('Mosaic 9x9 - R');
%    figure('Name','Mosaic 9x9 - G'); imagesc(buildMosaic(squeeze(LFdisp_RGB(:,:,:,:,2)))); axis image off; colormap gray; title('Mosaic 9x9 - G');
%    figure('Name','Mosaic 9x9 - B'); imagesc(buildMosaic(squeeze(LFdisp_RGB(:,:,:,:,3)))); axis image off; colormap gray; title('Mosaic 9x9 - B');
%    figure('Name','Mosaic 9x9 - Weight'); imagesc(buildMosaic(LFdisp_W)); axis image off; colormap gray; title('Mosaic 9x9 - Weight');
%end

%% ---------------- 3) MOSAIC 9x9 COLOR (RGB) ----------------
if showMosaicColor
    mosaicRGB = buildMosaicRGB(LFdisp_RGB);
    figure('Name','Mosaic 9x9 - RGB');
    imshow(mosaicRGB);
    title('Mosaic 9x9 - RGB');
end

%% ---------------- 4) CENTRAL VIEW (RGB + Weight) ----------------
uc = ceil(U/2);
vc = ceil(V/2);

centerRGB = squeeze(LFdisp_RGB(uc, vc, :, :, :)); % [Y X 3]
centerW   = squeeze(LFdisp_W(uc, vc, :, :));      % [Y X]

figure('Name','Central view - RGB'); imshow(centerRGB); title(sprintf('Central view RGB (u=%d,v=%d)', uc, vc));
figure('Name','Central view - Weight'); imagesc(centerW); axis image off; colormap gray; title(sprintf('Central view Weight (u=%d,v=%d)', uc, vc));

%% ---------------- 5) REFOCUS SWEEP (SHIFT-AND-SUM) ----------------
refocusStack = [];   % [Y X 3 nAlpha]
refocusWsum  = [];   % [Y X nAlpha]

if doRefocus
    fprintf("Refocus sweep: %d alphas from %.3f to %.3f\n", nAlpha, alphaMin, alphaMax);

    % Angular coordinates centered
    uCoord = (1:U) - uc;
    vCoord = (1:V) - vc;

    refocusStack = zeros(Y, X, 3, nAlpha);
    refocusWsum  = zeros(Y, X, nAlpha);

    for aIdx = 1:nAlpha
        a = alphaList(aIdx);

        accRGB = zeros(Y, X, 3);
        accW   = zeros(Y, X);

        for u = 1:U
            for v = 1:V
                img = squeeze(LF_RGB(u,v,:,:,:)); % double, not display-normalized
                w   = squeeze(LF_W(u,v,:,:));     % weights

                % Optional: weight shaping
                if useWeights
                    ww = max(w, 0) .^ weightGamma;
                else
                    ww = ones(size(w));
                end

                % Shifts: convention (dx from v, dy from u)
                dy = a * uCoord(u);
                dx = a * vCoord(v);

                imgS = imtranslate(img, [dx dy], 'linear', 'FillValues', 0);
                wS   = imtranslate(ww,  [dx dy], 'linear', 'FillValues', 0);

                accRGB = accRGB + imgS .* wS;
                accW   = accW   + wS;
            end
        end

        % Normalize by weight sum (avoid divide-by-zero)
        accW = max(accW, 1e-12);
        refocused = accRGB ./ accW;

        % Store (for fusion)
        refocusStack(:,:,:,aIdx) = refocused;
        refocusWsum(:,:,aIdx)    = accW;

        % Quick display every few steps
        if mod(aIdx, max(1, floor(nAlpha/8))) == 1 || aIdx == nAlpha
            refDisp = normalizeRGBForDisplay(refocused, pLow, pHigh);
            figure('Name', sprintf('Refocus alpha=%.3f', a));
            imshow(refDisp);
            title(sprintf('Refocus (weighted shift-sum) alpha=%.3f', a));
            drawnow;
        end
    end
end

%% ---------------- 6) ALL-IN-FOCUS (HARD) via ARGMAX ----------------
if doAllInFocus && ~isempty(refocusStack)
    fprintf("Building all-in-focus (HARD argmax) from refocus stack...\n");

    focusScore = zeros(Y, X, nAlpha);

    % Laplacian filter
    hLap = [0 1 0; 1 -4 1; 0 1 0];

    for aIdx = 1:nAlpha
        img = refocusStack(:,:,:,aIdx);

        % Luminance proxy
        Ylum = 0.2989*img(:,:,1) + 0.5870*img(:,:,2) + 0.1140*img(:,:,3);

        L = abs(imfilter(Ylum, hLap, 'replicate'));

        % Optional: smooth focus score locally for stability
        if focusWin > 1
            L = imboxfilt(L, focusWin);
        end

        focusScore(:,:,aIdx) = L;
    end

    % Winner alpha index per pixel
    [~, bestIdx] = max(focusScore, [], 3);   % [Y X], values in 1..nAlpha

    % Optional: stabilize bestIdx (reduces salt-and-pepper switching)
    % (median filter keeps discrete indices)
    bestIdx = medfilt2(bestIdx, [3 3], 'symmetric');

    % Build all-in-focus RGB by selecting refocusStack at bestIdx per pixel
    allInFocusHard = zeros(Y, X, 3);

    linPix = (1:(Y*X))';
    bestLin = bestIdx(:);

    for c = 1:3
        tmp = refocusStack(:,:,c,:);          % [Y X 1 nAlpha]
        tmp = reshape(tmp, Y*X, nAlpha);      % [Y*X, nAlpha]
        allInFocusHard(:,:,c) = reshape(tmp(sub2ind([Y*X, nAlpha], linPix, bestLin)), Y, X);
    end

    % Display (normalized for viewing)
    allInFocusDisp = normalizeRGBForDisplay(allInFocusHard, pLow, pHigh);
    figure('Name','All-in-focus (HARD argmax) fusion');
    imshow(allInFocusDisp);
    title('All-in-focus (HARD): per-pixel best alpha');

    % Show which alpha wins (map and also in physical alpha values)
    %figure('Name','Best alpha index map (HARD)');
    %imagesc(bestIdx); axis image off; colormap parula; colorbar;
    %title('Best alpha index per pixel');

    figure('Name','Best alpha value map (HARD)');
    imagesc(alphaList(bestIdx)); axis image off; colormap parula; colorbar;
    title('Best alpha value per pixel');
end

disp("Done.");

%% ===================== Helper functions =====================

function [RGBdisp, Wdisp] = normalizeForDisplay(RGB, W, pLow, pHigh)
    % RGB: [U V Y X 3], W: [U V Y X]
    RGBdisp = RGB;
    for c = 1:3
        tmp = RGB(:,:,:,:,c);
        lo = prctile(tmp(:), pLow);
        hi = prctile(tmp(:), pHigh);
        RGBdisp(:,:,:,:,c) = (tmp - lo) / max(eps, (hi - lo));
    end
    RGBdisp = max(0, min(1, RGBdisp));

    tmpW = W;
    lo = prctile(tmpW(:), pLow);
    hi = prctile(tmpW(:), pHigh);
    Wdisp = (tmpW - lo) / max(eps, (hi - lo));
    Wdisp = max(0, min(1, Wdisp));
end

function M = buildMosaic(LFuv)
    % LFuv: [U V Y X] -> mosaic [U*Y, V*X]
    sz = size(LFuv);
    U = sz(1); V = sz(2); Y = sz(3); X = sz(4);
    M = zeros(U*Y, V*X);
    for u = 1:U
        for v = 1:V
            tile = squeeze(LFuv(u,v,:,:));
            r0 = (u-1)*Y + 1;
            c0 = (v-1)*X + 1;
            M(r0:r0+Y-1, c0:c0+X-1) = tile;
        end
    end
end

function M = buildMosaicRGB(LFRGB)
    % LFRGB: [U V Y X 3] in [0,1]
    sz = size(LFRGB);
    U = sz(1); V = sz(2); Y = sz(3); X = sz(4);
    M = zeros(U*Y, V*X, 3);
    for u = 1:U
        for v = 1:V
            tile = squeeze(LFRGB(u,v,:,:,:)); % [Y X 3]
            r0 = (u-1)*Y + 1;
            c0 = (v-1)*X + 1;
            M(r0:r0+Y-1, c0:c0+X-1, :) = tile;
        end
    end
end

function RGBdisp = normalizeRGBForDisplay(RGB, pLow, pHigh)
    % RGB: [Y X 3] double
    RGBdisp = zeros(size(RGB));
    for c = 1:3
        tmp = RGB(:,:,c);
        lo = prctile(tmp(:), pLow);
        hi = prctile(tmp(:), pHigh);
        RGBdisp(:,:,c) = (tmp - lo) / max(eps, (hi - lo));
    end
    RGBdisp = max(0, min(1, RGBdisp));
end
