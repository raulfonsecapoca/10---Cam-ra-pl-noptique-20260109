%% Disparity / depth map from a 9x9 light field (simple multi-view block matching)
clear; close all; clc;

ImgName = 'IMG_0002';
matFile = ['Output\LF_from_' ImgName '.mat'];

maxDisp = 6;        % pixels (search range, try 4..12)
patchR  = 4;        % patch radius -> (2*patchR+1)^2, try 3..6
step    = 4;        % compute disparity every 'step' pixels (speed/quality trade)
useWeight = true;   % use LF(:,:,:,:,4) to mask unreliable areas (if present)
weightThresh = 0.05;

% --- Load LF ---
S = load(matFile);
LF = [];
names = fieldnames(S);
for k = 1:numel(names)
    v = S.(names{k});
    if isnumeric(v) && ndims(v)==5
        LF = double(v);
        fprintf("Using LF variable: %s\n", names{k});
        break;
    end
end
if isempty(LF), error("No 5D LF found."); end

[U,V,Y,X,C] = size(LF);
RGB = LF(:,:,:,:,1:3);

% --- Choose central view ---
uc = ceil(U/2); vc = ceil(V/2);

I0 = squeeze(RGB(uc,vc,:,:,:));  % [Y X 3]
I0g = rgb2gray01(I0);

% --- Pick a set of neighbor views (stereo baselines) ---
% Use a cross around the center to keep it simple
neighbors = [
    uc-2, vc;
    uc+2, vc;
    uc, vc-2;
    uc, vc+2;
    uc-2, vc-2;
    uc-2, vc+2;
    uc+2, vc-2;
    uc+2, vc+2
];
% keep only valid indices
neighbors = neighbors( ...
    neighbors(:,1)>=1 & neighbors(:,1)<=U & neighbors(:,2)>=1 & neighbors(:,2)<=V, :);

nN = size(neighbors,1);
fprintf("Using %d neighbor views.\n", nN);

% --- Optional reliability mask from weight channel ---
if useWeight && C>=4
    W0 = squeeze(LF(uc,vc,:,:,4));
    validMask = W0 > weightThresh;
else
    validMask = true(Y,X);
end

% --- Allocate disparity (dx,dy) and confidence ---
dxMap = nan(Y,X);
dyMap = nan(Y,X);
confMap = zeros(Y,X);

% --- Precompute padded reference for easy patch extraction ---
pad = patchR + maxDisp + 2;
I0p = padarray(I0g, [pad pad], 'replicate');

% --- Main loop on a grid (step) ---
ys = (1+patchR):(step):(Y-patchR);
xs = (1+patchR):(step):(X-patchR);

fprintf("Computing disparity on a %dx%d grid (step=%d)...\n", numel(ys), numel(xs), step);

for yi = 1:numel(ys)
    y = ys(yi);
    for xi = 1:numel(xs)
        x = xs(xi);

        if ~validMask(y,x)
            continue;
        end

        % reference patch
        pr = getPatch(I0p, y, x, pad, patchR);

        % accumulate cost over neighbors
        bestCost = inf;
        bestDx = 0; bestDy = 0;
        secondCost = inf;

        for dy = -maxDisp:maxDisp
            for dx = -maxDisp:maxDisp
                costSum = 0;

                for n = 1:nN
                    u = neighbors(n,1);
                    v = neighbors(n,2);

                    In = squeeze(RGB(u,v,:,:,:));
                    Ing = rgb2gray01(In);

                    % Patch in neighbor shifted by (dx,dy)
                    % NOTE: In LF, disparity direction depends on convention; this works as a practical match.
                    pn = getPatch(padarray(Ing,[pad pad],'replicate'), y+dy, x+dx, pad, patchR);

                    d = pr - pn;
                    costSum = costSum + mean(d(:).^2);  % SSD
                end

                if costSum < bestCost
                    secondCost = bestCost;
                    bestCost = costSum;
                    bestDx = dx;
                    bestDy = dy;
                elseif costSum < secondCost
                    secondCost = costSum;
                end
            end
        end

        % confidence: how much better is best vs second-best (bigger = more confident)
        conf = max(0, (secondCost - bestCost) / max(1e-12, secondCost));

        dxMap(y,x) = bestDx;
        dyMap(y,x) = bestDy;
        confMap(y,x) = conf;
    end
end

% --- Fill missing values + smooth ---
dxMap = inpaint_nans_simple(dxMap);
dyMap = inpaint_nans_simple(dyMap);

dxMap = imgaussfilt(dxMap, 1);
dyMap = imgaussfilt(dyMap, 1);
confMap = imgaussfilt(confMap, 1);

% --- Convert to "depth-like" map (relative): depth ~ 1 / disparity magnitude ---
dispMag = sqrt(dxMap.^2 + dyMap.^2);
depthRel = 1 ./ max(dispMag, 1e-3);

% --- Display ---
figure; imshow(normalize01(I0g)); title("Central view (grayscale)");

figure; imagesc(dxMap); axis image off; colorbar; title("Disparity dx (pixels)");
figure; imagesc(dyMap); axis image off; colorbar; title("Disparity dy (pixels)");
figure; imagesc(dispMag); axis image off; colorbar; title("Disparity magnitude (pixels)");
figure; imagesc(confMap); axis image off; colorbar; title("Confidence");

figure; imagesc(depthRel); axis image off; colorbar; title("Relative depth (1/|disp|)");


%% ---------- helper functions ----------
function g = rgb2gray01(I)
    % I: [Y X 3] arbitrary range -> grayscale double
    I = double(I);
    g = 0.2989*I(:,:,1) + 0.5870*I(:,:,2) + 0.1140*I(:,:,3);
end

function p = getPatch(Ipad, y, x, pad, r)
    % Extract patch centered at (y,x) from padded image
    yp = y + pad;
    xp = x + pad;
    p = Ipad(yp-r:yp+r, xp-r:xp+r);
end

function A = normalize01(A)
    A = double(A);
    lo = prctile(A(:), 1);
    hi = prctile(A(:), 99);
    A = (A - lo) / max(eps, (hi - lo));
    A = max(0, min(1, A));
end

function M = inpaint_nans_simple(M)
    % very simple fill: replace NaNs with nearest non-NaN (via scattered interpolant)
    [yy,xx] = ndgrid(1:size(M,1), 1:size(M,2));
    mask = ~isnan(M);
    F = scatteredInterpolant(xx(mask), yy(mask), M(mask), 'nearest', 'nearest');
    M(~mask) = F(xx(~mask), yy(~mask));
end
