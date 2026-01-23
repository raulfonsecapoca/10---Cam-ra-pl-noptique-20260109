%% Fast disparity from LF (grid + 1D search in dx and dy)
clear; close all; clc;

ImgName = 'IMG_0002';
matFile = ['Output\LF_from_' ImgName '.mat'];

maxDisp = 10;      % pixels, try 6..15
patchR  = 4;       % patch radius (9x9)
step    = 6;       % BIG speed lever: 4..10
useWeight = true;
weightThresh = 0.05;

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

uc = ceil(U/2); vc = ceil(V/2);

I0 = squeeze(RGB(uc,vc,:,:,:));
I0g = rgb2gray01(I0);

if useWeight && C>=4
    W0 = squeeze(LF(uc,vc,:,:,4));
    validMask = W0 > weightThresh;
else
    validMask = true(Y,X);
end

% Choose only 4 neighbor views (fast)
Lview = squeeze(RGB(uc, max(1,vc-2), :,:,:));  % left
Rview = squeeze(RGB(uc, min(V,vc+2), :,:,:));  % right
Uview = squeeze(RGB(max(1,uc-2), vc, :,:,:));  % up
Dview = squeeze(RGB(min(U,uc+2), vc, :,:,:));  % down

Lg = rgb2gray01(Lview);
Rg = rgb2gray01(Rview);
Ug = rgb2gray01(Uview);
Dg = rgb2gray01(Dview);

pad = patchR + maxDisp + 2;
I0p = padarray(I0g, [pad pad], 'replicate');
Lp  = padarray(Lg,  [pad pad], 'replicate');
Rp  = padarray(Rg,  [pad pad], 'replicate');
Up  = padarray(Ug,  [pad pad], 'replicate');
Dp  = padarray(Dg,  [pad pad], 'replicate');

ys = (1+patchR):(step):(Y-patchR);
xs = (1+patchR):(step):(X-patchR);

dxGrid = nan(numel(ys), numel(xs));
dyGrid = nan(numel(ys), numel(xs));
confGrid = zeros(numel(ys), numel(xs));

fprintf("Computing disparity on grid %dx%d (step=%d)...\n", numel(ys), numel(xs), step);

for iy = 1:numel(ys)
    y = ys(iy);
    for ix = 1:numel(xs)
        x = xs(ix);

        if ~validMask(y,x), continue; end

        pref = getPatch(I0p, y, x, pad, patchR);

        % ---- dx: match to left/right (1D search) ----
        bestCost = inf; bestDx = 0; second = inf;
        for dx = -maxDisp:maxDisp
            pL = getPatch(Lp, y, x+dx, pad, patchR);
            pR = getPatch(Rp, y, x-dx, pad, patchR); % opposite shift
            dL = pref - pL; dR = pref - pR;
            cost = mean(dL(:).^2) + mean(dR(:).^2);

            if cost < bestCost
                second = bestCost;
                bestCost = cost;
                bestDx = dx;
            elseif cost < second
                second = cost;
            end
        end
        confX = max(0, (second - bestCost) / max(1e-12, second));

        % ---- dy: match to up/down (1D search) ----
        bestCost = inf; bestDy = 0; second = inf;
        for dy = -maxDisp:maxDisp
            pU = getPatch(Up, y+dy, x, pad, patchR);
            pD = getPatch(Dp, y-dy, x, pad, patchR);
            dU = pref - pU; dD = pref - pD;
            cost = mean(dU(:).^2) + mean(dD(:).^2);

            if cost < bestCost
                second = bestCost;
                bestCost = cost;
                bestDy = dy;
            elseif cost < second
                second = cost;
            end
        end
        confY = max(0, (second - bestCost) / max(1e-12, second));

        dxGrid(iy,ix) = bestDx;
        dyGrid(iy,ix) = bestDy;
        confGrid(iy,ix) = 0.5*(confX+confY);
    end
end

% Interpolate grid to full resolution
[Xq,Yq] = meshgrid(1:X, 1:Y);
[Xg,Yg] = meshgrid(xs, ys);

dxMap = interp2(Xg, Yg, dxGrid, Xq, Yq, 'linear', 0);
dyMap = interp2(Xg, Yg, dyGrid, Xq, Yq, 'linear', 0);
confMap = interp2(Xg, Yg, confGrid, Xq, Yq, 'linear', 0);

% Smooth a bit
dxMap = imgaussfilt(dxMap, 1);
dyMap = imgaussfilt(dyMap, 1);
confMap = imgaussfilt(confMap, 1);

dispMag = hypot(dxMap, dyMap);
depthRel = 1 ./ max(dispMag, 1e-3);

figure; imshow(normalize01(I0g)); title("Central view (grayscale)");
figure; imagesc(dxMap); axis image off; colorbar; title("dx (fast)");
figure; imagesc(dyMap); axis image off; colorbar; title("dy (fast)");
figure; imagesc(dispMag); axis image off; colorbar; title("|disp| (fast)");
figure; imagesc(confMap); axis image off; colorbar; title("confidence (fast)");
figure; imagesc(depthRel); axis image off; colorbar; title("relative depth ~ 1/|disp| (fast)");

%% ---- helpers ----
function g = rgb2gray01(I)
I = double(I);
g = 0.2989*I(:,:,1) + 0.5870*I(:,:,2) + 0.1140*I(:,:,3);
end

function p = getPatch(Ipad, y, x, pad, r)
yp = y + pad; xp = x + pad;
p = Ipad(yp-r:yp+r, xp-r:xp+r);
end

function A = normalize01(A)
A = double(A);
lo = prctile(A(:), 1); hi = prctile(A(:), 99);
A = (A - lo) / max(eps, (hi - lo));
A = max(0, min(1, A));
end
