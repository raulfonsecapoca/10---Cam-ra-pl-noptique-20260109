% Choose best LensletGridModel by sampling centers and measuring contrast
clear; close all; clc;

rawFile = "IMG_0001__frame.raw";
imgMeta = LFReadMetadata("IMG_0001__frame_metadata.json");

W = imgMeta.image.width;
H = imgMeta.image.height;
bitPacking = sprintf("%dbit", imgMeta.image.rawDetails.pixelPacking.bitsPerPixel);

Iraw = LFReadRaw(rawFile, bitPacking, [W H]); % uint16 mosaic
I = double(Iraw);

% Load the two candidate grid models
M1 = LFReadMetadata("data.C.3__C__T1CALIB__MOD_0059.grid.json"); % focus 650
M2 = LFReadMetadata("data.C.3__C__T1CALIB__MOD_0061.grid.json"); % focus 600
models = {M1, M2};
names  = ["0059","0061"];

scores = zeros(2,1);

for k = 1:2
    G = models{k}.LensletGridModel;

    % Generate predicted lenslet centers (u,v) in pixel coords
    [xc, yc] = lensletCenters(G);

    % Keep centers inside image bounds (avoid border)
    ok = xc > 2 & xc < (W-1) & yc > 2 & yc < (H-1);
    xc = xc(ok); yc = yc(ok);

    % Sample a small patch around each center and compute a "sharpness" score
    % (centers should land in consistent bright spots)
    vals = zeros(numel(xc),1);
    for i = 1:numel(xc)
        x = round(xc(i)); y = round(yc(i));
        patch = I(y-1:y+1, x-1:x+1);
        vals(i) = std(patch(:));  % local contrast
    end

    scores(k) = median(vals); % robust aggregate
    fprintf("%s score = %.4f\n", names(k), scores(k));
end

[~,best] = max(scores);
fprintf("\nBest model: %s\n", names(best));

%% ---- helper: compute hex-row-major-ish centers from LensletGridModel ----
function [xc, yc] = lensletCenters(G)
    % G fields: HSpacing, VSpacing, HOffset, VOffset, Rot, UMax, VMax, FirstPosShiftRow
    % We build a grid and apply rotation.
    UMax = G.UMax;
    VMax = G.VMax;

    hs = G.HSpacing; vs = G.VSpacing;
    x0 = G.HOffset;  y0 = G.VOffset;
    rot = G.Rot;

    [uu,vv] = meshgrid(0:UMax-1, 0:VMax-1);

    % Hex shift: every other row shifts by half horizontal spacing
    shiftRow = mod(vv + (G.FirstPosShiftRow-1), 2); % 0/1
    x = x0 + uu*hs + shiftRow*(hs/2);
    y = y0 + vv*vs;

    % Apply rotation around origin (0,0). For small rot this is fine.
    xr = x*cos(rot) - y*sin(rot);
    yr = x*sin(rot) + y*cos(rot);

    xc = xr(:);
    yc = yr(:);
end
