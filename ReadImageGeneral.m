%% Read RAW + demosaic using metadata (LFToolbox + MATLAB)
clear; close all; clc;

% --- Files ---
metadataName = "IMG_0005__frame_metadata.json";
rawFile      = "IMG_0005__frame.raw";

% --- Read metadata ---
meta = LFReadMetadata(metadataName);

% 1) Extract RAW format info
W = meta.image.width;
H = meta.image.height;

bitsPerPixel = meta.image.rawDetails.pixelPacking.bitsPerPixel;   % e.g., 12
bitPacking = sprintf("%dbit", bitsPerPixel);                      % "12bit"

% 2) Bayer pattern from metadata
% meta says: mosaic.tile = "r,gr:gb,b" and upperLeftPixel = "b"
% MATLAB demosaic patterns are: 'rggb','bggr','grbg','gbrg'
bayerPattern = upperLeftToMatlabPattern(string(meta.image.rawDetails.mosaic.upperLeftPixel));

% 3) Black/white levels (use per-channel; we'll apply on mosaic positions)
black = meta.image.rawDetails.pixelFormat.black;
white = meta.image.rawDetails.pixelFormat.white;

% 4) Read RAW (Bayer mosaic) using LFToolbox
Iraw = LFReadRaw(rawFile, bitPacking, [W H]);  % uint16 mosaic

figure; imagesc(Iraw); axis image off; colormap gray;
title(sprintf("RAW Bayer mosaic (%s, %dx%d)", bitPacking, W, H));

% 5) Black subtraction + normalize to full uint16 range (apply per CFA position)
I16 = blackWhiteNormalizeMosaic(Iraw, black, white, bayerPattern);

% 6) Simple demosaic (linear) as in PlenCal
RGB = demosaic(I16, bayerPattern);

% 7) Optional: apply white balance gains from metadata (quick + simple)
if isfield(meta.image, "color") && isfield(meta.image.color, "whiteBalanceGain")
    wb = meta.image.color.whiteBalanceGain;
    RGB = applyWhiteBalance(RGB, wb);
end

% 8) Optional: apply CCM cameraRGB -> sRGB (from metadata)
if isfield(meta.image, "color") && isfield(meta.image.color, "ccmRgbToSrgbArray")
    ccmVec = meta.image.color.ccmRgbToSrgbArray;
    CCM = reshape(ccmVec, [3 3]).';   % JSON usually row-major
    RGB = applyCCM_uint16(RGB, CCM);
end

% 9) Display (imshow expects uint16 scaled to 0..65535 already)
figure; imshow(RGB);
title("Demosaic using metadata (black/white + WB + CCM)");

%% ===================== helpers =====================

function pat = upperLeftToMatlabPattern(upperLeft)
    % Map metadata upper-left pixel to MATLAB demosaic pattern
    % upperLeft: "r"|"gr"|"gb"|"b"
    % Most important: "b" -> 'bggr' (your case).
    switch lower(strtrim(upperLeft))
        case "r"
            pat = "rggb";
        case "b"
            pat = "bggr";
        case "gr"
            pat = "grbg";
        case "gb"
            pat = "gbrg";
        otherwise
            error("Unknown upperLeftPixel: %s", upperLeft);
    end
end

function I16 = blackWhiteNormalizeMosaic(Iraw, black, white, bayerPattern)
    % Applies per-channel black/white normalization directly on the mosaic
    % and expands to full uint16 [0..65535] for correct demosaic/display.

    I = double(Iraw);
    [H,W] = size(I);

    % Create masks for mosaic positions according to Bayer pattern
    [mR, mG1, mG2, mB] = bayerMasks(H, W, bayerPattern);

    % Subtract black per position
    I(mR)  = I(mR)  - double(black.r);
    I(mG1) = I(mG1) - double(black.gr);
    I(mG2) = I(mG2) - double(black.gb);
    I(mB)  = I(mB)  - double(black.b);

    % Clamp to [0, white-black] per position
    I(mR)  = min(max(I(mR),  0), double(white.r  - black.r));
    I(mG1) = min(max(I(mG1), 0), double(white.gr - black.gr));
    I(mG2) = min(max(I(mG2), 0), double(white.gb - black.gb));
    I(mB)  = min(max(I(mB),  0), double(white.b  - black.b));

    % Expand each position to full 16-bit range
    I(mR)  = I(mR)  * (65535 / max(1, double(white.r  - black.r)));
    I(mG1) = I(mG1) * (65535 / max(1, double(white.gr - black.gr)));
    I(mG2) = I(mG2) * (65535 / max(1, double(white.gb - black.gb)));
    I(mB)  = I(mB)  * (65535 / max(1, double(white.b  - black.b)));

    I16 = uint16(I);
end

function [mR, mG1, mG2, mB] = bayerMasks(H, W, pattern)
    % Returns logical masks for the 2x2 Bayer pattern positions.
    % G1 and G2 correspond to the two greens.

    rr = false(H,W); cc = false(H,W);
    rr(1:2:end,:) = true;     % odd rows
    cc(:,1:2:end) = true;     % odd cols

    switch lower(pattern)
        case "rggb"
            % [R G1; G2 B] with R at (1,1)
            mR  = rr & cc;
            mG1 = rr & ~cc;
            mG2 = ~rr & cc;
            mB  = ~rr & ~cc;

        case "bggr"
            % [B G1; G2 R] with B at (1,1)
            mB  = rr & cc;
            mG1 = rr & ~cc;
            mG2 = ~rr & cc;
            mR  = ~rr & ~cc;

        case "grbg"
            % [G1 R; B G2] with G at (1,1)
            mG1 = rr & cc;
            mR  = rr & ~cc;
            mB  = ~rr & cc;
            mG2 = ~rr & ~cc;

        case "gbrg"
            % [G1 B; R G2] with G at (1,1)
            mG1 = rr & cc;
            mB  = rr & ~cc;
            mR  = ~rr & cc;
            mG2 = ~rr & ~cc;

        otherwise
            error("Unsupported Bayer pattern: %s", pattern);
    end
end

function RGBout = applyWhiteBalance(RGB, wb)
    % Apply simple WB gains on uint16 RGB.
    RGBd = double(RGB);
    RGBd(:,:,1) = RGBd(:,:,1) * double(wb.r);
    % Use average of gr/gb as green gain
    gGain = 0.5*(double(wb.gr) + double(wb.gb));
    RGBd(:,:,2) = RGBd(:,:,2) * gGain;
    RGBd(:,:,3) = RGBd(:,:,3) * double(wb.b);

    RGBd = min(max(RGBd, 0), 65535);
    RGBout = uint16(RGBd);
end

function RGBout = applyCCM_uint16(RGB, CCM)
    % Apply 3x3 CCM to uint16 RGB (assumes values already in 0..65535).
    RGBd = double(RGB) / 65535;
    [H,W,~] = size(RGBd);
    X = reshape(RGBd, [], 3);
    Y = X * CCM.';               % apply rows
    RGBd2 = reshape(Y, H, W, 3);
    RGBd2 = min(max(RGBd2, 0), 1);
    RGBout = uint16(RGBd2 * 65535);
end
