Fname = 'IMG_0001__frame.raw';
BitPacking = '12bit';
ImgSize = [381,382];

BuffSize = prod(ImgSize) * 12/8;



Img = LFReadRaw( Fname, BitPacking);

%imshow(Img);
imagesc(double(Img))

patterns = {'rggb','bggr','grbg','gbrg'};

figure;
for k = 1:numel(patterns)
    RGB = demosaic(Img, patterns{k});
    subplot(2,2,k);
    imshow(RGB, []);
    title(patterns{k});
end





%% 
%% RAW 12-bit packed (big-endian) -> Bayer -> RGB (WB + CCM + Gamma)
clear; close all; clc;

rawFile = "IMG_0001__frame.raw";   % <-- your file
W = 3280; H = 3280;

% Metadata
black = 168;
white = 4095;

wb = struct('r', 1.04296875, 'gr', 1, 'gb', 1, 'b', 1.26953125);

% CCM camera RGB -> sRGB (row-major from metadata)
CCM = [ 3.1115827560424805,  -1.9393929243087769, -0.1721898615360260;
       -0.3629055917263031,   1.6408803462982178, -0.2779748141765595;
        0.0789670124650002,  -1.1558042764663696,  2.0768373012542725];

gammaExp = 0.41666001081466675;

% Bayer pattern from metadata:
% tile "r,gr:gb,b" with upperLeftPixel "b" => top-left is Blue.
% That corresponds to BGGR in MATLAB demosaic naming.
bayerPattern = 'bggr';

%% 1) Read packed 12-bit big-endian
I12 = readPacked12BE(rawFile, W, H);   % uint16 [H x W]

%% 2) Black/white normalization to [0,1] (still mosaiced)
I = double(I12);
I = max(0, I - black);
I = min(I, white - black);
I = I / (white - black);

figure; imagesc(I); axis image off; colormap gray;
title("Bayer mosaic (normalized)");

%% 3) Demosaic
% demosaic expects integer types usually; feed uint16 in full range
Iu = uint16(I * 65535);
RGB = demosaic(Iu, bayerPattern);   % uint16 RGB

% convert to double [0,1]
RGB = double(RGB) / 65535;

%% 4) White balance (apply per channel)
RGB(:,:,1) = RGB(:,:,1) * wb.r;
RGB(:,:,2) = RGB(:,:,2) * 0.5*(wb.gr + wb.gb);
RGB(:,:,3) = RGB(:,:,3) * wb.b;

%% 5) Apply CCM (camera RGB -> sRGB)
RGBccm = applyCCM(RGB, CCM);

%% 6) Clamp + Gamma for display
RGBccm = max(0, min(1, RGBccm));
RGBout = RGBccm .^ gammaExp;

figure; imshow(RGBout);
title("RGB (demosaic + WB + CCM + gamma)");

%% ------------ helpers ------------
function I = readPacked12BE(fname, W, H)
    % Reads 12-bit packed big-endian stream: 2 pixels in 3 bytes.
    % Returns uint16 image [H x W] with values 0..4095 (nominal).

    nPix = W * H;
    nBytes = nPix * 12 / 8;  % 1.5 bytes per pixel
    if mod(nBytes,1) ~= 0
        error("Byte count is not integer; check W/H.");
    end
    nBytes = uint64(nBytes);

    fid = fopen(fname, 'r');
    if fid < 0, error("Cannot open file."); end
    buf = fread(fid, nBytes, 'uint8=>uint8');
    fclose(fid);

    if numel(buf) ~= nBytes
        error("Read %d bytes, expected %d.", numel(buf), nBytes);
    end

    % Unpack: 3 bytes -> 2 pixels
    b0 = uint16(buf(1:3:end));
    b1 = uint16(buf(2:3:end));
    b2 = uint16(buf(3:3:end));

    % Common 12-bit packing (big-endian style):
    % p0 = (b0<<4) | (b1>>4)
    % p1 = ((b1 & 0x0F)<<8) | b2
    p0 = bitshift(b0, 4) + bitshift(b1, -4);
    p1 = bitshift(bitand(b1, 15), 8) + b2;

    % Interleave pixels
    out = zeros(2*numel(p0), 1, 'uint16');
    out(1:2:end) = p0;
    out(2:2:end) = p1;

    if numel(out) < nPix
        error("Unpack produced fewer pixels than expected.");
    end
    out = out(1:nPix);

    I = reshape(out, [W, H])'; % -> [H x W]
end

function RGB2 = applyCCM(RGB, CCM)
    % Apply 3x3 color correction matrix to an RGB image in [0,1]
    [H,W,~] = size(RGB);
    X = reshape(RGB, [], 3);
    Y = X * CCM.';              % apply rows to RGB vectors
    RGB2 = reshape(Y, H, W, 3);
end




