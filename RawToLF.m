%% RAW -> LF using a specific white set (MOD_0061), no database
clear; close all; clc;

ImgName = 'IMG_0005';

% --- Your capture ---
rawFile  = ['Dataset test\LFToolbox0.5_Samples\Images\F01\' ImgName '__frame.raw'];
metaFile = ['Dataset test\LFToolbox0.5_Samples\Images\F01\' ImgName '__frame_metadata.json'];



% --- White set you have ---
whiteGridFile = 'Dataset test\LFToolbox0.5_Samples\Cameras\A000424242\WhiteImages\data.C.3__C__T1CALIB__MOD_0061.grid.json';
whiteRawFile  = 'Dataset test\LFToolbox0.5_Samples\Cameras\A000424242\WhiteImages\data.C.3__C__T1CALIB__MOD_0061.RAW';
whiteMetaFile = 'Dataset test\LFToolbox0.5_Samples\Cameras\A000424242\WhiteImages\data.C.3__C__T1CALIB__MOD_0061.TXT';

% ---- Read LF metadata (your capture) ----
metaLF = LFReadMetadata(metaFile);

% ---- Read lenslet RAW (your capture) ----
W = metaLF.image.width;
H = metaLF.image.height;
LensletImage = LFReadRaw(rawFile, '12bit', [W H]);

% ---- Read white metadata (TXT: nested structure) ----
metaWhiteWhole = LFReadMetadata(whiteMetaFile);
FA = metaWhiteWhole.master.picture.frameArray;
if iscell(FA)
    frame1 = FA{1}.frame;
else
    frame1 = FA(1).frame;
end
metaWhite = frame1.metadata;

% ---- Read white RAW ----
Ww = metaWhite.image.width;
Hw = metaWhite.image.height;
WhiteImage = LFReadRaw(whiteRawFile, '12bit', [Ww Hw]);

% ---- Read grid model from processed white ----
gridData = LFReadMetadata(whiteGridFile);
LensletGridModel = gridData.LensletGridModel;

% ---- DecodeOptions (mirror what LFLytroDecodeImage sets for F01) ----
DecodeOptions = struct();
DecodeOptions.DemosaicOrder = 'bggr';

DecodeOptions.LevelLimits = [ ...
    metaLF.image.rawDetails.pixelFormat.black.gr, ...
    metaLF.image.rawDetails.pixelFormat.white.gr ];

DecodeOptions.ColourMatrix = reshape(metaLF.image.color.ccmRgbToSrgbArray, 3, 3);
DecodeOptions.ColourBalance = [ ...
    metaLF.image.color.whiteBalanceGain.r, ...
    metaLF.image.color.whiteBalanceGain.gb, ...
    metaLF.image.color.whiteBalanceGain.b ];

DecodeOptions.Gamma = metaLF.image.color.gamma^0.5;

% ---- Decode ----
fprintf('Decoding with LFDecodeLensletImageDirect...\n');
[LF, LFWeight, DecodeOptionsOut] = LFDecodeLensletImageDirect( ...
    LensletImage, WhiteImage, LensletGridModel, DecodeOptions );

LF(:,:,:,:,4) = LFWeight;

% ---- Display ----
if exist('LFDispTiles','file')
    LFDispTiles(LF);
else
    % fallback: show central subaperture (common quick sanity check)
    % (indexing order is [j,i,l,k,chan] in Dansereau's docs)
    midJ = round(size(LF,1)/2);
    midI = round(size(LF,2)/2);
    RGBc = squeeze(LF(midJ, midI, :, :, 1:3));
    figure; imshow(uint16(RGBc), []);
    title('Central subaperture (fallback)');
end

outputPath = ['Output\LF_from_' ImgName '.mat'];
save(outputPath, 'LF', 'DecodeOptionsOut', 'LensletGridModel', '-v7.3');
