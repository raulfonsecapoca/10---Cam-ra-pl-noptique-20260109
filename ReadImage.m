%% Minimal: Read RAW + show + simple demosaic (as in PlenCal)
clear; close all; clc;

metadataName = "IMG_0001__frame_metadata.json";
meta = LFReadMetadata(metadataName);
%fieldnames(meta);
%disp(meta);
%openvar meta;



% 2) Read RAW (Bayer mosaic)
rawFile = 'IMG_0001__frame.raw';  % <-- your file
bitPacking = "12bit";
imgSize = [3280 3280];            % Lytro F01 from your metadata

Iraw = LFReadRaw(rawFile, bitPacking, imgSize);  % uint16 mosaic

figure; imagesc(Iraw); axis image off; colormap gray;
title("RAW Bayer mosaic (uint16)");

black = 168; %black and white levels on metadata
white = 4095;

I = double(Iraw) - black;
I = max(0, min(I, white-black));

I16 = uint16(I * (65535/(white-black)));   % expande 12->16 bits
RGB = demosaic(I16, "bggr");

figure; imshow(RGB); title("Demosaic (12->16 expanded)");




%% 



DecodeOptions = struct();

DecodeOptions.WhiteProcDataFnameExtension = 'data.C.3__C__T1CALIB__MOD_0059.grid.json';

DecodeOptions.WhiteRawDataFnameExtension = 'data.C.3__C__T1CALIB__MOD_0059.RAW';

DecodeOptions.WhiteImageDatabasePath = '';

[LF, LFMetadata, WhiteImageMetadata, LensletGridModel, DecodeOptions] = ...
       LFLytroDecodeImage( rawFile, DecodeOptions );



