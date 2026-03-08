# Light Field Pipeline (Lytro) — RAW to LF + Visualization / Refocus / All-in-Focus

This project contains two MATLAB scripts for processing Lytro Light Field data:

1. **RAW to LF**: decodes a lenslet `.raw` image into a 5-D Light Field (LF).
2. **LF to Images**: loads the generated LF data and performs visualization and processing, including mosaics, central view extraction, refocus sweep, and all-in-focus fusion.

## Requirements

- MATLAB
- **Light Field Toolbox v0.5.3.0** by Donald Dansereau

Make sure the toolbox is installed and added to the MATLAB path before running the scripts.

## Dataset

The scripts are currently configured to use sample files stored in the `Dataset test` folder.

The pipeline uses the RAW image together with its metadata, white calibration files, and grid model data available in the dataset.

The input image is selected by changing the `ImgName` variable in both scripts.  
Available sample names are:

- `IMG_0001`
- `IMG_0002`
- `IMG_0003`
- `IMG_0004`
- `IMG_0005`
- `IMG_0006`

## Workflow

The scripts must be run in this order:

1. Run the **RAW to LF** script to generate the Light Field `.mat` file.
2. Run the **LF to Images** script to load that `.mat` file and produce the visualizations and processing results.

The first script saves an output file in the `Output` folder, for example:

`Output/LF_from_IMG_0005.mat`

This file is required by the second script.

## Notes

The scripts assume the folder structure and file paths defined in the code.  
You may need to adjust these paths depending on your local setup.