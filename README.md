# Light Field Pipeline (Lytro) — RAW → LF + Visualization / Refocus / All-in-Focus

This project contains two MATLAB scripts for processing **Lytro Light Field images**:
1. **Decoding** a lenslet `.raw` file into a **5-D Light Field (LF)** using a **specific white set** (`MOD_0061`), **without using a database**.
2. **Visualization and processing** of the LF, including sub-aperture mosaics, central view, **refocus sweep (shift-and-sum)** with weights, and **all-in-focus fusion**.

> **Required dependency:**
> **Light Field Toolbox Version 0.5.3.0 (18.2 MB)** by **Donald Dansereau**.

---

## Project structure (suggested)

project/
├─ Dataset test/
│ └─ LFToolbox0.5_Samples/
│ ├─ Images/
│ │ └─ F01/
│ │ ├─ IMG_0006__frame.raw
│ │ └─ IMG_0006__frame_metadata.json
│ └─ Cameras/
│ └─ A000424242/
│ └─ WhiteImages/
│ ├─ data.C.3__C__T1CALIB__MOD_0061.grid.json
│ ├─ data.C.3__C__T1CALIB__MOD_0061.RAW
│ └─ data.C.3__C__T1CALIB__MOD_0061.TXT
├─ Output/
│ └─ (generated .mat files)
├─ raw_to_lf_MOD_0061.m
└─ lf_visualize_refocus_allinfocus.m


Paths can be adjusted freely, as long as the scripts are updated accordingly.

---

## Installation Instructions

1. Download and extract **Light Field Toolbox v0.5.3.0**.
2. Add it to the MATLAB path:

```matlab
addpath(genpath('path/to/LFToolbox0.5'));
savepath;


