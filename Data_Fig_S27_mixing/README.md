# Mixing Visualization - Data & Code

MATLAB codes and processed data for green dye mixing analysis (flexible vs rigid nozzles).

## Structure

```
data_sharing/
├── data/
│   ├── metrics/       .mat files (area, intensity, time per video)
│   └── snapshots/     .png frames (R30mm, F20mm - front & side)
├── code/
│   ├── C01_front_area_intensity.m        Front view: area & intensity vs time
│   ├── C02_side_area_intensity.m         Side view: area & intensity vs time
│   ├── C03_nozzle_length_comparison.m    Area at tau/T=400 vs nozzle length
│   ├── C04_front_boundary_overlay.m      Front view: rigid vs flexible boundary overlay
│   ├── C05_side_boundary_overlay.m       Side view: rigid vs flexible boundary overlay
│   └── image_processing/                 Upstream algorithms (reference only)
│       ├── C01_detect_green_start.m      Detect green dye injection start frame
│       ├── C02_batch_green_detection.m   Batch process all videos
│       ├── C10_extract_snapshots.m       Extract time-stamped snapshots
│       └── C20_detect_green_boundary.m   Detect green region boundary & compute metrics
```

## Usage

1. Open MATLAB, `cd` into `code/`
2. Run `C01` through `C05` directly — all paths are relative to `code/`
3. `C01`-`C03` use `data/metrics/`, `C04`-`C05` use `data/snapshots/`

## Naming Convention

- **Folder/file prefix**: `F_` = front view, `S_` = side view
- **Nozzle type**: `F##mm` = flexible (## mm length), `R30mm` = rigid
- **Case**: `C01`, `C02`, ... = repeated trials

## Normalization

- Time: tau/T = t * U_jet / D_nozzle (U_jet = 616.22 px/s, D_nozzle = 13.5 px)
- Area: A / A_nozzle
- Nozzle length: tau_nozzle / T
