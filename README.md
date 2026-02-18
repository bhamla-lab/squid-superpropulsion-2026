# Squid-Inspired Superpropulsion

Data and code for: **"Squid-inspired superpropulsion"**

Choi, D., Singh, P., Bergerson, I., Kim, M., Park, J., Wallace, H.J., Zhang, K., Hsieh, S.Y., Asberry, A.T., Uyeno, T.A., Gilly, W.F., Park, H., Kang, D., Bose, C., & Bhamla, S.

## Repository Structure

```
├── data/
│   ├── fig1/                        # Main Figure 1: Squid funnel structure & jetting kinematics
│   ├── figS3/                       # Fig. S3: Chromatophore tracking measurements
│   ├── figS4/                       # Fig. S4: Funnel paralysis experiments
│   ├── figS6/                       # Fig. S6: PIV flow measurements
│   └── video/                       # Supplementary video data
├── code/
│   └── chromatophore_tracking/      # MATLAB pipeline for chromatophore-based deformation tracking
└── README.md
```

## Data

### Main Figure 1 — Funnel collagen architecture and phase-lagged recoil (`data/fig1/`)

| File | Panel | Description |
|------|-------|-------------|
| `Data_Fig01E_deformation_snapshots.txt` | 1E | Funnel and mantle deformation at three time points |
| `Data_Fig01G_funnel_mantle_width_pressure.xlsx` | 1G | Normalized funnel width, mantle width, and cavity pressure during a jet pulse (*D. opalescens*, condition $c_1$) |
| `Data_Fig01H_time_history.xlsx` | 1H | Time series of stimulus trigger, cavity pressure, mantle width, and funnel width |
| `Data_Fig01I_widening_across_conditions.xlsx` | 1I | Maximum funnel widening and minimum mantle width across conditions ($c_1$, $c_2$, $c_3$) |
| `Data_Fig01J_response_time_ratio.xlsx` | 1J | Funnel response-time ratio $\tau_\text{squid}/T_\text{squid}$ across conditions |

### Supplementary Figure S3 — Chromatophore tracking (`data/figS3/`)

| File | Panel | Description |
|------|-------|-------------|
| `Data_FigS3_GHIJ_chromatophore_tracking.xlsx` | S3 G–J | Chromatophore tracking results for immobilized *D. opalescens* |
| `Data_FigS3_K.xlsx` | S3 K | Additional chromatophore tracking metrics |
| `Data_FigS3_LMN.xlsx` | S3 L–N | PCA-based deformation analysis |
| `Data_FigS3_moving_squid.xlsx` | S3 | Additional data for freely moving *S. lessoniana* |

### Supplementary Figure S4 — Funnel paralysis (`data/figS4/`)

| File | Panel | Description |
|------|-------|-------------|
| `Data_FigS4_FG_funnel_paralysis.xlsx` | S4 F–G | Funnel deformation after nerve disruption |

### Supplementary Figure S6 — PIV flow measurements (`data/figS6/`)

| File | Panel | Description |
|------|-------|-------------|
| `Data_FigS6F_vortex_core.xlsx` | S6 F | Vortex core trajectories (rigid vs. compliant nozzle) |
| `Data_FigS6G_impulse_history.xlsx` | S6 G | Hydrodynamic impulse time history |
| `Data_FigS6H_entrainment.xlsx` | S6 H | Jet entrainment data |
| `Data_FigS6I_energy_history.xlsx` | S6 I | Jet kinetic energy time history |

### Supplementary Video data (`data/video/`)

| File | Description |
|------|-------------|
| `Data_Video_chromatophore_tracking.xlsx` | Chromatophore positions tracked from supplementary video |

## Code

### Chromatophore Tracking Pipeline (`code/chromatophore_tracking/`)

MATLAB scripts for quantifying mantle and funnel deformation from high-speed video of live squid by tracking chromatophore features (Fig. S3). Run sequentially:

| Step | Script | Description |
|------|--------|-------------|
| 1 | `step01_detect_black_dots.m` | Detect chromatophore features in a test frame |
| 2 | `step02_batch_convert_to_bw.m` | Batch convert images to binary |
| 3 | `step03_batch_detect_stationary_noise.m` | Identify and remove stationary artifacts |
| 4 | `step04_batch_track_dots.m` | Track chromatophore positions across frames |
| 5 | `step05_review_tracking.m` | Review and filter tracked trajectories |
| 6 | `step06_pca_rotation_analysis.m` | PCA-based width measurement of funnel and mantle |

**Requirements:** MATLAB with Image Processing Toolbox.

## Simulation Software

Fluid--structure interaction simulations used open-source tools:
- [OpenFOAM](https://www.openfoam.com/) — finite volume CFD solver
- [preCICE](https://precice.org/) — coupling library
- [CalculiX](http://www.calculix.de/) — finite element structural solver

## License

This data is provided for academic use accompanying the published manuscript.

## Contact

Saad Bhamla — saad.bhamla@colorado.edu
