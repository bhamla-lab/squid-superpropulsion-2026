# Squid-Inspired Superpropulsion: Data and Code Repository

Source data and analysis code for:

**"Squid-inspired superpropulsion"**
Choi, Singh, Bergerson, Kim, Park, Wallace, Zhang, Hsieh, Asberry, Uyeno, Gilly, Park, Kang, Bose & Bhamla

---

## Figure-to-Data Lookup

Every figure panel from the main text, Supplementary Information, and Supplementary Movie is listed below. `N/A` indicates panels that are schematics, photographs, or microscopy images with no associated source-data file.

### Main-Text Figures

| Figure | Panel | Description | Source Data | Code |
|--------|-------|-------------|-------------|------|
| **Fig. 1** | **A** | Squid mantle cavity and funnel photograph | N/A (image; credit: Michael Patrick O'Neill) | -- |
| | **B** | Mantle contraction and funnel expansion image | N/A (image) | -- |
| | **C** | Picrosirius Red staining, mantle (polarized light) | N/A (microscopy image) | -- |
| | **D** | Picrosirius Red staining, funnel (polarized light) | N/A (microscopy image) | -- |
| | **E** | Snapshots: initial state, peak funnel widening, recoil | [Data_Fig_01_E.txt](Data_Fig_01_E.txt) | -- |
| | **F** | Conceptual schematic of expansion--recoil sequence | N/A (schematic) | -- |
| | **G** | Normalized W_f\*, W_m\*, and pressure during jet pulse (c1) | [Data_Fig_01_G.xlsx](Data_Fig_01_G.xlsx) | [Code_Fig_S03_chromatophore_tracking_algorithm/](Code_Fig_S03_chromatophore_tracking_algorithm/) |
| | **H** | Time history: stimulus, pressure, W_m\*, W_f\* | [Data_Fig_01_H.xlsx](Data_Fig_01_H.xlsx) | [Code_Fig_S03_chromatophore_tracking_algorithm/](Code_Fig_S03_chromatophore_tracking_algorithm/) |
| | **I** | W_f,max\* and W_m,min\* for c1, c2, c3 | [Data_Fig_01_I.xlsx](Data_Fig_01_I.xlsx) | -- |
| | **J** | Funnel response-time ratio tau_squid/T_squid | [Data_Fig_01_J.xlsx](Data_Fig_01_J.xlsx) | -- |
| **Fig. 2** | **A** | PIV snapshots, rigid nozzle | N/A (derived from PIV; see Fig. S6) | -- |
| | **B** | PIV snapshots, compliant nozzle + nozzle deformation | N/A (derived from PIV; see Fig. S6) | -- |
| | **C** | 3D-CFD snapshots, rigid nozzle | N/A (CFD simulation output) | -- |
| | **D** | 3D-CFD snapshots, compliant nozzle | N/A (CFD simulation output) | -- |
| | **E** | Impulse ratio I_h_hat vs tau/T (exp + CFD + theory) | [Data_Fig_02_E.xlsx](Data_Fig_02_E.xlsx) | [Code_Fig_S11_superpropulsion_theory/](Code_Fig_S11_superpropulsion_theory/) |
| | **F** | Schematic: elastic strain power during store--release | N/A (schematic) | -- |
| | **G** | Spring--mass resonator analogy | N/A (schematic) | -- |
| | **H** | Time traces: elastic strain power and jet kinetic power | [Data_Fig_S06_I__energy_history_data.xlsx](Data_Fig_S06_I__energy_history_data.xlsx) | -- |
| | **I** | Measured tau/T vs theoretical estimates | [Data_Fig_02_I.xlsx](Data_Fig_02_I.xlsx) | [Code_Fig_S11_superpropulsion_theory/](Code_Fig_S11_superpropulsion_theory/) |
| | **J** | Schematic: stiff vs optimal vs compliant nozzle recoil | N/A (schematic) | -- |
| | **K** | Exit velocity ratio v_hat vs time (theory + exp) | [Data_Fig_S11_D,E,F.xlsx](Data_Fig_S11_D,E,F.xlsx) | [S01_velocity_history.m](Code_Fig_S11_superpropulsion_theory/S01_velocity_history.m) |
| | **L** | Impulse ratio I_h_hat vs time (theory + exp) | [Data_Fig_S11_D,E,F.xlsx](Data_Fig_S11_D,E,F.xlsx) | [S02_impulse_history.m](Code_Fig_S11_superpropulsion_theory/S02_impulse_history.m) |
| | **M** | Kinetic energy ratio E_k_hat vs time (theory + exp) | [Data_Fig_S11_D,E,F.xlsx](Data_Fig_S11_D,E,F.xlsx) | [S03_energy_history.m](Code_Fig_S11_superpropulsion_theory/S03_energy_history.m) |
| **Fig. 3** | **A** | Flying squid photograph | N/A (image; credit: Anthony Pierce) | -- |
| | **B** | Single-pulse aerial jet: 110% height gain | [Data_Fig_S15_B.xlsx](Data_Fig_S15_B.xlsx) | -- |
| | **C** | Input waveforms (single + repeated pulse) | [Data_Fig_S14_D.xlsx](Data_Fig_S14_D.xlsx) | -- |
| | **D** | Time sequence: jet development and nozzle deformation | N/A (image sequence) | -- |
| | **E** | Candle extinguishing range + averaged jet trajectories (inset) | [Data_Fig_03_E.xlsx](Data_Fig_03_E.xlsx), [Data_Fig_03_E,F/](Data_Fig_03_E,F/) | [Data_Fig_03_E,F/](Data_Fig_03_E,F/) |
| | **F** | Candle extinguishing comparison (rigid vs compliant) | [Data_Fig_03_E,F/](Data_Fig_03_E,F/) | [Data_Fig_03_E,F/](Data_Fig_03_E,F/) |
| | **G** | Nozzle set (heights 5--30 mm) photograph | N/A (photograph) | -- |
| | **H** | Squid-boat with compliant nozzle | N/A (photograph) | -- |
| | **I** | Boat translation, rigid nozzle | [Data_Fig_S16_D,E.../](Data_Fig_S16_D,E_largeBoat_tau_T_dependency/) | [chatsweep_poly.m](Data_Fig_S16_D,E_largeBoat_tau_T_dependency/chatsweep_poly.m) |
| | **J** | Boat translation, compliant nozzle | [Data_Fig_S16_D,E.../](Data_Fig_S16_D,E_largeBoat_tau_T_dependency/) | [chatsweep_poly.m](Data_Fig_S16_D,E_largeBoat_tau_T_dependency/chatsweep_poly.m) |
| | **K** | Squid inking photograph | N/A (image; credit: Blue Planet Archive) | -- |
| | **L** | Dye visualization, front view (compliant vs rigid) | [Data_Fig_S18_mixing/data/snapshots/](Data_Fig_S18_mixing/data/snapshots/) | [C04_front_boundary_overlay.m](Data_Fig_S18_mixing/code/C04_front_boundary_overlay.m) |
| | **M** | Dye visualization, side view (compliant vs rigid) | [Data_Fig_S18_mixing/data/snapshots/](Data_Fig_S18_mixing/data/snapshots/) | [C05_side_boundary_overlay.m](Data_Fig_S18_mixing/code/C05_side_boundary_overlay.m) |
| | **N** | Horizontal jet distance x/D and power P/P_0 vs tau/T | [Data_Fig_S15_E.xlsx](Data_Fig_S15_E.xlsx) | -- |
| | **O** | Boat speed and CoT vs tau/T | [Data_Fig_S16_D,E.../](Data_Fig_S16_D,E_largeBoat_tau_T_dependency/) | [chatsweep_poly.m](Data_Fig_S16_D,E_largeBoat_tau_T_dependency/chatsweep_poly.m) |
| | **P** | Dye dispersion area A/A_nozzle vs tau/T | [Data_Fig_S18_mixing/data/metrics/](Data_Fig_S18_mixing/data/metrics/) | [C03_tau_T_comparison.m](Data_Fig_S18_mixing/code/C03_tau_T_comparison.m) |
| | **Q** | Bar chart: mean horizontal jet distance +44.6% | Derived from panel N data | -- |
| | **R** | Bar chart: boat speed +41.1% | Derived from panel O data | -- |
| | **S** | Bar chart: acceleration +43.2% | Derived from panel O data | -- |
| | **T** | Bar chart: CoT -27.9% | Derived from panel O data | -- |
| | **U** | Bar chart: dye dispersion +40.5% | Derived from panel P data | -- |

### Supplementary Figures

| Figure | Panel | Description | Source Data | Code |
|--------|-------|-------------|-------------|------|
| **Fig. S1** | **A** | Monterey Bay collection site | N/A (photograph) | -- |
| | **B** | Circulating seawater tank, Hopkins Marine Station | N/A (photograph) | -- |
| | **C** | Gijang, Korea collection site | N/A (photograph) | -- |
| | **D** | Circulating seawater tank, Gijang Bay | N/A (photograph) | -- |
| **Fig. S2** | **A** | Fixed squid photograph | N/A (photograph) | -- |
| | **B** | Cross-section of funnel wall (Milligan's Trichrome) | N/A (microscopy image) | -- |
| | **C** | Funnel: H&E, Masson's, PSR brightfield + polarized | N/A (microscopy images) | -- |
| | **D** | Mantle: H&E, Masson's, PSR brightfield + polarized | N/A (microscopy images) | -- |
| **Fig. S3** | **A** | High-speed imaging setup schematic | N/A (schematic) | -- |
| | **B** | Raw high-speed image with chromatophores | N/A (image) | [Code_Fig_S03_chromatophore_tracking_algorithm/](Code_Fig_S03_chromatophore_tracking_algorithm/) |
| | **C** | Close-up: electrode and pressure sensor | N/A (image) | -- |
| | **D** | Image processing: intensity normalization | N/A (image) | [Code_Fig_S03_01_review_detect_black_dots.m](Code_Fig_S03_chromatophore_tracking_algorithm/Code_Fig_S03_01_review_detect_black_dots.m) |
| | **E** | Image processing: binarization | N/A (image) | [Code_Fig_S03_02_batch_convert_to_bw.m](Code_Fig_S03_chromatophore_tracking_algorithm/Code_Fig_S03_02_batch_convert_to_bw.m) |
| | **F** | Image processing: particle tracking with ROI | N/A (image) | [Code_Fig_S03_04_batch_track_dots.m](Code_Fig_S03_chromatophore_tracking_algorithm/Code_Fig_S03_04_batch_track_dots.m) |
| | **G** | Electrical stimulus timing | [Data_Fig_S03_G,H,I,J.xlsx](Data_Fig_S03_G,H,I,J.xlsx) | [Code_Fig_S03_chromatophore_tracking_algorithm/](Code_Fig_S03_chromatophore_tracking_algorithm/) |
| | **H** | Pressure sensor output | [Data_Fig_S03_G,H,I,J.xlsx](Data_Fig_S03_G,H,I,J.xlsx) | [Code_Fig_S03_chromatophore_tracking_algorithm/](Code_Fig_S03_chromatophore_tracking_algorithm/) |
| | **I** | Funnel width variation | [Data_Fig_S03_G,H,I,J.xlsx](Data_Fig_S03_G,H,I,J.xlsx) | [Code_Fig_S03_chromatophore_tracking_algorithm/](Code_Fig_S03_chromatophore_tracking_algorithm/) |
| | **J** | Mantle width variation | [Data_Fig_S03_G,H,I,J.xlsx](Data_Fig_S03_G,H,I,J.xlsx) | [Code_Fig_S03_chromatophore_tracking_algorithm/](Code_Fig_S03_chromatophore_tracking_algorithm/) |
| | **K** | Selected strong jets (N=30, n=8) | [Data_Fig_S03_K.xlsx](Data_Fig_S03_K.xlsx) | -- |
| | **L** | Averaged pressure profiles | [Data_Fig_S03_L,M,N.xlsx](Data_Fig_S03_L,M,N.xlsx) | -- |
| | **M** | Normalized funnel deformation | [Data_Fig_S03_L,M,N.xlsx](Data_Fig_S03_L,M,N.xlsx) | -- |
| | **N** | Normalized mantle deformation | [Data_Fig_S03_L,M,N.xlsx](Data_Fig_S03_L,M,N.xlsx) | -- |
| | -- | Additional data: freely moving *S. lessoniana* | [Data_Fig_S03_additional_data_for_moving_squid.xlsx](Data_Fig_S03_additional_data_for_moving_squid.xlsx) | -- |
| **Fig. S4** | **A** | Schematic of squid nervous system | N/A (schematic) | -- |
| | **B** | Nerve infundibulum, normal condition | N/A (microscope image) | -- |
| | **C** | Nerve infundibulum, paralyzed condition | N/A (microscope image) | -- |
| | **D** | Head-first reflex (SLEAP tracking) | N/A (tracking visualization) | -- |
| | **E** | Tail-first reflex (SLEAP tracking) | N/A (tracking visualization) | -- |
| | **F** | Directional velocity, intact squid (N=45) | [Data_Fig_S04_F,G.xlsx](Data_Fig_S04_F,G.xlsx) | -- |
| | **G** | Directional velocity, funnel-paralyzed squid (N=34) | [Data_Fig_S04_F,G.xlsx](Data_Fig_S04_F,G.xlsx) | -- |
| **Fig. S5** | **A** | Flange-clamped fabrication (thick nozzles) | N/A (schematic + photograph) | -- |
| | **B** | Flange-embedded fabrication (thin nozzles) | N/A (schematic + photograph) | -- |
| | **C** | Press-fit fabrication | N/A (schematic + photograph) | -- |
| **Fig. S6** | **A** | Single pulsed jet generator setup | N/A (schematic) | -- |
| | **B** | Overlaid raw images of jet and nozzle deformation | N/A (image) | -- |
| | **C** | Vorticity field, rigid nozzle (Eh = infinity) | N/A (PIV-derived visualization) | -- |
| | **D** | Vorticity field, Eh = 14.4 N/m | N/A (PIV-derived visualization) | -- |
| | **E** | Vorticity field, Eh = 7.0 N/m | N/A (PIV-derived visualization) | -- |
| | **F** | Vortex core travel distance vs t/T | [Data_Fig_S06_F_vortex_core_data.xlsx](Data_Fig_S06_F_vortex_core_data.xlsx) | -- |
| | **G** | Hydrodynamic impulse I_h vs t/T | [Data_Fig_S06_G_impulse_history_data.xlsx](Data_Fig_S06_G_impulse_history_data.xlsx) | -- |
| | **H** | Cumulative entrainment V/V_nozzle vs t/T | [Data_Fig_S06_H__entrainment_data.xlsx](Data_Fig_S06_H__entrainment_data.xlsx) | -- |
| | **I** | Energy partitioning: E_elastic and E_kinetic vs t/T | [Data_Fig_S06_I__energy_history_data.xlsx](Data_Fig_S06_I__energy_history_data.xlsx) | -- |
| **Fig. S7** | **A** | Thrust measurement setup | N/A (schematic) | -- |
| | **B** | Thrust force time history | [Data_Fig_S07_C.xlsx](Data_Fig_S07_C.xlsx) | -- |
| | **C** | Thrust force vs tau/T | [Data_Fig_S07_C.xlsx](Data_Fig_S07_C.xlsx) | -- |
| **Fig. S8** | -- | Schematic of 1D tube model | N/A (schematic) | -- |
| **Fig. S9** | **A** | Spring--damper--mass resonator analogy | N/A (schematic) | -- |
| | **B** | Electrical RLC resonator analogy | N/A (schematic) | -- |
| | **C** | Flexible nozzle system analogy | N/A (schematic) | -- |
| **Fig. S10** | **A** | Theoretical exit velocity: base vs collapse model | Theory output | [S01_velocity_history.m](Code_Fig_S11_superpropulsion_theory/S01_velocity_history.m) |
| | **B** | Velocity fields, rigid nozzle (t/T = 1.0, 2.0, 3.2) | N/A (PIV visualization) | -- |
| | **C** | Velocity fields, flexible nozzle (t/T = 1.0, 2.0, 3.2) | N/A (PIV visualization) | -- |
| | **D** | Theory vs exp: exit velocity, tau/T = 0.17 | [Data_Fig_S11_D,E,F.xlsx](Data_Fig_S11_D,E,F.xlsx) | [S01_velocity_history.m](Code_Fig_S11_superpropulsion_theory/S01_velocity_history.m) |
| | **E** | Theory vs exp: exit velocity, tau/T = 0.29 | [Data_Fig_S11_D,E,F.xlsx](Data_Fig_S11_D,E,F.xlsx) | [S01_velocity_history.m](Code_Fig_S11_superpropulsion_theory/S01_velocity_history.m) |
| | **F** | Theory vs exp: exit velocity, tau/T = 0.55 | [Data_Fig_S11_D,E,F.xlsx](Data_Fig_S11_D,E,F.xlsx) | [S01_velocity_history.m](Code_Fig_S11_superpropulsion_theory/S01_velocity_history.m) |
| **Fig. S11** | **A** | Exit velocity v_hat time histories (3 tau/T values) | [Data_Fig_S11_A,B,C.xlsx](Data_Fig_S11_A,B,C.xlsx) | [S01_velocity_history.m](Code_Fig_S11_superpropulsion_theory/S01_velocity_history.m) |
| | **B** | Impulse I_h_hat time histories (3 tau/T values) | [Data_Fig_S11_A,B,C.xlsx](Data_Fig_S11_A,B,C.xlsx) | [S02_impulse_history.m](Code_Fig_S11_superpropulsion_theory/S02_impulse_history.m) |
| | **C** | Kinetic energy E_k_hat time histories (3 tau/T values) | [Data_Fig_S11_A,B,C.xlsx](Data_Fig_S11_A,B,C.xlsx) | [S03_energy_history.m](Code_Fig_S11_superpropulsion_theory/S03_energy_history.m) |
| | **D** | Peak v_hat vs tau/T (exp + CFD + theory) | [Data_Fig_S11_D,E,F.xlsx](Data_Fig_S11_D,E,F.xlsx) | [S01_velocity_history.m](Code_Fig_S11_superpropulsion_theory/S01_velocity_history.m) |
| | **E** | Peak I_h_hat vs tau/T | [Data_Fig_S11_D,E,F.xlsx](Data_Fig_S11_D,E,F.xlsx) | [S02_impulse_history.m](Code_Fig_S11_superpropulsion_theory/S02_impulse_history.m) |
| | **F** | Peak E_k_hat vs tau/T | [Data_Fig_S11_D,E,F.xlsx](Data_Fig_S11_D,E,F.xlsx) | [S03_energy_history.m](Code_Fig_S11_superpropulsion_theory/S03_energy_history.m) |
| **Fig. S12** | **a** | Computational domain and boundary conditions | N/A (schematic) | -- |
| | **b** | Fluid domain mesh and nozzle mesh | N/A (mesh visualization) | -- |
| | **c** | Jet inlet velocity waveform | [Data_Fig_S12_C_velocity_comparison.xlsx](Data_Fig_S12_C_velocity_comparison.xlsx) | -- |
| | **d** | Normalized impulse vs dimensionless wave speed | [Data_Fig_S12_D_impulse_comparison.xlsx](Data_Fig_S12_D_impulse_comparison.xlsx) | -- |
| **Fig. S13** | -- | 3D nozzle deformation and Q-criterion iso-surfaces | N/A (CFD simulation output) | -- |
| **Fig. S14** | **A** | Diaphragm pump setup | N/A (schematic) | -- |
| | **B** | Sequential vortex ring ejection visualization | N/A (image sequence) | -- |
| | **C** | Vorticity field and velocity vectors at nozzle tip | N/A (PIV visualization) | -- |
| | **D** | Temporal velocity profiles (6 V and 12 V) | [Data_Fig_S14_D.xlsx](Data_Fig_S14_D.xlsx) | -- |
| | **E** | Pump voltage, power vs dominant frequency | [Data_Fig_S14_E.xlsx](Data_Fig_S14_E.xlsx) | -- |
| | **F** | Mean jet velocity and fluctuation vs frequency | [Data_Fig_S14_F.xlsx](Data_Fig_S14_F.xlsx) | -- |
| **Fig. S15** | **A** | Snapshots: max jet height for varying tau/T | N/A (image) | -- |
| | **B** | Normalized jet height H_hat vs tau/T (N=49, n=8) | [Data_Fig_S15_B.xlsx](Data_Fig_S15_B.xlsx) | -- |
| | **C** | Horizontal jet from rigid nozzle | N/A (image) | -- |
| | **D** | Horizontal jet from flexible nozzle | N/A (image) | -- |
| | **E** | Mean jet trajectories for various tau/T | [Data_Fig_S15_E.xlsx](Data_Fig_S15_E.xlsx) | -- |
| | **F** | Time series: nozzle deformation and liquid jet | N/A (image sequence) | -- |
| | **G** | Normalized jet distance vs input frequency (harmonics) | [Data_Fig_S15_G.xlsx](Data_Fig_S15_G.xlsx) | -- |
| | **H** | Top-down jet views at 20 Hz and 40 Hz | N/A (image) | -- |
| **Fig. S16** | **A** | Boat experimental setup | N/A (photograph) | -- |
| | **B** | Boat snapshots at 0 s and 10 s | N/A (image) | -- |
| | **C** | Horizontal position vs time (power in inset) | [Data_Fig_S16_D,E.../](Data_Fig_S16_D,E_largeBoat_tau_T_dependency/) | [chatsweep_poly.m](Data_Fig_S16_D,E_largeBoat_tau_T_dependency/chatsweep_poly.m) |
| | **D** | Velocity v_hat vs tau/T | [Data_Fig_S16_D,E.../](Data_Fig_S16_D,E_largeBoat_tau_T_dependency/) | [chatsweep_poly.m](Data_Fig_S16_D,E_largeBoat_tau_T_dependency/chatsweep_poly.m) |
| | **E** | Cost of transport CoT_hat vs tau/T | [Data_Fig_S16_D,E.../](Data_Fig_S16_D,E_largeBoat_tau_T_dependency/) | [chatsweep_poly.m](Data_Fig_S16_D,E_largeBoat_tau_T_dependency/chatsweep_poly.m) |
| | **F** | Velocity v_hat vs pump frequency | [Data_Fig_S16_F,G.../](Data_Fig_S16_F,G_largeBoat_frequency_dependency/) | [VoltageSweep_poly.m](Data_Fig_S16_F,G_largeBoat_frequency_dependency/VoltageSweep_poly.m) |
| | **G** | Acceleration a_hat vs pump frequency | [Data_Fig_S16_F,G.../](Data_Fig_S16_F,G_largeBoat_frequency_dependency/) | [VoltageSweep_poly.m](Data_Fig_S16_F,G_largeBoat_frequency_dependency/VoltageSweep_poly.m) |
| **Fig. S17** | **A** | Size comparison: small boat vs original | N/A (photograph) | -- |
| | **B** | Small nozzle on fingertip | N/A (photograph) | -- |
| | **C** | Small boat snapshots at 0 s and 10 s | [Data_Fig_S17_C_smallBoat/](Data_Fig_S17_C_smallBoat/) | [ChatComparison_poly.m](Data_Fig_S17_C_smallBoat/ChatComparison_poly.m) |
| **Fig. S18** | **A** | Mixing experimental setup | N/A (schematic) | -- |
| | **B** | Image processing pipeline (binarization + edge detection) | N/A (processing illustration) | [image_processing/](Data_Fig_S18_mixing/code/image_processing/) |
| | **C** | Dye cloud side view at t\*=274 | [data/snapshots/](Data_Fig_S18_mixing/data/snapshots/) | [C05_side_boundary_overlay.m](Data_Fig_S18_mixing/code/C05_side_boundary_overlay.m) |
| | **D** | Side view boundary overlay (t = 0--6 s) | [data/snapshots/](Data_Fig_S18_mixing/data/snapshots/) | [C05_side_boundary_overlay.m](Data_Fig_S18_mixing/code/C05_side_boundary_overlay.m) |
| | **E** | Dye cloud front view at t\*=274 | [data/snapshots/](Data_Fig_S18_mixing/data/snapshots/) | [C04_front_boundary_overlay.m](Data_Fig_S18_mixing/code/C04_front_boundary_overlay.m) |
| | **F** | Front view boundary overlay (t = 0--6 s) | [data/snapshots/](Data_Fig_S18_mixing/data/snapshots/) | [C04_front_boundary_overlay.m](Data_Fig_S18_mixing/code/C04_front_boundary_overlay.m) |
| | **G** | Side view: dye area vs normalized time (intensity inset) | [data/metrics/](Data_Fig_S18_mixing/data/metrics/) | [C02_side_area_intensity.m](Data_Fig_S18_mixing/code/C02_side_area_intensity.m) |
| | **H** | Front view: dye area vs normalized time | [data/metrics/](Data_Fig_S18_mixing/data/metrics/) | [C01_front_area_intensity.m](Data_Fig_S18_mixing/code/C01_front_area_intensity.m) |
| | **I** | Side view: dye area at t\*=400 vs tau/T | [data/metrics/](Data_Fig_S18_mixing/data/metrics/) | [C03_tau_T_comparison.m](Data_Fig_S18_mixing/code/C03_tau_T_comparison.m) |
| | **J** | Front view: dye area at t\*=400 vs tau/T | [data/metrics/](Data_Fig_S18_mixing/data/metrics/) | [C03_tau_T_comparison.m](Data_Fig_S18_mixing/code/C03_tau_T_comparison.m) |
| **Table S1** | -- | Summary of experimental and simulation conditions | [Data_Fig_S11_D,E,F.xlsx](Data_Fig_S11_D,E,F.xlsx) (partial) | -- |

### Supplementary Movie S1

| Timestamp | Description | Source Data |
|-----------|-------------|-------------|
| 00:05 | Squid mantle and funnel | N/A |
| 00:12 | Histology of flexible funnel | N/A |
| 00:16 | Squid escape jetting | N/A |
| 00:27 | Setup for deformation measurement | N/A |
| 00:33 | Chromatophore tracking results | [Data_Video_00m_33s_chromatophore_tracking.xlsx](Data_Video_00m_33s_chromatophore_tracking.xlsx) |
| 00:51 | Superpropulsion illustration | N/A |
| 00:58 | Flexible nozzle fabrication | N/A |
| 01:01 | PIV: faster vortex ring from flexible nozzle | See [Fig. S6 data](#supplementary-figures) |
| 01:20 | Wave-mediated nozzle deformation | See [Fig. 2 data](#main-text-figures) |
| 01:53 | Flying squid photograph | N/A |
| 01:56 | Aerial jet: 110% higher water column | See [Fig. S15 data](#supplementary-figures) |
| 02:05 | Aerial jet: greater horizontal distance | See [Fig. 3 E-F data](#main-text-figures) |
| 02:20 | Underwater propulsion | N/A |
| 02:26 | Squid-boat manufacturing | N/A |
| 02:31 | Flexible nozzle deformation on boat | N/A |
| 02:36 | Squid-boat: 41% speed increase | See [Fig. S16 data](#supplementary-figures) |
| 02:50 | Squid inking | N/A |
| 02:53 | Faster vortex ring disperses dye faster | See [Fig. S18 data](#supplementary-figures) |
| 03:03 | Dye dispersion: 41% wider cloud | See [Fig. S18 data](#supplementary-figures) |

---

## Repository Structure

```
.
├── README.md
│
│   ── Root-level data files ──
├── Data_Fig_01_E.txt
├── Data_Fig_01_G.xlsx
├── Data_Fig_01_H.xlsx
├── Data_Fig_01_I.xlsx
├── Data_Fig_01_J.xlsx
├── Data_Fig_01_N.xlsx
├── Data_Fig_01_Q.xlsx
├── Data_Fig_02_E.xlsx
├── Data_Fig_02_I.xlsx
├── Data_Fig_03_E.xlsx
├── Data_Fig_S03_G,H,I,J.xlsx
├── Data_Fig_S03_K.xlsx
├── Data_Fig_S03_L,M,N.xlsx
├── Data_Fig_S03_additional_data_for_moving_squid.xlsx
├── Data_Fig_S04_F,G.xlsx
├── Data_Fig_S06_F_vortex_core_data.xlsx
├── Data_Fig_S06_G_impulse_history_data.xlsx
├── Data_Fig_S06_H__entrainment_data.xlsx
├── Data_Fig_S06_I__energy_history_data.xlsx
├── Data_Fig_S07_C.xlsx
├── Data_Fig_S11_A,B,C.xlsx
├── Data_Fig_S11_D,E,F.xlsx
├── Data_Fig_S12_C_velocity_comparison.xlsx
├── Data_Fig_S12_D_impulse_comparison.xlsx
├── Data_Fig_S14_D.xlsx
├── Data_Fig_S14_E.xlsx
├── Data_Fig_S14_F.xlsx
├── Data_Fig_S15_B.xlsx
├── Data_Fig_S15_E.xlsx
├── Data_Fig_S15_G.xlsx
├── Data_Video_00m_33s_chromatophore_tracking.xlsx
│
│   ── Subdirectories ──
├── Code_Fig_S03_chromatophore_tracking_algorithm/
│   ├── Code_Fig_S03_01_review_detect_black_dots.m
│   ├── Code_Fig_S03_02_batch_convert_to_bw.m
│   ├── Code_Fig_S03_03_batch_detect_stationary_noise.m
│   ├── Code_Fig_S03_04_batch_track_dots.m
│   ├── Code_Fig_S03_05_review_tracking.m
│   └── Code_Fig_S03_06_pca_rotation_analysis.m
│
├── Code_Fig_S11_superpropulsion_theory/
│   ├── S01_velocity_history.m
│   ├── S02_impulse_history.m
│   ├── S03_energy_history.m
│   └── data/
│       ├── input_curve.mat
│       ├── velocity_exp.mat
│       ├── impulse_exp.mat
│       └── energy_exp.mat
│
├── Data_Fig_03_E,F/
│   ├── C01_average_frames.m
│   ├── C02_binarize_frame.m
│   ├── C02_polygon_roi.mat
│   ├── C03_batch_binarize.m
│   ├── C04_average_binary.m
│   ├── C05_overlay_result.m
│   └── sample_output/
│       ├── averaged_binary.png
│       └── overlay_result.png
│
├── Data_Fig_S16_D,E_largeBoat_tau_T_dependency/
│   ├── Trajectories/          35 CSV: {05,10,15,20,25,30}mm + Rigid, 5 trials each
│   ├── RawPower/              30 TXT: raw INA219 readings
│   ├── PowerReadingsCropped/  Cropped power per nozzle per trial
│   ├── chatsweep_poly.m
│   └── chatsweep_sgol.m
│
├── Data_Fig_S16_F,G_largeBoat_frequency_dependency/
│   ├── Trajectories/          54 CSV: Flex15mm + Rigid at 6 voltages, 3 trials each
│   ├── RawPowerData/          Raw power + FreqTest/ subfolder
│   ├── PowerReadingsCropped/  Cropped power readings
│   ├── VoltageSweep_poly.m
│   ├── VoltageSweep_sgol.m
│   └── velocity_analysis.m
│
├── Data_Fig_S17_C_smallBoat/
│   ├── 21Nov25TrajectoryData/
│   │   ├── FN.T{1,2,4}xypts.tsv    Flexible nozzle, 3 trials
│   │   └── RN.T{1,2,4}xypts.tsv    Rigid nozzle, 3 trials
│   ├── ChatComparison_poly.m
│   └── ChatComparison_sgol.m
│
└── Data_Fig_S18_mixing/
    ├── README.md
    ├── code/
    │   ├── C01_front_area_intensity.m
    │   ├── C02_side_area_intensity.m
    │   ├── C03_tau_T_comparison.m
    │   ├── C04_front_boundary_overlay.m
    │   ├── C05_side_boundary_overlay.m
    │   └── image_processing/
    │       ├── C01_detect_green_start.m
    │       ├── C02_batch_green_detection.m
    │       ├── C10_extract_snapshots.m
    │       └── C20_detect_green_boundary.m
    └── data/
        ├── metrics/     .mat files per view/nozzle/trial
        └── snapshots/   .png frames per experiment ID
```

### Naming conventions

- **Mixing data** (`Data_Fig_S18_mixing/`): `F_` = front view, `S_` = side view; `F##mm` = flexible nozzle (## mm height), `R30mm` = rigid; `C01`, `C02`, ... = trial replicates. See [Data_Fig_S18_mixing/README.md](Data_Fig_S18_mixing/README.md) for details.
- **Boat trajectories**: `[height]9VT[trial]Trajxypts.csv` for tau/T sweep; `NoWire_[nozzle]_[voltage]_T[trial]_CTrajxypts.csv` for frequency sweep.
- **Small boat**: `FN` = flexible nozzle, `RN` = rigid nozzle.

---

## Code Index

All analysis code is written in MATLAB.

| Code | Purpose | Figures |
|------|---------|---------|
| [Code_Fig_S03_chromatophore_tracking_algorithm/](Code_Fig_S03_chromatophore_tracking_algorithm/) | Chromatophore tracking: detection, binarization, noise removal, tracking, filtering, PCA width analysis | Fig. 1 G-J, Fig. S3 |
| [Code_Fig_S11_superpropulsion_theory/](Code_Fig_S11_superpropulsion_theory/) | Superpropulsion theory: lumped ODE (Eq. 1), collapse model, velocity/impulse/energy predictions | Fig. 2 E,I,K-M, Fig. S10-S11 |
| [Data_Fig_03_E,F/](Data_Fig_03_E,F/) | Aerial jet trajectory: frame averaging, binarization, binary averaging, overlay | Fig. 3 E-F |
| [Data_Fig_S16_D,E.../chatsweep_*.m](Data_Fig_S16_D,E_largeBoat_tau_T_dependency/) | Boat tau/T sweep: polynomial + Savitzky-Golay trajectory fitting, velocity/CoT | Fig. 3 I-J,O,R-T, Fig. S16 C-E |
| [Data_Fig_S16_F,G.../VoltageSweep_*.m](Data_Fig_S16_F,G_largeBoat_frequency_dependency/) | Boat frequency sweep: trajectory and velocity analysis | Fig. S16 F-G |
| [Data_Fig_S17_C_smallBoat/ChatComparison_*.m](Data_Fig_S17_C_smallBoat/) | Small boat trajectory comparison | Fig. S17 C |
| [Data_Fig_S18_mixing/code/](Data_Fig_S18_mixing/code/) | Mixing: dye detection, area/intensity quantification, boundary overlay | Fig. 3 L-M,P,U, Fig. S18 |

---

## Fluid-Structure Interaction Simulations

The 3D FSI simulations (Figs. 2C-D, S12, S13) were performed using:
- **OpenFOAM** (pimpleFoam) -- incompressible Navier-Stokes
- **CalculiX** -- linear-elastic structural solver
- **preCICE** -- fluid-solid coupling library

Parameters: Eh = 75--500 N/m, D = 15 mm, L = 41 mm, T_acc = 0.05 s, v_jet = 0.293 m/s. A representative case template is available from the corresponding author upon request.

---

## Key Parameters

| Symbol | Definition | Range |
|--------|-----------|-------|
| tau/T | Response-time ratio (nozzle wave-transit time / jet acceleration time) | 0.01--1.5 |
| tau | Nozzle response time = L/c, where c = sqrt(Eh/(rho\*D)) | 0.005--0.3 s |
| T | Jet acceleration time | 0.04--0.34 s |
| Eh | Structural stiffness (Young's modulus x wall thickness) | 7--10^7 N/m |
| D | Nozzle inner diameter | 2--20 mm |
| L | Nozzle length | 5--60 mm |

Optimal superpropulsion regime: **tau/T = 0.2--0.4**

---

## Software Requirements

- **MATLAB** (R2023b or later recommended) for all analysis scripts
- **DLTdv8** (MATLAB toolbox) for boat trajectory tracking
- **SLEAP** for squid behavioral tracking (Fig. S4)
- **OpenFOAM**, **CalculiX**, **preCICE** for FSI simulations

---

## Citation

If you use these data or code, please cite:

> Choi, D., Singh, P., Bergerson, I., Kim, M., Park, J., Wallace, H. J., Zhang, K., Hsieh, S. Y., Asberry, A. T., Uyeno, T. A., Gilly, W. F., Park, H., Kang, D., Bose, C. & Bhamla, S. Squid-inspired superpropulsion. (2026).

## Contact

Correspondence: Saad Bhamla (saad.bhamla@colorado.edu)
