# Squid-Inspired soft superpropulsion: Data and Code Repository

Source data and analysis code for:

**"Squid-inspired soft superpropulsion"**
Choi, Singh, Bergerson, Kim, Park, Wallace, Zhang, Hsieh, Asberry, Uyeno, Gilly, Park, Kang, Bose & Bhamla

---

## Figure-to-Data Lookup

Every figure panel from the main text, Supplementary Information, and Supplementary Movie is listed below. Figure numbers follow the revised manuscript. `N/A` indicates panels that are schematics, photographs, or microscopy images with no associated source-data file.

### Main-Text Figures

| Figure | Panel | Description <img width=400/> | Source Data | Code |
|--------|-------|-------------|-------------|------|
| **Fig. 1** | **A** | Squid mantle cavity and funnel photograph | N/A | -- |
| | **B** | Mantle contraction and funnel expansion image | N/A | -- |
| | **C** | Picrosirius Red staining, mantle (polarized light) | N/A | -- |
| | **D** | Picrosirius Red staining, funnel (polarized light) | N/A | -- |
| | **E** | Snapshots: initial state, peak funnel widening, recoil | [01_E.txt](Data_Fig_01_E.txt) | -- |
| | **F** | Conceptual schematic of expansion--recoil sequence | N/A | -- |
| | **G** | Normalized W_f\*, W_m\*, and pressure during jet pulse (c1) | [01_G.xlsx](Data_Fig_01_G.xlsx) | [S03 tracking/](Code_Fig_S03_chromatophore_tracking_algorithm/) |
| | **H** | Time history: stimulus, pressure, W_m\*, W_f\* | [01_H.xlsx](Data_Fig_01_H.xlsx) | [S03 tracking/](Code_Fig_S03_chromatophore_tracking_algorithm/) |
| | **I** | W_f,max\* and W_m,min\* for c1, c2, c3 (track jumps corrected, Fig. S5) | [01_I.xlsx](Data_Fig_01_I.xlsx) | -- |
| | **J** | Funnel response-time ratio tau_squid/T_squid (track jumps corrected) | [01_J.xlsx](Data_Fig_01_J.xlsx) | -- |
| **Fig. 2** | **A** | PIV snapshots, rigid nozzle | N/A | -- |
| | **B** | PIV snapshots, compliant nozzle + nozzle deformation | N/A | -- |
| | **C** | 3D-CFD snapshots, rigid nozzle | N/A | -- |
| | **D** | 3D-CFD snapshots, compliant nozzle | N/A | -- |
| | **E** | Impulse ratio I_h_hat vs tau/T (PIV + CFD + theory + load cell) | [02_E.xlsx](Data_Fig_02_E.xlsx) | [S18 theory/](Code_Fig_S18_superpropulsion_theory/) |
| | **F** | Schematic: elastic strain power during store--release | N/A | -- |
| | **G** | Spring--mass resonator analogy | N/A | -- |
| | **H** | Time traces: elastic strain power and jet kinetic power | [S11_I_energy.xlsx](Data_Fig_S11_I_energy_history_data.xlsx) | -- |
| | **I** | Measured tau/T vs theoretical estimates | [02_I.xlsx](Data_Fig_02_I.xlsx) | [S18 theory/](Code_Fig_S18_superpropulsion_theory/) |
| | **J** | Schematic: stiff vs optimal vs compliant nozzle recoil | N/A | -- |
| | **K** | Exit velocity ratio v_hat vs time (theory + exp) | [S18_A,B,C.xlsx](Data_Fig_S18_A,B,C.xlsx) | [S01_velocity.m](Code_Fig_S18_superpropulsion_theory/S01_velocity_history.m) |
| | **L** | Impulse ratio I_h_hat vs time (theory + exp) | [S18_A,B,C.xlsx](Data_Fig_S18_A,B,C.xlsx) | [S02_impulse.m](Code_Fig_S18_superpropulsion_theory/S02_impulse_history.m) |
| | **M** | Kinetic energy ratio E_k_hat vs time (theory + exp) | [S18_A,B,C.xlsx](Data_Fig_S18_A,B,C.xlsx) | [S03_energy.m](Code_Fig_S18_superpropulsion_theory/S03_energy_history.m) |
| **Fig. 3** | **A** | Flying squid photograph | N/A | -- |
| | **B** | Single-pulse aerial jet height | [S23_B.xlsx](Data_Fig_S23_B.xlsx) | -- |
| | **C** | Input waveforms (single + repeated pulse) | [S22_D.xlsx](Data_Fig_S22_D.xlsx) | -- |
| | **D** | Time sequence: jet development and nozzle deformation | N/A | -- |
| | **E** | Candle extinguishing range + averaged jet trajectories (inset) | [03_E.xlsx](Data_Fig_03_E.xlsx), [03_E,F/](Data_Fig_03_E,F/) | [03_E,F/](Data_Fig_03_E,F/) |
| | **F** | Candle extinguishing comparison (rigid vs compliant) | [03_E,F/](Data_Fig_03_E,F/) | [03_E,F/](Data_Fig_03_E,F/) |
| | **G** | Nozzle set (heights 5--30 mm) photograph | N/A | -- |
| | **H** | Squid-boat with compliant nozzle | N/A | -- |
| | **I** | Boat translation, rigid nozzle | [S25_D,E/](Data_Fig_S25_D,E_largeBoat_tau_T_dependency/) | [chatsweep_poly.m](Data_Fig_S25_D,E_largeBoat_tau_T_dependency/chatsweep_poly.m) |
| | **J** | Boat translation, compliant nozzle | [S25_D,E/](Data_Fig_S25_D,E_largeBoat_tau_T_dependency/) | [chatsweep_poly.m](Data_Fig_S25_D,E_largeBoat_tau_T_dependency/chatsweep_poly.m) |
| | **K** | Squid inking photograph | N/A | -- |
| | **L** | Dye visualization, front view (compliant vs rigid) | [S27 snapshots/](Data_Fig_S27_mixing/data/snapshots/) | [C04_front_boundary.m](Data_Fig_S27_mixing/code/C04_front_boundary_overlay.m) |
| | **M** | Dye visualization, side view (compliant vs rigid) | [S27 snapshots/](Data_Fig_S27_mixing/data/snapshots/) | [C05_side_boundary.m](Data_Fig_S27_mixing/code/C05_side_boundary_overlay.m) |
| | **N** | Horizontal jet distance x/D and power P/P_0 vs tau/T | [03_N.xlsx](Data_Fig_03_N.xlsx) | -- |
| | **O** | Boat speed and CoT vs tau/T | [03_O.xlsx](Data_Fig_03_O.xlsx) | -- |
| | **P** | Dye dispersion area A/A_nozzle vs tau/T | [S27 metrics/](Data_Fig_S27_mixing/data/metrics/) | [C03_tau_T.m](Data_Fig_S27_mixing/code/C03_tau_T_comparison.m) |
| | **Q** | Bar chart: mean horizontal jet distance +44.6% | [03_Q.xlsx](Data_Fig_03_Q.xlsx) | -- |
| | **R** | Bar chart: boat speed +42.6% | [03_R,S,T.xlsx](Data_Fig_03_R,S,T.xlsx) | -- |
| | **S** | Bar chart: acceleration +32.0% | [03_R,S,T.xlsx](Data_Fig_03_R,S,T.xlsx) | -- |
| | **T** | Bar chart: CoT -31.4% | [03_R,S,T.xlsx](Data_Fig_03_R,S,T.xlsx) | -- |
| | **U** | Bar chart: dye dispersion +40.5% | panel P | -- |

### Supplementary Figures

| Figure | Panel | Description <img width=400/> | Source Data | Code |
|--------|-------|-------------|-------------|------|
| **Fig. S1** | **A** | Monterey Bay collection site | N/A | -- |
| | **B** | Circulating seawater tank, Hopkins Marine Station | N/A | -- |
| | **C** | Gijang, Korea collection site | N/A | -- |
| | **D** | Circulating seawater tank, Gijang Bay | N/A | -- |
| **Fig. S2** | **A** | Fixed squid photograph | N/A | -- |
| | **B** | Cross-section of funnel wall (Milligan's Trichrome) | N/A | -- |
| | **C** | Funnel: H&E, Masson's, PSR brightfield + polarized | N/A | -- |
| | **D** | Mantle: H&E, Masson's, PSR brightfield + polarized | N/A | -- |
| **Fig. S3** | **A** | High-speed imaging setup schematic | N/A | -- |
| | **B** | Raw high-speed image with chromatophores | N/A | [S03 tracking/](Code_Fig_S03_chromatophore_tracking_algorithm/) |
| | **C** | Close-up: electrode and pressure sensor | N/A | -- |
| | **D** | Image processing: intensity normalization | N/A | [S03_01_detect.m](Code_Fig_S03_chromatophore_tracking_algorithm/Code_Fig_S03_01_review_detect_black_dots.m) |
| | **E** | Image processing: binarization | N/A | [S03_02_bw.m](Code_Fig_S03_chromatophore_tracking_algorithm/Code_Fig_S03_02_batch_convert_to_bw.m) |
| | **F** | Image processing: particle tracking with ROI | N/A | [S03_04_track.m](Code_Fig_S03_chromatophore_tracking_algorithm/Code_Fig_S03_04_batch_track_dots.m) |
| | **G** | Electrical stimulus timing | [S03_G,H,I,J.xlsx](Data_Fig_S03_G,H,I,J.xlsx) | [S03 tracking/](Code_Fig_S03_chromatophore_tracking_algorithm/) |
| | **H** | Pressure sensor output | [S03_G,H,I,J.xlsx](Data_Fig_S03_G,H,I,J.xlsx) | [S03 tracking/](Code_Fig_S03_chromatophore_tracking_algorithm/) |
| | **I** | Funnel width variation | [S03_G,H,I,J.xlsx](Data_Fig_S03_G,H,I,J.xlsx) | [S03 tracking/](Code_Fig_S03_chromatophore_tracking_algorithm/) |
| | **J** | Mantle width variation | [S03_G,H,I,J.xlsx](Data_Fig_S03_G,H,I,J.xlsx) | [S03 tracking/](Code_Fig_S03_chromatophore_tracking_algorithm/) |
| | **K** | Selected strong jets (N=30, n=8) | [S03_K.xlsx](Data_Fig_S03_K.xlsx) | -- |
| | **L** | Averaged pressure profiles | [S03_L,M,N.xlsx](Data_Fig_S03_L,M,N.xlsx) | -- |
| | **M** | Normalized funnel deformation | [S03_L,M,N.xlsx](Data_Fig_S03_L,M,N.xlsx) | -- |
| | **N** | Normalized mantle deformation | [S03_L,M,N.xlsx](Data_Fig_S03_L,M,N.xlsx) | -- |
| | -- | Additional data: freely moving *S. lessoniana* | [S03_additional.xlsx](Data_Fig_S03_additional_data_for_moving_squid.xlsx) | -- |
| **Fig. S4** | **A** | Membrane tensile setup | N/A | -- |
| | **B** | Frames at five grip displacements | N/A | -- |
| | **C** | DIC displacement field and tracked dots | N/A | -- |
| | **D** | Strain from DIC and from dot tracking vs grip displacement | [S04_D,E.xlsx](Data_Fig_S04_D,E.xlsx) | -- |
| | **E** | Dot displacement, dot tracking vs DIC | [S04_D,E.xlsx](Data_Fig_S04_D,E.xlsx) | -- |
| **Fig. S5** | **A** | Frames of a track jump | N/A | -- |
| | **B** | Funnel width, original vs jumps corrected | [S05_B,C.xlsx](Data_Fig_S05_B,C.xlsx) | -- |
| | **C** | Peak funnel and minimum mantle width per jet, original vs corrected | [S05_B,C.xlsx](Data_Fig_S05_B,C.xlsx) | -- |
| **Fig. S6** | **A** | Calibration and balloon setup, top views | N/A | -- |
| | **B** | Virtual cameras and views of the calibration target | N/A | -- |
| | **C** | Virtual cameras and views of the balloon | N/A | -- |
| | **D** | Area change of the dotted patch: 3D-DIC, 2D-DIC, dot tracking | [S06_D.xlsx](Data_Fig_S06_D.xlsx) | -- |
| | **E** | Frontal view before and after inflation | N/A | -- |
| | **F** | 3D-DIC area change field | N/A | -- |
| | **G** | 2D-DIC area change field | N/A | -- |
| | **H** | Dot-tracking area change field | N/A | -- |
| **Fig. S7** | **A** | Schematic of squid nervous system | N/A | -- |
| | **B** | Nerve infundibulum, intact condition | N/A | -- |
| | **C** | Nerve infundibulum, denervated condition | N/A | -- |
| | **D** | Head-first reflex (SLEAP tracking) | N/A | -- |
| | **E** | Tail-first reflex (SLEAP tracking) | N/A | -- |
| | **F** | Directional velocity, intact squid (N=45) | [S07_F,G.xlsx](Data_Fig_S07_F,G.xlsx) | -- |
| | **G** | Directional velocity, funnel-denervated squid (N=34) | [S07_F,G.xlsx](Data_Fig_S07_F,G.xlsx) | -- |
| **Fig. S8** | **A** | Flange-clamped fabrication (thick nozzles) | N/A | -- |
| | **B** | Flange-embedded fabrication (thin nozzles) | N/A | -- |
| | **C** | Press-fit fabrication | N/A | -- |
| **Fig. S9** | **A** | Sectioning and imaging of the nozzle wall | N/A | -- |
| | **B** | Wall thickness vs silicone thinner fraction | [S09_B.xlsx](Data_Fig_S09_B.xlsx) | -- |
| **Fig. S10** | **A** | Tensile stage | N/A | -- |
| | **B** | Speckle images and strain field | N/A | -- |
| | **C** | Tension vs strain for each thinner fraction | [S10_C.xlsx](Data_Fig_S10_C.xlsx) | -- |
| **Fig. S11** | **A** | Single pulsed jet generator setup | N/A | -- |
| | **B** | Overlaid raw images of jet and nozzle deformation | N/A | -- |
| | **C** | Vorticity field, rigid nozzle (Eh = infinity) | N/A | -- |
| | **D** | Vorticity field, Eh = 14.4 N/m | N/A | -- |
| | **E** | Vorticity field, Eh = 7.0 N/m | N/A | -- |
| | **F** | Vortex core travel distance vs t/T | [S11_F_vortex.xlsx](Data_Fig_S11_F_vortex_core_data.xlsx) | -- |
| | **G** | Hydrodynamic impulse I_h vs t/T | [S11_G_impulse.xlsx](Data_Fig_S11_G_impulse_history_data.xlsx) | -- |
| | **H** | Cumulative entrainment V/V_nozzle vs t/T | [S11_H_entrainment.xlsx](Data_Fig_S11_H_entrainment_data.xlsx) | -- |
| | **I** | Energy partitioning: E_elastic and E_kinetic vs t/T | [S11_I_energy.xlsx](Data_Fig_S11_I_energy_history_data.xlsx) | -- |
| **Fig. S12** | **A** | PIV impulse, rigid and flexible nozzle (T = 0.11 s) | [S11_G_impulse.xlsx](Data_Fig_S11_G_impulse_history_data.xlsx) | -- |
| | **B** | Reduced-order model for the same case | Theory output | [S18 theory/](Code_Fig_S18_superpropulsion_theory/) |
| | **C** | Impulse ratio vs tau/T, as reported and at t = 4T | [S12_C,D.xlsx](Data_Fig_S12_C,D.xlsx), [02_E.xlsx](Data_Fig_02_E.xlsx) | -- |
| | **D** | Value at t = 4T vs reported value | [S12_C,D.xlsx](Data_Fig_S12_C,D.xlsx) | -- |
| **Fig. S13** | **A** | Load-cell rig photograph | N/A | -- |
| | **B** | Rig schematic with drive and recording electronics | N/A | -- |
| | **C** | Rigid-nozzle shot, side-camera frames | N/A | -- |
| | **D** | Piston displacement of the filmed rigid shots | [S13_D.xlsx](Data_Fig_S13_D.xlsx) | -- |
| | **E** | Impact speed and push duration vs drive level | [S13_E.xlsx](Data_Fig_S13_E.xlsx) | -- |
| | **F** | Rigid nozzle, measured I_loadcell and prediction | [S13_F.xlsx](Data_Fig_S13_F.xlsx) | -- |
| **Fig. S14** | **A-C** | Force and impulse histories at 3.5, 5.0 and 7.0 V | [S14_A,B,C.xlsx](Data_Fig_S14_A,B,C.xlsx) | -- |
| | **D** | Corrected impulse vs tau/T | [S14_D.xlsx](Data_Fig_S14_D.xlsx) | -- |
| | **E** | Load-cell impulse ratio over the data of Fig. 2E | [S14_E.xlsx](Data_Fig_S14_E.xlsx), [02_E.xlsx](Data_Fig_02_E.xlsx) | -- |
| **Fig. S15** | -- | Schematic of 1D tube model | N/A | -- |
| **Fig. S16** | **A** | Spring--damper--mass resonator analogy | N/A | -- |
| | **B** | Electrical RLC resonator analogy | N/A | -- |
| | **C** | Flexible nozzle system analogy | N/A | -- |
| **Fig. S17** | **A** | Theoretical exit velocity: base vs collapse model | Theory output | [S01_velocity.m](Code_Fig_S18_superpropulsion_theory/S01_velocity_history.m) |
| | **B** | Velocity fields, rigid nozzle (t/T = 1.0, 2.0, 3.2) | N/A | -- |
| | **C** | Velocity fields, flexible nozzle (t/T = 1.0, 2.0, 3.2) | N/A | -- |
| | **D** | Theory vs exp: exit velocity, tau/T = 0.17 | [S18_A,B,C.xlsx](Data_Fig_S18_A,B,C.xlsx), inset: [S17_D_inset.xlsx](Data_Fig_S17_D_inset.xlsx) | [S01_velocity.m](Code_Fig_S18_superpropulsion_theory/S01_velocity_history.m) |
| | **E** | Theory vs exp: exit velocity, tau/T = 0.29 | [S18_A,B,C.xlsx](Data_Fig_S18_A,B,C.xlsx) | [S01_velocity.m](Code_Fig_S18_superpropulsion_theory/S01_velocity_history.m) |
| | **F** | Theory vs exp: exit velocity, tau/T = 0.55 | [S18_A,B,C.xlsx](Data_Fig_S18_A,B,C.xlsx) | [S01_velocity.m](Code_Fig_S18_superpropulsion_theory/S01_velocity_history.m) |
| **Fig. S18** | **A** | Exit velocity v_hat time histories (3 tau/T values) | [S18_A,B,C.xlsx](Data_Fig_S18_A,B,C.xlsx) | [S01_velocity.m](Code_Fig_S18_superpropulsion_theory/S01_velocity_history.m) |
| | **B** | Impulse I_h_hat time histories (3 tau/T values) | [S18_A,B,C.xlsx](Data_Fig_S18_A,B,C.xlsx) | [S02_impulse.m](Code_Fig_S18_superpropulsion_theory/S02_impulse_history.m) |
| | **C** | Kinetic energy E_k_hat time histories (3 tau/T values) | [S18_A,B,C.xlsx](Data_Fig_S18_A,B,C.xlsx) | [S03_energy.m](Code_Fig_S18_superpropulsion_theory/S03_energy_history.m) |
| | **D** | Peak v_hat vs tau/T (exp + CFD + theory) | [S18_D,E,F.xlsx](Data_Fig_S18_D,E,F.xlsx) | [S01_velocity.m](Code_Fig_S18_superpropulsion_theory/S01_velocity_history.m) |
| | **E** | I_h_hat vs tau/T | [S18_D,E,F.xlsx](Data_Fig_S18_D,E,F.xlsx) | [S02_impulse.m](Code_Fig_S18_superpropulsion_theory/S02_impulse_history.m) |
| | **F** | E_k_hat vs tau/T | [S18_D,E,F.xlsx](Data_Fig_S18_D,E,F.xlsx) | [S03_energy.m](Code_Fig_S18_superpropulsion_theory/S03_energy_history.m) |
| **Fig. S19** | **A** | Frames of the flexible nozzle | N/A | -- |
| | **B** | Cumulative volume through the plane y/D = -1 | [S11_H_entrainment.xlsx](Data_Fig_S11_H_entrainment_data.xlsx) | -- |
| | **C** | Hydrodynamic impulse vs time, with models | [S11_G_impulse.xlsx](Data_Fig_S11_G_impulse_history_data.xlsx) | [S18 theory/](Code_Fig_S18_superpropulsion_theory/) |
| | **D** | Jet kinetic energy and wall elastic energy vs time | [S11_I_energy.xlsx](Data_Fig_S11_I_energy_history_data.xlsx) | -- |
| | **E** | Impulse ratio vs kinetic-energy ratio | [S18_D,E,F.xlsx](Data_Fig_S18_D,E,F.xlsx) | -- |
| **Fig. S20** | **A** | Computational domain and boundary conditions | N/A | -- |
| | **B** | Fluid domain mesh and nozzle mesh | N/A | -- |
| | **C** | Jet inlet velocity waveform | [S20_C_velocity.xlsx](Data_Fig_S20_C_velocity_comparison.xlsx) | -- |
| | **D** | Normalized impulse vs dimensionless wave speed | [S20_D_impulse.xlsx](Data_Fig_S20_D_impulse_comparison.xlsx) | -- |
| **Fig. S21** | -- | 3D nozzle deformation and Q-criterion iso-surfaces | N/A | -- |
| **Fig. S22** | **A** | Diaphragm pump setup | N/A | -- |
| | **B** | Sequential vortex ring ejection visualization | N/A | -- |
| | **C** | Vorticity field and velocity vectors at nozzle tip | N/A | -- |
| | **D** | Temporal velocity profiles (6 V and 12 V) | [S22_D.xlsx](Data_Fig_S22_D.xlsx) | -- |
| | **E** | Pump voltage, power vs dominant frequency | [S22_E.xlsx](Data_Fig_S22_E.xlsx) | -- |
| | **F** | Mean jet velocity and fluctuation vs frequency | [S22_F.xlsx](Data_Fig_S22_F.xlsx) | -- |
| **Fig. S23** | **A** | Snapshots: max jet height for varying tau/T | N/A | -- |
| | **B** | Normalized jet height H_hat vs tau/T (N=49, n=8) | [S23_B.xlsx](Data_Fig_S23_B.xlsx) | -- |
| | **C** | Horizontal jet from rigid nozzle | N/A | -- |
| | **D** | Horizontal jet from flexible nozzle | N/A | -- |
| | **E** | Mean jet trajectories for various tau/T | [S23_E.xlsx](Data_Fig_S23_E.xlsx) | -- |
| | **F** | Time series: nozzle deformation and liquid jet | N/A | -- |
| | **G** | Normalized jet distance vs input frequency (harmonics) | [S23_G.xlsx](Data_Fig_S23_G.xlsx) | -- |
| | **H** | Top-down jet views at 20 Hz and 40 Hz | N/A | -- |
| **Fig. S24** | **A** | Downward jets, representative pulses | N/A | -- |
| | **B** | Front position z/D vs t\* | [S24_B,C.xlsx](Data_Fig_S24_B,C.xlsx) | -- |
| | **C** | Exit speed relative to the rigid nozzle at three launch angles | [S24_B,C.xlsx](Data_Fig_S24_B,C.xlsx) | -- |
| **Fig. S25** | **A** | Boat experimental setup | N/A | -- |
| | **B** | Boat snapshots at 0 s and 10 s | N/A | -- |
| | **C** | Horizontal position vs time (power in inset) | [S25_D,E/](Data_Fig_S25_D,E_largeBoat_tau_T_dependency/) | [chatsweep_poly.m](Data_Fig_S25_D,E_largeBoat_tau_T_dependency/chatsweep_poly.m) |
| | **D** | Velocity v_hat vs tau/T | [S25_D,E.xlsx](Data_Fig_S25_D,E.xlsx), [S25_D,E/](Data_Fig_S25_D,E_largeBoat_tau_T_dependency/) | [chatsweep_poly.m](Data_Fig_S25_D,E_largeBoat_tau_T_dependency/chatsweep_poly.m) |
| | **E** | Cost of transport CoT_hat vs tau/T | [S25_D,E.xlsx](Data_Fig_S25_D,E.xlsx), [S25_D,E/](Data_Fig_S25_D,E_largeBoat_tau_T_dependency/) | [chatsweep_poly.m](Data_Fig_S25_D,E_largeBoat_tau_T_dependency/chatsweep_poly.m) |
| | **F** | Velocity v_hat vs pump frequency | [S25_F,G/](Data_Fig_S25_F,G_largeBoat_frequency_dependency/) | [VoltageSweep_poly.m](Data_Fig_S25_F,G_largeBoat_frequency_dependency/VoltageSweep_poly.m) |
| | **G** | Acceleration a_hat vs pump frequency | [S25_F,G/](Data_Fig_S25_F,G_largeBoat_frequency_dependency/) | [VoltageSweep_poly.m](Data_Fig_S25_F,G_largeBoat_frequency_dependency/VoltageSweep_poly.m) |
| **Fig. S26** | **A** | Size comparison: small boat vs original | N/A | -- |
| | **B** | Small nozzle on fingertip | N/A | -- |
| | **C** | Small boat snapshots at 0 s and 10 s | [S26_C/](Data_Fig_S26_C_smallBoat/) | [ChatComparison_poly.m](Data_Fig_S26_C_smallBoat/ChatComparison_poly.m) |
| **Fig. S27** | **A** | Mixing experimental setup | N/A | -- |
| | **B** | Image processing pipeline (binarization + edge detection) | N/A | [image_processing/](Data_Fig_S27_mixing/code/image_processing/) |
| | **C** | Dye cloud side view at t\*=274 | [S27 snapshots/](Data_Fig_S27_mixing/data/snapshots/) | [C05_side_boundary.m](Data_Fig_S27_mixing/code/C05_side_boundary_overlay.m) |
| | **D** | Side view boundary overlay (t = 0--6 s) | [S27 snapshots/](Data_Fig_S27_mixing/data/snapshots/) | [C05_side_boundary.m](Data_Fig_S27_mixing/code/C05_side_boundary_overlay.m) |
| | **E** | Dye cloud front view at t\*=274 | [S27 snapshots/](Data_Fig_S27_mixing/data/snapshots/) | [C04_front_boundary.m](Data_Fig_S27_mixing/code/C04_front_boundary_overlay.m) |
| | **F** | Front view boundary overlay (t = 0--6 s) | [S27 snapshots/](Data_Fig_S27_mixing/data/snapshots/) | [C04_front_boundary.m](Data_Fig_S27_mixing/code/C04_front_boundary_overlay.m) |
| | **G** | Side view: dye area vs normalized time (intensity inset) | [S27 metrics/](Data_Fig_S27_mixing/data/metrics/) | [C02_side_area.m](Data_Fig_S27_mixing/code/C02_side_area_intensity.m) |
| | **H** | Front view: dye area vs normalized time | [S27 metrics/](Data_Fig_S27_mixing/data/metrics/) | [C01_front_area.m](Data_Fig_S27_mixing/code/C01_front_area_intensity.m) |
| | **I** | Side view: dye area at t\*=400 vs tau/T | [S27 metrics/](Data_Fig_S27_mixing/data/metrics/) | [C03_tau_T.m](Data_Fig_S27_mixing/code/C03_tau_T_comparison.m) |
| | **J** | Front view: dye area at t\*=400 vs tau/T | [S27 metrics/](Data_Fig_S27_mixing/data/metrics/) | [C03_tau_T.m](Data_Fig_S27_mixing/code/C03_tau_T_comparison.m) |
| **Table S1** | -- | Summary of experimental and simulation conditions | [S18_D,E,F.xlsx](Data_Fig_S18_D,E,F.xlsx) | -- |

### Supplementary Movie S1

| Timestamp | Description <img width=400/> | Source Data |
|-----------|-------------|-------------|
| 00:05 | Squid mantle and funnel | N/A |
| 00:12 | Histology of flexible funnel | N/A |
| 00:16 | Squid escape jetting | N/A |
| 00:27 | Setup for deformation measurement | N/A |
| 00:33 | Chromatophore tracking results | [Video_00m33s.xlsx](Data_Video_00m_33s_chromatophore_tracking.xlsx) |
| 00:51 | Superpropulsion illustration | N/A |
| 00:58 | Flexible nozzle fabrication | N/A |
| 01:01 | PIV: faster vortex ring from flexible nozzle | See [Fig. S11 data](#supplementary-figures) |
| 01:20 | Wave-mediated nozzle deformation | See [Fig. 2 data](#main-text-figures) |
| 01:31 | 3D FSI simulation and T/tau, theory vs measurement | See [Fig. 2 data](#main-text-figures) |
| 01:40 | Optimal condition for impulse generation | See [Fig. 2 data](#main-text-figures) |
| 01:46 | PIV of the optimally flexible nozzle jet | See [Fig. S11 data](#supplementary-figures) |
| 01:53 | Flying squid photograph | N/A |
| 01:56 | Aerial jet: higher water column | See [Fig. S23 data](#supplementary-figures) |
| 02:05 | Aerial jet: greater horizontal distance | See [Fig. 3 E-F data](#main-text-figures) |
| 02:20 | Underwater propulsion | N/A |
| 02:26 | Squid-boat manufacturing | N/A |
| 02:31 | Flexible nozzle deformation on boat | N/A |
| 02:36 | Squid-boat: 43% speed increase | See [Fig. S25 data](#supplementary-figures) |
| 02:50 | Squid inking | N/A |
| 02:53 | Faster vortex ring disperses dye faster | See [Fig. S27 data](#supplementary-figures) |
| 03:03 | Dye dispersion: wider cloud | See [Fig. S27 data](#supplementary-figures) |

---

## Code Index

All analysis code is written in MATLAB.

| Code | Purpose | Figures |
|------|---------|---------|
| [S03 tracking/](Code_Fig_S03_chromatophore_tracking_algorithm/) | Chromatophore tracking: detection, binarization, noise removal, tracking, filtering, PCA width analysis | Fig. 1 G-J, Fig. S3 |
| [S18 theory/](Code_Fig_S18_superpropulsion_theory/) | Superpropulsion theory: lumped ODE (Eq. 1), collapse model, velocity/impulse/energy predictions | Fig. 2 E,I,K-M, Fig. S17-S18 |
| [03_E,F/](Data_Fig_03_E,F/) | Aerial jet trajectory: frame averaging, binarization, binary averaging, overlay | Fig. 3 E-F |
| [S25_D,E/ chatsweep_*.m](Data_Fig_S25_D,E_largeBoat_tau_T_dependency/) | Boat tau/T sweep: polynomial + Savitzky-Golay trajectory fitting, velocity/CoT | Fig. 3 I-J, Fig. S25 C |
| [S25_F,G/ VoltageSweep_*.m](Data_Fig_S25_F,G_largeBoat_frequency_dependency/) | Boat frequency sweep: trajectory and velocity analysis | Fig. S25 F-G |
| [S26_C/ ChatComparison_*.m](Data_Fig_S26_C_smallBoat/) | Small boat trajectory comparison | Fig. S26 C |
| [S27 mixing code/](Data_Fig_S27_mixing/code/) | Mixing: dye detection, area/intensity quantification, boundary overlay | Fig. 3 L-M,P,U, Fig. S27 |

The plotted boat values of Fig. 3 O, R-T and Fig. S25 D-E (each position record fitted up to 2.4 body lengths) are in `Data_Fig_03_O.xlsx`, `Data_Fig_03_R,S,T.xlsx` and `Data_Fig_S25_D,E.xlsx`.

---

## Fluid-Structure Interaction Simulations

The 3D FSI simulations (Figs. 2C-D, S20, S21) were performed using:
- **OpenFOAM** (pimpleFoam) -- incompressible Navier-Stokes
- **CalculiX** -- linear-elastic structural solver
- **preCICE** -- fluid-solid coupling library

Parameters: Eh = 75--500 N/m, D = 15 mm, L = 41 mm, T_acc = 0.05 s, v_jet = 0.293 m/s. A representative case template is available from the corresponding author upon request.

---

## Software Requirements

- **MATLAB** (R2023b or later recommended) for all analysis scripts
- **DLTdv8** (MATLAB toolbox) for boat trajectory tracking
- **SLEAP** for squid behavioral tracking (Fig. S7)
- **OpenFOAM**, **CalculiX**, **preCICE** for FSI simulations
