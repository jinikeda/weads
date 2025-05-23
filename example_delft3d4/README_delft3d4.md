# File: README_delft3d4.md
# Delft3D–WEADS Coupling Framework

This project simulates eco-hydrological feedbacks using a grid-based coupling between **Delft3D4-FLOW** and a simplified WEADS ecological model. It is designed for multi-year simulations with yearly feedback of bathymetric and roughness updates based on MEM outputs.

---

## 📁 Directory Structure

```
basespace/
├── WEADS_delft3d4.py                  # Main driver script (argument interface)
├── run_delft3d4_loop.sh               # Bash loop for yearly simulation
├── clean_delft3d4.sh                  # Resets everything for a clean run
├── env_delft3d4.yml                   # Conda environment file
├── example_delft3d4/                  # Contains original Delft3D input files (restored before each run)
├── src_point/                         # Python module for preprocessing, MEM, and postprocessing
│   ├── basics.py
│   ├── preprocessing_delft3d4.py
│   ├── mem_delft3d4.py
│   ├── postprocessing_delft3d4.py
│   ├── tidaldatums_delft3d4.py
│   └── general_functions_delft3d4.py
└── src_delft3d4/                       # Additional scripts and MATLAB tools
    ├── plots_delft3d4.py
    ├── delft3d4_dep.py
    ├── delft3d4_rgh.py
    ├── waterlevel_delft3d4.m
    └── delft3d_matlab/                 # **Git submodule: required for MATLAB water level extraction**
```

---

## Initial Setup and Requirements

### MATLAB Submodule Setup (⚠️ Required)
This repository uses a Git submodule to include shared Delft3D MATLAB functions (e.g., `vs_use.m`) in:

```
src_delft3d4/delft3d_matlab/
```

To initialize the submodule after cloning:

```bash
git submodule update --init --recursive
```

> ✅ **This step is mandatory** before running the MATLAB script `waterlevel_delft3d4.m`

---

## 🧩 Input Files

Before running the model, all required **Delft3D input files** must be placed in the `example_delft3d4/` folder. These include:

- `45x45.grd` — grid file
- `45x45.enc` — enclosure
- `45x45.dep` — initial bathymetry
- `45x45.rgh` — roughness (Manning's n)
- `45x45.mdf`, `45x45.bnd`, `45x45.bct`, `45x45.obs`, `45x45.url` — other configuration files
- `config_dflow2d3d.xml` — Delft3D runtime driver

> ℹ️ Tip: You can use `/src_delft3d4/delft3d4_dep.py` and `/src_delft3d4/delft3d4_rgh.py` to generate structured `.dep` and `.rgh` files.

---

## ⚙️ Step-by-Step Instructions

Note: This project is intended to be run on the CEDS system, where Delft3D4 is installed at `/usr/local/delft3d4` and configured with MATLAB and Conda.

### 1. Create and Activate the Environment

```bash
conda env create -f env_delft3d4.yml
conda activate WEADS_delft3d4_env
```

> Only needed once. Later, use `conda activate WEADS_delft3d4_env` before any run.

---

### 2. Clean Up and Restore Initial Files

```bash
bash clean_delft3d4.sh
```

This will:
- Delete previous Delft3D and MEM outputs
- Remove `updated_*.dep`, `updated_*.rgh`, and `yearly_outputs/`
- Restore all original model inputs from `example_delft3d4/`

---

### 3. Run the Coupling Loop (Multi-Year)

```bash
bash run_delft3d4_loop.sh
```

Each iteration will:
- Run Delft3D
- Extract water levels via MATLAB
- Compute tidal metrics
- Run preprocessing → MEM → postprocessing
- Update `.dep` and `.rgh`
- Save all results in `yearly_outputs/year1/`, `year2/`, etc.

---

## 📊 Generate Summary Plots

```bash
python src_delft3d4/plots_delft3d4.py
```

Saved to:
```
yearly_outputs/summary_intertidal_plots.png
yearly_outputs/summary_subtidal_plots.png
yearly_outputs/summary_supratidal_plots.png
```

---

## 🗂️ Output Structure After Running

```
basespace/
├── updated_45x45.dep
├── updated_45x45.rgh
├── ecology.csv
├── tidal_metrics.csv
├── water_level_extracted.csv
├── grid_coords_extracted.csv
├── yearly_outputs/
│   ├── year1/
│   │   ├── ecology.csv
│   │   ├── updated_45x45.dep
│   │   ├── updated_45x45.rgh
│   │   ├── trim_45x45.nc
│   └── year2/
│       └── ...
```

---

## ✅ Final Checklist

- [x] Delft3D inputs exist in `example_delft3d4/`
- [x] Submodule initialized using `git submodule update --init --recursive`
- [x] Conda environment `WEADS_delft3d4_env` is created
- [x] Run `clean_delft3d4.sh` before each simulation
- [x] Run `run_delft3d4_loop.sh` to execute full coupling
- [x] Plots automatically generated after the loop
- [x] **Delft3D source code and binaries are installed and accessible on your system**

---

## 📬 Contact

Developed by: **Shabnam Mirheidarian**  
📧 smirhe1@lsu.edu
