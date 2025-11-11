# igUEL
Abaqus User ELement for Isogeometric Analysis.  
Developed to solve problems of classical and strain gradient elasticity.

## MATLAB Preprocessing Pipeline

### Quick Start
1. Navigate to `2D_PrePro_matlab/` directory
2. Edit `config.m` with your geometry file and parameters
3. Run `prepro_2D_main.m` with MATLAB

### Prerequisites (optionals)

**NURBS Toolbox (for geometry refining and plotting):**
#### Installation

Clone the repository into your MATLAB Add-Ons (or any preferred) folder; this version is updated comparing to the one from MATLAB Add-On Explorer:

```bash
git clone https://github.com/Institute-for-Risk-and-Reliability/NURBS-Toolbox.git
```

#### Setup in MATLAB

Add the toolbox to your MATLAB path (replace `path/to` with the actual location where you cloned the repo):

```matlab
addpath(genpath('path/to/NURBS-Toolbox'))
```

**IGES Support (only for iges geometries):**
#### Setup in MATLAB

- Install "IGES Toolbox" from MATLAB Add-On Explorer
- Navigate to addons folder in MATLAB command window 
```matlab
cd path/to/MATLAB Add-Ons/Toolboxes/IGES Toolbox/
```
- Compile MEX files: 
```matlab
makeIGESmex
```

### Supported File Formats
- **`.dat`** - Rhino IGA plugin format (native)
- **`.iges/.igs`** - IGES standard (requires IGES Toolbox)

### Configuration Parameters
Edit `config.m`:
- `geom_file` - Input geometry file path
- `inp_file` - Output Abaqus file path
- `deg_u, deg_v` - NURBS degrees for refinement
- `ref_u, ref_v` - Knot span subdivision factors
- `E, nu, ro, thk` - Material properties

### Output Files
- `inp_file` - Complete Abaqus input file in `outputs/` directory
- `KnotsWeights.dat` - Auxiliary NURBS data for Abaqus UEL in `outputs/` directory


## Solver 
Run Abaqus with UEL

## Python post-processing
Post-process results using Python scripts 
