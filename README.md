# igUEL

Abaqus User ELement (UEL) implementation for Isogeometric Analysis. 

Developed to solve two- and three-dimensional problems of classical and strain gradient elasticity.

The software consists of three modules:
<ol>
  <li>Module 1 (MATLAB): Preprocessing - geometry import, refinement, and Abaqus input file generation</li>
  <li>Module 2 (Fortran): Solution - UEL implementation for strain gradient elasticity</li>
  <li>Module 3 (Abaqus Python): Post-processing - result visualization in Abaqus</li>
</ol>

---

## Module 1: MATLAB Preprocessing Pipeline

### Quick Start

```bash
cd 2D_PrePro_matlab/    # or 3D_PrePro_matlab/ for 3D problems
matlab -batch "prepro_2D_main"
```
Or run `prepro_2D_main.m` interactively in MATLAB

### Optional Dependencies

#### NURBS Toolbox (required for geometry refinement and plotting):

1. Clone the following repository into your MATLAB Add-Ons folder:
```bash
git clone https://github.com/Institute-for-Risk-and-Reliability/NURBS-Toolbox.git
```

2. Add it to MATLAB path:
```matlab
addpath(genpath('path/to/NURBS-Toolbox'))
```

**Note**: This GitHub version is more up-to-date than the MATLAB Add-On Explorer version of the NURBS package.

#### IGES Toolbox (only required for `.iges/.igs` file import):

1. Install "IGES Toolbox" from MATLAB Add-On Explorer
2. Navigate to the toolbox folder and compile MEX files:
```matlab
cd('path/to/MATLAB Add-Ons/Toolboxes/IGES Toolbox/')
makeIGESmex
```

**Note**: C compiler is required for MEX compilation.

### Supported Geometry File Formats

- `.dat` - Simple text format (Rhinoceros IGA-plugin compatible)
  - Can be created manually or exported from Rhinoceros CAD software
  - Contains, by lines (for 2D or 3D, accordingly):
    - Header (any text line)
    - Polynomial Degrees (2 or 3 numbers with spaces between)
    - Number of Control Points (2 or 3 numbers with spaces between)
    - 2 or 3 lines with Knot Vectors
    - Lines with Control Points (1 line per CP): 4 values (x, y, z, weight). z=0 for 2D
    - One last line is ignored (heritage of IGA-plugin)

- `.iges` / `.igs` - IGES standard format
  - Requires IGES Toolbox
  - Must contain NURBS surface entity (type 128)
  - If multiple surfaces exist, only the first one is imported

**Note**: For 3D case, two files are required. They should contain the parametrically and geometrically identical initial 2D flat geometries and will be refined equally. 
The 3D structure is created between these flat geometries, which represent the "bottom" and "top" of the structure and should therefore be positioned accordingly in space.

### Configuration

Edit `config.m` to set file names, refinement, analysis and postprocessing options, and boundary conditions.

### Output Files

Generated in `analysis_input/` directory:
- `<inp_file>.inp` - Complete Abaqus input file containing:
  - Control point coordinates (*NODE)
  - User element definitions (*USER ELEMENT)
  - Element connectivity (*ELEMENT)
  - Material properties and other problem information (*UEL PROPERTY)
  - Boundary node sets (*NSET)
  - Analysis step (*STEP) with boundary conditions (*BOUNDARY)


- `Ks_Ws.dat` - NURBS data for UEL:
  - Line 1: Knot vector in u-direction
  - Line 2: Knot vector in v-direction
  - Line 3: Control point weights
  - For 3D, line 3 is for knot vector in w-direction, and weights are shifted to line 4

### Preprocessor components

- `geometryImport.m` auto-detects format and loads NURBS data
- `refineSurface.m` performs p- and h-refinement (maintains C^(p-1) continuity)
- `InputElements.m` generates connectivity arrays based on knot spans
- `INP_file_2D.m` generates creates Abaqus .inp file
- `prepro_2D_main.m` Main governing file, also generates `Ks_Ws.dat` file

---

## Module 2: Fortran UEL Solution

### Running the Solver

Execute Abaqus with the UEL:

```bash
#For 2D:
cd 2D_Solver/
abaqus job=<jobname> user=UEL_IGA_2D_SGE cpus=<ncpus> interactive
```

```bash
#For 3D:
cd 3D_Solver/
abaqus job=<jobname> user=UEL_IGA_3D_SGE cpus=<ncpus> interactive
```

### UEL Components

The UEL implementation includes the following subroutines:

- 2D UEL (`UEL_IGA_2D_SGE.f`):

  - `GaussLegendre.f` - Gauss quadrature points and weights
  - `Mat_Prop_2D.f` - Material property matrices (classical + gradient)
  - `NURBS_BF_2D.f` - NURBS basis functions and derivatives
  - `Strains_2D.f` - Strain and strain gradient computation
  - `Force_Stiff_2D.f` - Element force vector and stiffness matrix
  - `Output_2D.f` - Post-processing data generation

- 3D UEL (`UEL_IGA_3D_SGE.f`): Similar structure with 3D-specific implementations.

### Solution Output

Generated files in ./results/`:

Standard Abaqus output database `<jobname>.odb` could be ignored, however `<jobname>.dat` and other files still contain actual analysis logs.

- `results/` folder containing:
  - `Frames_description.dat`
  - `Nodes.dat`, `Nodes_Coords.dat`
  - `Elements.dat`, `Elements_Nodes.dat`
  - `U_Nodes.dat` - Displacement fields
  - `S_Nodes.dat` - Stress fields
  - `E_Nodes.dat` - Strain fields
  - `gSx_Nodes.dat`, `gSy_Nodes.dat` - Double stress (if `El_Output=1`)
  - `gEx_Nodes.dat`, `gEy_Nodes.dat` - Strain gradient (if `El_Output=1`)

---

## Module 3: Python Post-Processing

### Running Post-Processing

Generate a complete Abaqus ODB for visualization:

```bash
abaqus cae noGUI=PostProc_2D.py #2D Problems
```
```bash
abaqus cae noGUI=PostProc_3D.py #3D Problems
```
**Note:** The Python postprocessor generates extra PP files in `./results/`, stored to simplify troubleshooting.

### Post-Processing Output

- `results/2D.odb` or `results/3D.odb` - Complete Abaqus output database
  - Can be opened in Abaqus CAE for visualization as usually
  - Contains displacement, stress, and strain fields. Optionally: gradients of strains and stresses
  - Uses built-in CPE4 (2D) or C3D8 (3D) elements for visualization

---

## Examples

The `Examples/` folder contains:

- `2D_strip/` - Rectangular strip in shear (with analytical solution)
- `3D_square_prism/` - Square prism under body forces (with analytical solution)
- `3D_annular_prism/` - Annular prism in bending

Each example includes:

- Geometry files (`.dat`)
- Configuration file (`config.m`)
- Abaqus `<inputfile>.inp` file and auxiliary `Ks_Ws.dat`
- For some examples, there are manually edited versions of the `<inputfile>.inp` files. They contain some Node Sets and Boundary Conditions which are not included into the current MATLAB Preprocessor.
- `3D_square_prism/` also contains modified `UEL_IGA_3D_SGE.f`, which contains definitions of Body Forces as functions of coordinates. To run the simulation correctly, use this file instead of the default one in  `3D_SGE_igUEL/`.
- Folder `results/` with files generated during the solution and postprocessing. Final `.odb` is included only for 2D case, and excluded for 3D due to Github file size limitations.

---

## Complete 2D Workflow Example

```bash
# Step 1: Preprocessing (MATLAB)
cd 2D_PrePro_matlab/
# Edit config.m with your parameters
matlab -batch "prepro_2D_main"
# Copy resulting <inputfile>.inp and Ks_Ws.dat from 2D_PrePro_matlab/output to ../2D_SGE_igUEL/

# Step 2: Solution (Abaqus + Fortran UEL)
cd ../2D_SGE_igUEL/
abaqus job=<inputfile> user=UEL_IGA_2D_SGE cpus=4 interactive

# Step 3: Post-processing (Abaqus Python)
abaqus cae noGUI=PostProc_2D.py

# Step 4: Final result visualization (Abaqus CAE)
abaqus cae database=2D.odb
```

---

## Authors

- Sergei Khakalo (Aalto University) - sergei.khakalo@aalto.fi
- Viacheslav Balobanov (VTT) - vbalobanov@gmail.com

---

## Version History

- **v0.1.0** (2025) - Initial public release (beta)
- **v1.0.0** (2025) - Publication (SoftwareX) release
