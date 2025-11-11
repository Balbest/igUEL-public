% Configuration file for 3D preprocessing
% (c)2025 V. Balobanov

function configs = config()

% Input geometry files (two surfaces for 3D solid)
% The surfaces must be parametrically identical
configs.geom_file1 = "examples/Beam100-2-1_lower.dat"; % Lower surface
configs.geom_file2 = "examples/Beam100-2-1_upper.dat"; % Upper surface

% Output filename for Abaqus inp
configs.inp_file = "outputs/Beam100-2-1_100v4v2pqr3.inp";

% Do you want to refine the geometry? 
% Both surfaces are refined equally
configs.refinement = true;

% New polynomial degrees (p-refinement)
% If new degrees <= initial degrees, nothing changes
configs.deg_u = 3;
configs.deg_v = 3;
configs.deg_w = 3; % 3rd direction polynomial degree

% Number of knot spans used to divide each knot span (h-refinement)
configs.ref_u = 10;
configs.ref_v = 4;
configs.num_el_w = 2; % Number of elements in 3rd direction

% Do you want to plot the refined geometry? 
% The process can be time-consuming for large models
configs.plotKnots = false; % plot knot mesh (disabled by default for 3D)
configs.plotCPs   = false; % plot control points

% Number of integration points (based on maximum degree)
configs.intP = max([configs.deg_u, configs.deg_v, configs.deg_w]) + 1;

% Number of elements to divide each knot span in postprocessing
configs.numPP_u = 4;
configs.numPP_v = 4;
configs.numPP_w = 2;

% ELSET indicator for boundary conditions
configs.elset_indicator = 0;

% Type of problem: 1 - eigen value problem, 2 - static problem
configs.problem_type = 2; 

% Material properties
configs.E = 210000;  % Young's modulus [MPa]
configs.nu = 0.33;   % Poisson ratio
configs.ro = 1e-9;   % Mass density [kg/mm³] (used only for eigenvalue problems)
configs.g = 0;       % Strain gradient parameter (used only for SGE)

end