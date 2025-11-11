% Configuration file
% (c)2025 V. Balobanov

function configs = config()

% Input geometry file, supported formats: dat, iges
configs.geom_file = "examples/HypPar_Pure_geometry.dat";%"Scordelis_roof.igs";%

% Output filename for Abaqus inp
configs.inp_file = "configsputs/HypPar_Pure_geometry.inp";


% Do you want to refine the geometry? 
configs.refinement = true;
% New polynomial degree (p-refinement)
% If new degrees <= initial degrees, nothing changes
configs.deg_u = 6;
configs.deg_v = 6;

% Number of knot spans used to divide each knot span of the initial geometry
% (h-refiniment)
configs.ref_u = 3;
configs.ref_v = 5;

% Do you want to plot the refined geometry? 
% The process can be time-consuming for large models
configs.plotKnots = true; % plot knot mesh
configs.plotCPs   = true; % plot control points

% Number of integration points per knot span
configs.intP_u = configs.deg_u + 1;
configs.intP_v = configs.deg_v + 1;

% Number of elements to divide each knot span in postprocessing
% numPP_u*numPP_v
configs.numPP_u = 4;
configs.numPP_v = 4;

% type of problem, 1 - eigen value problem, 2 - static problem
configs.problem_type = 2; 

% material properties
configs.E = 210000; %Young's modulus
configs.nu = 0.33;  %Poisson ratio
configs.ro = 1E-9;  %mass density, used only for eigenvalue problems
configs.g = 0.0;    %strain gradient parameter, used only for SGE
configs.thk = 4.33; %shell thickness, used only for shells 
