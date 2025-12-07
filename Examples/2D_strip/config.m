% Configuration file for 2D_strip problem

function configs = config()

% Input geometry file, supported formats: dat, iges
configs.geom_file = "InitGeom_Rectangle_1x4.dat";

% Output filename for Abaqus inp
configs.inp_file = "Shear_strip.inp";

% Do you want to refine the geometry? 
configs.refinement = true;
% New polynomial degree (p-refinement)
% If new degrees <= initial degrees, nothing changes
configs.deg_u = 1;
configs.deg_v = 5;

% Number of knot spans used to divide each knot span of the initial geometry
% (h-refiniment)
configs.ref_u = 1;
configs.ref_v = 4;

% Do you want to plot the refined geometry? 
% The process can be time-consuming for large models
configs.plotKnots = true; % plot knot mesh
configs.plotCPs   = true; % plot control points

% Number of integration points per knot span in each direction
configs.intP_u = configs.deg_u + 1;
configs.intP_v = configs.deg_v + 1;

% Number of elements to divide each knot span in postprocessing
% numPP_u*numPP_v
configs.numPP_u = 1;
configs.numPP_v = 10;

% material properties
configs.E = 200000; %Young's modulus
configs.nu = 0.3;  %Poisson ratio
configs.g = 0.4;    %strain gradient parameter, used only for SGE

% Body force vector
configs.BF_vector = [0, 0]; % magnitudes in global x and y directions
% Other types of Body or Boundary Forces may require significant changes in 
% Fortran files, % see Examples/3D_cubic_prism/UEL_IGA_3D_SGE.f

% Boundary conditions 
% 0 = free boundary, 1 = zero displacement, 
% 2 = zero normal gradient of displacement
configs.BCs = [0,0,0,0]; 
% for boundaries in parametric space: 
% [left u=0, right u=1, bottom v=0, top v=1]
% Other types of BC require manual edits of the final .inp file

% Element output parameter (for the final simulation results): 
% 0 for classical displacements U, stresses S, strains E; 
% 1 for classical  U, S, E and gradient gS, gE 
% Note that El_Output = 1 may significantly increase simulation and postprocessing time
configs.El_Output = 1;
