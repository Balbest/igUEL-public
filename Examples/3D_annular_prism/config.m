% Configuration file for 3D_annular_prism problem

function configs = config()

% Input geometry files (two surfaces for 3D solid)
% The surfaces must be parametrically identical
configs.geom_file1 = "InitGeom_quarter_annulus_bot.dat"; % Lower surface
configs.geom_file2 = "InitGeom_quarter_annulus_top.dat"; % Upper surface

% Output filename for Abaqus inp
configs.inp_file = "annular_prism.inp";

% Do you want to refine the geometry? 
% Both surfaces are refined equally
configs.refinement = true;

% New polynomial degrees (p-refinement)
% If new degrees <= initial degrees, nothing changes
configs.deg_u = 3;
configs.deg_v = 3;
configs.deg_w = 3; 

% Number of knot spans used to divide each knot span (h-refinement)
configs.ref_u = 16;
configs.ref_v = 4;
configs.num_el_w = 4; 

% Do you want to plot the refined geometry? 
% The process can be time-consuming for large models
configs.plotKnots = true; % plot knot mesh (disabled by default for 3D)
configs.plotCPs   = true; % plot control points

% Number of integration points per knot span in each direction
configs.intP_u = configs.deg_u + 1;
configs.intP_v = configs.deg_v + 1;
configs.intP_w = configs.deg_w + 1;

% Number of elements to divide each knot span in postprocessing
configs.numPP_u = 4;
configs.numPP_v = 4;
configs.numPP_w = 4;

% Material properties
configs.E = 200000;  % Young's modulus [MPa]
configs.nu = 0.3;   % Poisson ratio
configs.g = 0.5;       % Strain gradient parameter (used only for SGE)

% Body force vector
configs.BF_vector = [0, 0, -100]; % magnitudes in global x,y,z directions
% Other types of Body or Boundary Forces may require significant changes in 
% Fortran files, % see Examples/3D_cubic_prism/UEL_IGA_3D_SGE.f

% Boundary conditions 
% 0 = free boundary, 1 = fixed displacement, 
% 2 = fixed normal gradient of displacement
configs.BCs = [1, 0, 0, 0, 0, 0]; 
% for boundaries in parametric space: 
% [left u=0, right u=1, bottom v=0, top v=1, front w=0, back w=1q]
% Other types of BC require manual edits of the final .inp file

% Element output parameter (for the final simulation results): 
% 0 for classical displacements U, stresses S, strains E; 
% 1 for classical  U, S, E and gradient gS, gE 
% Note that El_Output = 1 may significantly increase simulation and postprocessing time
configs.El_Output = 1;

end
