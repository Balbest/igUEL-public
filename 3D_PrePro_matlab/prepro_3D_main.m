% MATLAB 3D preprocessing pipeline for Abaqus IGA UEL
% Generate 3D solid geometry between two surfaces (1-patch domain)
% (c)2017,2025 V. Balobanov and S. Khakalo

clear;
configs = config(); % Read configuration parameters
addpath(genpath('../2D_PrePro_matlab/')) % Path to 2D functions used

%---------------------------Read geometry files---------------------------
% Read two surfaces for 3D solid generation
% Automatically detects file type (.dat, .iges, .igs) and imports NURBS data
% Output structure nurbs is compatible with NURBS toolbox.

nurbs1 = geometryImport(configs.geom_file1);
nurbs2 = geometryImport(configs.geom_file2);

%--------------------------------Refinement--------------------------------
% Degree elevation and h-refinement based on config settings
% Both surfaces are refined equally and must be parametrically identical
if configs.refinement
    nurbs1 = refineSurface(nurbs1, configs);
    nurbs2 = refineSurface(nurbs2, configs);
end

% Visualize refined geometry of the 1st surface
nTiles = configs.plotKnots + configs.plotCPs;
if nTiles > 0
    t = tiledlayout(1, nTiles);
    if configs.plotKnots
        nexttile
        nrbkntplot(nurbs1)
        title('Knot Mesh - Surface 1')
    end
    if configs.plotCPs
        nexttile
        nrbctrlplot(nurbs1)
        title('Control Points - Surface 1')
    end
end

%---------------------------Extract NURBS parameters----------------------
% Extract NURBS parameters from the 1st surface
p = nurbs1.order(1)-1;  xi = nurbs1.knots{1};
q = nurbs1.order(2)-1;  eta = nurbs1.knots{2};
NoCPsU = nurbs1.number(1);
NoCPsV = nurbs1.number(2);
NoCPs2D = NoCPsU * NoCPsV;

% Number of elements in 2D
NofEl_U = length(unique(xi)) - 1;
NofEl_V = length(unique(eta)) - 1;
NofEl_2D = NofEl_U * NofEl_V;

% Generate control point coordinates from the 1st surface
CP1 = nurbs1.coefs;
CP2 = nurbs2.coefs;

% Initialize coordinate arrays for the 1st surface
INP_CPs_Coords = zeros(4, NoCPs2D);
for j = 1:NoCPsV
    for i = 1:NoCPsU
        NoCP = NoCPsU*(j-1) + i;
        INP_CPs_Coords(1,NoCP) = NoCP;
        INP_CPs_Coords(2:4,NoCP) = CP1(1:3,i,j);
    end
end

% Store the 2nd surface coordinates in nurbs2 structure for later use
nurbs2.INP_CPs_Coords = zeros(4, NoCPs2D);
for j = 1:NoCPsV
    for i = 1:NoCPsU
        NoCP = NoCPsU*(j-1) + i;
        nurbs2.INP_CPs_Coords(1,NoCP) = NoCP;
        nurbs2.INP_CPs_Coords(2:4,NoCP) = CP2(1:3,i,j);
    end
end

%---------------------------3rd direction setup---------------------------
% Generate knot vector and control points in W-direction
[zeta, z_CPs] = nurbs3rdDir(configs.deg_w, configs.num_el_w);
NoCPsW = length(z_CPs);

% Create 3D weights vector
Vector_of_Weights = ones(1, NoCPs2D); % Assume unit weights for now
Weights = Vector_of_Weights;
for i = 2:NoCPsW-1
    Weights = [Weights Vector_of_Weights];
end
Weights = [Weights Vector_of_Weights]; % Last layer

%-------------------------Generate 3D control points----------------------
% Fill the volume between two surfaces with CPs
INP_CPs_Coords_3D = zeros(4, NoCPs2D * NoCPsW);

% 1st layer (surface 1)
INP_CPs_Coords_3D(:, 1:NoCPs2D) = INP_CPs_Coords;

% Intermediate layers
for i = 1:NoCPs2D
    for j = 2:NoCPsW-1
        idx = i + NoCPs2D*(j-1);
        INP_CPs_Coords_3D(1,idx) = idx;
        INP_CPs_Coords_3D(2:4,idx) = INP_CPs_Coords(2:4,i) + ...
            (nurbs2.INP_CPs_Coords(2:4,i) - INP_CPs_Coords(2:4,i)) * z_CPs(j);
    end

    % Last layer (surface 2)
    idx = i + NoCPs2D*(NoCPsW-1);
    INP_CPs_Coords_3D(1,idx) = idx;
    INP_CPs_Coords_3D(2:4,idx) = nurbs2.INP_CPs_Coords(2:4,i);
end

%------------------------------Elements------------------------------------
% Generate 2D element connectivity 1st
[INP_Elements_2D, KnMultU, KnMultV] = InputElements(nurbs1);

% Extend 2D elements to 3D solids
CPs_per_2D_elem = (p+1)*(q+1);
CPs_per_3D_elem = (p+1)*(q+1)*(configs.deg_w+1);
total_elements = NofEl_2D * configs.num_el_w;

INP_Elements = zeros(CPs_per_3D_elem + 1, total_elements);

% Build 3D elements by stacking 2D elements in W-direction
for el_w = 1:configs.num_el_w
    for el_2d = 1:NofEl_2D
        el_3d = (el_w-1)*NofEl_2D + el_2d;
        INP_Elements(1, el_3d) = el_3d;

        % Add connectivity for each W-layer of the element
        for layer = 1:configs.deg_w+1
            layer_offset = (layer-1) * NoCPs2D + (el_w-1) * NoCPs2D;
            start_idx = 2 + (layer-1) * CPs_per_2D_elem;
            end_idx = 1 + layer * CPs_per_2D_elem;

            INP_Elements(start_idx:end_idx, el_3d) = ...
                INP_Elements_2D(2:end, el_2d) + layer_offset;
        end
    end
end

%-------------------------Generate UEL Data File--------------------------
% Create outputs directory if it doesn't exist
if ~exist('outputs', 'dir')
    mkdir('outputs');
end

% Create KnotsWeights.dat for Fortran UEL with 3D format
fileID = fopen('outputs/KnotsWeights.dat','w');
formatSpec = '%.15f ';
fprintf(fileID,formatSpec,xi);
fprintf(fileID,'\n');
fprintf(fileID,formatSpec,eta);
fprintf(fileID,'\n');
fprintf(fileID,formatSpec,zeta);
fprintf(fileID,'\n');
fprintf(fileID,formatSpec,Weights);
fclose(fileID);

fprintf('outputs/KnotsWeights.dat file created successfully\n');

%-------------------------Generate Abaqus INP file------------------------
% Call the INP file generation function
INP_file_3D(nurbs1, nurbs2, INP_Elements, INP_CPs_Coords_3D, KnMultU, KnMultV, zeta, Weights, configs);

fprintf('3D preprocessing completed successfully\n');