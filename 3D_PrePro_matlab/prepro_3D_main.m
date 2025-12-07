% MATLAB 3D preprocessing pipeline for Abaqus IGA UEL
% Generate 3D solid geometry between two surfaces (1-patch domain)
% (c)2017,2025 V. Balobanov and S. Khakalo

clear;
configs = config(); % Read configuration parameters

% Path to 2D functions: geometryImport, refineSurface, inputElements
addpath(genpath('../2D_PrePro_matlab/')) 

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
        title('Grid of Control Points - Surface 1')
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

% Initialize coordinate arrays for the 1st surface
CP1 = zeros(4, NoCPs2D);
CP2 = zeros(4, NoCPs2D);
for j = 1:NoCPsV
    for i = 1:NoCPsU
        NoCP = NoCPsU*(j-1) + i;
        CP1(1,NoCP) = NoCP;
        % Convert CP coordinates back to Cartesian CS
        CP1(2:4,NoCP) = nurbs1.coefs(1:3,i,j) ./ nurbs1.coefs(4,i,j);
        CP2(2:4,NoCP) = nurbs2.coefs(1:3,i,j) ./ nurbs2.coefs(4,i,j);
    end
end

%---------------------------3rd direction setup---------------------------
% Generate knot vector and control points in W-direction
[zeta, z_CPs] = nurbs3rdDir(configs.deg_w, configs.num_el_w);
NoCPsW = length(z_CPs);

% Create 3D weights vector
Vector_of_Weights = squeeze(nurbs1.coefs(4, :, :));
Weights = Vector_of_Weights;
for i = 2:NoCPsW
    Weights = [Weights Vector_of_Weights];
end

%-------------------------Generate 3D control points----------------------
% Fill the volume between two surfaces with CPs
INP_CPs_Coords_3D = zeros(4, NoCPs2D * NoCPsW);

% 1st layer (surface 1)
INP_CPs_Coords_3D(:, 1:NoCPs2D) = CP1;

% Intermediate layers
for i = 1:NoCPs2D
    for j = 2:NoCPsW-1
        idx = i + NoCPs2D*(j-1);
        INP_CPs_Coords_3D(1,idx) = idx;
        INP_CPs_Coords_3D(2:4,idx) = CP1(2:4,i) + ...
            (CP2(2:4,i) - CP1(2:4,i)) * z_CPs(j);
    end

    % Last layer (surface 2)
    idx = i + NoCPs2D*(NoCPsW-1);
    INP_CPs_Coords_3D(1,idx) = idx;
    INP_CPs_Coords_3D(2:4,idx) = CP2(2:4,i);
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

%---------------------------Files generation-------------------------------
% Create outputs directory if it doesn't exist
if ~exist('analysis_input', 'dir')
    mkdir('analysis_input');
end

% Create Ks_Ws.dat
fileID = fopen('analysis_input/Ks_Ws.dat','w');
formatSpec = '%.16g ';
fprintf(fileID,formatSpec,xi);
fprintf(fileID,'\n');
fprintf(fileID,formatSpec,eta);
fprintf(fileID,'\n');
fprintf(fileID,formatSpec,zeta);
fprintf(fileID,'\n');
fprintf(fileID,formatSpec,Weights);
fclose(fileID);

fprintf('analysis_input/Ks_Ws.dat file created successfully\n');

% Call the INP file generation function
INP_file_3D(nurbs1, INP_Elements, INP_CPs_Coords_3D, zeta, configs);