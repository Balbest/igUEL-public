% MATLAB preprocessing pipeline for Abaqus IGA UEL.
% Read and refine 2D plane or shell nurbs geometry and 
% create input file for IGA Abaqus, 1-patch domain
% (c)2017,2025 Viacheslav Balobanov

clear;
configs = config(); % Read configuration parameters

%---------------------------Read geometry file-----------------------------
% Read geometry file
% Automatically detects file type (.dat, .iges, .igs) and imports NURBS data
% Output structure nurbs is compatible with NURBS toolbox.

nurbs = geometryImport(configs.geom_file); 

%--------------------------------Refinement--------------------------------
% Degree elevation and h-refinement based on config settings
if configs.refinement
    nurbs = refineSurface(nurbs,configs);
end

% Visualize refined geometry if requested

% Determine how many plots are needed
nTiles = configs.plotKnots + configs.plotCPs;
% Only create a layout if at least one plot is needed
if nTiles > 0
    t = tiledlayout(1, nTiles); 
    % Plot knot mesh if requested
    if configs.plotKnots
        nexttile
        nrbkntplot(nurbs)
        title('Knot Mesh')
    end
    % Plot control points if requested
    if configs.plotCPs
        nexttile
        nrbctrlplot(nurbs)
        title('Grid of Control Points')
    end
end
%---------------------------Files generation-------------------------------
% Create outputs directory if it doesn't exist
if ~exist('analysis_input', 'dir')
    mkdir('analysis_input');
end
% Generate IGA element connectivity and knot multiplicities
[INP_Elements,KnMultU,KnMultV] = InputElements(nurbs); 

% Generate Abaqus .inp with nodes, elements, materials, BCs
INP_file_2D (nurbs, INP_Elements, configs)

% Generate Ks_Ws.dat containing knot vectors and weights for UEL
fileID = fopen('analysis_input/Ks_Ws.dat', 'w');

fprintf(fileID, '%.16g ', nurbs.knots{1});
fprintf(fileID, '\n');
fprintf(fileID, '%.16g ', nurbs.knots{2});
fprintf(fileID, '\n');

% Extract and write CP weights
Weights = squeeze(nurbs.coefs(4, :, :));
fprintf(fileID, '%.16g ', Weights);
fclose(fileID);

fprintf('analysis_input/Ks_Ws.dat file created successfully\n');
