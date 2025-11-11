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
% Generate IGA element connectivity and knot multiplicities
[INP_Elements,KnMultU,KnMultV] = InputElements(nurbs); 

% Generate Abaqus .inp with nodes, elements, materials, BCs
INP_file_sh (nurbs, INP_Elements, configs)

% Calculate knot spans for UEL NURBS basis
NoElU = length(unique(nurbs.knots{1})) - 1; % Number of knot spans in U direction
NoElV = length(unique(nurbs.knots{2})) - 1; % Number of knot spans in V direction
KnSpansU = zeros(1, NoElU); 
KnSpansU(1) = KnMultU(1);
KnSpansV = zeros(1, NoElV); 
KnSpansV(1) = KnMultV(1);

for i = 2:NoElU
    KnSpansU(i) = KnSpansU(i-1) + KnMultU(i);
end
for i = 2:NoElV
    KnSpansV(i) = KnSpansV(i-1) + KnMultV(i);
end

% Generate KnotsWeights.dat containing knot vectors and weights for UEL
fileID = fopen('outputs/KnotsWeights.dat', 'w');

fprintf(fileID, '%.15f ', nurbs.knots{1});
fprintf(fileID, '\n');
fprintf(fileID, '%d ', KnSpansU);
fprintf(fileID, '\n');
fprintf(fileID, '%.15f ', nurbs.knots{2});
fprintf(fileID, '\n');
fprintf(fileID, '%d ', KnSpansV);
fprintf(fileID, '\n');

% Extract and write CP weights
Weights = squeeze(nurbs.coefs(4, :, :));
fprintf(fileID, '%.15f ', Weights);
fclose(fileID);

fprintf('outputs/KnotsWeights.dat file created successfully\n');
