% Combined geometry import function for .dat and .iges files
% Automatically detects file format and imports NURBS surface data
% (c) 2017,2025 V.Balobanov & S.Khakalo
% 
% Input:  geom_file - path to geometry file (.dat or .iges/.igs)
% Output: nurbs     - structure containing NURBS surface data
%         .U, .V    - knot vectors
%         .p, .q    - polynomial degrees
%         .CP       - control points array (4, NofCPs_U, NofCPs_V)
%
% .dat format file created by IGA-plugin for Rhinoceros
% 
% .iges import requires iges2matlab toolbox. 
% Install it as MATLAB Add-On and compile mex-files
% (c-compiler is needed for this purpose)

function nurbs = geometryImport(geom_file)

% Validate input file existence 
if ~exist(geom_file, 'file')
    error('Geometry file "%s" not found', geom_file);
end

% Extract file extension
[~, ~, ext] = fileparts(geom_file);
ext = lower(ext);

% Import based on file type
switch ext
    case '.dat'
        nurbs = import_dat(geom_file);
    case {'.iges', '.igs'}
        nurbs = import_iges(geom_file);
    otherwise
        error(['Unsupported geometry file format: %s. ' ...
               'Supported formats: .dat, .iges, .igs'], ext);
end

end

%--------------------------------------------------------------------------
% Import .dat file (Rhino IGA plugin format)
%--------------------------------------------------------------------------
function nurbs = import_dat(geom_file)

fileID = fopen(geom_file, 'r');
if fileID == -1
    error('Cannot open .dat file: %s', geom_file);
end

try
    % Skip first line (header)
    fgetl(fileID);
    
    % Read polynomial degrees
    tline = str2num(fgetl(fileID));
    if length(tline) ~= 2
        error('Invalid .dat file format: missing polynomial degrees');
    end
    p = tline(1);
    q = tline(2);
    
    % Read number of control points
    tline = str2num(fgetl(fileID));
    if length(tline) ~= 2
        error('Invalid .dat file format: missing number of control points');
    end
    NofCPs_U = tline(1);
    NofCPs_V = tline(2);
    NofCPs = NofCPs_U * NofCPs_V;
    
    % Read knot vectors
    xi = str2num(fgetl(fileID));
    eta = str2num(fgetl(fileID));
    
    if isempty(xi) || isempty(eta)
        error('Invalid .dat file format: empty knot vectors');
    end
    
    % Read control points
    CPs = zeros(NofCPs, 4);
    i = 0;
    while ~feof(fileID)
        tline = str2num(fgetl(fileID));
        if length(tline) == 4
            i = i + 1;
            CPs(i, :) = tline;
        end
    end
    
    % Validate control points count
    if i ~= NofCPs
        error('Incorrect number of control points: expected %d, got %d', NofCPs, i);
    end
    
    % Transpose for further processing
    CPs = CPs';
    
    % Convert to homogeneous coordinates
    for j = 1:size(CPs, 2)
        CPs(1, j) = CPs(1, j) * CPs(4, j);
        CPs(2, j) = CPs(2, j) * CPs(4, j);
        CPs(3, j) = CPs(3, j) * CPs(4, j);
    end
    
    % Reshape to 3D array for Nurbs toolbox format
    CPs = reshape(CPs, [4, NofCPs_U, NofCPs_V]);
    
    % Normalize knot vectors to [0,1] range
    U = (xi - xi(1)) ./ (xi(end) - xi(1));
    V = (eta - eta(1)) ./ (eta(end) - eta(1));

    % Create NURBS structure 
    nurbs = struct ('form', 'B-NURBS', ...
                    'dim', 4, ...
                    'number', [NofCPs_U,NofCPs_V], ...
                    'coefs', CPs, ...
                    'knots', {{U,V}}, ...
                    'order', [p+1,q+1]);
    
catch ME
    fclose(fileID);
    rethrow(ME);
end

fclose(fileID);

end

%--------------------------------------------------------------------------
% Import .iges file (IGES standard format)
%--------------------------------------------------------------------------
function nurbs = import_iges(geom_file)

% Check if iges2matlab toolbox is available
if ~exist('iges2matlab', 'file')
    error(['IGES import requires iges2matlab toolbox. ', ...
           'Install it as MATLAB Add-On and compile mex-files.']);
end

try
    % Load parameter data from IGES file
    [ParameterData,~,~,~,~] = iges2matlab(geom_file);
    % Find NURBS surface entities (type 128)
    Ptrs = findEntityIGES(ParameterData, 128);
    
    if isempty(Ptrs)
        error('No NURBS surfaces (entity type 128) found in IGES file');
    end
    
    if length(Ptrs) > 1
        warning(['Multiple NURBS surfaces found in IGES file. ' ...
                 'Using first surface only.']);
    end
    
    % Extract first surface
    Surf = ParameterData{Ptrs(1)};
    
    % Extract control points
    nurbs = Surf.nurbs;
    xi  = nurbs.knots{1};
    eta = nurbs.knots{2};
    % Normalize knot vectors to [0,1] range
    nurbs.knots{1} = (xi - xi(1)) ./ (xi(end) - xi(1));
    nurbs.knots{2} = (eta - eta(1)) ./ (eta(end) - eta(1));


catch ME
    error('Failed to import IGES file: %s', ME.message);
end

end

