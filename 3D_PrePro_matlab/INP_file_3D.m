% Create 3D solid .inp file for Abaqus IGA UEL (igUEL)
% (c)2017,2025 Sergei Khakalo, Viacheslav Balobanov
function INP_file_3D(nurbs1, INP_Elements, INP_CPs_Coords, zeta, configs)

%-------------------------Parameters definitions---------------------------
data_names = fieldnames (configs);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= configs.(data_names{iopt});']);
end

% Extract NURBS parameters from first surface
p = nurbs1.order(1)-1;
q = nurbs1.order(2)-1;
xi = nurbs1.knots{1};
eta = nurbs1.knots{2};

NoCPsU = nurbs1.number(1);
NoCPsV = nurbs1.number(2);
NoCPsW = length(zeta) - deg_w - 1;

% Element counts
NofEl_U = length (uniquetol(xi)) - 1;
NofEl_V = length (uniquetol(eta)) - 1;
NofEl_W = num_el_w;
NofEl_2D = NofEl_U * NofEl_V;

% 3D element properties
Num_of_Prop = 23;
Num_of_Coords = 3;
ActiveDOF = '1, 2, 3';
CPs_Per_EL = (p+1)*(q+1)*(deg_w+1);

NoKnots_u = length(xi);
NoKnots_v = length(eta);
Num_of_Knots_W = length(zeta);

NoCPs = size(INP_CPs_Coords, 2);

%-----------------------Header and CPs-coordinates-------------------------
heading = 'MATLAB generated Abaqus Input File for 3D igUEL subroutine ';
fileID = fopen(strcat('analysis_input/', inp_file),'w');
fprintf(fileID,'%s\n%s\n\n%s\n%s\n','*HEADING',heading,'**','*NODE,NSET=ALLNODES');

% Print the control points (number x y z)
formatSpec = '%d, %.16g, %.16g, %.16g\n';
fprintf(fileID,formatSpec,INP_CPs_Coords);

%-----------------------------Elements-------------------------------------
% Elements printing format for 3D
F_line = [];
formatSpec = '%d,\n'; % First line := No of Element
for i = 1:p+1
    F_line = [F_line '%d,']; % Format for the single line
end
F_line = [F_line '\n']; % Add end of line

for i = 1:(q+1)*(deg_w+1)
    formatSpec = [formatSpec F_line];
end
% Remove the last comma from the last row
formatSpec = [formatSpec(1:end-3) formatSpec(end-1:end)];

% Elements: UEL definition
fprintf(fileID,'%s\n%s%d%s%s%d%s%d\n%s\n',...
'**','*USER ELEMENT, NODES=',CPs_Per_EL,', TYPE=U3',...
', PROPERTIES=',Num_of_Prop,', COORDINATES=',Num_of_Coords,...
ActiveDOF);
fprintf(fileID,'%s\n%s%s\n','**','*ELEMENT, TYPE=U3',', ELSET=ELSET1');

% Print elements
fprintf(fileID,formatSpec,INP_Elements);

%-----------------------------PID properties-------------------------------
fprintf(fileID,...
['%s\n%s \n' ...
'%d, %d, %d, %d, %d, %d, %d, %d, \n' ...
'%d, %d, %d, %d, %d, %d, %d, %d, \n' ...
'%d, %.16g, %.16g, %.16g, %.16g, %.16g, %.16g \n'],...
'**','*UEL PROPERTY,ELSET=ELSET1',...
intP_u, intP_v, intP_w, NofEl_U, NofEl_V, NofEl_W, NoCPs, p, ...
q, deg_w, NoKnots_u, NoKnots_v, Num_of_Knots_W, numPP_u, numPP_v, numPP_w, ...
El_Output, BF_vector(1), BF_vector(2), BF_vector(3), E, nu, g);

%--------------------------Boundary surfaces NSETs-------------------------
% Calculate 2D layer size
NoCPs2D = NoCPsU * NoCPsV;

% Front and back surfaces (W-direction boundaries)
NSET_back1 = 1 : NoCPs2D;
NSET_front1 = NoCPs-NoCPs2D+1 : NoCPs;

NSET_back2 = NoCPs2D+1 : 2*NoCPs2D;
NSET_front2 = NoCPs-2*NoCPs2D+1 : NoCPs-NoCPs2D;

% Bottom and top surfaces (V-direction boundaries)
aux1 = 1 : NoCPsU;
aux2 = NoCPs2D-NoCPsU+1 : NoCPs2D;
NSET_bottom1 = aux1;
NSET_top1 = aux2;
for i = 2:NoCPsW
    aux1 = aux1 + NoCPs2D;
    NSET_bottom1 = [NSET_bottom1 aux1];
    aux2 = aux2 + NoCPs2D;
    NSET_top1 = [NSET_top1 aux2];
end
NSET_bottom2 = NSET_bottom1 + NoCPsU;
NSET_top2 = NSET_top1 - NoCPsU;

% Left and right surfaces (U-direction boundaries)
aux1 = 1:NoCPsU:NoCPs2D;
aux2 = NoCPsU:NoCPsU:NoCPs2D;
NSET_left1 = aux1;
NSET_right1 = aux2;
for i = 2:NoCPsW
    aux1 = aux1 + NoCPs2D;
    NSET_left1 = [NSET_left1 aux1];
    aux2 = aux2 + NoCPs2D;
    NSET_right1 = [NSET_right1 aux2];
end
NSET_left2 = NSET_left1 + 1;
NSET_right2 = NSET_right1 - 1;

% Writing of the NSETs
fprintf(fileID,'%s\n%s\n',...
'**','** --- Max number of nodes or elements in line is 10!');

% Write all boundary node sets
boundary_sets = {
    'NSET_left1', NSET_left1;
    'NSET_left2', NSET_left2;
    'NSET_right1', NSET_right1;
    'NSET_right2', NSET_right2;
    'NSET_bottom1', NSET_bottom1;
    'NSET_bottom2', NSET_bottom2;
    'NSET_top1', NSET_top1;
    'NSET_top2', NSET_top2;
    'NSET_back1', NSET_back1;
    'NSET_back2', NSET_back2;
    'NSET_front1', NSET_front1;
    'NSET_front2', NSET_front2
};

for i = 1:size(boundary_sets, 1)
    fprintf(fileID,'%s\n*NSET, NSET=%s\n','**',boundary_sets{i,1});
    fprintf(fileID,formatSpecNset(numel(boundary_sets{i,2})),boundary_sets{i,2});
end

%---------------------------------STEP-------------------------------------
% writing STEP and BCs
if numel(configs.BCs) ~= 6
    error('BCs vector should have exactly 4 elements, check config file!');
end

fprintf(fileID,'%s\n%s\n%s\n%s\n%s\n', ...
'**', ...
'*STEP,NAME=STEP1', ...
'*STATIC', ...
'**',...
'*BOUNDARY, TYPE=DISPLACEMENT');
writeBCs('left',configs.BCs(1),fileID)
writeBCs('right',configs.BCs(2),fileID)
writeBCs('bottom',configs.BCs(3),fileID)
writeBCs('top',configs.BCs(4),fileID)
writeBCs('front',configs.BCs(5),fileID)
writeBCs('back',configs.BCs(6),fileID)
fprintf(fileID, '*END STEP\n');

fclose(fileID); % End of .INP file

% Print success message
fprintf('%s file created successfully\n', inp_file);

end

%--------function define format for Set of Nodes printing------------------
function formatSpec = formatSpecNset(SetNumEl)
formatSpec = [];
for i = 1:SetNumEl-1
    if mod(i,10) ~= 0
        formatSpec = [formatSpec '%d, '];
    else
        formatSpec = [formatSpec '%d,\n'];
    end
end
formatSpec = [formatSpec '%d\n'];
end

%-------------Defining format for BCs printing-----------------------------
function writeBCs(Boundary, BC_type, fileID)
%BC_type=0 -- free boundary
if BC_type == 1 %zero displacement
    fprintf(fileID, 'NSET_%s1,1,3,0\n', Boundary);
elseif BC_type == 2 %zero normal gradient of displacement
    fprintf(fileID, 'NSET_%s1,1,3,0\n', Boundary);
    fprintf(fileID, 'NSET_%s2,1,3,0\n', Boundary);
elseif BC_type ~= 0 %Free boundary -- do nothing
    error('Invalid value for BC type, only 0, 1, or 2 allowed');
end
end
