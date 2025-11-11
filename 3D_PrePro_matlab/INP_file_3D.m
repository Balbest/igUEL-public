% Create 3D solid .inp file for Abaqus IGA UEL (igUEL)
% (c)2017,2025 Sergei Khakalo, Viacheslav Balobanov
function INP_file_3D(nurbs1, nurbs2, INP_Elements, INP_CPs_Coords, KnMultU, KnMultV, zeta, Weights, configs)

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
NofEl_U = length (unique(xi)) - 1;
NofEl_V = length (unique(eta)) - 1;
NofEl_W = num_el_w;
NofEl_2D = NofEl_U * NofEl_V;

% 3D element properties
Num_of_Prop = 19;
Num_of_Coords = 3;
Num_of_Var = 3;
ActiveDOF = '1, 2, 3';
CPs_Per_EL = (p+1)*(q+1)*(deg_w+1);
Num_of_IntPoints = intP * intP * intP; % Assuming same integration points in all directions

NoKnots_u = length(xi);
NoKnots_v = length(eta);
Num_of_Knots_W = length(zeta);

NoCPs = size(INP_CPs_Coords, 2);

%-----------------------Header and CPs-coordinates-------------------------
fileID = fopen(inp_file,'w');
fprintf(fileID,'%s\n\n%s\n%s\n','*HEADING','**','*NODE,NSET=ALLNODES');

% Print the control points (number x y z)
formatSpec = '%d, %10.17f, %10.17f, %10.17f\n';
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
fprintf(fileID,'%s\n%s%d%s%s%d%s%d%s%d%s%d\n%s\n',...
'**','*USER ELEMENT, NODES=',CPs_Per_EL,', TYPE=U1',...
', PROPERTIES=',Num_of_Prop,', COORDINATES=',Num_of_Coords,...
', INTEGRATION=',Num_of_IntPoints,', VARIABLES=',Num_of_Var,...
ActiveDOF);
fprintf(fileID,'%s\n%s%s\n','**','*ELEMENT, TYPE=U1',', ELSET=ELSET1');

% Print elements
fprintf(fileID,formatSpec,INP_Elements);

%-----------------------------PID properties-------------------------------
fprintf(fileID,...
'%s\n%s\n%.10f, %.10f, %.10f, %d, %d, %d, %d, %d,\n%d, %d, %d, %d, %d, %d, %d, %d,\n%d, %d, %.10f \n',...
'**','*UEL PROPERTY,ELSET=ELSET1',...
E, nu, ro, Num_of_IntPoints, NofEl_U, NofEl_V, NofEl_W, NoCPs,...
elset_indicator, p, q, deg_w, NoKnots_u, NoKnots_v, Num_of_Knots_W, numPP_u,...
numPP_v, numPP_w, g);

%--------------------------Boundary surfaces NSETs-------------------------
% Calculate 2D layer size
NoCPs2D = NoCPsU * NoCPsV;

% Bottom and top surfaces (W-direction boundaries)
NSET_hm1 = 1 : NoCPs2D;
NSET_hp1 = NoCPs-NoCPs2D+1 : NoCPs;

NSET_hm2 = NoCPs2D+1 : 2*NoCPs2D;
NSET_hp2 = NoCPs-2*NoCPs2D+1 : NoCPs-NoCPs2D;

% Front and back surfaces (V-direction boundaries)
aux1 = 1 : NoCPsU;
aux2 = NoCPs2D-NoCPsU+1 : NoCPs2D;
NSET_nm1 = aux1;
NSET_np1 = aux2;
for i = 2:NoCPsW
    aux1 = aux1 + NoCPs2D;
    NSET_nm1 = [NSET_nm1 aux1];
    aux2 = aux2 + NoCPs2D;
    NSET_np1 = [NSET_np1 aux2];
end
NSET_nm2 = NSET_nm1 + NoCPsU;
NSET_np2 = NSET_np1 - NoCPsU;

% Left and right surfaces (U-direction boundaries)
aux1 = 1:NoCPsU:NoCPs2D;
aux2 = NoCPsU:NoCPsU:NoCPs2D;
NSET_gm1 = aux1;
NSET_gp1 = aux2;
for i = 2:NoCPsW
    aux1 = aux1 + NoCPs2D;
    NSET_gm1 = [NSET_gm1 aux1];
    aux2 = aux2 + NoCPs2D;
    NSET_gp1 = [NSET_gp1 aux2];
end
NSET_gm2 = NSET_gm1 + 1;
NSET_gp2 = NSET_gp1 - 1;

% Writing of the NSETs
fprintf(fileID,'%s\n%s\n',...
'**','** --- Max number of nodes or elements in line is 10!!!!!');

% Write all boundary node sets
boundary_sets = {
    'NSET_gm1', NSET_gm1;
    'NSET_gm2', NSET_gm2;
    'NSET_gp1', NSET_gp1;
    'NSET_gp2', NSET_gp2;
    'NSET_nm1', NSET_nm1;
    'NSET_nm2', NSET_nm2;
    'NSET_np1', NSET_np1;
    'NSET_np2', NSET_np2;
    'NSET_hm1', NSET_hm1;
    'NSET_hm2', NSET_hm2;
    'NSET_hp1', NSET_hp1;
    'NSET_hp2', NSET_hp2
};

for i = 1:size(boundary_sets, 1)
    fprintf(fileID,'%s\n*NSET, NSET=%s\n','**',boundary_sets{i,1});
    fprintf(fileID,formatSpecNset(numel(boundary_sets{i,2})),boundary_sets{i,2});
end

%---------------------------------STEP-------------------------------------
if problem_type == 1
    fprintf(fileID,'%s\n%s\n%s\n%d%s\n%s\n%s\n','**',...
    '*STEP,NAME=STEP1,PERTURBATION','*FREQUENCY, eigensolver=lanczos, sim',...
    NoCPs*3,', ,','*BOUNDARY','*END STEP');
elseif problem_type == 2
    fprintf(fileID,'%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s',...
    '**','*STEP,NAME=STEP1','*STATIC','**',...
    '*BOUNDARY, TYPE=DISPLACEMENT',...
    'NSET_gm1,1,3,0','**NSET_gm2,1,3,0','**NSET_gp1,1,3,-1','**NSET_gp2,1,3,0',...
    '**NSET_nm1,1,3,0','**NSET_nm2,1,3,0','**NSET_np1,1,3,0','**NSET_np2,1,3,0',...
    '**NSET_hm1,1,3,0','**NSET_hm2,1,3,0','**NSET_hp1,1,3,0','**NSET_hp2,1,3,0',...
    '*END STEP');
end

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