% Create .INP file
% (c)2017,2025 Sergei Khakalo, Viacheslav Balobanov
function INP_file_2D (nurbs, INP_Elements, configs)
p = nurbs.order(1)-1; q = nurbs.order(2)-1; CP = nurbs.coefs;

%-------------------------Parameters definitions---------------------------
data_names = fieldnames (configs);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= configs.(data_names{iopt});']);
end

Num_of_Prop = 17;
Num_of_Coords = 2;
ActiveDOF='1, 2';

NoKnots_u = length(nurbs.knots{1});
NoKnots_v = length(nurbs.knots{2});

NofCPs_u = nurbs.number(1);
NofCPs_v = nurbs.number(2);
NofCPs = NofCPs_v*NofCPs_u;

NoEl_u = length (uniquetol(nurbs.knots{1})) - 1;
NoEl_v = length (uniquetol(nurbs.knots{2})) - 1;

CPs_Per_EL=(p+1)*(q+1);

%-----------------------Header and CPs-coordinates-------------------------
fileID = fopen(strcat('analysis_input/', inp_file),'w');
heading = 'MATLAB generated Abaqus Input File for 2D igUEL subroutine ';
fprintf(fileID,'%s\n%s\n\n%s\n%s\n','*HEADING',heading,'**','*NODE,NSET=ALLNODES');
formatSpec = '%d, %.16g, %.16g\n';
for j = 1:NofCPs_v
    for i = 1:NofCPs_u
        NoCP = NofCPs_u*(j-1) + i;
        % CP coordinates converted back to Cartesian CS, 3rd coordinate is
        % ignored for 2D
        fprintf(fileID, formatSpec, NoCP, CP(1:2,i,j)./CP(4,i,j));
    end
end

%-----------------------------Elements-------------------------------------
% Elements printing format
F_line = []; 
formatSpec = '%d,\n'; %first line := No of Element
for i = 1:p+1
    F_line = [F_line '%d,']; %format for the single line
end
F_line = [F_line '\n']; %add end of line

for i = 1:(q+1)
    formatSpec = [formatSpec F_line];
end
%remove the last comma from the last row
formatSpec = [formatSpec(1:end-3) formatSpec(end-1:end)];

% Elements writing
fprintf(fileID,'%s\n%s%d%s%s%d%s%d\n%s\n',...
'**', ...
'*USER ELEMENT, NODES=', CPs_Per_EL,', TYPE=U2',...
', PROPERTIES=',Num_of_Prop,', COORDINATES=',Num_of_Coords,...
ActiveDOF);
fprintf(fileID,'%s\n%s%s\n','**','*ELEMENT, TYPE=U2',', ELSET=ELSET1');

fprintf(fileID,formatSpec,INP_Elements);

%-----------------------------PID properties-------------------------------
fprintf(fileID,...
['%s\n%s\n',...
'%d, %d, %d, %d, %d, %d, %d, %d,\n'...
'%d, %d, %d, %d, %.16g, %.16g, %.16g, %.16g,\n'...
'%.16g\n'],...
'**','*UEL PROPERTY,ELSET=ELSET1',...
intP_u, intP_v, NoEl_u, NoEl_v, NofCPs, p, q, NoKnots_u, ...
NoKnots_v, numPP_u, numPP_v, El_Output, BF_vector(1), BF_vector(2), E, nu, ...
g);
%--------------------------------NSETS-------------------------------------
for i=1:NofCPs_u
    NSET_bottom1(i)=i;
    NSET_bottom2(i)=NofCPs_u+i;
    NSET_top1(i)=i+NofCPs_u*(NofCPs_v-1);
    NSET_top2(i)=i+NofCPs_u*(NofCPs_v-2);
end
for i=1:NofCPs_v
    NSET_left1(i)=1+(i-1)*NofCPs_u;
    NSET_left2(i)=2+(i-1)*NofCPs_u;
    NSET_right1(i)=NofCPs_u+(i-1)*NofCPs_u;
    NSET_right2(i)=NofCPs_u+(i-1)*NofCPs_u-1;
end

% writing NSETs
fprintf(fileID,'%s\n%s\n',...
'**','** NOTE: Max number of nodes or elements per line is 10');

fprintf(fileID,'%s\n%s\n','**','*NSET, NSET=NSET_bottom1');
fprintf(fileID,formatSpecNset(numel(NSET_bottom1)),NSET_bottom1);
fprintf(fileID,'%s\n%s\n','**','*NSET, NSET=NSET_bottom2');
fprintf(fileID,formatSpecNset(numel(NSET_bottom2)),NSET_bottom2);

fprintf(fileID,'%s\n%s\n','**','*NSET, NSET=NSET_top1');
fprintf(fileID,formatSpecNset(numel(NSET_top1)),NSET_top1);
fprintf(fileID,'%s\n%s\n','**','*NSET, NSET=NSET_top2');
fprintf(fileID,formatSpecNset(numel(NSET_top2)),NSET_top2);

fprintf(fileID,'%s\n%s\n','**','*NSET, NSET=NSET_left1');
fprintf(fileID,formatSpecNset(numel(NSET_left1)),NSET_left1);
fprintf(fileID,'%s\n%s\n','**','*NSET, NSET=NSET_left2');
fprintf(fileID,formatSpecNset(numel(NSET_left2)),NSET_left2);

fprintf(fileID,'%s\n%s\n','**','*NSET, NSET=NSET_right1');
fprintf(fileID,formatSpecNset(numel(NSET_right1)),NSET_right1);
fprintf(fileID,'%s\n%s\n','**','*NSET, NSET=NSET_right2');
fprintf(fileID,formatSpecNset(numel(NSET_right2)),NSET_right2);

% writing STEP and BCs
if numel(BCs) ~= 4
    error('BCs vector should have exactly 4 elements, check config file!');
end

fprintf(fileID,'%s\n%s\n%s\n%s\n%s\n', ...
'**', ...
'*STEP,NAME=STEP1', ...
'*STATIC', ...
'**',...
'*BOUNDARY, TYPE=DISPLACEMENT');
writeBCs('left',BCs(1),fileID)
writeBCs('right',BCs(2),fileID)
writeBCs('bottom',BCs(3),fileID)
writeBCs('top',BCs(4),fileID)

fprintf(fileID, '*END STEP\n');

%end of .INP file
fclose(fileID);
fprintf('analysis_input/%s file created successfully\n',inp_file);

end

%-------------Defining format for a Set of Nodes printing------------------
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