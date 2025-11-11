% Create .INP file
% (c)2017,2025 Sergei Khakalo, Viacheslav Balobanov
function INP_file_sh (nurbs, INP_Elements, configs)
p = nurbs.order(1)-1; q = nurbs.order(2)-1; CP = nurbs.coefs;

%-------------------------Parameters definitions---------------------------
data_names = fieldnames (configs);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= configs.(data_names{iopt});']);
end

Num_of_Prop = 17;
Num_of_Coords = 3;
Num_of_IntPoints = intP_u*intP_v;
Num_of_Var = 3;

if Num_of_Coords == 1
   ActiveDOF='1';
elseif Num_of_Coords == 2
   ActiveDOF='1, 2';
elseif Num_of_Coords == 3
   ActiveDOF='1, 2, 3';
end

NoKnots_u = length(nurbs.knots{1});
NoKnots_v = length(nurbs.knots{2});

NofCPs_u = nurbs.number(1);
NofCPs_v = nurbs.number(2);
NofCPs = NofCPs_v*NofCPs_u;

NoEl_u = length (nurbs.knots{1}) - 1;
NoEl_v = length (nurbs.knots{2}) - 1;

CPs_Per_EL=(p+1)*(q+1);

%-----------------------Header and CPs-coordinates-------------------------
fileID = fopen(inp_file,'w');
fprintf(fileID,'%s\n\n%s\n%s\n','*HEADING','**','*NODE,NSET=ALLNODES');
formatSpec = '%d, %.20f, %.20f, %.20f\n';
for j = 1:NofCPs_v
    for i = 1:NofCPs_u
        NoCP = NofCPs_u*(j-1) + i;
        fprintf(fileID, formatSpec, NoCP, CP(1:3,i,j));
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
fprintf(fileID,'%s\n%s%d%s%s%d%s%d%s%d%s%d\n%s\n',...
'**','*USER ELEMENT, NODES=',CPs_Per_EL,', TYPE=U1',...
', PROPERTIES=',Num_of_Prop,', COORDINATES=',Num_of_Coords,...
', INTEGRATION=',Num_of_IntPoints,', VARIABLES=',Num_of_Var,...
ActiveDOF);
fprintf(fileID,'%s\n%s%s\n','**','*ELEMENT, TYPE=U1',', ELSET=ELSET1');

fprintf(fileID,formatSpec,INP_Elements);

%-----------------------------PID properties-------------------------------
fprintf(fileID,...
['%s\n%s\n',...
'%f, %f, %f, %f, %d, %d, %d, %d,\n',...
'%d, %d, %d, %d, %d, %d, %d,  0,\n',...
'%f\n'],...
'**','*UEL PROPERTY,ELSET=ELSET1',...
thk, E, nu, ro, intP_u, intP_v, NoEl_u, NoEl_v,...
NofCPs, p, q, NoKnots_u, NoKnots_v, numPP_u, numPP_v, ...
g);

%--------------------------------NSETS-------------------------------------
for i=1:NofCPs_u
    NSET_nm1(i)=i;
    NSET_nm2(i)=NofCPs_u+i;
    NSET_np1(i)=i+NofCPs_u*(NofCPs_v-1);
    NSET_np2(i)=i+NofCPs_u*(NofCPs_v-2);
end
for i=1:NofCPs_v
    NSET_gm1(i)=1+(i-1)*NofCPs_u;
    NSET_gm2(i)=2+(i-1)*NofCPs_u;
    NSET_gp1(i)=NofCPs_u+(i-1)*NofCPs_u;
    NSET_gp2(i)=NofCPs_u+(i-1)*NofCPs_u-1;
end
%
%
%for i=1:NoElU
%    ELSET_nm1(i)=i;
%    ELSET_np1(i)=i+NoElU*(NoElV-1);
%end
%for i=1:NoElV
%    ELSET_gm1(i)=1+(i-1)*NoElU;
%    ELSET_gp1(i)=NoElU+(i-1)*NoElU;
%end

% writing of the NSETs
fprintf(fileID,'%s\n%s\n',...
'**','** NOTE: Max number of nodes or elements per line is 10');

fprintf(fileID,'%s\n%s\n','**','*NSET, NSET=NSET_bottom1');
fprintf(fileID,formatSpecNset(numel(NSET_nm1)),NSET_nm1);
fprintf(fileID,'%s\n%s\n','**','*NSET, NSET=NSET_bottom2');
fprintf(fileID,formatSpecNset(numel(NSET_nm2)),NSET_nm2);

fprintf(fileID,'%s\n%s\n','**','*NSET, NSET=NSET_top1');
fprintf(fileID,formatSpecNset(numel(NSET_np1)),NSET_np1);
fprintf(fileID,'%s\n%s\n','**','*NSET, NSET=NSET_top2');
fprintf(fileID,formatSpecNset(numel(NSET_np2)),NSET_np2);

fprintf(fileID,'%s\n%s\n','**','*NSET, NSET=NSET_left1');
fprintf(fileID,formatSpecNset(numel(NSET_gm1)),NSET_gm1);
fprintf(fileID,'%s\n%s\n','**','*NSET, NSET=NSET_left2');
fprintf(fileID,formatSpecNset(numel(NSET_gm2)),NSET_gm2);

fprintf(fileID,'%s\n%s\n','**','*NSET, NSET=NSET_right1');
fprintf(fileID,formatSpecNset(numel(NSET_gp1)),NSET_gp1);
fprintf(fileID,'%s\n%s\n','**','*NSET, NSET=NSET_right2');
fprintf(fileID,formatSpecNset(numel(NSET_gp2)),NSET_gp2);

if problem_type == 1
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(fileID,'%s\n%s\n%s\n%d%s\n%s\n%s\n','**',...
'*STEP,NAME=STEP1,PERTURBATION','*FREQUENCY, eigensolver=lanczos, sim',...
NofCPs*2,', ,','*BOUNDARY','*END STEP');
fclose(fileID);
%
%
elseif problem_type == 2
%
%
fprintf(fileID,'%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s',...
'**','*STEP,NAME=STEP1','*STATIC','**',...
'*BOUNDARY, TYPE=DISPLACEMENT',...
'**NSET_bottom1,1,1,0',...
'**NSET_bottom1,3,3,0',...
'**NSET_top1,1,1,0',...
'**NSET_top1,3,3,0',...
'NSET_left1,1,3,0',...
'NSET_left2,1,3,0',...
'**NSET_right1,2,2,0',...
'**NSET_right2,2,2,0',...
'*DLOAD',...
'ELSET1, U7, -1',...
'*END STEP');
fclose(fileID);
end

%-End of .INP file
fprintf('%s file created successfully\n',inp_file);

end

%-------------function define format for Set of Nodes printing-------------
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---Creating of equation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%fprintf('%s%d\n', ...
%        'U direction. Equation Set Number',i);
% CP_u_1=NofCPs_u-1;
% CP_u_2=NofCPs_u*NofCPs_v-1;
% %
% %fprintf('%s%d\n', ...
% %        'V direction. Equation Set Number',i);
% CP_v_1=1;
% CP_v_2=NofCPs_u*(NofCPs_v-1)+1;
% %
% fprintf(fileID,'%s\n%s\n%s\n%s\n','**','** EQUATIONS','**','*EQUATION');
% %for j=CP_u_1:CP_u_2
% %    fprintf(fileID,'%s,\n%d, %s, %s, %d, %s, %s\n','2',j,'1','1.',...
% %                   j+NofCPs_u,'1','-1.');
% %end
% %
% for j=fix((CP_v_1-1)/NofCPs_u)+1:fix((CP_v_2-1)/NofCPs_u)+1
%     fprintf(fileID,'%s,\n%d, %s, %s, %d, %s, %s\n','2',...
%                    CP_v_1+(j-1)*NofCPs_u,'2','1.',...
%                    CP_v_1+(j-1)*NofCPs_u+1,'2','-1.');
% end
% %
% for j=fix((CP_u_1-1)/NofCPs_u)+1:fix((CP_u_2-1)/NofCPs_u)+1
%     fprintf(fileID,'%s,\n%d, %s, %s, %d, %s, %s\n','2',...
%                    CP_u_1+(j-1)*NofCPs_u,'1','1.',...
%                    CP_u_1+(j-1)*NofCPs_u+1,'1','-1.');
% end
% fprintf(fileID,'%s\n','**************************************');
%
% --- For Annulus only
%
%
%Num_low=input('Number of the lower CP = ');
%
%fprintf(fileID,'%s\n%s\n%s\n%s\n','**','** EQUATIONS','**','*EQUATION');
%for j=1:NofCPs_v
%    jj=1+NofCPs_u*(j-1);
%    fprintf(fileID,'%s,\n%d, %s, %s, %d, %s, %s\n','2',NofCPs_u*j,'1','1.',...
%                   jj,'1','-1.');
%    fprintf(fileID,'%s,\n%d, %s, %s, %d, %s, %s\n','2',NofCPs_u*j,'2','1.',...
%                   jj,'2','-1.');
%end
%
%fprintf(fileID,'%s\n%s\n%s\n','**','**','**');
%
%for j=1:NofCPs_v
%    jj=1+NofCPs_u*(j-1);
%    fprintf(fileID,'%s,\n%d, %s, %s, %d, %s, %s, %d, %s, %s\n','3',jj,'1','1.',...
%                   jj+1,'1','-0.5',NofCPs_u*j-1,'1','-0.5');
%    fprintf(fileID,'%s,\n%d, %s, %s, %d, %s, %s, %d, %s, %s\n','3',jj,'2','1.',...
%                   jj+1,'2','-0.5',NofCPs_u*j-1,'2','-0.5');
%end
%
%fprintf(fileID,'%s\n%s\n%s\n','**','**','**');
%
%for j=1:NofCPs_v
%    jj=Num_low+NofCPs_u*(j-1);
%    fprintf(fileID,'%s,\n%d, %s, %s, %d, %s, %s, %d, %s, %s\n','3',jj,'1','1.',...
%                   jj-1,'1','-0.5',jj+1,'1','-0.5');
%    fprintf(fileID,'%s,\n%d, %s, %s, %d, %s, %s, %d, %s, %s\n','3',jj,'2','1.',...
%                   jj-1,'2','-0.5',jj+1,'2','-0.5');
%end
%
%fprintf(fileID,'%s\n','**');
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%