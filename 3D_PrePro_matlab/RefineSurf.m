% (c) V. Balobanov, 2017, Aalto, Dpt. of Civil Engineering
% uses J. Kiendl's and D. Proserpio's functions
% Refine a surface nurbs structure
%
% IMPUTS:
% nurbs - initial geometry
% deg_u,deg_v - new degree (if it lower than the current, save the old one)
% refu,refv - each of the initial elements will be divided by this number of els
% plotting - boolean: plot or not?

function nurbs = RefineSurf(nurbs,deg_u,deg_v,refu,refv,plotting)

addpath('nurbs_functions')
% Extract the fields from the data structure of the first layer into local variables
data_names = fieldnames (nurbs);
for iopt  = 1:numel (data_names)
	eval ([data_names{iopt} '= nurbs.(data_names{iopt});']);
end

% reshape array CPs to Josef's format
% 1dim: num of CP in u-dir, 2dim: num of CP in v-dir, 3dim: num of coord
CP(:,1,:) = CPs_Coords;
CP = permute(CP,[3 2 1]);
CP = reshape (CP,[NoCPsU NoCPsV 4]);

nurbsaux.U  = xi;
nurbsaux.V  = eta;
nurbsaux.p  = p;
nurbsaux.q  = q;
nurbsaux.CP = CP;

if p < deg_u && q >= deg_v
    nurbsaux = degree_elevate_surf(nurbsaux,deg_u-p,0);
elseif p >= deg_u && q < deg_v
    nurbsaux = degree_elevate_surf(nurbsaux,0,deg_v-q);
elseif p < deg_u && q < deg_v
    nurbsaux = degree_elevate_surf(nurbsaux,deg_u-p,deg_v-q);
end
Ru = refinement_vec(nurbsaux.U,refu);
Rv = refinement_vec(nurbsaux.V,refv);
nurbsaux = knot_refine_surf(nurbsaux,Ru,Rv);

% plot refined geometry - optional
if plotting
    addpath('plot_functions');
    plotNURBS_surf(nurbsaux.p, nurbsaux.q, nurbsaux.U, nurbsaux.V, nurbsaux.CP)
end

nurbs.xi  = nurbsaux.U;
nurbs.eta = nurbsaux.V;
nurbs.p   = nurbsaux.p;
nurbs.q   = nurbsaux.q;
nurbs.NoCPsU = size(nurbsaux.CP,1);
nurbs.NoCPsV = size(nurbsaux.CP,2);
nurbs.NofEl_U = length (unique(nurbs.xi)) - 1;
nurbs.NofEl_V = length (unique(nurbs.eta)) - 1;

%return to the initial format
nurbsaux.CP = reshape (nurbsaux.CP,[], 4); 
nurbsaux.CP = nurbsaux.CP';

nurbs.INP_CPs_Coords = []; 
nurbs.Vector_of_Weights = [];
nurbs.INP_CPs_Coords(1,:) = 1:size(nurbsaux.CP,2);
nurbs.INP_CPs_Coords(2:4,:) = nurbsaux.CP(1:3,:);
nurbs.Vector_of_Weights = nurbsaux.CP(4,:); 
