% Copyright (C) 2017, 2025 Viacheslav Balobanov
%
% IMPUTS:
% nurbs - initial geometry structure in NURBS toolbox format
% config - configuration structure containing information about refinement

function nurbs = refineSurface(nurbs,config)

% Degree elevation and checks
p_old = nurbs.order(1) - 1;
q_old = nurbs.order(2) - 1;
p_new = config.deg_u;
q_new = config.deg_v;

if p_old < p_new && q_old >= q_new
    nurbs = nrbdegelev(nurbs, [p_new-p_old,0]);
    disp('No degree elevation for v direction performed');
elseif p_old >= p_new && q_old < q_new
    nurbs = nrbdegelev(nurbs, [0,q_new-q_old]);
    disp('No degree elevation for u direction performed');
elseif p_old < p_new && q_old < q_new
    nurbs = nrbdegelev(nurbs, [p_new-p_old,q_new-q_old]);
else
    disp('No degree elevation for u and v directions performed')
end

% Knot refinement

new_knots_u = knots_refine(nurbs.knots{1}, config.ref_u);
new_knots_v = knots_refine(nurbs.knots{2}, config.ref_v);
nurbs = nrbkntins(nurbs, {new_knots_u new_knots_v});


function new_knots = knots_refine(kntvec,refN)
kntvec_uniq = unique(kntvec);
new_knots = zeros(1,refN*(length(kntvec_uniq)-1));
for i=1:length(kntvec_uniq)-1
    new_knots_i = linspace(kntvec_uniq(i),kntvec_uniq(i+1),refN+2);
    new_knots_i = new_knots_i(2:end-1);
    new_knots(((i-1)*refN+1):(i*refN)) = new_knots_i;
end

