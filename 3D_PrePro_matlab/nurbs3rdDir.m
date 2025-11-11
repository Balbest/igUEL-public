%creates 3rd direction nurbs-structure for IGA_3d_ABAQUS
%nurbs package is needed
%inputs: degree and numEl is final properties of output structure
function [knots, CPs] = nurbs3rdDir(r, nsub)

knots = [0 0 1 1];
CPs = [0 1; 0 0; 0 0];

[CPs, knots] = bspdegelev(1, CPs, knots, r-1);
if nsub > 1 
    [~, ~, Nknots] = kntrefine(knots, nsub-1, r, r-1);
    [CPs, knots] = bspkntins(r, CPs, knots, Nknots);
end

CPs = CPs(1,:);
end

% z_CPs is a given normalized vector containing the distribution the points along the volume direction
%(depends on number of elements along 3rd drection, RPol and regularity)
%CPs = [0 0.1 0.3 0.5 0.7 0.9 1.0];
%5 els: [0 0.04 0.12 0.24 0.4 0.6 0.76 0.88 0.96 1];
%for the cube: [0 0.3125 0.9375 1.875 3.125 4.375 5.625 6.875 8.125 9.0625 9.6875 10]/10; 


%zeta = [0 0 0 0 0 0 0.5 1 1 1 1 1 1]; %knot vector [0 0 0 0 0 0 0.2 0.4 0.6 0.8 1 1 1 1 1 1];