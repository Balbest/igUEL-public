% function creates "connectivity" array INP_elements
% containing control points for each element
% (c)2017,2025 S.Khakalo, V.Balobanov

function [INP_Elements, KnMultU, KnMultV] = InputElements(nurbs)

p = nurbs.order(1)-1;  xi = nurbs.knots{1}; 
q = nurbs.order(2)-1; eta = nurbs.knots{2};
NoCPsU = nurbs.number(1); 

% Number of elements
NofElU = length (unique(xi)) - 1;
NofElV = length (unique(eta)) - 1;
NofEl = NofElU*NofElV;

KnMultU = zeros(1,NofElU+1); % vector of knots repetitions
if NofElU >1
    Kn_Ind_Num = 2+p;
   for i = 2:NofElU+1
    Kn_Ind_Num = Kn_Ind_Num+KnMultU(i-1);
       for j = Kn_Ind_Num:Kn_Ind_Num+p-2
         if xi(j) == xi(j+1) && xi(Kn_Ind_Num)-xi(j+1)==0
         KnMultU(i)=KnMultU(i)+1;
         end
       end
       KnMultU(i)=KnMultU(i)+1;
   end
end
KnMultU(1)=1;
KnMultU(NofElU+1)=1;

KnMultV=zeros(1,NofElV+1); % vector of knots repetitions
if NofElV >1
    Kn_Ind_Num = 2+q;
   for i = 2:NofElV+1
    Kn_Ind_Num = Kn_Ind_Num+KnMultV(i-1);
       for j = Kn_Ind_Num:Kn_Ind_Num+q-2
         if eta(j) == eta(j+1) && eta(Kn_Ind_Num)-eta(j+1)==0
         KnMultV(i)=KnMultV(i)+1;
         end
       end
       KnMultV(i)=KnMultV(i)+1;
   end
end
KnMultV(1)=1;
KnMultV(NofElV+1)=1;

% Element construction (valid for repeated inner knots)
% INP_Elements = zeros((p+1)*(q+1)+1, NofEl);
for j=1:NofEl
    INP_Elements(1,j)=j+0;
end

for i=1:q+1
    for j=1:p+1
        f=j+(i-1)*(p+1);
        i_v=0;
        for m=1:NofElV
            i_u=0;
            i_v=i_v+KnMultV(m);
            for k=1:NofElU
                r=k+(m-1)*NofElU;
                i_u=i_u+KnMultU(k);

                INP_Elements(f+1,r)=0+i_u+(i_v-1)*NoCPsU+...
                                  (j-1)+(i-1)*NoCPsU;
            end
        end
    end
end

% Multiplications of the knots for output
KnMultU(1) = p + 1;
KnMultU(NofElU+1) = p + 1;

KnMultV(1) = q + 1;
KnMultV(NofElV+1) = q + 1;

end