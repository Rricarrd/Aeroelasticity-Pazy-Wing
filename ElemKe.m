function [Ke] = ElemKe(Li,p)
% Parameters
EI = p.EI;
GJ = p.GJ;

% Element stiffness matrix [Ke]
K1 = (EI/Li^3) * [12 6*Li 0 -12 6*Li 0;
                 6*Li 4*Li^2 0 -6*Li 2*Li^2 0;
                 0 0 0 0 0 0;
                 -12 -6*Li 0 12 -6*Li 0;
                 6*Li 2*Li^2 0 -6*Li 4*Li^2 0;
                 0 0 0 0 0 0];

K2 = (GJ/Li) * [0 0 0 0 0 0;
               0 0 0 0 0 0;
               0 0 1 0 0 -1;
               0 0 0 0 0 0;
               0 0 0 0 0 0;
               0 0 -1 0 0 1];

Ke = K1 + K2;
end

