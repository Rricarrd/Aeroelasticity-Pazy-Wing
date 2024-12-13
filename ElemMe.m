function [Me] = ElemMe(Li,p)
% Parameters
mu = p.mu;
rS = p.rS;
r0 = p.r0;
d = p.d;


%%%%%%%%%%%%%%%%%% REPASSAR MATRIUS (CHAT GPT ES CA EQUIVOCAR)
% Element mass matrix [Me]
M1 = (mu*Li/420) * [156 22*Li 0 54 -13*Li 0;
                   22*Li 4*Li^2 0 13*Li -3*Li^2 0;
                   0 0 0 0 0 0;
                   54 13*Li 0 156 -22*Li 0;
                   -13*Li -3*Li^2 0 -22*Li 4*Li^2 0;
                   0 0 0 0 0 0];

M2 = (mu*r0^2/(30*Li)) * [36 3*Li 0 -36 -3*Li 0;
                          3*Li Li^2 0 -3*Li -Li^2 0;
                          0 0 0 0 0 0;
                          -36 -3*Li 0 36 3*Li 0;
                          -3*Li -Li^2 0 3*Li Li^2 0;
                          0 0 0 0 0 0];

M3 = (mu*rS^2*Li/6) * [0 0 0 0 0 0;
                       0 2 0 0 0 1;
                       0 0 0 0 0 0;
                       0 0 0 0 0 0;
                       0 0 0 0 2 0;
                       0 1 0 0 0 2];

M4 = (mu*d*Li/60) * [0 21 0 0 9 0;
                    21 3*Li 0 0 0 2*Li;
                    0 0 0 0 0 0;
                    0 0 0 0 21 3*Li;
                    9 0 0 21 -2*Li -3*Li;
                    0 2*Li 0 -3*Li 0 0];

Me = M1 + M2 + M3 + M4;


end

