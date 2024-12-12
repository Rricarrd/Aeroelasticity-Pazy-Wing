function [y,Tn,Tr,Nnode,Nelem,Ndof] = DiscretizeWing(p)

% Parameters
Rs = p.Rs
Lw = p.Lw
Rn = p.Rn
nj = p.nj
ng = p.ng

% Rib positions
yj = 0:Rs:Lw
yj(15) = 550


% Nodal coordinates
y = 0
for i = 1:length(yj)-1
    yi = linspace(yj(i),yj(i+1),Rn);
    y = [y,yi(2:end)];
end

% Element connectivities
Tn = [1:length(y); 2:length(y)+1]'; % node conectivities
Tr = (Rn:Rn:length(y))';
Tr(end) = length(y)+1; % rib conectivities

% Node number, element and DOFs
Nnode = length(y);
Nelem = Nnode-1;
Ndof = Nnode*nj*2;

end

