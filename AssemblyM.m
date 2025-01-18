function [M] = AssemblyM(y,Tn,Tr,p)

% Parameters
Ndof = p.Ndof;
Nelem = p.Nelem;
nj = p.nj;
en = p.en;

% Preallocation
M = zeros(Ndof,Ndof);

for e = 1:Nelem % For each element
    
    Li = y(Tn(e,2))-y(Tn(e,1)); % Element size
    Me = ElemMe(Li,p); % Mass element matrix 6x6 (6 = 2 nodes x 3 DOFs)
    
    I = zeros(1,en*nj); % Array of positions corresponding to the global matrix
    for i = 1:en     % For each node of the element
        for j = 1:nj     % For each DOF of the node
            I(nj*(i-1)+j) = IndexDOF(p,Tn(e,i),j); % Filling the I positions
        end
    end
    
    M(I,I) = M(I,I) + Me; % Assembling the element matrix into the desired global positions
    
end

% Introduction of the ribs mass
Mrib = computeRibMass(p);

for e = 1:size(Tr,1)

    I_ribs = zeros(1,nj);
    k = 0;
    for j = 1:nj
        k = k + 1;
        I_ribs(1,k) = IndexDOF(p,Tr(e),j);
    end

    M(I_ribs,I_ribs) = M(I_ribs,I_ribs) + Mrib;
end

end

