function [A] = AssemblyA(Ae,y,Tn,p)

% Parameters
Ndof = p.Ndof;
Nelem = p.Nelem;
nj = p.nj;
en = p.en;

% Preallocation
A = zeros(Ndof,Ndof);

for e = 1:Nelem % For each element
    
    Li = y(Tn(e,2))-y(Tn(e,1)); % Element size
    Ai = Ae*Li; % Aerodynamic element matrix 6x6 (6 = 2 nodes x 3 DOFs) (Still missing Uing and Ck (if needed))
    
    I = zeros(1,en*nj); % Array of positions corresponding to the global matrix

    for i = 1:en     % For each node of the element
        for j = 1:nj     % For each DOF of the node
            I(nj*(i-1)+j) = IndexDOF(p,Tn(e,i),j); % Filling the I positions
        end
    end
    
    A(I,I) = A(I,I) + Ai; % Assembling the element matrix into the desired global positions
    
end

if p.sparse
    A = sparse(A);
end

end
