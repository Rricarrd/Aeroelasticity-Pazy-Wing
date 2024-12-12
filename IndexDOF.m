function I = IndexDOF(p,elem,dof)
I = p.nj*(elem-1)+ dof;
end

