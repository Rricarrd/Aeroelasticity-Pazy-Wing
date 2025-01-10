function [U_,p_,Vp_] = FlutterPMethod(Umax,Nm,p,K,M,A0,A1c,A1nc)

c = p.c;
U_ = linspace(0.001,Umax,100); % Velocities vector

for i = 1:length(U_)
    
    % Computing the effective matrices
    Ceff = 0.5*c*U_(i)*(A1c-A1nc);
    Keff = K - U_(i)^2*A0;
    
    % Extending the system
    Ndof = size(K,1);
    A = [Keff, zeros(Ndof,Ndof); zeros(Ndof,Ndof), eye(Ndof,Ndof)];
    B = [-Ceff, -M; eye(Ndof,Ndof), zeros(Ndof,Ndof)];
    
    % Solving the problem
    [Vp,Dp] = eigs(A,B,Nm,'sm'); 
    
    % Sorting the values so that to avoid the jumps jumps
    p = diag(Dp);
    [~,indsorted] = sort(real(p),'descend');
    p = p(indsorted);
    Vp = Vp(:,indsorted);
    Vp_(:,:,i) = Vp;
    p_(:,i) = p;
end

% Finding flutter speed for p
Ufp = 0;
i = 1;
while Ufp == 0
    for j = 1:length(p_(:,i))
        if real(p_(j,i))>0
            Ufp = U_(i);
        end
    end
    i = i + 1;
end
disp(["The flutter speed for the p method is: ", num2str(Ufp), "m/s"])



end

