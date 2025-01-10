function[U_, gam_,w_] = FlutterPKMethod(Umax,Nm,p,K,M,A0,A1c,A1nc)
% Numerical variables
tol = 1e-6; 
max_iter = 100;

% Parameters
U_ = linspace(0.001,Umax,100); % Velocities vector
c = p.c

% Defining the Theodorsen's function:
C = @(k) 1 - 0.165/(1-1i*0.045/k) - 0.335/(1-1i*0.3/k);

% Preallocation
Ndof = size(K,1); 
gam_ = zeros(Ndof, length(U_));
w_ = zeros(Ndof, length(U_));
Vpk_= zeros(Ndof, length(U_));

% First guess for frequencies
[~,Ds] = eig(K,M);
w_(:,1) = sqrt(diag(Ds));

% Loop through velocites
for i = 1:length(U_)

    if i>1
        w_(:,i) = w_(:,i-1);
        gam_(:,i) = gam_(:,i-1);
    end
    
    for j = 1:Nm
        conv = 1;
        iter = 0;
        
        while conv > tol && iter < max_iter
            iter = iter+1;
            
            % Get reduced frequency
            k = w_(j,i)*c/(2*U_(i));
            
            % Compute effective matrices
            Keff = K - U_(i)^2*(C(k)*A0 - 1i*k*(C(k)*A1c-A1nc));
            
            % Extend system matrices
            A = [Keff, zeros(Ndof,Ndof); zeros(Ndof,Ndof), eye(Ndof,Ndof)];
            B = [zeros(Ndof), -M; eye(Ndof,Ndof), zeros(Ndof,Ndof)];

            % Solve eigenvalues
            [Vpk, Dpk] = eig(A,B,Nm,'sm');
            p = diag(Dpk);
            
            % Look for closest node to initial guess (w_ and gam_)
            [conv,jmin] = min(abs(real(p)-gam_(j,i)) + abs(imag(p)-w_(j,i)));
            
            % Store the newer values as initial guesses or previous values
            w_(j,i) = imag(p(jmin));
            gam_(j,i) = real(p(jmin));
            Vpk_(:,j,i) = Vpk(1:Ndof,jmin);
            
            
        
        end
        
    end
    
end

% Finding flutter speed for pk
Ufpk = 0;
i = 1;
while Ufpk == 0
    for j = 1:length(gam_(:,i))
        if gam_(j,i)>0
            Ufpk = U_(i);
        end
    end
    i = i + 1;
end
disp(["The flutter speed for the pk method is: ", num2str(Ufpk), "m/s"])


end

