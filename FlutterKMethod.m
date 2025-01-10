function[Uk_,gk_,wk_] = FlutterKMethod(Umax,Nm,p,K,M,A0,A1c,A1nc)
% k-method

% Parameters
U_ = linspace(0.001,Umax,100); % Velocities vector
c = p.c;


% Defining the Theodorsen's function:
C = @(k) 1 - 0.165/(1-1i*0.045/k) - 0.335/(1-1i*0.3/k);

% Inverse k vector (k^-1 = (2*Uinf)/(omega*c))
invk_ = linspace(0.001,5,100);

% Preallocation
Ndof = size(K,1); 
l_ = zeros(Nm,length(invk_));
Vk_ = zeros(Ndof,Nm,length(invk_));

% Sweep through inverse of k
for i = 1:length(invk_)
    
    % Get reduced frequency
    k = 1/invk_(i);
    
    % Compute effective mass expression
    Meff = M + c^2/(4*k^2)*(C(k)*A0 - 1i*k*(C(k)*A1c - A1nc));
    
    % Solve the eigenvalues problem
    [Vk,Dk] = eigs(Meff,K,Nm,'sm'); % We use the opposite (usually it is eigs(K,M) but now it the opposite: see equations of page T2.2-8 --> we have the eigenvalue (lambda) multiplying the K matrix instead of the Meff matrix)
   
    % Sort (the following algorithm analyses the previous iteration and saves the colsest one to the actual iteration, and we will elimiate this value)
    if i > 1
        % it could be added to the p-method to elimiate the jupms
        l = diag(Dk);
        tosort = 1:Nm;
        for j = 1:Nm
            [~,jmin] = min(abs(real(l(tosort))-real(l_(j,i-1))) + abs(imag(l(tosort))-imag(l_(j,i-1)))); % to find from all the previous eigenvalues which is the closest one (difference between the real part and the imaginary part)
            l_(j,i) = l(tosort(jmin));
            Vk_(:,j,i) = Vk(:,tosort(jmin));
            tosort(jmin) = [];
        end
    else
        l_(:,1) = diag(Dk);
    end
end
wk_ = sqrt(1./real(l_));
gk_ = imag(l_)./real(l_);
Uk_ = wk_*c.*invk_/2;


% Finding flutter speed for k
Ufk = 0;
i = 1;
while Ufk == 0
    for j = 1:length(gk_(:,i))
        if real(gk_(j,i))>0
            Ufk = Uk_(i);
        end
    end
    i = i + 1;
end
disp(["The flutter speed for the k method is: ", num2str(Ufk), "m/s"])


end

