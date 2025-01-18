function [Ud] = Divergence(p,Knum,Anum,k)




% Eigenproblem
[~,d] = eigs(Knum,Anum,k,'sm');

% Velocities at the diagonal of the eigenvalues
Ud = sqrt(diag(d)/(0.5*p.rho));

% Classifying the velocities
for i = 1:size(d,1)
    if Ud(i) < 0                % No negative values
        Ud(i) = 0;
    elseif isreal(Ud(i)) == 0   % No complex values
        Ud(i) = 0;
    elseif Ud(i) > 10^10        % No infinite values
        Ud(i) = 0;
    end
end

% Getting the minimum positive divergence speed different from 0
if min(Ud(Ud>0)) ~= 0
    Ud = min(Ud(Ud>0));
else
    Ud = 0;
end

disp(["The divergence speed is: ",num2str(Ud),"m/s"])
end

