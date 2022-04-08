function [Dup,Ddown,D2] = FDcoeff_GL(x)

% creates FD derivation matrices
% x: discrete points of coordinate
% nota bene: this implementation implies that psi is strictly zero 
%            everywhere outside the numerical domain.

n=length(x);
dx=x(2)-x(1);
if sum(diff(diff(x)))>1000*eps
    disp(sum(diff(diff(x))));
    error('Derivation matrices require uniform dx.');
end

Dup  =zeros(n);
Ddown=zeros(n);
D2   =zeros(n);

Dup   = 1.5*eye(n) - 2*diag(ones(n-1,1),-1) + 0.5*diag(ones(n-2,1),-2);
Dup   = Dup/dx;
Ddown = -1.5*eye(n) + 2*diag(ones(n-1,1),1) - 0.5*diag(ones(n-2,1),2);
Ddown = Ddown/dx;
D2    = -2*eye(n) + diag(ones(n-1,1),-1) + diag(ones(n-1,1),1);
D2    = D2/(dx^2);

end