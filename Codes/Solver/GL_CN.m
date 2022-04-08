function [t,q] = GL_CN(dt,A,q0,f)

nt = size(f,2);
nx = length(q0);

q = zeros(nx,nt);
q(:,1)=q0;

M    = eye(size(A))-0.5*dt*A;
% Minv = inv(M);
Am = eye(size(A))+0.5*dt*A;

%Use sparse matrices
Am = sparse(Am);
M  = sparse(M);

for it = 2:nt
    b= Am*q(:,it-1) + 0.5*dt*(f(:,it)+f(:,it-1));
    q(:,it)=M\b;
end
q=q(:,1:it);
t = (0:size(q,2)-1)*dt;

