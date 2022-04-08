function [x,A] = GL_lin_init(U,gamma,mu,xmin,xmax,dx)

% build x
x=(xmin:dx:xmax).';

% build derivative matrices
[Dup, ~, D2] = FDcoeff_GL(x); 
% Dup: upwinding for convection
% D2: second derivative for diffusion

% Direct GL equation
mu = mu(x);
A = diag(mu) - U*Dup + gamma*D2;
end
