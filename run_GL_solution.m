clear all
addpath(genpath('Codes'));

% Equation parameters, wave-packets (Cavalieri et al.) article
U = 10;
gamma = (1-1i)/10;
A = 0.6; %0.6, 1, 1.25
mu = @(x) A*(1-x/10);
xmin = 0;
xmax = 30;

% General parameters
dx   = 0.1;  
dt   = 0.01; 
Tmax = 500; 

% Rounding values
nt   = round(Tmax/dt); 
xmax = dx*round(xmax/dx); % fixes rounding errors
xmin = dx*round(xmin/dx);

% Solution
[q,f,x,t] = GL_Integration(U,gamma,mu,dx,dt,nt,xmin,xmax);

% Linear matrix for resolvent
GL_LinearMatrix(U,gamma,mu,xmin,xmax,dx);

% Plot fields
plot_solution(q,f,x,t)