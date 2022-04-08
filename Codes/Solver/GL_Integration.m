function [q,f,x,t] = GL_Integration(U,gamma,mu,dx,dt,nt,xmin,xmax)

% build x and A operator
[x,A] = GL_lin_init(U,gamma,mu,xmin,xmax,dx); % initialize the linear operator
%Boundary Conditions
A = A(2:end-1,2:end-1); 
x = x(2:end-1);
nx = length(x);

% Original values
q0 = zeros(nx,1); %Initial Condition
f = randn(nx,nt)+randn(nx,nt)*1i;  %random white-noise force

%Smooth forcing signal using a low pass filter 
m = [1 1 0 0];
fs = [0 0.6 0.6 1];  
myFIR = fir2(30,fs,m);
f = filter(myFIR,1,f,[],2);

[t,q]= GL_CN(dt,A,q0,f); 

% Time to cross the domain
T = (xmax-xmin)/U;

% Window necessary to achieve T
wT = T/dt;

save(strcat("GL_solution.mat"),'q','f','x','t','T','wT','-v7.3');
end