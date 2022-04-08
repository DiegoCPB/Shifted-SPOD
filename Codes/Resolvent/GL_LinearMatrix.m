function GL_LinearMatrix(U,gamma,mu,xmin,xmax,dx)
    
[x,L] = GL_lin_init(U,gamma,mu,xmin,xmax,dx);
%Boundary Conditions
L = L(2:end-1,2:end-1); 
x = x(2:end-1);

save(strcat("GL_linearMatrix.mat"),'L','x','-v7.3');
end

