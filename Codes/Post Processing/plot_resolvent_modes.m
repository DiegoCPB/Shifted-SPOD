function plot_resolvent_modes(LinearMatrix_path,index_mode,omega)

[U,~,~,x] = GL_solveResolvent(LinearMatrix_path,omega);

U = U(:,index_mode);

% Padding
x = [0;x;30];
U = [0;U;0];

figure;
plot(x,real(U),x,imag(U),x,abs(U)); 
title(['Response mode ' int2str(index_mode), ', \omega = ', num2str(omega), ' rad/s']);
xlabel('x');
ylabel('q');
legend('Real part','Imag. part','Abs. value');

end