function plot_SPOD_Modes(SPOD_path,index_mode,index_omega)
    data = matfile(SPOD_path);
    x = data.x;
    mode = data.Psi(index_omega,:,index_mode)';
    omega = data.St;
    omega = 2*pi*omega(index_omega);
    
    %adjust phase
    mode = mode/exp(1i*angle(mode(2)));
    
    % Padding
    x = [0;x;30];
    mode = [0;mode;0];
    
    figure;
    plot(x,real(mode),x,imag(mode),x,abs(mode)); 
    title(['SPOD mode ', int2str(index_mode), ', \omega = ', num2str(omega), ' rad/s']);
    xlabel('x');legend('Real part','Imag. part','Abs. value');
end