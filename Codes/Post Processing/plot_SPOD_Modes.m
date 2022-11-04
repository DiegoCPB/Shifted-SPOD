function plot_SPOD_Modes(SPOD_path,index_mode,omega)
    data = matfile(SPOD_path);
    x = data.x;
    St = data.St;
    [~,index_omega] = min(abs(2*pi*St-omega));

    mode = data.Psi(index_omega,:,index_mode).';
    omega = St;
    omega = 2*pi*omega(index_omega);
    
    % Padding
    x = [0;x;30];
    mode = [0;mode;0];
    
    figure;
    plot(x,real(mode),x,imag(mode),x,abs(mode)); 
    title(['SPOD mode ', int2str(index_mode), ', \omega = ', num2str(omega), ' rad/s']);
    xlabel('x');legend('Real part','Imag. part','Abs. value');
end