clear all
addpath(genpath('Codes'));

%% Input Parameters
% I/O
input_file = "GL_solution.mat";

%SPOD
Nfft = 100;
Olap = 75;
U = 10;
xp = 15;

%Plot
mode2plot = 1;
omega2plot = 2*pi;

%% Calculation
data = load(input_file);
t = data.t;
x = data.x;
Q = data.q;
W = weightVector(x); % SPOD weights
dt = t(2)-t(1);

disp(['Calculating Shifted SPOD with window = ', num2str(Nfft), ' and overlap = ', num2str(Olap), '%'])
[Psi, Lambda, Qhat, St, Nb] = shifted_spod(Q,W,x,xp,U,dt,Nfft,Olap);
outfile = sprintf('ShiftedSPOD_nb%d_w%d_o%d.mat',Nb,Nfft,Olap);
save(outfile,'Psi','Lambda','Qhat','St','x','t','Nb','Nfft','Olap','-v7.3');

%% PLot mode
[~,idx_omega] = min(abs(2*pi*St-omega2plot));
plot_SPOD_Modes(outfile,mode2plot,idx_omega)
