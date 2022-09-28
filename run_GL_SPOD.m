clear all
addpath(genpath('Codes'));

%% Input Parameters
% I/O
input_file = "GL_solution.mat";

%SPOD
Nfft = 100;
Olap = 75;

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

outfile = sprintf('SPOD_w%d_o%d.mat',Nfft,Olap);
if ~isfile(outfile)
    disp(['Calculating SPOD with window = ', num2str(Nfft), ' and overlap = ', num2str(Olap), '%'])
    [Psi, Lambda, Qhat, St, Nb] = spod(Q,W,dt,Nfft,Olap);
    save(outfile,'Psi','Lambda','Qhat','St','x','t','Nb','Nfft','Olap','-v7.3');
end

%% PLot mode
plot_SPOD_Modes(outfile,mode2plot,omega2plot)
