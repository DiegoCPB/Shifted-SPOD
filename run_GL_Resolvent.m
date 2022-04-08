clear all
addpath(genpath('Codes'));

% I/O
LinearMatrix_file = 'GL_linearMatrix.mat';
mode2plot = 1;
omega2plot = 2*pi;


plot_resolvent_modes(LinearMatrix_file,mode2plot,omega2plot)