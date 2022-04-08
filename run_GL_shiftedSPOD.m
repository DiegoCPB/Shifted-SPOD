clear all

%% Input Parameters
% I/O
input_folder = "~/SPOD_convergence/Dados/U12_A100";
input_file = "GL_whiteForcing_FIR.mat";
output_folder = "~/SPOD_convergence/SPOD_test/";
U = 12;
xp = 15;

%SPOD
dw = 100; % Window step
min_w = 1000; % Minimal window, must be divisible by 4.
max_w = 1000; % Max window must be multiple of min window.
overlap_vals = [75]; %[50 75]; % DO NOT CHANGE

%% Other parameters
addpath(genpath('Codes'));
% Dumping folders
SPOD_folder = strcat(output_folder);

%% Calculation
% Check
if mod(min_w,4) ~= 0
    throw(MException('MATLAB:NotDivisible', 'min_w must be divisible by 4'));
end

if mod(max_w,min_w) ~= 0
    throw(MException('MATLAB:NotDivisible', 'max_w must be multiple of min_w'));
end

data = load(input_folder+'/'+input_file);
t = data.t;
x = data.x;
Q = data.q;

dt = t(2)-t(1);
n_snapshots = length(t);
W = weightVector(x); % SPOD weights

% Saving SPOD modes to file
if ~exist(SPOD_folder, 'dir')
    mkdir(SPOD_folder)
end

disp(' ');
for overlap = overlap_vals
    for window = min_w:dw:max_w
        disp(['Calculating SPOD with window = ', num2str(window), ' and overlap = ', num2str(overlap), '%'])
        % Usual SPOD
%         [Psi, Lambda, Theta, Qhat, St, n_blocks] = spod_snapshot(Q,W,dt,window,overlap);

        % Shifted Snapshot SPOD v5
%         [Psi, Lambda, Theta, Qhat, St, n_blocks] = shifted_spod_snapshot_v5(Q,W,x,xp,U,dt,window,overlap);
        
        % Shifted SPOD
        [Psi, Lambda, Theta, Qhat, St, n_blocks] = shifted_spod(Q,W,x,xp,U,dt,window,overlap);
       
        save_filename = sprintf('%s/SPOD_nb%d_w%d_o%d.mat',SPOD_folder,n_blocks,window,overlap);
        save(save_filename,'Psi','Lambda','Theta','Qhat','St','x','t','n_blocks','window','overlap','-v7.3');
        
    end
end
