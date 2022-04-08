function [filepaths, betas, t_vec, grid, baseFlow2D] = dataMatrices(x_lims,t0,sim_path,save_folder)
% This script is used to create the SPOD 3D data matrix

% INPUTS:
% x_lims : [x0, x1], first and last positions in x direction
% sim_path : Path to the folder containing the simulation files

% OUTPUTS:
% filenames : Names of the saved files containing the data matrices for each beta 
% betas : Values for beta

% Initial parameters
[string_files,t_vec] = getFileNames(t0,sim_path);

s = strcat(sim_path,'/bflow.u');
fprintf('\nReading base flow file: %s\n',s);
[bf.vel,~,~,~]  = readFile(s);

filepaths = {};
for ii = 1:length(t_vec)
    input_filepath = strcat(sim_path,'/',string_files{ii});
    % Load each flowfield
    if ii == 1
        [vel,grid.x,grid.y,grid.z] = readFile(input_filepath,string_files{ii});
        
        % Remove grid values outside x_lims
        [~,ind_x0] = min(abs(grid.x-x_lims(1)));
        [~,ind_x1] = min(abs(grid.x-x_lims(2)));
        grid.x = grid.x(ind_x0:ind_x1);
        
        % Cutting base flow domain
        baseFlow2D.u = bf.vel.u(ind_x0:ind_x1,:);
        baseFlow2D.v = bf.vel.v(ind_x0:ind_x1,:);
        baseFlow2D.w = bf.vel.w(ind_x0:ind_x1,:);
        
        % Size of domain 
        n_z = length(grid.z);
        f_z = 1/(grid.z(2)-grid.z(1)); % sampling frequency for the z vector
        
        % Generate betas and filenames to save
        [filepaths,betas] = generateFilepathsAndBetas(n_z,f_z,save_folder,'DataMatrix');
        n_betas = length(betas);
        
        % Verify if the data is already available, as it's extraction is quite time consuming 
        if all(isfile(filepaths))
            fprintf('\nAll data files already exist inside %s\n\n', save_folder)
            break
        end
    else
        [vel,~,~,~] = readFile(input_filepath,string_files{ii});
    end

    % Remove velocity values outside x_lims
    vel.u = vel.u(ind_x0:ind_x1,:,:);
    vel.v = vel.v(ind_x0:ind_x1,:,:);
    vel.w = vel.w(ind_x0:ind_x1,:,:);
    
    % Fourier transform along the z direction, which has periodic boundary
    % conditions
    vel.u = fft(vel.u,[],3);
    vel.v = fft(vel.v,[],3);
    vel.w = fft(vel.w,[],3);
    
    for jj = 1:n_betas
        u_ = vel.u(:,:,jj);
        v_ = vel.v(:,:,jj);
        w_ = vel.w(:,:,jj);
        var = [u_(:); v_(:); w_(:)];
        if ii == 1
            Q = var;
            if ~exist(save_folder, 'dir')
               mkdir(save_folder)
            end
            % The -v7.3 flag enables partial loading of the .mat file
            save(filepaths{jj},'Q','-v7.3');
        else
            mat = matfile(filepaths{jj},'Writable',true);
            s = size(mat, 'Q');
            mat.Q(:,s(2)+1) = var;
        end
    end
end
end

function [vel,x,y,z] = readFile(file_path,filename)

if isfile(file_path)
    [vel,x,y,z,~,~,~,~] = read_flowfield(file_path);
else
    throw(MException('File:read','No file %s in directory.',filename));
end

end