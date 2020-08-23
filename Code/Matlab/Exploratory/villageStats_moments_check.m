%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check D matrix values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd /home/sakina/Github/MFdiffusion_replication

%% qN = qP (model 1)
load('data_model_1_mom_1 trimesters .mat')

% Calculate the village avarage for each moment for each simulation
villmnD1 = cell(length(qN),1);

for i=1:length(qN)
    villmnD1{i} = mean(D{i});
end

% Convert to a matrix
villmnDMat1 = permute(reshape([villmnD1{:}], [], size(villmnD1, 1), size(villmnD1, 2)), [2 1]);

% Calculate mean
mean(villmnDMat1, 1)

% For sanity check 
mean(villmnD1{1}, 1)

%% qN != qP (model 2)
load('data_model_3_mom_1 trimesters .mat')

% Calculate the village avarage for each moment for each simulation
villmnD2 = cell(length(qN),length(qP));

for i=1:length(qN)
        for j=1:length(qP)
           villmnD2{i,j} = mean(D{i,j});
        end
end

% Turn villmnD as a 3 dimensional matrix (31X39X5)
villmnDMat2 = permute(reshape([villmnD2{:}], [], size(villmnD2, 1), size(villmnD2, 2)), [2 3 1]);

% Calculate the average for each moment
mean(squeeze(mean(villmnDMat2)), 1)
