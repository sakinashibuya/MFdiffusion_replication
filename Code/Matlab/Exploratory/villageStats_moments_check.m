%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check D matrix values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd /home/sakina/Github/MFdiffusion_replication

%% qN = qP (model 1)
load('data_model_1_mom_1 trimesters .mat')

% Calculate the village avarage for each moment for each simulation
villmnD1 = cell(length(qN),1);
villmnE1 = cell(length(qN),1);
villmnS1 = cell(length(qN),1);

for i=1:length(qN)
    villmnD1{i} = mean(D{i});
    villmnE1{i} = mean(EmpiricalMoments{i});
    villmnS1{i} = mean(MeanSimulatedMoments{i});
end

% Convert to a matrix
villmnDMat1 = permute(reshape([villmnD1{:}], [], size(villmnD1, 1), size(villmnD1, 2)), [2 1]);
villmnEMat1 = permute(reshape([villmnE1{:}], [], size(villmnE1, 1), size(villmnE1, 2)), [2 1]);
villmnSMat1 = permute(reshape([villmnS1{:}], [], size(villmnS1, 1), size(villmnS1, 2)), [2 1]);

% Calculate mean
mean(villmnDMat1, 1)
mean(villmnEMat1, 1)
mean(villmnSMat1, 1)


%% qN != qP (model 2)
load('data_model_3_mom_1 trimesters .mat')

% Calculate the village avarage for each moment for each simulation
villmnD2 = cell(length(qN),length(qP));
villmnE2 = cell(length(qN),length(qP));
villmnS2 = cell(length(qN),length(qP));

for i=1:length(qN)
        for j=1:length(qP)
           villmnD2{i,j} = mean(D{i,j});
           villmnE2{i,j} = mean(EmpiricalMoments{i,j});
           villmnS2{i,j} = mean(MeanSimulatedMoments{i,j});
        end
end

% Turn villmnD as a 3 dimensional matrix (31X39X5)
villmnDMat2 = permute(reshape([villmnD2{:}], [], size(villmnD2, 1), size(villmnD2, 2)), [2 3 1]);
villmnEMat2 = permute(reshape([villmnE2{:}], [], size(villmnE2, 1), size(villmnE2, 2)), [2 3 1]);
villmnSMat2 = permute(reshape([villmnS2{:}], [], size(villmnS2, 1), size(villmnS2, 2)), [2 3 1]);

% Calculate the average for each moment
mean(squeeze(mean(villmnDMat2)), 1)
mean(squeeze(mean(villmnEMat2)), 1)
mean(squeeze(mean(villmnSMat2)), 1)
