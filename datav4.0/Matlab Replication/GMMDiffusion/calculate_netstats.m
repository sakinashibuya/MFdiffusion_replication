%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main_models_1_3

% Outline

% 0. Model type
% 1. Select villages to consider
% 2. Select moments
% 3. Select time vector and number of repetitions per trial
% 4. Select parameter grid
% 5. Load data
% 6. Logistic fit to get coefficients for covariates and the constant
% 7. RUNNING THE MODEL
% 8. RUNNING THE AGGREGATOR

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
tic;
% Adjusting so this runs off dropbox
cd ..
location = pwd;
addpath(genpath(location));

%% Parameters
% 0. Model type
modelType = 3; % modelType 1 if q, modelType 3 if qN and qP, qN != qP

% 1. Select villages to consider
vills = [1:4,6,9, 12, 15, 19:21, 23:25, 29, 31:33, 36, 39, 42, 43, 45:48, 50:52, 55, 57, 59:60, 62, 64:65, 67:68, 70:73, 75];
G = length(vills); % Number of graphs

% 2. Select moments
version = 1;

switch version
    case 1
        m = 5;
    case 2
        m = 3;
    case 3
        m = 3;
    case 4
        m = 3;
end

% 3. Select time vector and number of repetitions per trial
S = 75; % Number of simulations
timeVector = 'trimesters'
TMonths = [31 35 15 35 13 2 32 5 31 35 31 29 19 22 25 25 23 23 24 25 26 24 17 16 17 13 19 20 20 19 22 14 12 15 10 19 18 18 19 19 19 17 17]; % Months

switch timeVector
    case 'months'
        T = TMonths + 1
    case 'quarters'
        T = ceil(TMonths./3) + 1 %Quarters have 3 months in them
    case 'trimesters'
        T = ceil(TMonths./4) + 1 %Trimesters have 4 months in them
end
assert(G == numel(T))


% 4. Select parameter grid
if modelType == 1,
    qN = [(0:0.001:0.01), (0.05:0.05:1)];
elseif modelType == 3,
    qN = [(0:0.001:0.01), (0.05:0.05:1)];
    qP = [(0:0.005:0.1), (0.15:0.05:1)];
end


% 5. Load data
%% Pre-allocation
X = cell(G,1);
inGiant = cell(G,1);
leaders = cell(G,1);
TakeUp = cell(G,1);
TakingLeaders = cell(G,1);
ZLeaders = cell(G,1);
Covars = [];
Outcome = [];
EmpRate = zeros(G,1);
rdist = cell(G,1);
dist = cell(G,1);
hermits = cell(G,1);
W = cell(G,1);
Z = cell(G,1);


% Load the household connection adjacency matrix.
X = load(['India Networks/adjacencymatrix.mat']);
X = X.X;

%% Construct data
counter = 0;
for vilnum = vills
    counter = counter + 1;
    
    % Load the Leader data
    templeaders = load(['./India Networks/HHhasALeader' num2str(vilnum) '.csv']);
    leaders{counter} = templeaders(:,2);
    
    % Load the Take-Up data
    TakeUp{counter} = load(['./India Networks/MF' num2str(vilnum) '.csv']);
    EmpRate(counter) = mean(TakeUp{counter}(~leaders{counter}));
    
    % Load the giant component data
    inGiant{counter} = load(['./India Networks/inGiant' num2str(vilnum) '.csv']);
    
    % Generate hermits
    d = sum(X{counter},2);
    hermits{counter}=(d==0);
    
    % Load the Covariates data
    W{counter} = load(['./India Networks/hhcovariates' num2str(vilnum) '.csv']);
    
    % Which covariates to use - for instance we want to add a PCA
    Z{counter} = [W{counter}(:,1:6)]; % for instance, take the first covariate only
    
    % prune other stats
    leaders{counter} = leaders{counter}(logical(inGiant{counter}));
    TakeUp{counter} = TakeUp{counter}(logical(inGiant{counter}));
    Z{counter} = Z{counter}(logical(inGiant{counter}),:);
    
    TakingLeaders{counter} = TakeUp{counter}(logical(leaders{counter}));
    ZLeaders{counter} = Z{counter}(logical(leaders{counter}),:);
    Outcome = [Outcome; TakingLeaders{counter}];
    Covars = [Covars; ZLeaders{counter}];
    
    % Second neighbors
    Sec{counter} = (X{counter}^2>0);
    for i=1:length(X{counter})
        Sec{counter}(i,i) = 0;
    end
    Sec{counter}=(Sec{counter}-X{counter}>0);
end


counter = 0;
for vilnum = vills
    counter = counter + 1;

    N = size(X{counter},1);
    
    [R,D] = breadthdistRAL(X{counter}, leaders{counter});
    
    minDistFromLeaders = min(D(:,logical(leaders{counter}))')';
    avgDistFromLeaders = mean(D(:,logical(leaders{counter}))')';
    
    if sum(TakeUp{counter}.*leaders{counter})>0
        minDistInfectedLeaders = min(D(:,logical(TakeUp{counter}.*leaders{counter}))')';
    else
        minDistInfectedLeaders = 0;
    end
    if sum((1-TakeUp{counter}).*leaders{counter})>0
        minDistNonInfectedLeaders = min(D(:,logical((1-TakeUp{counter}).*leaders{counter}))')';
    else
        minDistNonInfectedLeaders = 0;
    end
    
    netstats(counter).minDistFromLeaders = minDistFromLeaders;
    netstats(counter).avgDistFromLeaders = avgDistFromLeaders;
    netstats(counter).minDistInfectedLeaders = minDistInfectedLeaders;
    netstats(counter).minDistNonInfectedLeaders = minDistNonInfectedLeaders;
    
    [size(netstats(counter).minDistFromLeaders) counter];
    [size(netstats(counter).minDistInfectedLeaders) counter];
    [size(netstats(counter).minDistNonInfectedLeaders) counter];
    [size(+(minDistNonInfectedLeaders==1)) counter];
    [size(+(minDistInfectedLeaders==1)) counter];
    [size(+(minDistNonInfectedLeaders==1)) counter];
    
    netstats(counter).neighborOfInfected = ((+(minDistInfectedLeaders==1) - +(minDistNonInfectedLeaders==1))>0);
    netstats(counter).neighborOfNonInfected = ((+(minDistInfectedLeaders==1) - +(minDistNonInfectedLeaders==1))<0);

    netstats(counter).degree = sum(X{counter},2);
    netstats(counter).num_edges = sum(X{counter}(:));
    netstats(counter).leaderneighborhood = ((leaders{counter}*ones(1,N) + ones(N,1)*leaders{counter}')>0);
    netstats(counter).num_leaders = sum(leaders{counter});
    netstats(counter).num_leader_edges = (leaders{counter}'*X{counter}*leaders{counter});

    D_list{counter} = D;
end

save("netstats.mat", "netstats")
save("D.mat", "D_list")



