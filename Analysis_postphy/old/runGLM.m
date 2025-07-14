%% Load data

addpath('D:\GitHub\glmnet_matlab');
nTrials = rawData.nTrials; % number of trials
unitOfTime = 'ms';
binSize = 50; % TODO some continuous observations might need up/down-sampling if binSize is not 1!?
tvList = fieldnames(rawData.trial);



%% initial run, single variable model
expt = buildGLM.initExperiment(unitOfTime, binSize, [], rawData.param);
% expt = buildGLM.registerTiming(expt, 'stimon', 'stim onset');
% expt = buildGLM.registerTiming(expt, 'stimoff', 'stim offset');
% expt = buildGLM.registerTiming(expt, 'outoff', 'stim offset');
expt.trial = rawData.trial;

y = getBinnedSpikeTrain2(expt, Pool,n);
% y = full(y);
y = y.';
for f = 1:length(tvList)
    expt = buildGLM.registerValue(expt,tvList{f},[]);
end
%% 
binfun = expt.binfun;
stimHandle = {};
stimHandle{1} = @(trial, expt) trial.stim5 * basisFactory.boxcarStim(binfun(trial.stimon), binfun(trial.stimoff), binfun(trial.duration));
stimHandle{2} = @(trial, expt) trial.stim10 * basisFactory.boxcarStim(binfun(trial.stimon), binfun(trial.stimoff), binfun(trial.duration));
stimHandle{3} = @(trial, expt) trial.r1Hit * basisFactory.boxcarStim(binfun(trial.stimon), binfun(trial.stimoff), binfun(trial.duration));
stimHandle{4} = @(trial, expt) trial.r2Hit * basisFactory.boxcarStim(binfun(trial.stimon), binfun(trial.stimoff), binfun(trial.duration));
stimHandle{5} = @(trial, expt) trial.rtHit * basisFactory.boxcarStim(binfun(trial.stimon), binfun(trial.stimoff), binfun(trial.duration));
stimHandle{6} = @(trial, expt) trial.r1Miss * basisFactory.boxcarStim(binfun(trial.stimon), binfun(trial.stimoff), binfun(trial.duration));
stimHandle{7} = @(trial, expt) trial.r2Miss * basisFactory.boxcarStim(binfun(trial.stimon), binfun(trial.stimoff), binfun(trial.duration));
stimHandle{8} = @(trial, expt) trial.rtMiss * basisFactory.boxcarStim(binfun(trial.stimon), binfun(trial.stimoff), binfun(trial.duration));
stimHandle{9} = @(trial, expt) trial.r1CR * basisFactory.boxcarStim(binfun(trial.stimon), binfun(trial.stimoff), binfun(trial.duration));
stimHandle{10} = @(trial, expt) trial.r2CR * basisFactory.boxcarStim(binfun(trial.stimon), binfun(trial.stimoff), binfun(trial.duration));
stimHandle{11} = @(trial, expt) trial.rtCR * basisFactory.boxcarStim(binfun(trial.stimon), binfun(trial.stimoff), binfun(trial.duration));
stimHandle{12} = @(trial, expt) trial.r1FA * basisFactory.boxcarStim(binfun(trial.stimon), binfun(trial.stimoff), binfun(trial.duration));
stimHandle{13} = @(trial, expt) trial.r2FA * basisFactory.boxcarStim(binfun(trial.stimon), binfun(trial.stimoff), binfun(trial.duration));
stimHandle{14} = @(trial, expt) trial.rtFA * basisFactory.boxcarStim(binfun(trial.stimon), binfun(trial.stimoff), binfun(trial.duration));




%%


error = zeros(1,14);
for f = 1:14 %length(tvList)
    dspec = buildGLM.initDesignSpec(expt);
    binfun = expt.binfun;
    if f <3
        bs = basisFactory.makeSmoothTemporalBasis('raised cosine', 1500, 30, binfun);
    else
        bs = basisFactory.makeSmoothTemporalBasis('raised cosine', 7000, 70, binfun);
    end
    dspec = buildGLM.addCovariate(dspec, tvList{f}, tvList{f}, stimHandle{f}, bs);
    trialIndices = 1:(nTrials-1); % use all trials except the last one
    dm = buildGLM.compileSparseDesignMatrix(dspec, trialIndices);
    dm = buildGLM.removeConstantCols(dm);
    researchLambda=1;
%     endTrialIndices = cumsum(binfun([expt.trial(trialIndices).duration]));
%     X = dm.X(1:endTrialIndices(3),:);
%     X = dm.X;
%     mv = max(abs(X), [], 1); mv(isnan(mv)) = 1;
%     X = bsxfun(@times, X, 1 ./ mv);
%     figure(742); clf; imagesc(X.');
    disp('making cv specific lambda search since full lambda value could not converge in this kfold')
    mdl =  fitglm(dm.X,full(y(1:end-1)),'linear','Distribution','poisson');
    error(f) = mdl.SSE;
    
end
%     cvfit = cvglmnet(dm.X,y(1:end-1),'poisson',options,[],[],[],false); % do not do parallel processing as you are already inside a parfor loop
%     bestLambda=find(cvfit.lambda==cvfit.lambda_min); % I am using minimum lambda (best cross validated performance), alternatively use llse)
%     w= cvfit.glmnet_fit.beta(:,bestLambda);
%     












%% Specify the fields to load
expt = buildGLM.initExperiment(unitOfTime, binSize, [], rawData.param);
expt = buildGLM.registerSpikeTrain(expt, 'sptrain', 'Our Neuron'); % Spike train!!!
expt = buildGLM.registerValue(expt, 'stim5', '5kHz'); % information on the trial, but not associated with time
expt = buildGLM.registerValue(expt, 'stim10', '10kHz'); % information on the trial, but not associated with time
expt = buildGLM.registerTiming(expt, 'stimon', 'stim onset');
expt = buildGLM.registerTiming(expt, 'stimoff', 'stim offset');
expt = buildGLM.registerTiming(expt, 'outoff', 'stim offset');

expt = buildGLM.registerValue(expt, 'r1Hit', 'r1 Hit'); % information on the trial, but not associated with time
expt = buildGLM.registerValue(expt, 'r2Hit', 'r2 Hit'); % information on the trial, but not associated with time
expt = buildGLM.registerValue(expt, 'rtHit', 'r2 Hit'); % information on the trial, but not associated with time

expt = buildGLM.registerValue(expt, 'r1Miss', 'r1 Hit'); % information on the trial, but not associated with time
expt = buildGLM.registerValue(expt, 'r2Miss', 'r2 Hit'); % information on the trial, but not associated with time
expt = buildGLM.registerValue(expt, 'rtMiss', 'r2 Hit'); % information on the trial, but not associated with time

expt = buildGLM.registerValue(expt, 'r1CR', 'r1 Hit'); % information on the trial, but not associated with time
expt = buildGLM.registerValue(expt, 'r2CR', 'r2 Hit'); % information on the trial, but not associated with time
expt = buildGLM.registerValue(expt, 'rtCR', 'r2 Hit'); % information on the trial, but not associated with time

expt = buildGLM.registerValue(expt, 'r1FA', 'r1 Hit'); % information on the trial, but not associated with time
expt = buildGLM.registerValue(expt, 'r2FA', 'r2 Hit'); % information on the trial, but not associated with time
expt = buildGLM.registerValue(expt, 'rtFA', 'r2 Hit'); % information on the trial, but not associated with time

%% Convert the raw data into the experiment structure
expt.trial = rawData.trial;


%% Build 'designSpec' which specifies how to generate the design matrix
% Each covariate to include in the model and analysis is specified.
dspec = buildGLM.initDesignSpec(expt);
binfun = expt.binfun;
% bs = basisFactory.makeSmoothTemporalBasis('boxcar', 100, 10, binfun);
% bs.B = 0.1 * bs.B;


%% stim
% a box car that depends on the coh value
bs = basisFactory.makeSmoothTemporalBasis('raised cosine', 1500, 15, binfun);
stimHandle = @(trial, expt) trial.stim5 * basisFactory.boxcarStim(binfun(trial.stimon), binfun(trial.stimoff), binfun(trial.duration));

dspec = buildGLM.addCovariate(dspec, 'stimKer5', '5kHz', stimHandle, bs);

stimHandle = @(trial, expt) trial.stim10 * basisFactory.boxcarStim(binfun(trial.stimon), binfun(trial.stimoff), binfun(trial.duration));

dspec = buildGLM.addCovariate(dspec, 'stimKer10', '10kHz', stimHandle, bs);

%% outcome

bs = basisFactory.makeSmoothTemporalBasis('raised cosine', 7000, 7, binfun);
stimHandle = @(trial, expt) trial.r1Hit * basisFactory.boxcarStim(binfun(trial.stimon), binfun(trial.outoff), binfun(trial.duration));
dspec = buildGLM.addCovariate(dspec, 'r1HitKer', 'r1Hit', stimHandle, bs);

stimHandle = @(trial, expt) trial.r2Hit * basisFactory.boxcarStim(binfun(trial.stimon), binfun(trial.outoff), binfun(trial.duration));
dspec = buildGLM.addCovariate(dspec, 'r2HitKer', 'r2Hit', stimHandle, bs);

stimHandle = @(trial, expt) trial.rtHit * basisFactory.boxcarStim(binfun(trial.stimon), binfun(trial.outoff), binfun(trial.duration));
dspec = buildGLM.addCovariate(dspec, 'rtHitKer', 'rtHit', stimHandle, bs);

stimHandle = @(trial, expt) trial.r1Miss * basisFactory.boxcarStim(binfun(trial.stimon), binfun(trial.outoff), binfun(trial.duration));
dspec = buildGLM.addCovariate(dspec, 'r1MissKer', 'r1Miss', stimHandle, bs);

stimHandle = @(trial, expt) trial.r2Miss * basisFactory.boxcarStim(binfun(trial.stimon), binfun(trial.outoff), binfun(trial.duration));
dspec = buildGLM.addCovariate(dspec, 'r2MissKer', 'r2Miss', stimHandle, bs);

stimHandle = @(trial, expt) trial.rtMiss * basisFactory.boxcarStim(binfun(trial.stimon), binfun(trial.outoff), binfun(trial.duration));
dspec = buildGLM.addCovariate(dspec, 'rtMissKer', 'rtMiss', stimHandle, bs);

stimHandle = @(trial, expt) trial.r1CR * basisFactory.boxcarStim(binfun(trial.stimon), binfun(trial.outoff), binfun(trial.duration));
dspec = buildGLM.addCovariate(dspec, 'r1CRKer', 'r1CR', stimHandle, bs);

stimHandle = @(trial, expt) trial.r2CR * basisFactory.boxcarStim(binfun(trial.stimon), binfun(trial.outoff), binfun(trial.duration));
dspec = buildGLM.addCovariate(dspec, 'r2CRKer', 'r2CR', stimHandle, bs);

stimHandle = @(trial, expt) trial.rtCR * basisFactory.boxcarStim(binfun(trial.stimon), binfun(trial.outoff), binfun(trial.duration));
dspec = buildGLM.addCovariate(dspec, 'rtCRKer', 'rtCR', stimHandle, bs);

stimHandle = @(trial, expt) trial.r1FA * basisFactory.boxcarStim(binfun(trial.stimon), binfun(trial.outoff), binfun(trial.duration));
dspec = buildGLM.addCovariate(dspec, 'r1FAKer', 'r1FA', stimHandle, bs);

stimHandle = @(trial, expt) trial.r2FA * basisFactory.boxcarStim(binfun(trial.stimon), binfun(trial.outoff), binfun(trial.duration));
dspec = buildGLM.addCovariate(dspec, 'r2FAKer', 'r2FA', stimHandle, bs);

stimHandle = @(trial, expt) trial.rtFA * basisFactory.boxcarStim(binfun(trial.stimon), binfun(trial.outoff), binfun(trial.duration));
dspec = buildGLM.addCovariate(dspec, 'rtfAKer', 'rtFA', stimHandle, bs);

%% Compile the data into 'DesignMatrix' structure
trialIndices = 1:10; %(nTrials-1); % use all trials except the last one

dm = buildGLM.compileSparseDesignMatrix(dspec, trialIndices);



%% Visualize the design matrix
endTrialIndices = cumsum(binfun([expt.trial(trialIndices).duration]));
% X = dm.X(1:endTrialIndices(3),:);
X = dm.X;
mv = max(abs(X), [], 1); mv(isnan(mv)) = 1;
X = bsxfun(@times, X, 1 ./ mv);
figure(742); clf; imagesc(X.');
% buildGLM.visualizeDesignMatrix(dm, 1); % optionally plot the first trial
