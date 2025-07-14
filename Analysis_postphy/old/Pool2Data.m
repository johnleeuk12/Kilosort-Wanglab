% function trial = Pool2Data(n)
% 


nTrials = length(Pool(n).eventtimes);

% preallocate structure in memory
trial = struct();
trial(nTrials).duration = 0; % preallocate


for kTrial = 1:nTrials-1
    trial_time2 = Pool(n).eventtimes(kTrial+1);
    trial_time1 = Pool(n).eventtimes(kTrial);
    duration = trial_time2-trial_time1;
    trial(kTrial).duration = duration*1e3;
    
    trial(kTrial).stim5 = ismember(Pool(n).xb.trial_type(kTrial,4),[1,3,6,8]);
    trial(kTrial).stim10 = ismember(Pool(n).xb.trial_type(kTrial,4),[2,4,5,7]);
    
    trial(kTrial).r1Hit = (Pool(n).xb.hit_code(kTrial,1) == 1)*ismember(Pool(n).xb.trial_type(kTrial,4),[1,2]);
    trial(kTrial).r2Hit = (Pool(n).xb.hit_code(kTrial,1) == 1)*ismember(Pool(n).xb.trial_type(kTrial,4),[3,4]);
    trial(kTrial).rtHit = (Pool(n).xb.hit_code(kTrial,1) == 1)*ismember(Pool(n).xb.trial_type(kTrial,4),[7,8]);
    
    trial(kTrial).r1Miss = (Pool(n).xb.hit_code(kTrial,2) == 1)*ismember(Pool(n).xb.trial_type(kTrial,4),[1,2]);
    trial(kTrial).r2Miss = (Pool(n).xb.hit_code(kTrial,2) == 1)*ismember(Pool(n).xb.trial_type(kTrial,4),[3,4]);
    trial(kTrial).rtMiss = (Pool(n).xb.hit_code(kTrial,2) == 1)*ismember(Pool(n).xb.trial_type(kTrial,4),[7,8]);
    
    trial(kTrial).r1CR = (Pool(n).xb.hit_code(kTrial,3) == 1)*ismember(Pool(n).xb.trial_type(kTrial,4),[1,2]);
    trial(kTrial).r2CR = (Pool(n).xb.hit_code(kTrial,3) == 1)*ismember(Pool(n).xb.trial_type(kTrial,4),[3,4]);
    trial(kTrial).rtCR = (Pool(n).xb.hit_code(kTrial,3) == 1)*ismember(Pool(n).xb.trial_type(kTrial,4),[7,8]);
    
    trial(kTrial).r1FA = (Pool(n).xb.hit_code(kTrial,4) == 1)*ismember(Pool(n).xb.trial_type(kTrial,4),[1,2]);
    trial(kTrial).r2FA = (Pool(n).xb.hit_code(kTrial,4) == 1)*ismember(Pool(n).xb.trial_type(kTrial,4),[3,4]);
    trial(kTrial).rtFA = (Pool(n).xb.hit_code(kTrial,4) == 1)*ismember(Pool(n).xb.trial_type(kTrial,4),[7,8]);
    
    
    trial(kTrial).stimon = 1;
    trial(kTrial).stimoff = 50;
    trial(kTrial).outoff = 50;
end



rawData = {};
rawData.trial = trial;
rawData.nTrials = nTrials;
rawData.param = [];
