
function [raster2, rate, Lick2, Lick_rate] = gather_raster_ephys2(PreStim, PostStim, Pool)

%{
from variable Pool, gather spk data in raster and rate format.
Gather Lick into lick raster and lick rate format.

Edit log
==========================================================================
2025/02/11 JHL
data_new variable now reads from column 5,6 which depicts TTR instead of
conditioning stages 


%}


% PreStim and PostStim should be in s. 
fprintf('calculating rate and raster \n')

% initialize variables

rate = {};
raster = {};
Lick_rate = {};
Lick_raster = {};
Lick_raster.stim = {};
Lick_raster.rep = {};
Lick_raster.spikes = {};
licks_pooled = {};



raster.stim = {};
raster.rep = {};
raster.spikes = {};
spikes_pooled = {};


rate.PSTH = {};
rate.behav = {};
Lick_rate.PSTH = {};

raster2 ={};
Lick2 = {};


tic
for p = 1:length(Pool)
    
    if mod(p,10) ==1
        fprintf(['%4d /' num2str(length(Pool)) ' time : %6.2f sec \n'],p,toc')
    end

    data_new = Pool(p).xb.trial_type(:,[5,6]); % modified from gather_raster_ephys
    hit_code = zeros(1,length(data_new));
    for h = 1:length(data_new)
        if h <201
            hit_code(1,h) = find(Pool(p).xb.hit_code(h,:) == 1);
%         elseif 200 <h & h<261
%             hit_code(1,h) = 0;
        elseif h >260 && h<461
             hit_code(1,h) = find(Pool(p).xb.hit_code(h-60,:) == 1);
        end
    end
    raster.stim{p} = [];
    raster.rep{p} = [];
    raster.spikes{p} = [];
    spikes_pooled{p} = [];
    
    Lick_raster.stim{p} = [];
    Lick_raster.rep{p} = [];
    Lick_raster.spikes{p} = [];
    licks_pooled{p} = [];

    
    
    
    TotalReps= size(Pool(p).xb.trial_type,1);
    rep_limit= 520;
    if TotalReps > rep_limit
        TotalReps = rep_limit;
    end
    % for raster
    for rep = 1:TotalReps
        spikes1 = Pool(p).spiketimes(find(Pool(p).spiketimes>=Pool(p).eventtimes(rep)-PreStim & ...
            Pool(p).spiketimes<=Pool(p).eventtimes(rep)+ PostStim)).';
        spikes1 = spikes1 - Pool(p).eventtimes(rep);
        spikes_pooled{p} = [spikes_pooled{p} spikes1];
        raster.stim{p} = [raster.stim{p} data_new(rep,1)*ones(size(spikes1))];
        raster.rep{p} = [raster.rep{p} data_new(rep,2)*ones(size(spikes1))];
        raster.spikes{p} = [raster.spikes{p} spikes1];
        
        licks1 = Pool(p).licktimes(find(Pool(p).licktimes>=Pool(p).eventtimes(rep)-PreStim & ...
            Pool(p).licktimes<=Pool(p).eventtimes(rep)+ PostStim)).';
        licks1 = licks1 - Pool(p).eventtimes(rep);
        licks_pooled{p} = [licks_pooled{p} spikes1];
        Lick_raster.stim{p} = [Lick_raster.stim{p} data_new(rep,1)*ones(size(licks1))];
        Lick_raster.rep{p} = [Lick_raster.rep{p} data_new(rep,2)*ones(size(licks1))];
        Lick_raster.spikes{p} = [Lick_raster.spikes{p} licks1];
        
    
    
    %for rate
        
        TrialLength = PreStim + PostStim;
        rep_rate = zeros(1,(round(TrialLength)*1e3));
        spikes4rate = spikes1 + PreStim;
        for st = spikes4rate
            if st <= 0
                st1 = 1;
                if ceil(st1*1e3) <= length(rep_rate)
                    rep_rate(1,ceil(st1*1e3)) = rep_rate(1,ceil(st1*1e3))+1;
                end
            
            elseif ceil(st*1e3) <= length(rep_rate)
                rep_rate(1,ceil(st*1e3)) = rep_rate(1,ceil(st*1e3))+1;
            end
            
        end
        
        lick_rep_rate = zeros(1,(round(TrialLength)*1e3));
        licks4rate = licks1 + PreStim;
        for st = licks4rate
            if st <= 0
                st1 = 1;
                if ceil(st1*1e3) <= length(lick_rep_rate)
                    lick_rep_rate(1,ceil(st1*1e3)) = lick_rep_rate(1,ceil(st1*1e3))+1;
                end
            
            elseif ceil(st*1e3) <= length(lick_rep_rate)
                lick_rep_rate(1,ceil(st*1e3)) = lick_rep_rate(1,ceil(st*1e3))+1;
            end
            
        end
        
        
%         rate_total = [rate_total ; rate*1000];
        % calculate PSTH
        rate.behav{p,data_new(rep,1)}(data_new(rep,2)) = hit_code(rep);
        rate.PSTH{p,data_new(rep,1)}(data_new(rep,2),:) = rep_rate*1e3;
        raster.nrep{p} = [];
        for st  = 1:size(rate.PSTH(p,:),2)
            raster.nrep{p} = [raster.nrep{p},size(rate.PSTH{p,st},1)];
        end
        
        Lick_rate.PSTH{p,data_new(rep,1)}(data_new(rep,2),:) = lick_rep_rate*1e3;
        
        
        
    end
    raster2{p} = {};
    Lick2{p} = {};
    stim_list = unique(data_new(:,1));
    for st  = 1:length(stim_list)
        ind = raster.stim{p} == stim_list(st);
        raster2{p}(stim_list(st)).rep = raster.rep{p}(1,ind);
        raster2{p}(stim_list(st)).spikes = raster.spikes{p}(1,ind);
        raster2{p}(stim_list(st)).nreps = raster.nrep{p}(stim_list(st));
        
        ind2 = Lick_raster.stim{p} ==stim_list(st);
        Lick2{p}(stim_list(st)).rep = Lick_raster.rep{p}(1,ind2);
        Lick2{p}(stim_list(st)).licks = Lick_raster.spikes{p}(1,ind2);
        Lick2{p}(stim_list(st)).nreps = raster.nrep{p}(stim_list(st));
    end
    
    
    
end

% % nreps =30;
% % p = 1;
% figure(p)
% stim_dur = 0.5;
% delay = 1.0;
% reward = 2.0;
% TotalReps= size(Pool(p).xb.trial_type,1);
% 
% nreps1 = cumsum(raster.nrep{p});
%         rectangle('Position',[0 0,stim_dur nreps1(end)],'FaceColor',[0.9,0.9,0.9],'EdgeColor','none')
%         
%         rectangle('Position',[(stim_dur + delay) 0,(stim_dur + delay+reward) nreps1(end)],'FaceColor',[0.9,0.9,0.9],'EdgeColor','none')
%         hold on
% nreps = nreps1([raster.stim{p}]);
% plot(raster.spikes{p},nreps+raster.rep{p}-nreps1(1),'k.','MarkerSize',15);
% %         yticks([1:nreps*2:nreps*stim_num]+5)
%         stim_name = {};
%         for st = 1:stim_num
%             stim_name{st} = num2str(stim_set(st));
%         end
%         yticklabels(stim_name(1:2:end));
%         xlabel('time (s)')
%         axis([-PreStim max(stim_dur) + PostStim 0 nreps*stim_num+1])
%         set(gcf, 'Position', [800 800 1000 400]); 


