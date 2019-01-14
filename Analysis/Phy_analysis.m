function Phy_analysis()

addpath(genpath('C:\Users\John.Lee\Documents\GitHub\copy\KiloSort')) % path to kilosort folder
addpath(genpath('C:\Users\John.Lee\Documents\GitHub\npy-matlab')) % path to npy-matlab scripts
addpath('C:\DATA\OpenEphys\M44D');
addpath(genpath('C:\Users\John.Lee\Documents\GitHub\analysis-tools'))
fpath    = 'C:\DATA\OpenEphys\M44D\2019-01-10_13-56-05'; % where on disk do you want the simulation? ideally and SSD...
fs = 30000;

spike_times_all = readNPY(fullfile(fpath, 'spike_times.npy'));
spike_times_all = double(spike_times_all)/fs;
spike_clusters = readNPY(fullfile(fpath, 'spike_clusters.npy'));

[Cids Cgroups] = readClusterGroupsCSV(fullfile(fpath, 'cluster_group.tsv'));
% for Cgroups:
% - 0 = noise
% - 1 = mua
% - 2 = good
% - 3 = unsorted

x = M44D1084;
% in the case where single recording in two different sessions
% x2 = M44D1076;


[data, timestamps, info] = load_open_ephys_data(fullfile(fpath, 'all_channels.events'));


%% Extracting spiketimes

SU = find(Cgroups == 2);    
MU = find(Cgroups == 1);


PreStim = x.pre_stimulus_record_time*1e-3; %s
PostStim = x.post_stimulus_record_time*1e-3; %s
StimDur = x.stimulus_ch1(1,5)*1e-3;

stim_info = x.data(find(x.data(:,3) == 1 & x.data(:,4) == -1),:); 

% stim_info2 = x2.data(find(x2.data(:,3) == 1 & x2.data(:,4) == -1),:); 
% stim_info2(:,1) = stim_info2(:,1)+ max(stim_info(:,1));
% stim_info2(:,2) = stim_info2(:,2)+ 20;
%
% stim_info = [stim_info ; x2.data(find(x2.data(:,3) == 1 & x2.data(:,4) == -1),:)]; 
%
data_new = stim_info;
% data_new = [stim_info ; stim_info2];

% stim_info_tags = x.data_tags;
start_stim_times = timestamps(find(info.eventId ==1));
end_stim_times = timestamps(find(info.eventId ==0));


%if kwe
filename = 'experiment1.kwe';
rec_start_time = 211384/fs;
event_times = h5read('experiment1.kwe','/event_types/TTL/events/time_samples');
event_times = event_times/fs-rec_start_time;
start_stim_times = double(event_times(1:2:end));
end_stim_times = double(event_times(2:2:end));

nreps = x.stimulus_ch1(1,4);
nStim = max(x.stimulus_ch1(:,1));

%
% nStim = nStim + max(x2.stimulus_ch1(:,1));
%
TotalReps = nStim*nreps;

raster.stim = {};
raster.rep = {};
raster.spikes = {};
spikes_pooled = {};
spike_times = {};
rate_stim = {};


for id = 1:length(SU)
    raster.stim{id} = [];
    raster.rep{id} = [];
    raster.spikes{id} = [];
    spikes_pooled{id} = [];  
    rate_stim{id} = [];
    spike_times{id} = spike_times_all(find(spike_clusters == Cids(SU(id))));
    for rep = 1:TotalReps
        %for raster
        spikes1 = spike_times{id}(find(spike_times{id}>=start_stim_times(rep)-PreStim & ...
            spike_times{id}<=end_stim_times(rep)+ PostStim)).';
        spikes1 = spikes1 - start_stim_times(rep);
        spikes_pooled{id} = [spikes_pooled{id} spikes1];
        raster.stim{id} = [raster.stim{id} data_new(rep,1)*ones(size(spikes1))];
        raster.rep{id} = [raster.rep{id} data_new(rep,2)*ones(size(spikes1))];
        raster.spikes{id} = [raster.spikes{id} spikes1];
        
        %for rate
        spikes2 = spike_times{id}(find(spike_times{id}>=start_stim_times(rep) & ...
            spike_times{id}<=end_stim_times(rep))).';
        rate_stim{id}(data_new(rep,1),data_new(rep,2)) = length(spikes2);
        spikes3 = spike_times{id}(find(spike_times{id}<=start_stim_times(rep) & ...
            spike_times{id}>=start_stim_times(rep)-PreStim)).';
        rate_pre{id}(data_new(rep,1),data_new(rep,2)) = length(spikes3);
    end
end


%extract MU activity

rasterMU.stim = [];
rasterMU.rep = [];
rasterMU.spikes = [];
spike_timesMU = [];
rate_stimMU =[];
for id = 1:length(MU)
    spike_timesMU = spike_times_all(find(spike_clusters == Cids(MU(id))));
    for rep = 1:TotalReps
        spikes = spike_timesMU(find(spike_timesMU >= start_stim_times(rep)-PreStim & ...
            spike_timesMU <= end_stim_times(rep)+ PostStim)).';
        spikes = spikes - start_stim_times(rep);
        rasterMU.stim = [rasterMU.stim data_new(rep,1)*ones(size(spikes))];
        rasterMU.rep = [rasterMU.rep data_new(rep,2)*ones(size(spikes))];
        rasterMU.spikes = [rasterMU.spikes spikes];
        spikes2 = spike_timesMU(find(spike_timesMU >= start_stim_times(rep) & ...
            spike_timesMU <= end_stim_times(rep))).';
        rate_stimMU(data_new(rep,1),data_new(rep,2)) = length(spikes2);
    end
end

stim_label  = x.stimulus_ch1(:,8);

SUrate = {};
figure
for id = 1:length(SU)
    
    SUrate{id}.mean = mean(rate_stim{id},2)-mean(rate_pre{id},2);
    
    
    SUrate{id}.error = std(rate_stim{id},1,2)/400;
    
    %     subplot(3,4,id)
    %     errorbar(2:2:50,SUrate{id}.mean(2:2:50),SUrate{id}.error(2:2:50))
    %
    %     hold on
    %     errorbar(1:2:49,SUrate{id}.mean(1:2:49),SUrate{id}.error(1:2:49))
    %
    
    xs = 1:37;
    h =1.0;
    for i = 1:37
        ys1(i) = gaussian_kern_reg(xs(i),xs,SUrate{id}.mean.',h);
        %         ys2(i) = gaussian_kern_reg(xs(i),xs,SUrate{id}.mean(2:2:50).',h);
    end
        subplot(4,5,id)
    plot(stim_label,ys1,'LineWidth',2);
    
    %     hold on
    %     plot(2:2:50,ys2,'LineWidth',2);
    legend on
    %     axis([0 50 0 2])
    %     xticklabels({'8' '11' '14' '17' '20' '23' '26' '29' '32'})
    xlabel('Hz')
    ylabel('Average Spike count')
    drawnow()
end

% 8:3:32
MUrate = {};
MUrate.mean = mean(rate_stimMU,2);
test = MUrate.mean(1:2:49);
test2=  MUrate.mean(2:2:50);
figure
plot(1:2:49,test)
hold on 
plot(2:2:50,test2)

MUrate.error = std(rate_stimMU,1);
figure
plot(MUrate.mean);
edges = [-PreStim:1e-3:StimDur + PostStim];
[N,edges] = histcounts(spikes_pooled{1, 1},edges);
N = N*1e3;
xs = 1:600;

 h = 10;
 for i = 1:length(edges)-1
     ys(i) = gaussian_kern_reg(xs(i),xs,N,h);
 end

        plot(ys)

% figure
% ylabel('Repetition rate (Hz)')
for id = 1:length(SU)
%     subplot(ceil(sqrt(length(SU))),ceil(sqrt(length(SU))),id)
    figure
    hold on
    title(['cluster nb ' num2str(id)])
    area([0 StimDur StimDur 0],[0 TotalReps+5 0 TotalReps+5],'LineStyle','none','FaceColor',[.85 .85 1]);
    plot(raster.spikes{id},nreps*(raster.stim{id}-1)+raster.rep{id},'k.','MarkerSize',9);
    xlabel('time (s)')
    ylabel('reps')
    axis([-PreStim StimDur+ PostStim 0 TotalReps+1])
end


figure
% ylabel('Repetition rate (Hz)')


hold on
title('Multi Units ')
area([0 StimDur StimDur 0],[0 TotalReps+5 0 TotalReps+5],'LineStyle','none','FaceColor',[.85 .85 1]);
plot(rasterMU.spikes,nreps*(rasterMU.stim-1)+rasterMU.rep,'k.','MarkerSize',9);
xlabel('time (s)')
ylabel('reps')
axis([-PreStim StimDur+ PostStim 0 TotalReps+1])


