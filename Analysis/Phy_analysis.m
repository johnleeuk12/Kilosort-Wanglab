function Phy_analysis()

addpath(genpath('C:\Users\John.Lee\Documents\GitHub\copy\KiloSort')) % path to kilosort folder
addpath(genpath('C:\Users\John.Lee\Documents\GitHub\npy-matlab')) % path to npy-matlab scripts
addpath('C:\DATA\OpenEphys\M44D');
addpath(genpath('C:\Users\John.Lee\Documents\GitHub\analysis-tools'))
fpath    = 'C:\DATA\OpenEphys\M44D\2018-11-21_15-31-22'; % where on disk do you want the simulation? ideally and SSD...
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

x = M44D1092;



[data, timestamps, info] = load_open_ephys_data(fullfile(fpath, 'all_channels.events'));


%% Extracting spiketimes

SU = find(Cgroups == 2);    
MU = find(Cgroups == 1);


PreStim = 0.3; %s
PostStim = 0.3; %s
StimDur = x.stimulus_ch1(1,5)*1e-3;

stim_info = x.data(find(x.data(:,3) == 1 & x.data(:,4) == -1),:); 
data_new = stim_info;
stim_info_tags = x.data_tags;
start_stim_times = timestamps(find(info.eventId ==1));
end_stim_times = timestamps(find(info.eventId ==0));
nreps = x.stimulus_ch1(4);
TotalReps = max(x.stimulus_ch1(1))*nreps;

raster.stim = {};
raster.rep = {};
raster.spikes = {};
spikes_pooled = {};
spike_times = {};

for id = 1:length(SU)
    raster.stim{id} = [];
    raster.rep{id} = [];
    raster.spikes{id} = [];
    spikes_pooled{id} = [];   
    spike_times{id} = spike_times_all(find(spike_clusters == Cids(SU(id))));
    for rep = 1:TotalReps
        spikes1 = spike_times{id}(find(spike_times{id}>=start_stim_times(rep)-PreStim & ...
            spike_times{id}<=end_stim_times(rep)+ PostStim)).';
        spikes1 = spikes1 - start_stim_times(rep);
        spikes_pooled{id} = [spikes_pooled{id} spikes1];
        raster.stim{id} = [raster.stim{id} data_new(rep,1)*ones(size(spikes1))];
        raster.rep{id} = [raster.rep{id} data_new(rep,2)*ones(size(spikes1))];
        raster.spikes{id} = [raster.spikes{id} spikes1];
    end
end


%extract MU activity

rasterMU.stim = [];
rasterMU.rep = [];
rasterMU.spikes = [];
spike_timesMU = [];

for id = 1:length(MU)
    spike_timesMU = spike_times_all(find(spike_clusters == Cids(MU(id))));
    for rep = 1:TotalReps
        spikes = spike_timesMU(find(spike_timesMU >= start_stim_times(rep)-PreStim & ...
            spike_timesMU <= end_stim_times(rep)+ PostStim)).';
        spikes = spikes - start_stim_times(rep);
        rasterMU.stim = [rasterMU.stim data_new(rep,1)*ones(size(spikes))];
        rasterMU.rep = [rasterMU.rep data_new(rep,2)*ones(size(spikes))];
        rasterMU.spikes = [rasterMU.spikes spikes];
    end
end



        

figure
% ylabel('Repetition rate (Hz)')
for id = 1:length(SU)
    subplot(ceil(sqrt(length(SU))),ceil(sqrt(length(SU))),id)
    hold on
    title(['cluster nb ' num2str(id)])
    area([0 StimDur StimDur 0],[0 TotalReps+5 0 TotalReps+5],'LineStyle','none','FaceColor',[.85 .85 1]);
    plot(raster.spikes{id},nreps*(raster.stim{id}-1)+raster.rep{id},'k.','MarkerSize',9);
    xlabel('time (s)')
    ylabel('reps')
    axis([-0.5 1.5 0 TotalReps+1])
end


figure
% ylabel('Repetition rate (Hz)')


hold on
title('Multi Units ')
area([0 StimDur StimDur 0],[0 TotalReps+5 0 TotalReps+5],'LineStyle','none','FaceColor',[.85 .85 1]);
plot(rasterMU.spikes,nreps*(rasterMU.stim-1)+rasterMU.rep,'k.','MarkerSize',9);
xlabel('time (s)')
ylabel('reps')
axis([-0.3 1.3 0 TotalReps+1])


