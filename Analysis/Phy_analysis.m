 

addpath(genpath('C:\Users\John.Lee\Documents\GitHub\copy\KiloSort')) % path to kilosort folder
addpath(genpath('C:\Users\John.Lee\Documents\GitHub\npy-matlab')) % path to npy-matlab scripts
addpath('C:\DATA\OpenEphys\M44D');
addpath(genpath('C:\Users\John.Lee\Documents\GitHub\analysis-tools'))
fpath    = 'C:\DATA\OpenEphys\M44D\2019-01-10_13-56-05'; % where on disk do you want the simulation? ideally and SSD...
fs = 30000;

x = M44D1204; %// add option to choose file


kwe = 0; %// option for debugging, use when events are saved in .kwe format

%// Loading spike times and cluster information, output from Kilosort
spike_times_all = readNPY(fullfile(fpath, 'spike_times.npy'));
% spike_times_all = double(spike_times_all)/fs;
spike_times_all = double(spike_times_all);
spike_clusters = readNPY(fullfile(fpath, 'spike_clusters.npy'));
[Cids Cgroups] = readClusterGroupsCSV(fullfile(fpath, 'cluster_group.tsv'));

%/*
% for Cgroups:
% - 0 = noise
% - 1 = mua
% - 2 = good
% - 3 = unsorted
%*/
Cchannels = csvread(fullfile(fpath, 'channel_data.csv'),1,0);
channels = unique(Cchannels(:,3));
%/* first row = id, 3rd row = ch number
create this file using phy, so that you have matching info on file number, channel number
%*/

%// in the case where single recording in two different sessions
% x2 = M44D1076;

%%// Extracting raw waveforms 

%// This allows us to obtain raw/filtered waveforms for offline amplitude/ PCA analysis.

file_type = '126';
parfor ch = 1:length(channels)
    disp(['loading data from channel ' num2str(ch))
    [data1,timestamps{ch},info{ch}] = load_open_ephys_data_faster([fpath filesep file_type '_CH' num2str(channels(ch)) '.continuous' ]);    
    datach{ch} = data1;
end

%// Preprocessing data
[b,a] = butter(4, [0.0244 0.6104]);


total_sample = length(datach{1});
filtData = zeros(length(channels),total_sample);
parfor ch = 1:length(channels)
    disp(['filtering CH_' num2str(ch)])
    filtData(ch,:) = filtfilt(b,a,datach{ch});
    
end

%// common median referencing
CommonMedian = median(filtData);
st_dev = zeros(1,64);
parfor ch = 1:64
    filtData(ch,:) = filtData(ch,:)-CommonMedian;
    st_dev(ch) = median(abs(filtData(ch,:))/0.6745);
    disp(ch)
end


%// extracting waveforms from processed data
waveforms = {};
for tid = 1:length(Cchannels(:,1))
    id = Cchannels(tid,1);
    spike_times{tid} = spike_times_all(find(spike_clusters == id)); %*fs;
    waveforms.raw{tid} = zeros(length(spike_times{tid}),60);
    for t = 1:length(spike_times{tid})
        waveforms.raw{tid}(t,:) = filtData(find(channels == Cchannels(find(Cchannels(:,1)==id),3)), ...
            round(spike_times{tid}(t,1)-19):round(spike_times{tid}(t,1)+40,1));
    end
    waveforms.mean{tid} = mean(waveforms{tid},1);
    waveforms.std{tid} = std(waveforms{tid},1,1);
end


%// PCA if needed 
[Wi,score,latent] = pca(waveforms{14} );
scatter(score(:,1),score(:,2));

X = [score(:,1) score(:,2)]; 
idx = kmeans(X,3);


figure
gscatter(score(:,1),score(:,2),idx)
%         
% 
% figure
% for n=  1:80
%     plot(waveforms{4}(n,:))
%     hold on
% end
%         


%%
if kwe ~=1
    file_type = '116';
    %ch =1;
    %[data1,timestamps1,info] = load_open_ephys_data_faster([fpath filesep file_type '_CH' num2str(ch) '.continuous' ]); 
    % 
    [data, timestamps, info] = load_open_ephys_data(fullfile(fpath, 'all_channels.events'));
    rec_start_time = timestamps{ch}(1);
    start_stim_times = timestamps(find(info.eventId ==1))-rec_start_time;
    end_stim_times = timestamps(find(info.eventId ==0))-rec_start_time;
end


%%// Extracting spiketimes

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

%if kwe
if kwe ==1
    filename = 'experiment1.kwe';
    rec_start_time = 211384/fs;
    event_times = h5read('experiment1.kwe','/event_types/TTL/events/time_samples');
    event_times = double(event_times)/fs -rec_start_time;
    start_stim_times = double(event_times(1:2:end));
    end_stim_times = double(event_times(2:2:end));
end
nreps = x.stimulus_ch1(1,4);
nStim = max(x.stimulus_ch1(:,1));

%
% nStim = nStim + max(x2.stimulus_ch1(:,1));
%
TotalReps = nStim*nreps;
false_start = length(start_stim_times)-TotalReps;
start_stim_times = start_stim_times(false_start+1:end);
end_stim_times = end_stim_times(false_start+1:end);


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
            spike_times{id}>=start_stim_times(rep)-StimDur)).' ;
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

stim_label  = log(x.stimulus_ch1(:,8));

SUrate = {};
% figure
for id = 1:length(SU)
    figure
    
    set(gcf, 'Position', get(gcf,'Position').*[1 1 0 0] + [0 -600 1000 800]);
    subplot(2,2,1)
    SUrate{id}.mean = mean(rate_stim{id},2)-mean(rate_pre{id},2);
    
    
    SUrate{id}.error = std(rate_stim{id},1,2)/400;
    
    %     subplot(3,4,id)
    %     errorbar(2:2:50,SUrate{id}.mean(2:2:50),SUrate{id}.error(2:2:50))
    %
    %     hold on
    %     errorbar(1:2:49,SUrate{id}.mean(1:2:49),SUrate{id}.error(1:2:49))
    %
    
    xs = 1:25;
    h =1.0;
    for i = 1:25
        ys1(i) = gaussian_kern_reg(xs(i),xs,SUrate{id}.mean(1:2:50).',h);
        ys2(i) = gaussian_kern_reg(xs(i),xs,SUrate{id}.mean(2:2:50).',h);
    end
%     subplot(3,3,id)
%     plot(stim_label(1:2:end),ys1,'LineWidth',2);
    
    hold on
    plot(stim_label(2:2:end),ys2,'LineWidth',2);
    %     legend on
    %         axis([0 50 0 2])
    %         xticklabels({'8' '11' '14' '17' '20' '23' '26' '29' '32'})
    xlabel('Hz')
    ylabel('Average Spike count')
    title(['cluster ' num2str(Cids(SU(id)))])
    
        subplot(2,2,[2 4])
    area([0 StimDur StimDur 0],[0 TotalReps+5 0 TotalReps+5],'LineStyle','none','FaceColor',[.85 .85 1]);
    hold on
    plot(raster.spikes{id},nreps*(raster.stim{id}-1)+raster.rep{id},'k.','MarkerSize',9);
%     pause
    xlabel('time (s)')
    ylabel('reps')
    axis([-PreStim StimDur+ PostStim 0 TotalReps+1])
    hold off
    drawnow()
end
% 
% for id = 1:length(SU)
% %     subplot(ceil(sqrt(length(SU))),ceil(sqrt(length(SU))),id)
% 
% end
% 
% % 8:3:32
% MUrate = {};
% MUrate.mean = mean(rate_stimMU,2);
% test = MUrate.mean(1:2:49);
% test2=  MUrate.mean(2:2:50);
% figure
% plot(1:2:49,test)
% hold on 
% plot(2:2:50,test2)
% 
% MUrate.error = std(rate_stimMU,1);
% figure
% plot(MUrate.mean);
% edges = [-PreStim:1e-3:StimDur + PostStim];
% [N,edges] = histcounts(spikes_pooled{1, 1},edges);
% N = N*1e3;
% xs = 1:600;
% 
%  h = 10;
%  for i = 1:length(edges)-1
%      ys(i) = gaussian_kern_reg(xs(i),xs,N,h);
%  end
% 
%         plot(ys)
% 
% % figure
% % ylabel('Repetition rate (Hz)')
% 
% 
% 
% figure
% % ylabel('Repetition rate (Hz)')
% 
% 
% hold on
% title('Multi Units ')
% area([0 StimDur StimDur 0],[0 TotalReps+5 0 TotalReps+5],'LineStyle','none','FaceColor',[.85 .85 1]);
% plot(rasterMU.spikes,nreps*(rasterMU.stim-1)+rasterMU.rep,'k.','MarkerSize',9);
% xlabel('time (s)')
% ylabel('reps')
% axis([-PreStim StimDur+ PostStim 0 TotalReps+1])


