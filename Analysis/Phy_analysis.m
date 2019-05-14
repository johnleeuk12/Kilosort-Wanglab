function Phy_analysis()


PC_name = '426_Analysis';
animal_name = 'M94W';
session_name = '2019-04-05_14-05-45';

fpath = directories(PC_name,animal_name,session_name);

fs = 30000;

addpath(genpath('C:\DATA\xbz'));




x = M94W0239; %// add option to choose file
Nb_ch = 64;

kwe = 0; %// option for debugging, use when events are saved in .kwe format

%// Loading spike times and cluster information, output from Kilosort
spike_times_all = readNPY(fullfile(fpath, 'spike_times.npy'));
spike_times_all = double(spike_times_all);
spike_clusters = readNPY(fullfile(fpath, 'spike_clusters.npy'));
[Cids Cgroups] = readClusterGroupsCSV(fullfile(fpath, 'cluster_group.tsv'));
spike_templates = readNPY(fullfile(fpath, 'spike_templates.npy'));
Cchannels = [];
templates = readNPY(fullfile(fpath, 'templates.npy'));
%// Find channel information for each cluster
for Cid = Cids
    if Cid <= max(spike_templates)
        [~, Ch] = max(mean(abs(templates(Cid+1,:,:))));
    else
        Cid2 = spike_templates(find(spike_clusters == Cid,1));
        [~, Ch] = max(mean(abs(templates(Cid2+1,:,:))));
    end
    Cchannels = [Cchannels Ch];
    
end

Cchannels = Cchannels.';

%/*
% for Cgroups:
% - 0 = noise
% - 1 = mua
% - 2 = good
% - 3 = unsorted
make_HarrisChannelMap(fpath);


chanMap = [27 32 21 3 25 30 19 5 23 17 24 7 20 29 26 9 ...
    22 31 28 11 16 13 1 15 18 12 14 10 8 4 6 2 ...
    64 60 62 58 56 52 54 48 49 53 51 55 57 61 59 63 ...
    38 42 40 45 43 35 33 47 36 50 34 44 46 39 41 37];



channels  = chanMap(Cchannels);

data_all = [];
timestamps1 = {};
info1 = {};
file_type = '100';
tic
fprintf('Time %3.0fs. Loading data... \n', toc);
parfor ch = 1:Nb_ch
%     disp(['loading data from channel ' num2str(ch)])
    [data1,timestamps1{ch},info1{ch}] = load_open_ephys_data_faster([fpath filesep file_type '_CH' num2str(ch) '.continuous' ]);    
    data_all(ch,:) = data1;
end

%/ the channel number (from 1 to 64) in the data is arranged by channel 
%site id and not by channel mapping.

%// filtering data
fprintf('Time %3.0fs. filtering data... \n', toc);

gpu_data = gpuArray(int16(data_all));
filtData = single([]);
% [b,a] = butter(4, [0.0244 0.6104],'bandpass');
[b,a] = butter(4, 500/fs,'low');
Nbuff = 32*32;
Nbatch = length(data_all(1,:))/Nbuff;
for batch = 1:Nbatch
    ipoint = (batch-1)*Nbuff +1;
    datr = filter(b,a,single(gpu_data(:,ipoint:ipoint+Nbuff-1)).');
    datr = flipud(datr);
    datr = filter(b,a,datr);
    datr = flipud(datr);
    filtData(:,ipoint:ipoint+Nbuff-1) = gather_try(int16(datr.'));
end


chanMap2 = reshape(chanMap,4,16).';
dataCSD = zeros(30000,4,16);
for nx = 1:4
    for ny = 1:16
        dataCSD(:,nx,ny) = filtData(chanMap2(ny,nx),1:30000);
    end
end
%  
filtData3 = filtfilt(b,a,data_all.').';
filtData2 = zeros(64,length(filtData));
for ch = 1:64
    filtData2(ch,:) = filtData3(chanMap(ch),:);
end
%  
figure
for ch = 1:64
    plot(filtData2(ch,1:30000)+ch*200)
    hold on
end
csdout = CSD(filtData2(:,1:1e5).',30000,2e-5,'inverse',1e-4);
% filtData2 = [];
% for ch = 1:64
% filtData2(:,ch) = filtfilt(b,a,data_all(ch,:));
% end
% filtData2 = filtData2.';
% 
% // common median referencing
fprintf('Time %3.0fs. Common median referencing... \n', toc);
CommonMedian = median(filtData);
st_dev = zeros(1,Nb_ch);
for ch = 1:Nb_ch
    filtData(ch,:) = filtData(ch,:)-CommonMedian;
    st_dev(ch) = median(abs(filtData(ch,:))/0.6745);
end

fprintf('Time %3.0fs. Preprocessing complete! \n', toc);

filtData = double(filtData);
% 
% plot raw data
% figure
% for ch = 1:64
%     plot(filtData(ch,300000:400000)-chanMap(ch)*250)
%     drawnow
%     hold on
% end


fprintf('Time %3.0fs. Extracting waveforms... \n', toc);
 %// extracting waveforms from processed data
waveforms.raw = {};
waveforms.mean = {};
waveforms.std = {};
for tid = 1:length(Cids)
    id = Cids(tid);
    spike_times{tid} = spike_times_all(find(spike_clusters == id)); %*fs;
    waveforms.raw{tid} = zeros(length(spike_times{tid}),60);
    for t = 1:length(spike_times{tid})
        waveforms.raw{tid}(t,:) = filtData(channels(tid), ...
            round(spike_times{tid}(t,1)-19):round(spike_times{tid}(t,1)+40,1)); %
    end
    waveforms.mean{tid} = mean(waveforms.raw{tid},1);
    waveforms.std{tid} = std(waveforms.raw{tid},[],1);
    peak_to_peak = max(waveforms.mean{tid})-min(waveforms.mean{tid});
    noise = mean(waveforms.std{tid}(1:10));
    waveforms.SNR{tid} = 20.*log10(peak_to_peak./noise);
end

for ch = 1:64
    tid = 17;
    waveforms_test{ch} = zeros(length(spike_times{tid}),60);
    for t = 1:length(spike_times{tid})
        waveforms_test{ch}(t,:) = filtData(ch, ...
            round(spike_times{tid}(t,1)-19):round(spike_times{tid}(t,1)+40,1)); %
    end
    waveforms_mean{ch} = mean(waveforms_test{ch},1);
end


% figure
% for ch = 1:64
% % for n = 1:100
% %     plot(waveforms_test{ch}(n,:))
% %     hold on
% % end
% subplot(16,4,(ch))
% plot(waveforms_mean{chanMap(ch)},'LineWidth',2)
% axis([0 60 -200 100])
% 
% hold on
% drawnow
% title(num2str(chanMap(ch)))
% end
%     


% figure
% bt = 1;
% for tid = [17,18,19,27,28,29]
%     subplot(3,2,bt)
%     
%     for nb = 1:100
%     plot(waveforms.raw{tid}(nb,:))
%     hold on
%     end
%     plot(waveforms.mean{tid},'Linewidth',2)
%     bt = bt+1;
%     title(num2str(channels(tid)))
% end

% SNR
% function [ SNR, meanSnippet, stdSnippet ] = AnalyzeSpikeSnippet( Snippets )
% %ANALYZESPIKESNIPPET Computes the SNR and average spike snippet.
% %   Snippets must be passed in as a m x n array where m is the number of
% %   snippets and n is the number of points in each snippet.

% SNR = zeros(1,length(Cchannels(:,1)));
% for tid = 1:length(Cchannels(:,1))  
%     Snippets = waveforms.raw{tid};
%     peak_to_peak = max(meanSnippet) - min(meanSnippet);
%     noise = mean(stdSnippet(1:10));
%     SNR(ch) = 20.*log10(peak_to_peak./noise);
% end
% Output.SNR = SNR.';



%%
if kwe ~=1
    file_type = '100';
    ch =1;
    rec_start_time = timestamps1{ch}(1);
    %[data1,timestamps1,info] = load_open_ephys_data_faster([fpath filesep file_type '_CH' num2str(ch) '.continuous' ]); 
    % 
    [data, timestamps, info] = load_open_ephys_data(fullfile(fpath, 'all_channels.events'));
    
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
nreps = 10;

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
spike_timesSU = {};
rate_stim = {};


for id = 1:length(SU)
    raster.stim{id} = [];
    raster.rep{id} = [];
    raster.spikes{id} = [];
    spikes_pooled{id} = [];  
    rate_stim{id} = [];
    spike_timesSU{id} = spike_times_all(find(spike_clusters == Cids(SU(id))))/fs;
    for rep = 1:TotalReps
        %for raster
        spikes1 = spike_timesSU{id}(find(spike_timesSU{id}>=start_stim_times(rep)-PreStim & ...
            spike_timesSU{id}<=end_stim_times(rep)+ PostStim)).';
        spikes1 = spikes1 - start_stim_times(rep);
        spikes_pooled{id} = [spikes_pooled{id} spikes1];
        raster.stim{id} = [raster.stim{id} data_new(rep,1)*ones(size(spikes1))];
        raster.rep{id} = [raster.rep{id} data_new(rep,2)*ones(size(spikes1))];
        raster.spikes{id} = [raster.spikes{id} spikes1];
        
        %for rate
        spikes2 = spike_timesSU{id}(find(spike_timesSU{id}>=start_stim_times(rep) & ...
            spike_timesSU{id}<=end_stim_times(rep))).';
        rate_stim{id}(data_new(rep,1),data_new(rep,2)) = length(spikes2);
        spikes3 = spike_timesSU{id}(find(spike_timesSU{id}<=start_stim_times(rep) & ...
            spike_timesSU{id}>=start_stim_times(rep)-StimDur)).' ;
        rate_pre{id}(data_new(rep,1),data_new(rep,2)) = length(spikes3);
    end
end


%extract MU activity

% rasterMU.stim = [];
% rasterMU.rep = [];
% rasterMU.spikes = [];
% spike_timesMU = [];
% rate_stimMU =[];
% for id = 1:length(MU)
%     spike_timesMU = spike_times_all(find(spike_clusters == Cids(MU(id))))/fs;
%     for rep = 1:TotalReps
%         spikes = spike_timesMU(find(spike_timesMU >= start_stim_times(rep)-PreStim & ...
%             spike_timesMU <= end_stim_times(rep)+ PostStim)).';
%         spikes = spikes - start_stim_times(rep);
%         rasterMU.stim = [rasterMU.stim data_new(rep,1)*ones(size(spikes))];
%         rasterMU.rep = [rasterMU.rep data_new(rep,2)*ones(size(spikes))];
%         rasterMU.spikes = [rasterMU.spikes spikes];
%         spikes2 = spike_timesMU(find(spike_timesMU >= start_stim_times(rep) & ...
%             spike_timesMU <= end_stim_times(rep))).';
%         rate_stimMU(data_new(rep,1),data_new(rep,2)) = length(spikes2);
%     end
% end

stim_label  = x.stimulus_ch1(:,8);

SUrate = {};
% figure
for id = 1:length(SU)
    figure
    suptitle(['cluster ' num2str(Cids(SU(id))) ' channel ' num2str(Cchannels(id)) ' SNR: ' num2str(waveforms.SNR{SU(id)})])
    set(gcf, 'Position', get(gcf,'Position').*[1 1 0 0] + [0 -600 1000 800]);
    
    %// this should change depending on stim_set. make it modular?
    subplot(2,2,1)
    SUrate{id}.mean = mean(rate_stim{id},2); %-mean(rate_pre{id},2);
    SUrate{id}.error = std(rate_stim{id},1,2);
    SUrate{id}.spont = mean(mean(rate_pre{id}));
    
    xs = 1:length(stim_label);
    h =1.0;
    for i = 1:length(stim_label)
        ys1(i) = gaussian_kern_reg(xs(i),xs,SUrate{id}.mean.',h);
        %         ys2(i) = gaussian_kern_reg(xs(i),xs,SUrate{id}.mean(2:2:50).',h);
    end
    %     subplot(3,3,id)
    errorbar(stim_label,ys1/StimDur,SUrate{id}.error,'LineWidth',2);
    hold on
    
    plot(stim_label,ones(1,length(stim_label))*SUrate{id}.spont/StimDur,'--k');
    %     plot(stim_label(2:2:end),ys2,'LineWidth',2);
    %     legend on
    %         axis([0 50 0 2])
    %         xticklabels({'8' '11' '14' '17' '20' '23' '26' '29' '32'})
    xlabel('Hz')
    ylabel('firing rate (spikes/s)')
    title('tuning curve')
    
    %waveform
    subplot(2,2,3)
    time_step = [1:60]/fs*1e3;
%     Nb_wave = min(size(waveforms.raw{find(Cchannels ==Cchannels(id))},1),100);
        Nb_wave = min(size(waveforms.raw{SU(id)},1),100);

    for n = 1:Nb_wave    
        plot(time_step,waveforms.raw{SU(id)}(n,:))
        hold on
    end
        plot(time_step,waveforms.mean{SU(id)},'Linewidth',2,'Color','k')
    title('waveform')
    xlabel('time (ms)')
    ylabel('mV');
hold off
    %Raster plot
    subplot(2,2,[2 4])
    area([0 StimDur StimDur 0],[0 TotalReps+5 0 TotalReps+5],'LineStyle','none','FaceColor',[.85 .85 1]);
    hold on
    plot(raster.spikes{id},nreps*(raster.stim{id}-1)+raster.rep{id},'k.','MarkerSize',9);
    %     pause
    xlabel('time (s)')
    ylabel('reps')
    axis([-PreStim StimDur+ PostStim 0 TotalReps+1])
    hold off
    title('rasterplot')
    
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
% for id = 1:length(SU)
%     subplot(2,4,id)
%     time_step = [1:60]/fs*1e3;
%     for n = 1:100    
%         plot(time_step,waveforms.raw{find(Cchannels ==Cchannels(id))}(n,:))
%         hold on
%     end
%         plot(time_step,waveforms.mean{find(Cchannels ==Cchannels(id))},'Linewidth',2,'Color','k')
%     title([' channel ' num2str(Cchannels(id))])
%     xlabel('time (ms)')
%     ylabel('mV');
% hold off
% end



% 
% 
% figure
% % ylabel('Repetition rate (Hz)')
% 
% 
% hold on
% figure
% title('Multi Units ')
% area([0 StimDur StimDur 0],[0 TotalReps+5 0 TotalReps+5],'LineStyle','none','FaceColor',[.85 .85 1]);
% hold on 
% plot(rasterMU.spikes,nreps*(rasterMU.stim-1)+rasterMU.rep,'k.','MarkerSize',9);
% xlabel('time (s)')
% ylabel('reps')
% axis([-PreStim StimDur+ PostStim 0 TotalReps+1])


