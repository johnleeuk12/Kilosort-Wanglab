function out = Spike_Analysis_app(fpath,x,raw_id, spike_id)

addpath(genpath('C:\Users\skoehler\Documents\GitHub\KiloSort')) % path to kilosort folder
addpath(genpath('C:\Users\skoehler\Documents\GitHub\npy-matlab')) % path to npy-matlab scripts
addpath('D:\Data\Experiments\M44D');
addpath(genpath('C:\Users\skoehler\Documents\GitHub\analysis-tools'))
% fpath    = 'D:\Data\Experiments\M44D\2019-01-10_13-56-05'; % where on disk do you want the simulation? ideally and SSD...
% fpath = 'D:\Data\Experiments\M94W\2019-03-21_17-06-34';
%'14-53-43'
fs = 30000;
N = 64;
% Alldata = {};
spikes = {};
timestamps = {};
info = {};

out = {}; %main output



% %test variables
% fpath=  'D:\Data\Experiments\M44D\2019-01-10_13-56-05';
% x  = eval('M44D1083');
% spike_id = '133';
% raw_id = '126';

% %test variables
% fpath=  'D:\Data\Experiments\M94W\2019-04-02_14-00-11';
addpath('C:\Data\Experiments\M94W');
% 
% x  = eval('M94W0020');
% spike_id = '109';
% raw_id = '100';

%% extracting spike times from openephys spike detection
% parfor i = 1:N
for i = 1:N
    if isfile(fullfile(fpath, sprintf(['SEp' spike_id '.0n' '%d.spikes'],i-1)))
        [data1, timestamps3, info1] = load_open_ephys_data_faster(fullfile(fpath,sprintf(['SEp' spike_id '.0n' '%d.spikes'],i-1)));
        mb_not_noise = [];
        for n = 1:size(data1,1)
            if min(data1(n,:))<-50 && max(data1(n,:))<1000 && min(data1(n,:))>-1500
                mb_not_noise = [mb_not_noise n];
            end
        end
        spikes{i} = data1(mb_not_noise,:);
        %         spikes{i} = data1;
        timestamps{i} = timestamps3(mb_not_noise,:);
        info{i} = info1;
        
    else
        spikes{i} = [];
        timestamps{i} = [];
        info{i} = [];
    end
end



%% extracting event times from kwd (should be replaced with openephys)
% the only info necessary is relative start-time of each events 
% use for testing and troubleshooting

% filename = 'experiment1.kwe';
% rec_start_time = 211384/fs;
% %     rec_start_time = timestamps{ch}(1);
% 
% event_times = h5read('experiment1.kwe','/event_types/TTL/events/time_samples');
% event_times = event_times/fs-rec_start_time;
% start_stim_times = double(event_times(1:2:end));
% end_stim_times = double(event_times(2:2:end));



%% extract event times with Openephys format
    file_type = raw_id;
    ch =1;
    [data1,timestamps1,info] = load_open_ephys_data_faster([fpath filesep file_type '_CH' num2str(ch) '.continuous' ]);
    rec_start_time = timestamps1(1);

    [data, timestamps2, info] = load_open_ephys_data(fullfile(fpath, 'all_channels.events'));

    start_stim_times = timestamps2(find(info.eventId ==1))-rec_start_time;
    end_stim_times = timestamps2(find(info.eventId ==0))-rec_start_time;

%% extracting stim and trial info

% x = M44D1083;
stim_info = x.data(find(x.data(:,3) == 1 & x.data(:,4) == -1),:);
nreps = x.stimulus_ch1(1,4);
nStim = max(x.stimulus_ch1(:,1));
TotalReps = nStim*nreps;
StimDur = x.stimulus_ch1(:,5)*1e-3;
PreStim = x.pre_stimulus_record_time*1e-3;
PostStim = x.post_stimulus_record_time*1e-3;

out.PreStim = PreStim;
out.PostStim = PostStim;
out.TotalReps = TotalReps;
out.StimDur = StimDur;
out.nreps = nreps;
out.nStim = nStim;

%% extracting spikes

raster.stim = {};
raster.rep = {};
raster.spikes = {};
% spikes_pooled = {};
% spike_times = {};
rate_stim = {};
rate_pre = {};
%
% % timestamps{52} = timestamps{52}(test);

for id = 1:N
    raster.stim{id} = [];
    raster.rep{id} = [];
    raster.spikes{id} = [];
    rate_stim{id} = [];
    rate_pre{id} = [];
    timestamps{id} = timestamps{id}-rec_start_time;
    for rep  = 1:TotalReps
        if length(timestamps{id}) > 10

            spikes1 = timestamps{id}(find(timestamps{id}>=start_stim_times(rep)-PreStim & ...
                timestamps{id}<=end_stim_times(rep)+ PostStim)).';
            spikes1 = spikes1 - start_stim_times(rep);

            raster.spikes{id} = [raster.spikes{id} spikes1];
            raster.stim{id} = [raster.stim{id} stim_info(rep,1)*ones(size(spikes1))];
            raster.rep{id} = [raster.rep{id} stim_info(rep,2)*ones(size(spikes1))];

            spikes_stim = timestamps{id}(find(timestamps{id}>=start_stim_times(rep) & ...
                timestamps{id}<=end_stim_times(rep))).';
            spikes_pre = timestamps{id}(find(timestamps{id}<=start_stim_times(rep) & ...
                timestamps{id}>=start_stim_times(rep)-PreStim)).';

            rate_stim{id}(stim_info(rep,1),stim_info(rep,2)) = length(spikes_stim);
            rate_pre{id}(stim_info(rep,1),stim_info(rep,2)) = length(spikes_pre);
        end
    end
end


out.raster = raster;
out.spikes = spikes;
out.rate.rate_stim = rate_stim;
out.rate.rate_pre = rate_pre;

%
%% Figures
% % raster plot
% %
% figure
% nn =0;
% % ylabel('Repetition rate (Hz)')
% for id = 1:N
%     if ~isempty(raster.stim{id}) % I'm only plotting channels that contain spikes during stim presentation
%         %     else
%         nn = nn+1;
%         subplot(3,3,nn)
%         hold on
%         title(['cluster nb ' num2str(id)])
%         area([0 StimDur StimDur 0],[0 TotalReps+5 0 TotalReps+5],'LineStyle','none','FaceColor',[.85 .85 1]);
%         plot(raster.spikes{id},nreps*(raster.stim{id}-1)+raster.rep{id},'k.','MarkerSize',9);
%         xlabel('time (s)')
%         ylabel('reps')
%         axis([-PreStim StimDur+ PostStim 0 TotalReps+1])
%     end
% end
%
% %plotting waveforms
% figure
% for ch = 1:63
%     subplot(8,8,ch)
%     title(['cluster nb ' num2str(ch)])
%     for n = 1:min([100,size(spikes{ch},1)])
%         plot(spikes{ch}(n,:))
%         hold on
%     end
%     drawnow
% end
%
%
% stim_tags = x.stimulus_ch1(:,8);
% plot PSTH
% figure
% nn = 0;
% for ch = 1:63
%     if ~isempty(rate_stim{ch})
%         nn = nn+1;
%         subplot(3,3,nn)
%         xs = 1:length(stim_tags);
%         h =1.0;
%         for i = 1:length(stim_tags)
%             ys(i) = gaussian_kern_reg(xs(i),xs,((mean(rate_stim{ch},2)-mean(rate_pre{ch},2))*5).',h);
%         end
%         plot(log(stim_tags),ys,'LineWidth',2)
%         title([num2str(ch)])
% 
%     end
% end
% 
% 
% 
% figure
% for n = 1: length(test)
%     nn= test(n);
%     plot(spikes{52}(nn,:))
%     hold on
% end
%
%
%
%
% figure
% % % PCA and clustering
% [Wi,score,latent] = pca(spikes{52} );
% scatter(score(:,1),score(:,2));
%
% X = [score(:,1) score(:,2)];
% idx = kmeans(X,3);
%
%
% figure
% gscatter(score(:,1),score(:,2),idx)
% %
%
% test = find(idx ==2);













