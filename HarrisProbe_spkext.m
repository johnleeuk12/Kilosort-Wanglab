function HarrisProbe_spkext()

close all




% doc McsHDF5
% time unit = 10^-6 seconds 
% Voltage unit = 10^-9 V
% First trigger is trigger of recording ( needs to be verified)
addpath(genpath('C:\Users\Seth\Documents\GitHub\copy\KiloSort')) % path to kilosort folder
addpath(genpath('C:\Users\Seth\Documents\GitHub\npy-matlab')) % path to npy-matlab scripts
addpath(genpath('C:\Users\Seth\Documents\GitHub\analysis-tools-master'))
animal = 'M94W';
% filenb = '2018-12-05_13-11-28';
session_name = '2019-04-05_14-05-45';
filepath = ['C:\DATA\OpenEphys' filesep animal filesep session_name];
file_type = '100';
if ~exist(filepath, 'dir'); mkdir(filepath); end
% filepath2 = ':\Data\Experiments\M44D'; %file path for .m files
addpath(genpath(filepath))
% data1 = McsHDF5.McsData([filepath filesep animal filenb '.h5']);

% for ch = 1:64
% data{ch} = {};
% timestamps{ch} = {};
% info{ch} = {};
% end


% file_type = '116';

parfor ch = 1:64    
    [data1,timestamps{ch},info{ch}] = load_open_ephys_data_faster([filepath filesep file_type '_CH' num2str(ch) '.continuous' ]);   
    datach{ch} = data1;
    disp(ch)
end

% All_Ch.data = datacch;
% All_Ch.timestamps = timestamps;
% All_Ch.info = info;


% butterworth bandpass filter
[b,a] = butter(4, [0.0244 0.6104]);

total_sample = length(datach{1});
filtData = zeros(64,total_sample);
parfor ch = 1:64
    disp(['filtering CH_' num2str(ch)])
    filtData(ch,:) = filtfilt(b,a,datach{ch});
end

%Common median referencing

CommonMedian = median(filtData);
st_dev = zeros(1,64);
parfor ch = 1:64
    filtData(ch,:) = filtData(ch,:)-CommonMedian;
    st_dev(ch) = median(abs(filtData(ch,:))/0.6745);
    disp(ch)
end

threshold = st_dev*-5; %modifiable. 

%% wavelet extraction

disp('extracting wavelets...')
Output.wavelets = {};
Output.spiketime = {};
for ch = 1:64
    tic
    disp(ch)
    [peak,Output.spiketime{ch}] = findpeaks(-filtData(ch,:),'MinPeakHeight',-threshold(ch));
    Output.spiketime{ch}(find(peak>1e8)) =[];
    test = 1;
    Output.spiketime{ch}(find(Output.spiketime{ch}>length(CommonMedian)-50)) = [];
    Output.spiketime{ch}(find(Output.spiketime{ch}<30)) = [];
    Output.wavelets{ch} = zeros(length(Output.spiketime{ch}),61);
%     spikenb1= 1;
%     spikenb2 = 1;
    for spikenb = 1:length(Output.spiketime{ch}) 
        Output.wavelets{ch}(spikenb,:) = filtData(ch,Output.spiketime{ch}(spikenb)-19:Output.spiketime{ch}(spikenb)+41);
%         if max(abs(Output.wavelets{ch}(spikenb2,:)))<2*1e8
%             spikenb2 = spikenb2+1;
%         end
%         spikenb1 = spikenb1+1;
    end
    disp(['total number of spikes = ' num2str(spikenb)])
    toc
end



% function [ SNR, meanSnippet, stdSnippet ] = AnalyzeSpikeSnippet( Snippets )
% %ANALYZESPIKESNIPPET Computes the SNR and average spike snippet.
% %   Snippets must be passed in as a m x n array where m is the number of
% %   snippets and n is the number of points in each snippet.
% 
SNR = zeros(1,64);
for ch = 1:64    
    Snippets = Output.wavelets{ch};
    meanSnippet = mean(Snippets,1);
    stdSnippet = std(Snippets,[],1);
    peak_to_peak = max(meanSnippet) - min(meanSnippet);
    noise = mean(stdSnippet(1:10));
    SNR(ch) = 20.*log10(peak_to_peak./noise);
end
Output.SNR = SNR.';


% channelmap = [5 19 29 3 20 64 17 6 62 18 61 37 54 44 53 41; ...
%     14 12 24 16 28 27 1 30 58 33 34 45 50 40 57 49; ...
%     7 9 22 10 26 32 31 8 36 63 46 59 42 56 43 51; ...
%     23 13 15 25 11 21 2 4 60 47 48 35 38 52 55 39];
% 
% channelmap = [23 7 14 5 13 9 12 19 15 22 24 29 25 10 16 3 11 26 ...
%     28 20 21 32 27 64 2 31 1 17 4 8 30 6 60 36 58 62 47 63 33 18 ...
%     48 46 34 61 35 59 45 37 38 42 50 54 52 56 40 44 55 43 57 53 ...
%     39 51 39 31];

% channelmap = [30 20 29 32 31 12 5 8 25 21 22 19 16 15 4 9 27 23 17 18 13 3 7 6 26 24 28 2 1 14 11 10];

% channelmap = [28 26 25 19 20 33 18 41 49 57 64 56 5 10 3 13 ...    
%     31 29 22 30 32 36 37 45 53 61 60 52 15 4 63 8 ...
%     24 23 17 21 27 48 35 43 51 59 62 54 11 2 1 7 ...
%     46 44 42 40 38 34 39 47 55 16 58 50 9 12 14 6];

chanMap = [23 8 14 6 14 10 11 19 16 22 24 29 25 9 15 4 12 26 ...
    28 20 21 32 27 63 1 31 2 17 3 7 30 5 59 36 57 61 ...
    47 64 33 18 48 46 34 62 35 60 45 37 38 42 49 53 ...
    51 55 40 44 56 43 58 54 39 52 50 41]; 


figureon = 0;

if figureon == 1
    figure('units','normalized','outerposition',[0 0 1 1])
    
    for ch = 1:64
%         if ismember(ch,[1 3:9 11:31])
            subplot(4,16,find(channelmap == ch))
            if size(Output.wavelets{1,ch},2) < 40
                for ind = 1:size(Output.wavelets{1,ch})
                    plot(Output.wavelets{ch}(ind,:))
                    axis([-inf inf -100 100])   
                    hold on
                end
            else
                for ind = randsample([1:size(Output.wavelets{1,ch},2)],40)
                    plot(Output.wavelets{ch}(ind,:))
                    hold on
                end
            end
            plot(mean(Output.wavelets{ch}(:,:)),'-k','Linewidth',3)
            axis([-inf inf -100 100]) 
            title(num2str(ch))
            drawnow
%         end
    end
    
    saveas(gcf,[animal filenb '_waveform3.png'])
end




All_wavelets = [];
for ch = 1:64
    All_wavelets = [All_wavelets ; Output.wavelets{ch}];
end

[Wi,score,latent] = pca(All_wavelets);
% 

disp('Saving files...')
% save([animal filenb '_proc.mat'],'Output','-v7.3')

save(fullfile(filepath, '\PCspikes_94W.mat'),'Wi')

test

% [coeff,score,latent] = pca(Output.wavelets);
% %% preprocessing data
% disp('Pre processing data...')
% % Output = {};
% 
% % butterworth bandpass filter
% 
% 
% figureon = 0; % trigger for figures
% %% extract data information
% infoname = [animal filenb];
% dataInfo = eval(infoname);
% % dataInfo = M94W0194;
% nbreps = dataInfo.stimulus_ch1(1,4);
% StimDur = dataInfo.stimulus_ch1(1,5);
% 
% 
% 
% stim_index = dataInfo.data(find(dataInfo.data(:,3) == 1 & dataInfo.data(:,4) == -1),1);
% stim_nb = dataInfo.stimulus_ch1(:,1);
% stim_name = num2str(dataInfo.stimulus_ch1(:,10));
% 
% rep_index = dataInfo.data(find(dataInfo.data(:,3) == 1 & dataInfo.data(:,4) == -1),2);
% 
% if size(data1.Recording{1, 1}.EventStream,2) == 1
%     event_times = data1.Recording{1, 1}.EventStream{1, 1}.Events{1,1}(1,:).';
% else
%     event_times = data1.Recording{1, 1}.EventStream{1, 2}.Events{1,1}(1,:).';
% end
% 
% start_trigger = event_times(1);
% event_times = event_times(2:end);
% timeStamp = data1.Recording{1, 1}.AnalogStream{1, 2}.ChannelDataTimeStamps;
% 
% %% preprocessing data
% disp('Pre processing data...')
% Output = {};
% % butterworth bandpass filter
% Output.rawData = data1.Recording{1, 1}.AnalogStream{1, 2}.ChannelData;
% [b,a] = butter(4, [0.0244 0.6104]);
% Output.filtData = zeros(32, length(data1.Recording{1, 1}.AnalogStream{1, 2}.ChannelData));
% 
% for ch = 1:32
%     Output.filtData(ch,:) = filtfilt(b,a,data1.Recording{1, 1}.AnalogStream{1, 2}.ChannelData(ch,:));
% end
% 
% % remove data until first trigger
% 
% Output.filtData = Output.filtData(:,find(timeStamp == start_trigger):end);
% Output.timeStampNew = timeStamp(find(timeStamp == start_trigger):end);
% 
% %Common median referencing
% 
% CommonMedian = median(Output.filtData);
% st_dev = zeros(1,32);
% for ch = 1:32
%     Output.filtData(ch,:) = Output.filtData(ch,:)-CommonMedian;
%     st_dev(ch) = median(abs(Output.filtData(ch,:))/0.6745);
% end
% 
% 
% % st_dev = median(abs(Output.filtData),0,2);
% Output.threshold = st_dev*-5; %modifiable. 
% % time between events should be 600ms, but with sampling more often than
% % not becomes 600.05ms. This accumulates. 
% % So, we take the event start time, and take 600ms of data. 
% 
% nb_rep = max(unique(rep_index));
% for ch = 1:32
%     disp(ch)
%     Output.bytrial{ch} = zeros(length(stim_index),round(600/0.05));
%     for ind = 1:length(stim_index)
%         timeInd = find(Output.timeStampNew == event_times(ind));
%         Output.bytrial{ch}((stim_index(ind)-1)*nb_rep+rep_index(ind),:) = Output.filtData(ch,timeInd:timeInd+round(600/0.05)-1);
%     end
% end
% 
% % test = std(Output.bytrial{4}*10^-9);
% 
% 
% % plotting raw data
% figure
% for ch = [1 3:9 11:31]
%     plot(Output.bytrial{ch}(1:100000)*10^-9 + ch*0.1)
%     hold on 
% end
% drawnow
%     
% % for ind = 1:380
% %     test(ind,:) = test(ind,:) + ind;
% % end
% % 
% % figure
% % 
% % for ind = 1:380
% %     plot(test(ind,:))
% %     hold on
% % end
% 
% %% wavelet extraction
% 
% disp('extracting wavelets...')
% Output.wavelets = {};
% Output.spiketime = {};
% for ch = [1 3:9 11:31]
%     disp(ch)
%     [peak,Output.spiketime{ch}] = findpeaks(-Output.filtData(ch,:),'MinPeakHeight',-Output.threshold(ch));
%     Output.spiketime{ch}(find(peak>1e8)) =[];
%     test = 1;
%     Output.spiketime{ch}(find(Output.spiketime{ch}>length(CommonMedian)-50)) = [];
%     Output.spiketime{ch}(find(Output.spiketime{ch}<30)) = [];
%     Output.wavelets{ch} = zeros(length(Output.spiketime{ch}),60);
% %     spikenb1= 1;
% %     spikenb2 = 1;
%     for spikenb = 1:length(Output.spiketime{ch}) 
%         Output.wavelets{ch}(spikenb,:) = Output.filtData(ch,Output.spiketime{ch}(spikenb)-19:Output.spiketime{ch}(spikenb)+40);
% %         if max(abs(Output.wavelets{ch}(spikenb2,:)))<2*1e8
% %             spikenb2 = spikenb2+1;
% %         end
% %         spikenb1 = spikenb1+1;
%     end
% end
% 
% % function [ SNR, meanSnippet, stdSnippet ] = AnalyzeSpikeSnippet( Snippets )
% % %ANALYZESPIKESNIPPET Computes the SNR and average spike snippet.
% % %   Snippets must be passed in as a m x n array where m is the number of
% % %   snippets and n is the number of points in each snippet.
% % 
% SNR = zeros(1,32);
% for ch = [1 3:9 11:31]
%     
%     Snippets = Output.wavelets{ch};
%     meanSnippet = mean(Snippets,1);
%     stdSnippet = std(Snippets,[],1);
%     peak_to_peak = max(meanSnippet) - min(meanSnippet);
%     noise = mean(stdSnippet(1:10));
%     SNR(ch) = 20.*log10(peak_to_peak./noise);
% end
% Output.SNR = SNR.';
% 
% 
% 
% 
% channelmap = [30 20 29 32 31 12 5 8 25 21 22 19 16 15 4 9 27 23 17 18 13 3 7 6 26 24 28 2 1 14 11 10];
% 
% if figureon == 1
%     figure('units','normalized','outerposition',[0 0 1 1])
%     
%     for ch = 1:32
%         if ismember(ch,[1 3:9 11:31])
%             subplot(4,8,find(channelmap == ch))
%             if size(Output.wavelets{1,ch},2)<1000
%                 for ind = 1:size(Output.wavelets{1,ch})
%                     plot(Output.wavelets{ch}(ind,:))
%                     hold on
%                 end
%             else
%                 for ind = randsample([1:size(Output.wavelets{1,ch},2)],1000)
%                     plot(Output.wavelets{ch}(ind,:))
%                     hold on
%                 end
%             end
%             plot(mean(Output.wavelets{ch}(:,:)),'-k','Linewidth',3)
%             title(num2str(ch))
%             drawnow
%         end
%     end
%     
%     saveas(gcf,[animal filenb '_waveform.png'])
% end
% 
% %% rasterplot
% %     figure ('units','normalized','outerposition',[0 0 1 1])
% %     axis([0 12000 0 382])
% %     rectangle('Position',[4000 0 2000 381],'FaceColor',[.9 .9 .9],'LineStyle','None')
% %     if isempty(Output.spiketime2{ch,trial},0)
% %     plot(Output.spiketime2{ch,trial},
% 
% 
% 
% 
% 
% warning('off','all')
% trialmax= max(stim_nb)*max(unique(rep_index));
% for ch = [1 3:9 11:31]
%     for trial = 1:trialmax
%         [peak,Output.spiketime2{ch,trial}] = findpeaks(-Output.bytrial{ch}(trial,:),'MinPeakHeight',-Output.threshold(ch));
%         Output.spiketime{ch,trial}(find(peak>1e8)) =[];
%     end
%     test= 1;
% end
% warning('on','all')
% 
% edges = [0:200:12000];
% 
% for ch = [1 3:9 11:31]
%     for stim = 1:length(stim_nb)
%         [Output.density{ch,stim},edges] = histcounts(horzcat(Output.spiketime2{ch,(stim-1)*nb_rep+1:(stim)*nb_rep}),edges);
%     end 
%     Output.heatraster{ch} = vertcat(Output.density{ch,:});
% end
% 
% if figureon ==1
%     figure('units','normalized','outerposition',[0 0 1 1])
%     for ch = 1:32
%         if ismember(ch,[1 3:9 11:31])
%             subplot(4,8,find(channelmap == ch))
%             
%             imagesc(Output.heatraster{ch})
%             colormap(flipud(gray))
%             caxis([0 15])
%             rectangle('Position',[20 0 30 max(stim_nb)+1],'FaceColor','None','Linewidth',2)
%             title(num2str(ch))
%             yticks([1:length(stim_nb)]);
%             yticklabels(stim_name);
%             drawnow
%         end
%     end
%     saveas(gcf,[animal filenb '_raster.png'])
% end
% 
% save([animal filenb '_proc.mat'],'Output','-v7.3')

