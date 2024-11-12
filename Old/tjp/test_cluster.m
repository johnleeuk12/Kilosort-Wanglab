


list  = {'2019-07-25_13-42-03','2019-07-25_13-50-02','2019-07-25_13-52-32','2019-07-25_13-57-28','2019-07-25_14-07-17', ...
    '2019-07-25_14-16-18','2019-07-25_14-25-19','2019-07-25_14-27-22'};

addpath('C:\Users\Seth\Documents\GitHub\Kilosort-Wanglab\Analysis')
x.PC_name = '426_Analysis';
x.animal_name = 'M12E';
x.session_name = list{1};
x.file_type = '100';
x.xbz_file_name = 'M12E0428';
x.fpath = directories(x.PC_name,x.animal_name,x.session_name);
x.fs = 30000;
x.Nb_ch = 64;
x.chanMap = [27 32 21 3 25 30 19 5 23 17 24 7 20 29 26 9 ...
    22 31 28 11 16 13 1 15 18 12 14 10 8 4 6 2 ...
    64 60 62 58 56 52 54 48 49 53 51 55 57 61 59 63 ...
    38 42 40 45 43 35 33 47 36 50 34 44 46 39 41 37];

x.conData = 'C:\Data\OpenEphys\M12E\Test_concat\test_binary.dat';
x.fpath2 = 'C:\Data\OpenEphys\M12E\Test_concat\';




%%
[data1,timestamps1{1},info1{1}] = load_open_ephys_data_faster([x.fpath filesep x.file_type '_CH' num2str(1) '.continuous' ]);

test = (timestamps1{1}(2)-timestamps1{1}(1))*30000;
% length(data1)-test;
[data, timestamps, info] = load_open_ephys_data(fullfile(x.fpath, 'all_channels.events'));


% fid = x.conData;
fid = fopen(x.conData,'r');
Nbuff = 1024*32;

Nbatch = length(data1)/(Nbuff);
Nbatch = floor(Nbatch);
DATA = single([]);
tic
[b1,a1] = butter(4, [0.0244 0.6104],'bandpass');
% Nbuff = 1024*32;
for batch = 1:Nbatch
    ipoint = (batch-1)*Nbuff +1;
%     disp(batch);
    buff = fread(fid,[64,Nbuff], '*int16');
    
    dataRAW = gpuArray(buff);
    dataRAW = dataRAW';
    dataRAW = single(dataRAW);
    dataRAW = dataRAW(:, x.chanMap);
    
    % subtract the mean from each channel
    dataRAW = dataRAW - mean(dataRAW, 1);
    
    datr = filter(b1, a1, dataRAW);
    datr = flipud(datr);
    datr = filter(b1, a1, datr);
    datr = flipud(datr);
    datr = datr - median(datr, 2);
    
    DATA(ipoint:ipoint+Nbuff-1,:)= gather_try(int16(datr));
    
end

toc
fclose(fid);

seg_1 = 12737536;

DATA = double(DATA);

% addpath(fullfile('C:\Data\OpenEphys', filesep, obj.params.animal_name));
spike_times_all = readNPY(fullfile(x.fpath2, 'spike_times.npy'));
spike_clusters = readNPY(fullfile(x.fpath2, 'spike_clusters.npy'));
[Cids, Cgroups] = readClusterGroupsCSV(fullfile(x.fpath2, 'cluster_group.tsv'));

SU = find(Cgroups == 2);
id = 1;
spike_timesSU{id} = spike_times_all(find(spike_clusters == Cids(SU(id)))); %/obj.params.fs;
waveforms{id} = zeros(64,60,length(spike_timesSU{id}));

temp_spike_times = spike_timesSU{id}(find(spike_timesSU{id} < length(filtData)+seg_1 & seg_1< spike_timesSU{id}));
temp_spike_times = temp_spike_times - seg_1;
for ch =1 : 64
    for s = 1:length(temp_spike_times)
        waveforms{id}(ch,:,s) = DATA(temp_spike_times(s,1)-19:temp_spike_times(s,1)+40,ch);
    end
end

figure
for ch = 1:2:64
    subplot(8,8,ch)
    for p = 1:100
        plot(waveforms{id}(ch,:,p))
        hold on
    end
    drawnow
end
tid = 1;
obj.waveforms.raw = {};
%             for tid = 1:length(Cids)
%                 id = Cids(tid);
%                 spike_times{tid} = spike_times_all(find(spike_clusters == id)); %*fs;
                obj.waveforms.raw{tid} = zeros(length(temp_spike_times),64,60);
%                 for t = 1:length(spike_timesSU{tid})
                    for t = 1:length(temp_spike_times)
                    obj.waveforms.raw{tid}(t,:,:) = filtData(:, ...
                        round(temp_spike_times(t,1)-19):round(temp_spike_times(t,1)+40)); %
                end
                obj.waveforms.mean{tid} = mean(obj.waveforms.raw{tid},1);
                obj.waveforms.std{tid} = std(obj.waveforms.raw{tid},[],1);
                peak_to_peak = max(obj.waveforms.mean{tid})-min(obj.waveforms.mean{tid});
                noise = mean(obj.waveforms.std{tid}(1:10));
                obj.waveforms.SNR(tid) = 20.*log10(peak_to_peak./noise);
            end



temp = []
for ch = 1:64
    temp(1,:) = obj.waveforms.mean{1}(1,ch,:);
    plot(temp)
    hold on
end


