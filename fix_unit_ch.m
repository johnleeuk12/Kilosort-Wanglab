clear all
addpath(genpath(fullfile('D:\Data\M12E\Units')));
addpath('C:\Users\Seth\Documents\GitHub\Kilosort-Wanglab\Analysis')
list_data = load('M12E_unit_list.mat');

x.PC_name = '426_Analysis';
x.animal_name = 'M12E';
x.session_name = '2019-07-25_14-07-17';
x.file_type = '100';
x.xbz_file_name = 'M12E0431';
fpath = directories(x.PC_name,x.animal_name,x.session_name);
x.chanMap = [27 32 21 3 25 30 19 5 23 17 24 7 20 29 26 9 ...
    22 31 28 11 16 13 1 15 18 12 14 10 8 4 6 2 ...
    64 60 62 58 56 52 54 48 49 53 51 55 57 61 59 63 ...
    38 42 40 45 43 35 33 47 36 50 34 44 46 39 41 37];





save_dir = 'D:\Data\M12E\Units';
% spike_times_all = readNPY(fullfile(obj.params.fpath, 'spike_times.npy'));
% spike_times_all = double(spike_times_all);
spike_clusters = readNPY(fullfile(fpath, 'spike_clusters.npy'));
[Cids, Cgroups] = readClusterGroupsCSV(fullfile(fpath, 'cluster_group.tsv'));
spike_templates = readNPY(fullfile(fpath, 'spike_templates.npy'));
new_channels = readNPY(fullfile(fpath, 'channel_map.npy'));
new_channels = new_channels+1;
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

CchanMap = zeros(1,length(new_channels));
for c = 1:length(new_channels)
    CchanMap(c) = find(x.chanMap == new_channels(c));
end
CchanMap = CchanMap.';
Cchannels = Cchannels.';

for id = 268:306
    unit_file_name = 'M12Eu000';
    unit_file_name = [unit_file_name(1:end-size(num2str(id),2)) num2str(id) '.mat'];
    x = load(unit_file_name);
    x.s_unit.ch = CchanMap(x.s_unit.ch);
    s_unit = x.s_unit;
    save(fullfile(save_dir,unit_file_name),'s_unit')
end