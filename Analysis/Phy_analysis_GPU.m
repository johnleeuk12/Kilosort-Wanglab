% rmpath(genpath('C:\Users\Seth\Documents\GitHub\KiloSort'))
% addpath(genpath('C:\Users\Seth\Documents\GitHub\KiloSort2')) % path to kilosort folder
% addpath(genpath('C:\Users\Seth\Documents\GitHub\npy-matlab')) % path to npy-matlab scripts
% addpath(genpath('C:\Users\Seth\Documents\GitHub\analysis-tools-master')) % path to npy-matlab scripts


PC_name = '426_John';
animal_name = 'M94W';
session_name = '2019-03-22_14-15-50';

fpath = directories(PC_name,animal_name,session_name);

fs = 30000;

addpath(genpath('C:\DATA\xbz'));



%%
x = M94W0239; %// add option to choose file
Nb_ch = 64;

kwe = 0; %// option for debugging, use when events are saved in .kwe format

%// Loading spike times and cluster information, output from Kilosort
spike_times_all = readNPY(fullfile(fpath, 'spike_times.npy'));
% spike_times_all = double(spike_times_all)/fs;
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



%%

data_all = [];
timestamps = {};
info = {};
file_type = '100';
parfor ch = 1:Nb_ch
    disp(['loading data from channel ' num2str(ch)])
    [data1,timestamps{ch},info{ch}] = load_open_ephys_data_faster([fpath filesep file_type '_CH' num2str(ch) '.continuous' ]);    
    data_all(ch,:) = data1;
end

for ch = 1:Nb_ch
    