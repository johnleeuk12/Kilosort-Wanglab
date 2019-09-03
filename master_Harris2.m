function master_Harris2()
% This function is to make master_Harris.m compatible with Kilosort2.

%Currently Kilosort2 only supports GPU processing
useGPU = 1; %else 1  % do you have a GPU? Kilosorting 1000sec of 32chan simulated data takes 55 seconds on gtx 1080 + M2 SSD.


% session_name ='2019-07-26_17-38-31';
session_name = 'H6T3S1_concat';
Animal_name = 'M12E';
fpath    = fullfile('C:\DATA\OpenEphys', filesep, Animal_name, filesep, session_name); % where on disk do you want the simulation? ideally and SSD...
file_type = '100';

if ~exist(fpath, 'dir'); mkdir(fpath); end

rootZ = fpath;


% rmpath(genpath('C:\Users\Seth\Documents\GitHub\KiloSort'))
addpath(genpath('C:\Users\Seth\Documents\GitHub\KiloSort2')) % path to kilosort folder
addpath(genpath('C:\Users\Seth\Documents\GitHub\npy-matlab')) % path to npy-matlab scripts

pathToYourConfigFile = 'C:\Users\Seth\Documents\GitHub\Kilosort-Wanglab'; % take from Github folder and put it somewhere else (together with the master_file)
run(fullfile(pathToYourConfigFile, 'Harris_config2.m'))

ops.trange = [0 Inf]; % time range to sort
ops.ephys_type = file_type;

% % concatanate file
% list = {'2019-07-27_14-24-54','2019-07-27_14-26-39', '2019-07-27_14-28-36', '2019-07-27_14-30-40','2019-07-27_14-33-57','2019-07-27_14-38-11',   ...
% '2019-07-27_14-40-28', '2019-07-27_14-50-19', '2019-07-27_14-59-45', '2019-07-27_15-07-30','2019-07-27_15-10-49','2019-07-27_15-18-16', ...
% '2019-07-27_15-21-46','2019-07-27_15-26-37','2019-07-27_15-30-37','2019-07-27_15-39-42'};
% 
% list= {'2019-07-26_16-31-51'
% '2019-07-26_16-34-50'
% '2019-07-26_16-42-16'
% '2019-07-26_16-47-30'
% '2019-07-26_16-56-03'
% '2019-07-26_17-05-34'
% '2019-07-26_17-08-48'
% '2019-07-26_17-15-29'
% '2019-07-26_17-21-53'
% '2019-07-26_17-26-29'
% '2019-07-26_17-31-04'
% '2019-07-26_17-32-54'
% '2019-07-26_17-34-41'
% '2019-07-26_17-36-30'
% '2019-07-26_17-38-31'
% };
% 
%  
% fname = 'test_binary.dat';
% filepath    = fullfile('C:\DATA\OpenEphys', filesep, Animal_name, filesep, list, filesep, fname);
% buff= {};
% % Big_buff = [];
% tic
% fileID = fopen(fullfile(fpath,filesep,fname),'w');
% 
% for t = 1:length(list)
%     
%     fid         = fopen(filepath{t}, 'r');
%     buff{t} = fread(fid, '*int16');
%     fclose(fid);
%     disp(t)
%     fwrite(fileID,buff{t},'*int16');
%     
%     toc
% end
% fclose(fileID);





% k = kilosort
%% this block runs all the steps of the algorithm
fprintf('Looking for data inside %s \n', rootZ)

% is there a channel map file in this folder?
fs = dir(fullfile(rootZ, 'chan*.mat'));
if ~isempty(fs)
    ops.chanMap = fullfile(rootZ, fs(1).name);
else
    make_HarrisChannelMap(fpath)
end



tic; % start timer
%
if ops.GPU     
    gpuDevice(1); % initialize GPU (will erase any existing GPU arrays)
end

% Openephys to binary converter
ops.ephys_type = file_type;

disp('converting data...')
tic
% test = 1;
if strcmp(ops.datatype , 'openEphys')
   ops = convertOpenEphysToRawBInary(ops);  % convert data, only for OpenEphys
end
disp('converting data... Done')
toc

% preprocess data to create temp_wh.dat
rez = preprocessDataSub(ops);

% time-reordering as a function of drift
rez = clusterSingleBatches(rez);
save(fullfile(rootZ, 'rez.mat'), 'rez', '-v7.3');

% main tracking and template matching algorithm
rez = learnAndSolve8b(rez);

% final merges
rez = find_merges(rez, 1);

% final splits by SVD
rez = splitAllClusters(rez, 1);

% final splits by amplitudes
rez = splitAllClusters(rez, 0);

% decide on cutoff
rez = set_cutoff(rez);

fprintf('found %d good units \n', sum(rez.good>0))

% write to Phy
fprintf('Saving results to Phy  \n')
rezToPhy(rez, rootZ);