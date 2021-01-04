function master_Harris2()
% This function is to make master_Harris.m compatible with Kilosort2.

%Currently Kilosort2 only supports GPU processing
useGPU = 1; %else 1  % do you have a GPU? Kilosorting 1000sec of 32chan simulated data takes 55 seconds on gtx 1080 + M2 SSD.


% session_name ='2020-01-08_13-50-00';
session_name = 'H1T3S2_concat';
Animal_name = 'M60F';
fpath    = fullfile('C:\DATA\OpenEphys', filesep, Animal_name, filesep, session_name); % where on disk do you want the simulation? ideally and SSD...
file_type = '100';

if ~exist(fpath, 'dir'); mkdir(fpath); end

rootZ = fpath;


rmpath(genpath('C:\Users\Seth\Documents\GitHub\KiloSort2'))

addpath(genpath('C:\Users\Seth\Documents\GitHub\KiloSort2_old')) % path to kilosort folder
addpath(genpath('C:\Users\Seth\Documents\GitHub\npy-matlab')) % path to npy-matlab scripts

pathToYourConfigFile = 'C:\Users\Seth\Documents\GitHub\Kilosort-Wanglab'; % take from Github folder and put it somewhere else (together with the master_file)
run(fullfile(pathToYourConfigFile, 'Harris_config2.m'))

ops.trange = [0 Inf]; % time range to sort
ops.ephys_type = file_type;


% edit for neuropixels (change preprocess sub too)
% rootZ = 'D:\Data\JHU_Johns_data';
% run(fullfile('C:\Users\Seth\Documents\GitHub\Kilosort2\configFiles\configFile384.m'));

% % rootH = 'C:\Data\OpenEphys\M12E\Test_concat';
% ops.fproc       = fullfile(rootZ, 'temp_wh.dat'); % proc file on a fast SSD
% % ops.chanMap = fullfile(pathToYourConfigFile, 'neuropixPhase3A_kilosortChanMap.mat');
% 
% 
% 
% 
% ops.NchanTOT    = 385;
% % ops.fbinary             = 'test_binary.dat';
% ops.root = rootZ;
% fs          = [dir(fullfile(rootZ, '*.bin')) dir(fullfile(rootZ, '*.dat'))];
% ops.fbinary = fullfile(rootZ, fs(1).name);
% %

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
    gpuDevice(1) % initialize GPU (will erase any existing GPU arrays)
end

% Openephys to binary converter
ops.ephys_type = file_type;

% disp('converting data...')
% tic
% % test = 1;
% % if strcmp(ops.datatype , 'openEphys')
% %    ops = convertOpenEphysToRawBInary(ops);  % convert data, only for OpenEphys
% % end
% disp('converting data... Done')
% toc

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