function master_Harris3()
% 11/24/2020
% This function is to make master_Harris.m compatible with Kilosort2.5

%Currently Kilosort2 only supports GPU processing
useGPU = 1; %else 1  % do you have a GPU? Kilosorting 1000sec of 32chan simulated data takes 55 seconds on gtx 1080 + M2 SSD.

% session_name ='2020-01-08_13-50-00';
session_name = 'H2T3S1_concat';
concat_OE(session_name);

Animal_name = 'M160E';
fpath    = fullfile('C:\DATA\OpenEphys', filesep, Animal_name, filesep, session_name); % where on disk do you want the simulation? ideally and SSD...
file_type = '100';

if ~exist(fpath, 'dir'); mkdir(fpath); end

rootZ = fpath;


% rmpath(genpath('C:\Users\Seth\Documents\GitHub\KiloSort'))
addpath(genpath('C:\Users\Seth\Documents\GitHub\KiloSort2')) % path to kilosort folder
% rmpath('C:\Users\Seth\Documents\GitHub\KiloSort2');
% addpath(genpath('C:\Users\Seth\Documents\GitHub\KiloSort-main'))
addpath(genpath('C:\Users\Seth\Documents\GitHub\npy-matlab')) % path to npy-matlab scripts

pathToYourConfigFile = 'C:\Users\Seth\Documents\GitHub\Kilosort-Wanglab'; % take from Github folder and put it somewhere else (together with the master_file)
run(fullfile(pathToYourConfigFile, 'Harris_config2.m'))

ops.trange = [0 Inf]; % time range to sort
ops.ephys_type = file_type;


%% this block runs all the steps of the algorithm
fprintf('Looking for data inside %s \n', rootZ)

% main parameter changes from Kilosort2 to v2.5
ops.sig        = 20;  % spatial smoothness constant for registration
ops.fshigh     = 300; % high-pass more aggresively
ops.nblocks    = 5; % blocks for registration. 0 turns it off, 1 does rigid registration. Replaces "datashift" option. 

% main parameter changes from Kilosort2.5 to v3.0
% ops.Th       = [15 10];


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





% find the binary file
% fs          = [dir(fullfile(rootZ, '*.bin')) dir(fullfile(rootZ, '*.dat'))];
% ops.fbinary = fullfile(rootZ, fs(1).name);

% preprocess data to create temp_wh.dat
rez = preprocessDataSub(ops);
%
% NEW STEP TO DO DATA REGISTRATION
rez = datashift2(rez, 1); % last input is for shifting data

% ORDER OF BATCHES IS NOW RANDOM, controlled by random number generator
iseed = 1;
                 
% main tracking and template matching algorithm
rez = learnAndSolve8b(rez, iseed);

% OPTIONAL: remove double-counted spikes - solves issue in which individual spikes are assigned to multiple templates.
% See issue 29: https://github.com/MouseLand/Kilosort/issues/29
rez = remove_ks2_duplicate_spikes(rez);

% final merges
rez = find_merges(rez, 1);

% final splits by SVD
rez = splitAllClusters(rez, 1);

% decide on cutoff
rez = set_cutoff(rez);
% eliminate widely spread waveforms (likely noise)
rez.good = get_good_units(rez);

fprintf('found %d good units \n', sum(rez.good>0))

% write to Phy
fprintf('Saving results to Phy  \n')
rezToPhy(rez, rootZ);

% %% if you want to save the results to a Matlab file...
% 
% % discard features in final rez file (too slow to save)
% rez.cProj = [];
% rez.cProjPC = [];
% 
% % final time sorting of spikes, for apps that use st3 directly
% [~, isort]   = sortrows(rez.st3);
% rez.st3      = rez.st3(isort, :);
% 
% % Ensure all GPU arrays are transferred to CPU side before saving to .mat
% rez_fields = fieldnames(rez);
% for i = 1:numel(rez_fields)
%     field_name = rez_fields{i};
%     if(isa(rez.(field_name), 'gpuArray'))
%         rez.(field_name) = gather(rez.(field_name));
%     end
% end
% 
% % save final results as rez2
% fprintf('Saving final results in rez2  \n')
% fname = fullfile(rootZ, 'rez2.mat');
% save(fname, 'rez', '-v7.3');