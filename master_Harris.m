function master_Harris()


% default options are in parenthesis after the comment
useGPU = 1; %else 1  % do you have a GPU? Kilosorting 1000sec of 32chan simulated data takes 55 seconds on gtx 1080 + M2 SSD.
ops.ephys_type = file_type;


session_name = '2019-04-04_15-17-23';
animal_name = 'M94W';
file_type = '100';
fpath    = fullfile('C:\DATA\OpenEphys\',animal_name, filesep , session_name); % where on disk do you want the simulation? ideally and SSD...
if ~exist(fpath, 'dir'); mkdir(fpath); end
addpath(genpath([Gitpath 'KiloSort'])) % path to kilosort folder
addpath(genpath([Gitpath 'npy-matlab'])) % path to npy-matlab scripts
pathToYourConfigFile = [Gitpath 'Kilosort-Wanglab']; % take from Github folder and put it somewhere else (together with the master_file)

make_HarrisChannelMap(fpath)

% addpath(genpath('C:\Users\John.Lee\Documents\GitHub\copy\KiloSort')) % path to kilosort folder
addpath(genpath('C:\Users\Seth\Documents\GitHub\copy\KiloSort')) % path to kilosort folder

% addpath(genpath('C:\Users\John.Lee\Documents\GitHub\npy-matlab')) % path to npy-matlab scripts
addpath(genpath('C:\Users\Seth\Documents\GitHub\npy-matlab')) % path to npy-matlab scripts

% pathToYourConfigFile = 'C:\Users\John.Lee\Documents\GitHub\Kilosort-Wanglab'; % take from Github folder and put it somewhere else (together with the master_file)
pathToYourConfigFile = 'C:\Users\Seth\Documents\GitHub\Kilosort-Wanglab'; % take from Github folder and put it somewhere else (together with the master_file)

run(fullfile(pathToYourConfigFile, 'Harris_config.m'))

%
if ops.GPU     
    gpuDevice(1); % initialize GPU (will erase any existing GPU arrays)
end
% run convertOpenEphysToRawBInary if haven't already
% test = 1;
if strcmp(ops.datatype , 'openEphys') && ~exist(fullfile(ops.root,ops.fbinary),'file')
   tic
   disp('converting data...')
   ops = convertOpenEphysToRawBInary(ops);  % convert data, only for OpenEphys
   disp('converting data... Done')
   toc
else
    fprintf('Binary file already created for session %1s\n',session_name)
end
%
if ~exist(fullfile(ops.root,'params.py'),'file')
    tic
    [rez, DATA, uproj] = preprocessData(ops); % preprocess data and extract spikes for initialization
    rez                = fitTemplates(rez, DATA, uproj);  % fit templates iteratively
    rez                = fullMPMU(rez, DATA);% extract final spike times (overlapping extraction)


    % load(fullfile(ops.root,  'rez.mat'));
    % AutoMerge. rez2Phy will use for clusters the new 5th column of st3 if you run this)
    rez = merge_posthoc2(rez);
    % save matlab results file
    save(fullfile(ops.root,  'rez.mat'), 'rez', '-v7.3');

    % save python results file for Phy
    rezToPhy(rez, ops.root);

    toc
    fprintf('KiloSort finished for session %1s\n',session_name)
    % remove temporary file
    delete(ops.fproc);
else
    fprintf('KiloSort already processed session %1s\n',session_name)
end
%%
