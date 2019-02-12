function master_Harris(session_name, file_type)


% default options are in parenthesis after the comment
useGPU = 1; %else 1  % do you have a GPU? Kilosorting 1000sec of 32chan simulated data takes 55 seconds on gtx 1080 + M2 SSD.
ops.ephys_type = file_type;

% file paths
Gitpath = 'C:\Users\Seth\Documents\GitHub\';
datapath = '\\datacenterchx.bme.jhu.edu\Project_TNT\Data\Experiments\M44D\';
fpath    = fullfile(datapath, session_name); % where on disk do you want the simulation? ideally and SSD...
if ~exist(fpath, 'dir'); mkdir(fpath); end
addpath(genpath([Gitpath 'KiloSort'])) % path to kilosort folder
addpath(genpath([Gitpath 'npy-matlab'])) % path to npy-matlab scripts
pathToYourConfigFile = [Gitpath 'Kilosort-Wanglab']; % take from Github folder and put it somewhere else (together with the master_file)

% run HarrisProbe_spkext if haven't already, then run config file 
if ~exist(fullfile(fpath,'PCspikes3.mat'),'file')
    HarrisProbe_spkext(fpath,file_type)
    run(fullfile(pathToYourConfigFile, 'Harris_config.m'))
else
    fprintf('Spike templates already extracted for session %1s\n',session_name)
    run(fullfile(pathToYourConfigFile, 'Harris_config.m'))
end
% run make_HarrisChannelMap if haven't already
if ~exist(fullfile(fpath, 'chanMap.mat'),'file')
    make_HarrisChannelMap(fpath)
else
    fprintf('chanMap file already created for session %1s\n',session_name)
end

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
