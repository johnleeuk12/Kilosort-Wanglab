% default options are in parenthesis after the comment
useGPU = 0; %else 1  % do you have a GPU? Kilosorting 1000sec of 32chan simulated data takes 55 seconds on gtx 1080 + M2 SSD.



fpath    = 'D:\DATA\OpenEphys\M44D\2018-11-20_13-13-43'; % where on disk do you want the simulation? ideally and SSD...
if ~exist(fpath, 'dir'); mkdir(fpath); end

make_HarrisChannelMap(fpath)

addpath(genpath('C:\Users\John.Lee\Documents\GitHub\copy\KiloSort')) % path to kilosort folder
addpath(genpath('C:\Users\John.Lee\Documents\GitHub\npy-matlab')) % path to npy-matlab scripts

pathToYourConfigFile = 'C:\Users\John.Lee\Documents\GitHub\Harrisprobe'; % take from Github folder and put it somewhere else (together with the master_file)
run(fullfile(pathToYourConfigFile, 'Harris_config.m'))

tic; % start timer
%
if ops.GPU     
    gpuDevice(1); % initialize GPU (will erase any existing GPU arrays)
end

disp('converting data...')
tic
test = 1;
if strcmp(ops.datatype , 'openEphys')
   ops = convertOpenEphysToRawBInary(ops);  % convert data, only for OpenEphys
end
disp('converting data... Done')
toc
%
[rez, DATA, uproj] = preprocessData(ops); % preprocess data and extract spikes for initialization
rez                = fitTemplates(rez, DATA, uproj);  % fit templates iteratively
rez                = fullMPMU(rez, DATA);% extract final spike times (overlapping extraction)


% load(fullfile(ops.root,  'rez.mat'));
% AutoMerge. rez2Phy will use for clusters the new 5th column of st3 if you run this)
%     rez = merge_posthoc2(rez);

% save matlab results file
save(fullfile(ops.root,  'rez.mat'), 'rez', '-v7.3');

% save python results file for Phy
rezToPhy(rez, ops.root);

toc
% remove temporary file
delete(ops.fproc);
%%