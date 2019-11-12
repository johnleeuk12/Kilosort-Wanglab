function concat_OE()


%% Configs. must change
useGPU = 1; %else 1  % do you have a GPU? Kilosorting 1000sec of 32chan simulated data takes 55 seconds on gtx 1080 + M2 SSD.

session_list = {
'2019-07-26_16-31-51'
'2019-07-26_16-34-50'
'2019-07-26_16-42-16'
'2019-07-26_16-47-30'
'2019-07-26_16-56-03'
'2019-07-26_17-05-34'
'2019-07-26_17-08-48'
'2019-07-26_17-15-29'
'2019-07-26_17-21-53'
'2019-07-26_17-26-29'
'2019-07-26_17-31-04'
'2019-07-26_17-32-54'
'2019-07-26_17-34-41'
'2019-07-26_17-36-30'
'2019-07-26_17-38-31'

};
concat_folder_name = 'H6T3S1_concat';


Animal_name = 'M12E';
addpath(genpath('C:\Users\Seth\Documents\GitHub\KiloSort2')) % path to kilosort folder
addpath(genpath('C:\Users\Seth\Documents\GitHub\npy-matlab')) % path to npy-matlab scripts
file_type = '100';
pathToYourConfigFile = 'C:\Users\Seth\Documents\GitHub\Kilosort-Wanglab'; % take from Github folder and put it somewhere else (together with the master_file)



for s = 1:length(session_list)
    session_name =session_list{s};
    % session_name = 'H6T3S1_concat';
    fpath    = fullfile('C:\DATA\OpenEphys', filesep, Animal_name, filesep, session_name); % where on disk do you want the simulation? ideally and SSD...
    rootZ = fpath;
    
    if ~exist(fpath, 'dir'); mkdir(fpath); end
    
    % rmpath(genpath('C:\Users\Seth\Documents\GitHub\KiloSort'))
    addpath(genpath('C:\Users\Seth\Documents\GitHub\KiloSort2')) % path to kilosort folder
    addpath(genpath('C:\Users\Seth\Documents\GitHub\npy-matlab')) % path to npy-matlab scripts
    
    pathToYourConfigFile = 'C:\Users\Seth\Documents\GitHub\Kilosort-Wanglab'; % take from Github folder and put it somewhere else (together with the master_file)
    run(fullfile(pathToYourConfigFile, 'Harris_config2.m'))
    
    ops.trange = [0 Inf]; % time range to sort
    ops.ephys_type = file_type;
    ops.datatype = 'openEphys';

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
    
    
    
end








%% concatenating files

fname = 'test_binary.dat';
filepath    = fullfile('C:\DATA\OpenEphys', filesep, Animal_name, filesep, session_list, filesep, fname);

fpath    = fullfile('C:\DATA\OpenEphys', filesep, Animal_name, filesep, concat_folder_name); % where on disk do you want the simulation? ideally and SSD...

if ~exist(fpath, 'dir'); mkdir(fpath); end

    

run(fullfile(pathToYourConfigFile, 'Harris_config2.m'))
ops.datatype = 'dat';

buff= {};
% Big_buff = [];
tic
fileID = fopen(fullfile(fpath,filesep,fname),'w');

for t = 1:length(session_list)

    fid         = fopen(filepath{t}, 'r');
    buff{t} = fread(fid, '*int16');
    fclose(fid);
    disp(t)
    fwrite(fileID,buff{t},'*int16');

    toc
end
fclose(fileID);
