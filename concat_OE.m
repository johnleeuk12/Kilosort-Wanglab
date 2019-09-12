function concat_OE()


%% Configs. must change
useGPU = 1; %else 1  % do you have a GPU? Kilosorting 1000sec of 32chan simulated data takes 55 seconds on gtx 1080 + M2 SSD.

session_list = {'2019-07-27_16-07-19'
'2019-07-27_16-09-07'
'2019-07-27_16-10-59'
'2019-07-27_16-12-48'
'2019-07-27_16-14-56'
'2019-07-27_16-25-04'
'2019-07-27_16-33-40'
'2019-07-27_16-43-42'
'2019-07-27_16-48-35'
'2019-07-27_16-57-15'
'2019-07-27_17-05-31'
'2019-07-27_17-20-48'
'2019-07-27_17-24-17'
'2019-07-27_17-28-41'
'2019-07-27_17-34-34'
'2019-07-27_17-37-32'
'2019-07-27_17-41-36'






};
concat_folder_name = 'H6T4S2_concat';


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
