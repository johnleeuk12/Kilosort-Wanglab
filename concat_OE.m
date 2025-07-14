function concat_OE(concat_folder_name)


%% Configs. must change
useGPU = 1; %else 1  % do you have a GPU? Kilosorting 1000sec of 32chan simulated data takes 55 seconds on gtx 1080 + M2 SSD.

exp_name = 'Reversal Learning';
Animal_name = 'MPM04';
session_list = {
'2025-05-30_15-17-46'
'2025-05-30_16-52-26'

};
% concat_folder_name = 'H11T6S1_concat';


% addpath(genpath('C:\Users\Seth\Documents\GitHub\KiloSort2')) % path to kilosort folder
% addpath(genpath('C:\Users\Seth\Documents\GitHub\npy-matlab')) % path to npy-matlab scripts

addpath(genpath('D:\GitHub\Kilosort2')) % path to kilosort folder

addpath(genpath('D:\\GitHub\npy-matlab')) % path to npy-matlab scripts
% file_type = '100';
% pathToYourConfigFile = 'D:\GitHub\Kilosort-Wanglab'; % take from Github folder and put it somewhere else (together with the master_file)
% 

%% init parameters

NchanTOT = 64*3 + 8;
Nchan = 64;
ntbuff = 64;
NT= 32*1024 + ntbuff;
NTbiff = NT; %+ 3*ntbuff; 


%% concatenating files
concat_folder_name = 'concat_0529';
fname1 = 'continuous.dat';
fname2 = 'continuous.dat';
filepath    = fullfile('D:\DATA',filesep,exp_name, filesep, Animal_name, filesep, session_list, filesep, fname1);

fpath    = fullfile('D:\DATA',filesep,exp_name, filesep, Animal_name, filesep, concat_folder_name); % where on disk do you want the simulation? ideally and SSD...

if ~exist(fpath, 'dir'); mkdir(fpath); end

    

% run(fullfile(pathToYourConfigFile, 'Harris_config2.m'))
% ops.datatype = 'dat';

buff= {};
% Big_buff = [];
tic
fileID = fopen(fullfile(fpath,filesep,fname2),'w');

for t = 1:length(session_list)
    fid = fopen(filepath{t}, 'r');
    bytes       = get_file_size(filepath{t}); % size in bytes of raw binary
    nTimepoints = floor(bytes/NchanTOT/2); % number of total timepoints
    Nbatch      = ceil(nTimepoints /NT); 
    for ibatch = 1:Nbatch
        if mod(ibatch,100) ==1
        fprintf(['%4d /' num2str(Nbatch) ' time : %6.2f sec \n'],ibatch,toc')
        end
        % offset = max(0,2*NchanTOT*(NT * (2-1))); % number of samples to start reading at.
        offset = 2*NchanTOT*NTbiff*(ibatch-1);
        fseek(fid, offset, 'bof'); % fseek to batch start in raw file
        buff = fread(fid, [NchanTOT NTbiff], '*int16');
        fwrite(fileID,buff,'*int16');
    end

    fclose(fid);

    toc
end
fclose(fileID);



%%
% 
% for s = 1:length(session_list)
%     session_name =session_list{s};
%     % session_name = 'H6T3S1_concat';
%     fpath    = fullfile('D:\DATA', filesep,exp_name,filesep, Animal_name, filesep, session_name); % where on disk do you want the simulation? ideally and SSD...
%     rootZ = fpath;
%     
%     try exist(fpath, 'dir');
%         %         ; mkdir(fpath); end
%     catch
%         disp('Error. No folder found')
%         disp([fpath ' does not exist'])
%         break
% 
%     end
% 
%     % rmpath(genpath('C:\Users\Seth\Documents\GitHub\KiloSort'))
% %     addpath(genpath('C:\Users\Seth\Documents\GitHub\KiloSort2')) % path to kilosort folder
% %     addpath(genpath('C:\Users\Seth\Documents\GitHub\npy-matlab')) % path to npy-matlab scripts
%     
% %     pathToYourConfigFile = 'C:\Users\Seth\Documents\GitHub\Kilosort-Wanglab'; % take from Github folder and put it somewhere else (together with the master_file)
% %     run(fullfile(pathToYourConfigFile, 'Harris_config3.m'))
% %     
% %     ops.trange = [0 Inf]; % time range to sort
% %     ops.ephys_type = file_type;
% %     ops.datatype = 'openEphys';
% % 
% %     % this block runs all the steps of the algorithm
% %     fprintf('Looking for data inside %s \n', rootZ)
% %     
% %     % is there a channel map file in this folder?
% %     fs = dir(fullfile(rootZ, 'chan*.mat'));
% %     if ~isempty(fs)
% %         ops.chanMap = fullfile(rootZ, fs(1).name);
% %     else
% %         make_HarrisChannelMap(fpath)
% %     end
%     
%     
%     
%     tic; % start timer
%     %
% %     if ops.GPU
% %         gpuDevice(1); % initialize GPU (will erase any existing GPU arrays)
% %     end
%     
%     % Openephys to binary converter
%     ops.ephys_type = file_type;
%     
%     disp('converting data...')
%     tic
%     % test = 1;
%     if strcmp(ops.datatype , 'openEphys')
%         ops = convertOpenEphysToRawBInary(ops);  % convert data, only for OpenEphys
%     end
%     disp('converting data... Done')
%     toc
%     
%     
%     
% end
% 




