function [fpath, savepath] = directories(PC_name,animal_name,session_name)
% Select PC name : Chamber_T, 426_Analysis 426_John

switch PC_name
    case '426_John'
    start_path = 'C:\Users\John.Lee\';
    kilo_path = 'copy\KiloSort';
    analysis_path  ='\analysis-tools';
    

    case {'Chamber_T','426_Analysis'}
    start_path = 'C:\Users\Seth\';
    kilo_path = 'Kilosort2';
    analysis_path = '\analysis-tools-master';
end

addpath(genpath(fullfile(start_path, 'Documents\GitHub\', kilo_path))); % path to kilosort folder
addpath(genpath(fullfile(start_path, 'Documents\GitHub\npy-matlab'))); % path to npy-matlab scripts
addpath(genpath(fullfile(start_path,'Documents\GitHub', analysis_path)));

addpath(fullfile('D:\Data\Experiments\', animal_name));
fpath = fullfile('D:\Data\Experiments\', animal_name, filesep);
savepath = fullfile('D:\Data\Units\',animal_name,filesep);
%, session_name);% where on disk do you want the simulation? ideally and SSD...


% addpath(fullfile('D:\DATA\Experiments\', animal_name));
% fpath = fullfile('D:\DATA\Experiments\', animal_name, filesep, session_name);% where on disk do you want the simulation? ideally and SSD...