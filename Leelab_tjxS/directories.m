function [fpath, savepath] = directories(PC_name)
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
    
    case 'LeeLab_John'
        start_path = 'D:\';
        kilo_path = 'Kilosort2';
        analysis_path  ='\analysis-tools';
        
    case 'LeeLab_analysis'
        start_path = 'K:\';
        kilo_path = 'Kilosort2';
        analysis_path  ='\analysis-tools';

end

addpath(genpath(fullfile(start_path, 'GitHub\', kilo_path))); % path to kilosort folder
addpath(genpath(fullfile(start_path, 'GitHub\npy-matlab'))); % path to npy-matlab scripts
addpath(genpath(fullfile(start_path,'GitHub', analysis_path)));
addpath(genpath(fullfile(start_path,'GitHub\Kilosort-Wanglab\Analysis')));
addpath(genpath(fullfile(start_path,'GitHub\Kilosort-Wanglab\TDTMatlabSDK')));
addpath(genpath(fullfile(start_path,'GitHub\NPMK')));




fpath = fullfile(start_path,'Data');
savepath = fullfile(start_path,'Data\Units');
if ~isfolder(savepath)
    mkdir(savepath);
end
