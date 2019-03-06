% File path to Kilosort-Wanglab folder
GitPath = 'C:\Users\skoehler\Documents\GitHub\';
ourGitpath = fullfile(GitPath, 'Kilosort-Wanglab');
loadOEpath = fullfile(GitPath,'analysis-tools');
addpath(ourGitpath,loadOEpath)

% File path to Log 
logpath =  'C:\Data\Data Logs\Neurophysiology Log - Multi Electrode Array_KiloSort.xlsx';
LogTable = readtable(logpath,'Sheet','ExpSessions');
sessions = LogTable.Session;
OE_File = string(LogTable.OE_File);
Date    = datetime(LogTable.Date,'Format','yyyy-MM-dd');

% File path to data
datapath = '\\datacenterchx.bme.jhu.edu\Project_TNT\Data\Experiments\M44D\';

Nsessions = length(sessions);
for sesh = 1:Nsessions
    curr_date = char(Date(sesh));
    Header.curr_sesh = sessions(sesh);
    OE_folder = [curr_date '_' char(OE_File(sesh))];
%     OE_events_filepath = [datapath animal filesep OE_folder filesep 'all_channels.events'];
    if curr_date > datetime(2018,12,3) %after this date, a new profile began being used on Open Ephys, which changed this number in the file name
        file_type = '116';
    else
        file_type = '100';
    end
    
    master_Harris(OE_folder,file_type,GitPath,datapath)
end