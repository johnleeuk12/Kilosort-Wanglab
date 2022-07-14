function sr = TDTtoDAT(fpath,fname)
addpath('D:\GitHub\TDTMatlabSDK');
% curr_path=pwd; % path where .tsq, .tnt, .tin, .tev, .Tdx, .Tbk file exist
% curr_path = 'D:\Data\example_TDT recording data';


TDT_data = TDTbin2mat(fpath,'STORE','SU_1');% load TDT data in matlab
%%
data=TDT_data.streams.SU_1.data; % single unit data after filtering and referencing

sr=TDT_data.streams.SU_1.fs; % sampling rate
clearvars TDT_data % remove unused data from memory

%%
newFilename = [curr_path filesep fname '.dat'];


FIDw = fopen(newFilename, 'w');

% Writing data into file
disp('Writing the converted data into the new .dat file...');
fwrite(FIDw, data, 'int16');
fclose(FIDw);