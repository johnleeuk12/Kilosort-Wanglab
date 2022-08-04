addpath('D:\GitHub\TDTMatlabSDK');
fpath = fullfile('D:\Data\example_TDT recording data');
TDTdata = TDTbin2mat(fpath,'TYPE',{'epocs'});