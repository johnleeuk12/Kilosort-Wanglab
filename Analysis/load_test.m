addpath('C:\Users\Seth\Documents\GitHub\analysis-tools-master\')



tic
fpath = 'D:\Data\Experiments\M60F\2020-12-07_16-18-38';
addpath(fpath);
fname = 'test_binary.dat';
filepath = fullfile(fpath,filesep,fname);

fid = fopen(filepath,'r');
buff = fread(fid,'*int16');
fclose(fid);
toc


file_type = '100';

tic
fprintf('Time %3.0fs. Loading data... \n', toc);
parfor ch = 1:64
    %     disp(['loading data from channel ' num2str(ch)])
    [data1,timestamps1{ch},info1{ch}] = load_open_ephys_data_faster([fpath filesep file_type '_CH' num2str(ch) '.continuous' ]);
    data_all(ch,:) = data1;
end
obj.params.rec_start_time = timestamps1{1}(1);
obj.raw_data = data_all;
fprintf('Time %3.0fs. Loading data... Done \n', toc);