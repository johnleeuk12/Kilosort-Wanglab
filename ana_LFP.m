
addpath(genpath('D:\GitHub\Kilosort2')) % path to kilosort folder

addpath(genpath('D:\\GitHub\npy-matlab')) % path to npy-matlab scripts
addpath(genpath('D:\\GitHub\spikes')) % path to npy-matlab scripts



fpath = fullfile('D:\DATA\WG2\2023-12-06_11-47-37');
eventpath = 'events\OE_FPGA_Acquisition_Board-114.Rhythm Data\TTL';
timepath = 'OE_DAQ';
fs = 30000;
% savedir = 'D:\DATA\Units\RL';
% animal_id = 'RLn605';
% region = 'IC';
% date = '2023-10-17_17-49-49';

%% 

gwfparams.dataDir = fpath;    % KiloSort/Phy output folder           
filepath = dir(fullfile(fpath, '*.dat'));

if strcmp(filepath(1).name,'temp_wh.dat')
    filepath = fullfile(fpath,filepath(2).name);
else
    filepath = fullfile(fpath,filepath(1).name);
end

tic
fprintf('Time %3.0fs. loading data... \n', toc);

fid = fopen(filepath,'r');
gwfparams.buff = fread(fid,'*int16');
fclose(fid);
toc

fprintf('Time %3.0fs. loading data... complete! \n', toc);


% fprintf('Time %3.0fs. Extracting waveforms... \n', toc);

gwfparams.buff = reshape(gwfparams.buff,72,[]);

lick = gwfparams.buff(65,:);
save(fullfile(fpath,filesep,'lick.mat'),'lick')

gwfparams.buff = gwfparams.buff(1:64,:);


%% time stamps


event.state = readNPY(fullfile(fpath, filesep, eventpath, filesep, 'states.npy'));
event.time = readNPY(fullfile(fpath, filesep, eventpath, filesep, 'timestamps.npy'));
timestamps = event.time(event.state == 1);

Y = readNPY(fullfile(fpath,filesep,timepath, filesep,'timestamps.npy'));
timestamps = timestamps-Y(1);
clear Y

timestamps_SR = timestamps*fs;

%% 

nb_tr = length(timestamps);
nb_tr = nb_tr-1;
R = zeros(64,5*fs,nb_tr);
L = zeros(1,5*fs,nb_tr);

for tr = 1:nb_tr
    R(:,:,tr) = gwfparams.buff(:,int32(timestamps_SR(tr,1)-fs):int32(timestamps_SR(tr,1)+4*fs-1));
    L(:,:,tr) = lick(:,int32(timestamps_SR(tr,1)-fs):int32(timestamps_SR(tr,1)+4*fs-1));
%     R(:,:,tr) = R(:,:,tr)- median(R(:,:,tr),1);
end
    
R = imresize(R,[64,5000]);
R = R(1:end,700:2000,:);
R = R-median(R,1);
R_mean = mean(R,3);
R_mean = R_mean-mean(R_mean(:,1:200),2);
% R_mean = imgaussfilt(R_mean,[1,50]);
R_std = std(R,[],3);
R_std = R_std/sqrt(64);

figure
% % % plot(R_mean(12,:))
% 
ch = 7;
% Rp = R(1:12,750:2000,:);
pre =200;
shadedErrorBar(-pre+1:size(R_mean,2)-pre,R_mean(ch,:),R_std(ch,:),'lineProps','b')
hold on 
% plot(-pre+1:size(R_mean,2)-pre,

% Rp_mean = imgaussfilt(mean(mean(R,3),1),[1,50]);
% Rp_std = mean(std(R,[],3),1)/sqrt(12*nb_tr);
% shadedErrorBar(1:5000,Rp_mean,Rp_std,'lineProps','b')
% hold on
% for tr = 1:nb_tr
%     plot(-pre+1:size(R_mean,2)-pre,R(ch,:,tr))
%     hold on
% end

xlim([-200,1000])
% hold off
% hold on
%     pause

% 
% figure
% imagesc(imgaussfilt(R_mean(:,:),[4,70]))
% caxis([-150,100])
%% Stats

figure
% % % plot(R_mean(12,:))
% 
ch = 7;
% Rp = R(1:12,750:2000,:);
pre =200;
SD = mean(R_std(ch,1:200),2);
shadedErrorBar(-pre+1:size(R_mean,2)-pre,R_mean(ch,:),R_std(ch,:),'lineProps','b')
hold on 

for t = 1: 1200
    if R_mean(ch,t) > 1.5*SD
        scatter(t-pre,600,15,'k','filled')
    end
end
hold off





