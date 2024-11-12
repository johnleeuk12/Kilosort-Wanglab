





addpath(genpath('D:\GitHub\Kilosort2')) % path to kilosort folder

addpath(genpath('D:\\GitHub\npy-matlab')) % path to npy-matlab scripts
addpath(genpath('D:\\GitHub\spikes')) % path to npy-matlab scripts



fpath = fullfile('D:\DATA\Reversal Learning\TSn003\2024-10-30_13-50-34\probeC');
event_id = 0;
eventpath = 'events\OE_FPGA_Acquisition_Board-114.Rhythm Data\TTL';
timepath = 'OE_DAQ';
savedir = 'D:\DATA\Units\TS_2';
animal_id = 'TSn003';
region = 'V1';
date = '2024-10-30_13-50-34';



%% loading spikes and waveforms 

sp = loadKSdir(fpath);

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


%             buff = buff(1:obj.params.Nb_ch_real);
            
fprintf('Time %3.0fs. loading data... complete! \n', toc);


fprintf('Time %3.0fs. Extracting waveforms... \n', toc);

switch event_id
    case 1
        gwfparams.nCh = 72;                      % Number of channels that were streamed to disk in .dat file
        gwfparams.buff = reshape(gwfparams.buff,72,[]);
        lick = gwfparams.buff(65,:);
        save(fullfile(fpath,filesep,'lick.mat'),'lick')
    case 0
        gwfparams.buff = reshape(gwfparams.buff,64,[]);
        gwfparams.nCh = 64;                      % Number of channels that were streamed to disk in .dat file
end

gwfparams.fileName = filepath;         % .dat file containing the raw
gwfparams.dataType = 'int16';            % Data type of .dat file (this should be BP filtered)
gwfparams.wfWin = [-20 40];              % Number of samples before and after spiketime to include in waveform
gwfparams.nWf = 2000;                    % Number of waveforms per unit to pull out

wf = {};
for c = 1:length(sp.cgs)
    temp = [];
    gwfparams.spikeTimes = ceil(sp.st(sp.clu==sp.cids(c))*30000); % Vector of cluster spike times (in samples) same length as .spikeClusters
    gwfparams.spikeClusters = sp.clu(sp.clu==sp.cids(c));
    wf{c} = getWaveForms(gwfparams);
    temp(:,:) =  wf{c}.waveFormsMean(1,:,:);
    [~,wf{c}.bestch] = max(mean(abs(temp),2));

end

clear gwfparams
% % 
% test = [];
% 
% test(:,:) =  wf{2}.waveFormsMean(1,:,:);
% 
% test2 = mean(abs(test),2);
% % figure
% for ch = 1:64
%     plot(test(ch,:))
%     hold on 
% end

fprintf('Time %3.0fs. Extracting waveforms... complete! \n', toc);


%% extracting event time file

switch event_id
    case 1
        event.state = readNPY(fullfile(fpath, filesep, eventpath, filesep, 'states.npy'));
        event.time = readNPY(fullfile(fpath, filesep, eventpath, filesep, 'timestamps.npy'));
        timestamps = event.time(event.state == 1);

        Y = readNPY(fullfile(fpath,filesep,timepath, filesep,'timestamps.npy'));
        timestamps = timestamps-Y(1);
        clear Y
    case 0
        timestamps = [];
end
%% save units


save_dir = fullfile(savedir,filesep,animal_id);
if ~isfolder(save_dir)
    mkdir(save_dir);
end
unit_list_name = [animal_id '_unit_list.mat'];
if exist(fullfile(save_dir,filesep,unit_list_name))
    load(fullfile(save_dir,unit_list_name));
    test_list = find(strcmp(unit_list.data(:,5),date));
    if  ~isempty(test_list)
        if strcmp(unit_list.data{test_list(end),7},region) 
            error('Units already saved')
        end
    end
    
else
    unit_list = {};
    unit_list.tags{1} = 'id';
    unit_list.tags{2} = 'SNR';
    unit_list.tags{3} = 'Amplitude';
    unit_list.tags{4} = 'Cid';
    unit_list.tags{5} = 'filename';
    unit_list.tags{6} = 'animal_id';
    unit_list.tags{7} = 'experiment';
    
    
    unit_list.data = {};
    save(fullfile(save_dir,unit_list_name),'unit_list');
end

% tic
for c = 1:length(sp.cgs)
    if sp.cgs(c) ==2
        s_unit = {};
        
        % Get waveform information
        s_unit.best_ch = wf{c}.bestch;
%         %                 s_unit.ch = CchanMap(Cchannels(SU(id)));
%         s_unit.templates = obj.templates{id};
        s_unit.cid = sp.cids(c);
%         s_unit.SU = SU;
%         s_unit.SU_good = obj.SU_good;
%         s_unit.SNR = obj.waveforms.SNR(id);
        s_unit.waveforms = wf{c};
        %                 s_unit.spiketimes = spike_times_all(find(spike_clusters == Cids(SU(id))))/obj.params.fs;
        s_unit.spiketimes = sp.st(sp.clu==sp.cids(c));
        %                     s_unit.xbz_file_name = obj.params.xbz_file_name;
        s_unit.event_time = timestamps;
%         s_unit.amplitude = max(abs(obj.waveforms.mean{id}));
        %                     s_unit.start_times = obj.start_times;
        %                     s_unit.end_times =obj.end_times;
        unitname= [animal_id 'u00001']; %03/03/20 changed from u001 to u00001
        if isempty(unit_list.data)
            save(fullfile(save_dir,unitname),'s_unit')
            unit_list.data{1,1} = 1;
%             unit_list.data{1,2} = s_unit.SNR;
%             unit_list.data{1,3} = s_unit.amplitude;
            unit_list.data{1,4} = s_unit.cid;
            unit_list.data{1,5} = date;
            unit_list.data{size((unit_list.data),1),6} = animal_id;
            unit_list.data{size((unit_list.data),1),7} = region;
        else
            unitname = [unitname(1:end-length(num2str(size((unit_list.data),1)+1))) num2str(size((unit_list.data),1)+1)];
            save(fullfile(save_dir,unitname),'s_unit')
            unit_list.data{size((unit_list.data),1)+1,1} = size((unit_list.data),1)+1;
%             unit_list.data{size((unit_list.data),1),2} = s_unit.SNR;
%             unit_list.data{size((unit_list.data),1),3} = s_unit.amplitude;
            unit_list.data{size((unit_list.data),1),4} = s_unit.cid;
            unit_list.data{size((unit_list.data),1),5} = date;
            unit_list.data{size((unit_list.data),1),6} = animal_id;
            unit_list.data{size((unit_list.data),1),7} = region;
        end
        save(fullfile(save_dir,unit_list_name),'unit_list')
        fprintf('Time %3.0fs. unit %3.0f/ %3.0f.  \n', toc,c,length(sp.cgs));

    end
    
end

fprintf('%3.0f/ %3.0f Units saved  \n', sum(sp.cgs == 2),length(sp.cgs));

disp('Units saved')



