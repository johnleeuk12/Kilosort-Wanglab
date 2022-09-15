
%{
    tjx is the first stage in the reclustering step in the pipeline. 
the inputs to tjx are cluster ids and spike times from Phy and Kilosort,
and the output are spike times of units and other relevent unit information
such as waveforms, best channel, and information on when and where it was
recorded. the output is sent to tjd, which allows users to merge and split 
units, as well as decide which to define as single and multiunits.

    tjx runs in the order as shown in run_tjx.m. However, tjx can also
extract layer information through CSD analysis. 

    DO NOT MODIFY THE CODE without permission from JHL

%%%%%%%%%%%%%% EDIT LOG %%%%%%%%%%%%%%%%%%%%%%%%
11/27/2020 JHL
    turned xblaster stim data (parameter x) into a class property obj.x
    
11/29/2020 Saul
    Copied part of the Single_units() code into extract_units() without
    modifying original functionality

12/3/2020 JHL
    modified code for save_units() so that code isnt specific to M12E
    modified save directory for units

07/27/2022 JHL
    Changing code to fit with LeeLab, including TDT data support.
    changes also include epoch data

%}




classdef tjxS < handle
    
    properties(Constant)

    end
    
    properties
        params;
        raw_data = {};
        filt_data = {};
        waveforms = {};
        analysis_code = [];
        analysis_type = {};
        plot_type = {};
        spike_times = {};
        manual = 0; % variable to indicate whether the clusters were sorted manually
        template = [];
        seg_length = [];
        templates = {};
        SU_good = [];
        start_times = [];
        end_times = [];
        clusters = []; %NEW ENTRY
        idds = []; %NEW
        SUU = []; %NEW
        Cidds = []; %NEW
        good_ch = []; %NEW
        x = {};
        Eventdata = {};
    end
    
    
    methods(Static)
        
    end
    methods

        %%
        function obj = tjxS(pf,a,b,c) %New entry, good one
            obj.params = pf(a,b,c);
        end
        %%
        function obj = initialize(obj)
            global spike_times_all spike_clusters Cchannels Cids Cgroups new_channels
            
            % load sorted spikes for entire session
            try
                spike_times_all = readNPY(fullfile(obj.params.fpath2, 'spike_times_original.npy'));
            catch
                spike_times_all = readNPY(fullfile(obj.params.fpath2, 'spike_times.npy'));
            end
            spike_times_all = double(spike_times_all);
            spike_clusters = readNPY(fullfile(obj.params.fpath2, 'spike_clusters.npy'));
            [Cids, Cgroups] = readClusterGroupsCSV(fullfile(obj.params.fpath2, 'cluster_group.tsv'));
            %             spike_templates = readNPY(fullfile(obj.params.fpath2, 'spike_templates.npy'));
            new_channels = readNPY(fullfile(obj.params.fpath2, 'channel_map.npy'));
            new_channels = new_channels+1;
            Cchannels = [];
            obj.Cidds = [obj.Cidds Cids];
            
        end
        %% Loading event files 
        function obj = load_data_OE(obj)
            
            data_all = [];
            timestamps1 = {};
            info1 = {};
            file_type = '100';
            ch = 1;
            [~,timestamps1{ch},~] = load_open_ephys_data_faster([obj.params.fpath obj.params.session_name filesep obj.params.file_type '_CH' num2str(ch) '.continuous' ]);

            obj.params.rec_start_time = timestamps1{1}(1);

        end
        
        function obj = load_data_TDT(obj)
            obj.Eventdata = TDTbin2mat(obj.params.fpath2,'TYPE',{'epocs'});
            
        end
        
        function obj = load_data_NEV(obj)
            fname = erase(obj.params.fname, '.ns5.dat');
            fname = [fname,'.nev'];
            
            temp = openNEV(fullfile(obj.params.fpath2, fname),'nosave');
            obj.Eventdata = temp.Data.SerialDigitalIO;
        end
        %%
        function obj = extract_units(obj)
            tic
            global spike_times_all spike_clusters Cchannels Cids Cgroups new_channels SU
            %// loading data
            fprintf('Time %3.0fs. loading data... \n', toc);
            %             fpath = 'D:\Data\Experiments\M60F\2020-12-07_16-18-38';
            %             addpath(fpath);
%             fname = 'test_binary.dat';
%             filepath = fullfile(obj.params.fpath, obj.params.session_name, filesep, fname);
            %             filepath = fullfile(fpath,filesep,fname);
            
            filepath = dir(fullfile(obj.params.fpath2, '*.dat'));
            
            if strcmp(filepath(1).name,'temp_wh.dat')
                filepath = fullfile(obj.params.fpath2,filepath(2).name);
            else
                filepath = fullfile(obj.params.fpath2,filepath(1).name);
            end
%             filepath = fullfile(obj.params.fpath2,obj.params.fname);
            fid = fopen(filepath,'r');
            buff = fread(fid,'*int16');
            fclose(fid);
            toc
            
                       
            buff = reshape(buff,obj.params.Nb_ch,[]);
            
            fprintf('Time %3.0fs. loading data... complete! \n', toc);
            
            %// preprocessing
            fprintf('Time %3.0fs. filtering data... \n', toc);
            [b,a] = butter(4, [0.0244 0.6104],'bandpass');
            %// extracting waveforms from processed data
            obj.waveforms.raw = {};
            obj.waveforms.mean = {};
            obj.waveforms.std = {};
            obj.waveforms.SNR = [];
            obj.templates = {};
            %             spike_times_all = {};
            SU = find(Cgroups == 2);
            
            spike_times_SU = {};
            for id = 1:length(SU)
                spike_times_SU{id} = spike_times_all(find(spike_clusters == Cids(SU(id)))); %/obj.params.fs;
%                 if obj.params.session_id > 1
%                     spike_times_SU{id} = spike_times_SU{id}(find(spike_times_SU{id} < (obj.seg_length(obj.params.session_id)-40) & (obj.seg_length(obj.params.session_id-1)+20)< spike_times_SU{id}));
%                     spike_times_SU{id} = spike_times_SU{id} - obj.seg_length(obj.params.session_id-1);
%                 else
%                     spike_times_SU{id} = spike_times_SU{id}(find(spike_times_SU{id} < obj.seg_length(obj.params.session_id)-60));
%                 end
                obj.waveforms.raw{id} = zeros(32,60,length(spike_times_SU{id}));
            end
            
            obj.SUU = [obj.SUU SU]; %NEW LINE
            
            obj.SU_good = ones(length(SU),1);
            
            if strcmp(obj.params.file_type,'BR') %blackrock systems
                voltage_amp = 0.1;
            elseif strcmp(obj.params.file_type,'TDT') %TDT systems
                voltage_amp = 1;
            else
                voltage_amp = 0.195;
            end
            
            
            
            Nbuff = 32*1024;
            Nbatch = length(buff(1,:))/Nbuff;
            for batch = 1:Nbatch
                ipoint = (batch-1)*Nbuff +1;
                % filtering
                if batch<Nbatch
                    datr = filter(b,a,single(gpuArray(buff(:,ipoint:ipoint+Nbuff-1+60))).');
                else
                    datr = filter(b,a,single(gpuArray(buff(:,ipoint:ipoint+Nbuff-1))).');
                end
                datr = flipud(datr);
                datr = filter(b,a,datr);
                datr = flipud(datr);
                % common median referencing
                datr = datr-median(datr,2);
                %                 filtData(:,ipoint:ipoint+Nbuff-1) = gather_try(int16(datr.'));
                if mod(batch,100) ==0
                    lineLength = fprintf('Time %3.0fs. batch %3.0f/ %3.0f.  \n', toc,batch,Nbatch);
                end
                %                 datr = datr.';
                %                 datr = datr(obj.params.chanMap,:);
                filtData = gather_try(double(datr.'));
                filtData = filtData(obj.params.chanMap,:);
                for id = 1:length(SU)
                    
                    
                    %                     obj.waveforms.raw{id} = zeros(64,60,length(spike_times_SU{id}));
                    %                     obj.waveforms.raw{id} = [];
                    %                     p =0;
                    for t = 1:length(spike_times_SU{id})
                        if batch<Nbatch
                            if spike_times_SU{id}(t,1)-19 > ipoint && spike_times_SU{id}(t,1)+40 < ipoint+Nbuff-1+60
                                %                                 p = p+1;
                                %                                 disp(p)
                                obj.waveforms.raw{id}(:,:,t) = voltage_amp*filtData(:, ...
                                    round(spike_times_SU{id}(t,1)-ipoint-19):round(spike_times_SU{id}(t,1)-ipoint+40));
                            end
                        else
                            if spike_times_SU{id}(t,1)-19 > ipoint && spike_times_SU{id}(t,1)+40 < ipoint+Nbuff-1
                                obj.waveforms.raw{id}(:,:,t) = voltage_amp*filtData(:, ...
                                    round(spike_times_SU{id}(t,1)-ipoint-19):round(spike_times_SU{id}(t,1)-ipoint+40));
                            end
                        end
                        %
                        
                    end
                    
                    
                end
            end
            
            
            fprintf('Time %3.0fs. Preprocessing complete! \n', toc);
            clear buff;
            
            
            
            %// removing bad clusters and extracting templates;
            for id = 1:length(SU)
                if size(obj.waveforms.raw{id},3) < 10 % if number of spikes is less than 10, remove cluster
                    obj.SU_good(id) = 0;
                else
                    % extracting templates via SVD
                    obj.templates{id}.data = zeros(obj.params.Nb_ch,60);
                    %                  obj.templates{id}.best_ch = [];
                    for ch = 1:obj.params.Nb_ch
                        temp = [];
                        temp(:,:) = obj.waveforms.raw{id}(ch,:,:);
                        [Y,Sig,X] = svd(temp,'econ');
                        %                 sig = diag(Sig);%figure; semilogy(sig(sig>1),'kx-')
                        k = 1:3;
                        P = Y(:,k)*Sig(k,k)*X(:,k)';
                        obj.templates{id}.data(ch,:) = mean(P,2).';
                    end
                    [~, obj.templates{id}.best_ch] = max(mean(abs(obj.templates{id}.data(:,10:50)),2));
                    %                 plot(obj.templates{id}(2,:))
                    
                    obj.waveforms.mean{id}(1,:) = mean(obj.waveforms.raw{id}(obj.templates{id}.best_ch,:,:),3);
                    obj.waveforms.std{id}(1,:) = std(obj.waveforms.raw{id}(obj.templates{id}.best_ch,:,:),[],3);
                    peak_to_peak = max(obj.waveforms.mean{id})-min(obj.waveforms.mean{id});
                    noise = mean(obj.waveforms.std{id}(1:10));
                    obj.waveforms.SNR(id) = 20.*log10(peak_to_peak./noise);
                    
                end
                obj.spike_times{id} =  spike_times_SU{id}./obj.params.fs;
            end
            
            %NEW CODE BEGINS (copied part from Single_units code)
            for id = 1:length(SU)
                
                if obj.SU_good(id) ==1
                    obj.good_ch=[obj.good_ch obj.templates{id}.best_ch];
                    obj.clusters=[obj.clusters Cids(SU(id))];
                    obj.idds = [obj.idds id];
                end
                
            end
            %NEW CODE ENDS
            
            fprintf('Time %3.0fs. Extracting waveforms... Done \n', toc)
            
            
        end
        %%
        
        function obj = save_Units(obj, save_dir,project_name,fname,animalname,experiment)
            global   Cchannels Cids  SU
            
            
            %save_dir = 'D:\Data\M12E\Units_training';
            %save_dir = 'C:\Users\Seth\Desktop\Saul\Dummy_Save_Folder';
            %obj.start_times = start_times; %NEW LINE
            %obj.end_times = end_times; %NEW LINE
            
            %             save_dir = a; %NEW
            %             save_dir = obj.params.savepath;
            if ~isfolder(save_dir)
                mkdir(save_dir);
            end
            unit_list_name = [project_name '_unit_list.mat'];
            if exist(fullfile(save_dir,filesep,unit_list_name))
                load(fullfile(save_dir,unit_list_name));
                if  ~isempty(find(strcmp(unit_list.data(:,5),fname)))
                    error('Units already saved')
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
            tic
            for id = 1:length(SU)
                if obj.SU_good(id) ==1
                    s_unit = {};
                    s_unit.id = 1;
                    s_unit.depth = [];
                    %                 s_unit.ch = CchanMap(Cchannels(SU(id)));
                    s_unit.templates = obj.templates{id};
                    s_unit.cid = Cids(SU(id));
                    s_unit.SU = SU;
                    s_unit.SU_good = obj.SU_good;
                    s_unit.SNR = obj.waveforms.SNR(id);
                    s_unit.waveforms = obj.waveforms.raw(id);
                    %                 s_unit.spiketimes = spike_times_all(find(spike_clusters == Cids(SU(id))))/obj.params.fs;
                    s_unit.spiketimes = obj.spike_times{id};
                    %                     s_unit.xbz_file_name = obj.params.xbz_file_name;
                    s_unit.stim_data = obj.Eventdata;
                    s_unit.amplitude = max(abs(obj.waveforms.mean{id}));
                    %                     s_unit.start_times = obj.start_times;
                    %                     s_unit.end_times =obj.end_times;
                    unitname= [project_name 'u00001']; %03/03/20 changed from u001 to u00001
                    if isempty(unit_list.data)
                        save(fullfile(save_dir,unitname),'s_unit')
                        unit_list.data{1,1} = 1;
                        unit_list.data{1,2} = s_unit.SNR;
                        unit_list.data{1,3} = s_unit.amplitude;
                        %                         unit_list.data{1,4} = s_unit.xbz_file_name;
                        unit_list.data{1,4} = s_unit.cid;
                        unit_list.data{1,5} = fname;
                        unit_list.data{size((unit_list.data),1),6} = animalname;
                        unit_list.data{size((unit_list.data),1),7} = experiment;
                    else
                        unitname = [unitname(1:end-length(num2str(size((unit_list.data),1)+1))) num2str(size((unit_list.data),1)+1)];
                        save(fullfile(save_dir,unitname),'s_unit')
                        unit_list.data{size((unit_list.data),1)+1,1} = size((unit_list.data),1)+1;
                        unit_list.data{size((unit_list.data),1),2} = s_unit.SNR;
                        unit_list.data{size((unit_list.data),1),3} = s_unit.amplitude;
                        %                         unit_list.data{size((unit_list.data),1),4} = s_unit.xbz_file_name;
                        unit_list.data{size((unit_list.data),1),4} = s_unit.cid;
                        unit_list.data{size((unit_list.data),1),5} = fname;
                        unit_list.data{size((unit_list.data),1),6} = animalname;
                        unit_list.data{size((unit_list.data),1),7} = experiment;
                    end
                    save(fullfile(save_dir,unit_list_name),'unit_list')
                end
                fprintf('Time %3.0fs. unit %3.0f/ %3.0f.  \n', toc,id,length(SU));
                
            end
            disp('Units saved')
            
        end
        
        %%
        
        function [obj, out] = Single_units(obj)
            
            %{
            This function combines spike time information with stimulus
            information. It will read from xbz files that are outputs of
            xb3 during recording sessions.
            %}
            
            
            global SU  Cids
            
            % .events files have timestamps of when each stimuli was
            % played. This needs to be synchronized with the timestamps
            % determined by xb3.
            
            [~, timestamps, info] = load_open_ephys_data(fullfile(obj.params.fpath, obj.params.session_name, 'all_channels.events'));
            
            %             start_stim_times = timestamps(find(info.eventId ==1))-29.9349;
            start_stim_times = timestamps(find(info.eventId ==1))-obj.params.rec_start_time; %$-2.2870; 
            end_stim_times = timestamps(find(info.eventId ==0))-obj.params.rec_start_time; %-2.2870;
%             start_stim_times = start_stim_times(2:end);
%             end_stim_times = end_stim_times(3:end);
            obj.start_times = start_stim_times;
            obj.end_times = end_stim_times;
            obj.x = eval(obj.params.xbz_file_name);

            PreStim = obj.x.pre_stimulus_record_time*1e-3; %s
            PostStim = obj.x.post_stimulus_record_time*1e-3; %s
            StimDur = obj.x.stimulus_ch1(:,5)*1e-3;
            
            data_new = obj.x.data(find(obj.x.data(:,3) == 1 & obj.x.data(:,4) == -1),:);
            
            
%             data_new = stim_info;
            
            nreps = obj.x.stimulus_ch1(1,4);
            nStim = max(obj.x.stimulus_ch1(:,1));
            %             nreps = 10;
            
            %
            % nStim = nStim + max(x2.stimulus_ch1(:,1));
            %
            TotalReps = nStim*nreps;
            false_start = length(start_stim_times)-TotalReps;
            start_stim_times = start_stim_times(false_start+1:end);
            end_stim_times = end_stim_times(false_start+1:end); %(false_start+1:end);
            
%             start_stim_times = start_stim_times(1:end-false_start);
%             end_stim_times = end_stim_times(1:end-false_start); %(false_start+1:end);
            
            
            raster.stim = {};
            raster.rep = {};
            raster.spikes = {};
            spikes_pooled = {};
            spike_timesSU = {};
            rate.stim = {};
            rate.pre = {};
            %             for id = 1:length(SU)
            %             obj.spike_times{id} = obj.spike_times{id}.*30000;
            %             end
            for id = 1:length(SU)
            
                if obj.SU_good(id) ==1
                    %obj.good_ch=[obj.good_ch obj.templates{id}.best_ch]; %NEW ENTRY
                    %obj.clusters=[obj.clusters Cids(SU(id))]; %NEW ENTRY, GOOD ONE
                    %obj.idds = [obj.idds id]; %NEW ENTRY, GOOD ONE
                    
                    %                     obj.spike_times{id} = obj.spike_times{id}./30000;
                    
                    raster.stim{id} = [];
                    raster.rep{id} = [];
                    raster.spikes{id} = [];
                    spikes_pooled{id} = [];
                    %                 rate_stim{id} = [];
                    %                 if obj.manual == 0
                    %                     spike_timesSU{id} = spike_times_all(find(spike_clusters == Cids(SU(id))))/obj.params.fs;
                    %                 else
                    %                     spike_timesSU{id} = obj.spike_times{id};
                    %                 end
                    for rep = 1:TotalReps
                        %for raster
                        %                     if id ==4
                        %                         spike_timesSU{id} = a.spike_times{4};
                        %                     end
                        spikes1 = obj.spike_times{id}(find(obj.spike_times{id}>=start_stim_times(rep)-PreStim & ...
                            obj.spike_times{id}<=end_stim_times(rep)+ PostStim)).';
                        spikes1 = spikes1 - start_stim_times(rep);
                        spikes_pooled{id} = [spikes_pooled{id} spikes1];
                        raster.stim{id} = [raster.stim{id} data_new(rep,1)*ones(size(spikes1))];
                        raster.rep{id} = [raster.rep{id} data_new(rep,2)*ones(size(spikes1))];
                        raster.spikes{id} = [raster.spikes{id} spikes1];
                        
                        %for rate
                        spikes2 = obj.spike_times{id}(find(obj.spike_times{id}>=start_stim_times(rep) & ...
                            obj.spike_times{id}<=end_stim_times(rep))).';
                        rate.stim{id}(data_new(rep,1),data_new(rep,2)) = length(spikes2)/StimDur(data_new(rep,1));
                        spikes3 = obj.spike_times{id}(find(obj.spike_times{id}<=start_stim_times(rep) & ...
                            obj.spike_times{id}>=start_stim_times(rep)-PreStim)).' ;
%                         spikes4 = obj.spike_times{id}(find(obj.spike_times{id}<=end_stim_times(rep)+ PostStim & ...
%                             obj.spike_times{id}>=end_stim_times(rep))).' ;
                        rate.pre{id}(data_new(rep,1),data_new(rep,2)) = length(spikes3)/PreStim;
%                         rate.post{id}(data_new(rep,1),data_new(rep,2)) = length(spikes4)/PostStim
                    end
                    SUrate{id}.mean = mean(rate.stim{id},2); %-mean(rate_pre{id},2);
                    SUrate{id}.error = std(rate.stim{id},1,2)/sqrt(nreps);
                    SUrate{id}.spont = mean(mean(rate.pre{id}));
%                     SUrate{id}.post = mean(rate.stim{id},2);
                end
            end
            
            out.SUrate = SUrate;
            out.spikes1 = spikes1;
            out.spikes2 = spikes2;
            out.spikes3 = spikes3;
            out.rate = rate;
            out.raster = raster;
            out.StimDur = StimDur;
            out.nreps = nreps;
            out.TotalReps = TotalReps;
            out.PreStim = PreStim;
            out.PostStim = PostStim;
            

        end
        
        
    end
end




















