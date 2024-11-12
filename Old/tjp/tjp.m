classdef tjp < handle
    
    properties(Constant)
    end
    
    properties
        params = eval('parameters');
        raw_data = {};
        filt_data = {};
        waveforms = {};
        analysis_code = [];
        analysis_type = {};
        plot_type = {};
        spike_times = {};
        manual = 0; % variable to indicate whether the clusters were sorted manually
        template = [];
    end
    
    
    methods(Static)
        
        
        
    end
    
    methods
        function obj = initialize(obj)
            global spike_times_all spike_clusters Cchannels Cids Cgroups new_channels CchanMap
            addpath(fullfile('C:\Data\OpenEphys', filesep, obj.params.animal_name)); 
            spike_times_all = readNPY(fullfile(obj.params.fpath, 'spike_times.npy'));
            spike_times_all = double(spike_times_all);
            spike_clusters = readNPY(fullfile(obj.params.fpath, 'spike_clusters.npy'));
            [Cids, Cgroups] = readClusterGroupsCSV(fullfile(obj.params.fpath, 'cluster_group.tsv'));
            spike_templates = readNPY(fullfile(obj.params.fpath, 'spike_templates.npy'));
            new_channels = readNPY(fullfile(obj.params.fpath, 'channel_map.npy'));
            new_channels = new_channels+1;
            Cchannels = [];
            templates = readNPY(fullfile(obj.params.fpath, 'templates.npy'));
            %// Find channel information for each cluster
            for Cid = Cids
                if Cid <= max(spike_templates)
                    [~, Ch] = max(mean(abs(templates(Cid+1,:,:))));
                else
                    Cid2 = spike_templates(find(spike_clusters == Cid,1));
                    [~, Ch] = max(mean(abs(templates(Cid2+1,:,:))));
                end
                Cchannels = [Cchannels Ch];
                
            end
            
            CchanMap = zeros(length(new_channels));
            for c = 1:length(new_channels)
                CchanMap(c) = find(obj.params.chanMap == new_channels(c));
            end
            
            Cchannels = Cchannels.';
            
        end
        
        function obj = load_data_OE(obj)
            
            data_all = [];
            timestamps1 = {};
            info1 = {};
            file_type = '100';
            tic
            fprintf('Time %3.0fs. Loading data... \n', toc);
            parfor ch = 1:obj.params.Nb_ch
                %     disp(['loading data from channel ' num2str(ch)])
                [data1,timestamps1{ch},info1{ch}] = load_open_ephys_data_faster([obj.params.fpath filesep obj.params.file_type '_CH' num2str(ch) '.continuous' ]);
                data_all(ch,:) = data1;
            end
            obj.params.rec_start_time = timestamps1{1}(1);
            obj.raw_data = data_all;
            fprintf('Time %3.0fs. Loading data... Done \n', toc);
        end
        
        function obj = load_LFP(obj)
            tic
            %// filtering data
            fprintf('Time %3.0fs. filtering data... \n', toc);
            
            gpu_data = gpuArray(int16(obj.raw_data));
            filtData = single([]);
            % [b,a] = butter(4, [0.0244 0.6104],'bandpass');
            [b,a] = butter(4, 0.05,'low');
            Nbuff = 32*32;
            Nbatch = length(obj.raw_data(1,:))/Nbuff;
            for batch = 1:Nbatch
                ipoint = (batch-1)*Nbuff +1;
                datr = filter(b,a,single(gpu_data(:,ipoint:ipoint+Nbuff-1)).');
                datr = flipud(datr);
                datr = filter(b,a,datr);
                datr = flipud(datr);
                filtData(:,ipoint:ipoint+Nbuff-1) = gather_try(int16(datr.'));
            end
            
            
            filtData2 = zeros(64,length(filtData));
            for ch = 1:64
                filtData2(ch,:) = filtData(obj.params.chanMap(ch),:);
            end
            
            obj.filt_data = filtData2;
            fprintf('Time %3.0fs. filtering data... Done \n', toc);
            %         if plot_var == 1
            %             figure
            %             % for ch = 1:64
            %             %     filtData2(ch,:) = filtData3(chanMap(ch),:);
            %             % end
            %         end
        end
                
        function obj = CSD_analysis(obj)
%         csdout =  CSD(obj.filt_data(1:4:end,:).',30000,2e-5,'inverse',5e-4);
        
        [data, timestamps, info] = load_open_ephys_data(fullfile(obj.params.fpath, 'all_channels.events'));
        
        start_stim_times = (timestamps(find(info.eventId ==1))-obj.params.rec_start_time)*obj.params.fs;
        end_stim_times = (timestamps(find(info.eventId ==0))-obj.params.rec_start_time)*obj.params.fs;
        
        
        
        x = eval(obj.params.xbz_file_name);        
        PreStim = x.pre_stimulus_record_time*1e-3; %s
        PostStim = x.post_stimulus_record_time*1e-3; %s
        StimDur = x.stimulus_ch1(:,5)*1e-3;        
        stim_info = x.data(find(x.data(:,3) == 1 & x.data(:,4) == -1),:);          
        data_new = stim_info;
        Stim_label = x.stimulus_ch1(:,8);    
        nreps = x.stimulus_ch1(1,4);
        nStim = max(x.stimulus_ch1(:,1));
        TotalReps = nStim*nreps;
        false_start = length(start_stim_times)-TotalReps;
        start_stim_times = start_stim_times(false_start+1:end);
        end_stim_times = end_stim_times(false_start+1:end);
        LFP = {};
        for rep = 1:TotalReps
            LFP{data_new(rep,1)}(:,:,data_new(rep,2)) = obj.filt_data(1:4:end, ...
               round( start_stim_times(rep)-PreStim*obj.params.fs):round(start_stim_times(rep)+(StimDur(data_new(rep,2))+PostStim)*obj.params.fs)); 
        end
        for st = 1:length(LFP)
            LFP_mean{st} = mean(LFP{st},3);
            csdout{st} = CSD(LFP_mean{st}.',30000,2E-5,0,'inverse',5E-4);
            drawnow
            C{st} = smoothdata(csdout{st},1,'gaussian',10);
        end
        
        for st =1:length(LFP)
            figure(st)
            subplot(1,2,2)
            M = max(max(abs(C{st}))); % abosolute maximum CSD, for the colormap scale
            clims = [-M M]; % gives the upper and lower limit for the colormap
            ylabel('Electrode');
            xlabel('Time (ms)');
            title('CSD (\color{red}sink, \color{blue}source\color{black})');
            imagesc(C{st},clims); % CSD as heatmap
            colormap(flipud(jet)); % blue = source; red = sink
            colorbar('SouthOutside');
            
            subplot(1,2,1)
            for ch = 1:16
                plot(LFP_mean{st}(ch,:)+ch*50)
                hold on
            end
            drawnow
            title(['Stim =' num2str(Stim_label(st))])
        end
        
        
            
            
            
        end
               
        function obj = extract_units(obj)
            tic
            global spike_times_all spike_clusters Cchannels Cids Cgroups new_channels
            %// filtering data
            fprintf('Time %3.0fs. filtering data... \n', toc);
            
            gpu_data = gpuArray(int16(obj.raw_data));
            filtData = single([]);
            [b,a] = butter(4, [0.0244 0.6104],'bandpass');
            % [b,a] = butter(4, 500/fs,'low');
            Nbuff = 32*32;
            Nbatch = length(obj.raw_data(1,:))/Nbuff;
            for batch = 1:Nbatch
                ipoint = (batch-1)*Nbuff +1;
                datr = filter(b,a,single(gpu_data(:,ipoint:ipoint+Nbuff-1)).');
                datr = flipud(datr);
                datr = filter(b,a,datr);
                datr = flipud(datr);
                filtData(:,ipoint:ipoint+Nbuff-1) = gather_try(int16(datr.'));
            end
            
            %// common median referencing
            fprintf('Time %3.0fs. Common median referencing... \n', toc);
            CommonMedian = median(filtData);
            st_dev = zeros(1,obj.params.Nb_ch);
            for ch = 1:obj.params.Nb_ch
                filtData(ch,:) = filtData(ch,:)-CommonMedian;
                st_dev(ch) = median(abs(filtData(ch,:))/0.6745);
            end
            
            fprintf('Time %3.0fs. Preprocessing complete! \n', toc);
            
            filtData = double(filtData);
            
            
            fprintf('Time %3.0fs. Extracting waveforms... \n', toc);
            %// extracting waveforms from processed data
            obj.waveforms.raw = {};
            obj.waveforms.mean = {};
            obj.waveforms.std = {};
            obj.waveforms.SNR = [];
            
            channels  = new_channels(Cchannels);
            
            for tid = 1:length(Cids)
                id = Cids(tid);
                spike_times{tid} = spike_times_all(find(spike_clusters == id)); %*fs;
                obj.waveforms.raw{tid} = zeros(length(spike_times{tid}),60);
                for t = 1:length(spike_times{tid})
                    obj.waveforms.raw{tid}(t,:) = filtData(channels(tid), ...
                        round(spike_times{tid}(t,1)-19):round(spike_times{tid}(t,1)+40,1)); %
                end
                obj.waveforms.mean{tid} = mean(obj.waveforms.raw{tid},1);
                obj.waveforms.std{tid} = std(obj.waveforms.raw{tid},[],1);
                peak_to_peak = max(obj.waveforms.mean{tid})-min(obj.waveforms.mean{tid});
                noise = mean(obj.waveforms.std{tid}(1:10));
                obj.waveforms.SNR(tid) = 20.*log10(peak_to_peak./noise);
            end
            
%             obj.waveforms = waveforms;
            fprintf('Time %3.0fs. Extracting waveforms... Done \n', toc)
            
        end
        
        
        function obj = save_Units(obj)
            global spike_times_all spike_clusters Cchannels Cids Cgroups SU CchanMap
            
%             [data, timestamps, info] = load_open_ephys_data(fullfile(obj.params.fpath, 'all_channels.events'));
%             
%             start_stim_times = timestamps(find(info.eventId ==1))-obj.params.rec_start_time;
%             end_stim_times = timestamps(find(info.eventId ==0))-obj.params.rec_start_time;
            
            
            save_dir = 'D:\Data\M12E\Units';
            if exist(fullfile(save_dir,filesep,'M12E_unit_list.mat'))
                load(fullfile(save_dir,'M12E_unit_list.mat'));
                
                if  ~isempty(find(strcmp(M12E_unit_list.data(:,4),obj.params.xbz_file_name)))
                    error('Units already saved')
                end
            else
                M12E_unit_list = {};
                M12E_unit_list.tags{1} = 'id';
                M12E_unit_list.tags{2} = 'SNR';
                M12E_unit_list.tags{3} = 'Amplitude';
                M12E_unit_list.tags{4} = 'xbz_filename';
                M12E_unit_list.data = {};
                save(fullfile(save_dir,'M12E_unit_list'),'M12E_unit_list');
            end
            
            for id = 1:length(SU)
                s_unit = {};
                s_unit.id = 1;
                s_unit.depth = [];
                s_unit.ch = CchanMap(Cchannels(SU(id)));
                s_unit.SNR = obj.waveforms.SNR(SU(id));
                s_unit.waveforms = obj.waveforms.raw(SU(id));
%                 s_unit.spiketimes = spike_times_all(find(spike_clusters == Cids(SU(id))))/obj.params.fs;
                s_unit.spiketimes = obj.spike_times{id};
                s_unit.xbz_file_name = obj.params.xbz_file_name;
                s_unit.amplitude = max(abs(obj.waveforms.mean{SU(id)}));
                unitname= 'M12Eu001';
                if isempty(M12E_unit_list.data)
                    save(fullfile(save_dir,unitname),'s_unit')
                    M12E_unit_list.data{1,1} = 1;
                    M12E_unit_list.data{1,2} = s_unit.SNR;
                    M12E_unit_list.data{1,3} = s_unit.amplitude;
                    M12E_unit_list.data{1,4} = s_unit.xbz_file_name;
                else
                    unitname = [unitname(1:end-length(num2str(size((M12E_unit_list.data),1)+1))) num2str(size((M12E_unit_list.data),1)+1)];
                    save(fullfile(save_dir,unitname),'s_unit')
                    M12E_unit_list.data{size((M12E_unit_list.data),1)+1,1} = size((M12E_unit_list.data),1)+1;
                    M12E_unit_list.data{size((M12E_unit_list.data),1),2} = s_unit.SNR;
                    M12E_unit_list.data{size((M12E_unit_list.data),1),3} = s_unit.amplitude;
                    M12E_unit_list.data{size((M12E_unit_list.data),1),4} = s_unit.xbz_file_name;
                    
                end
                save(fullfile(save_dir,'M12E_unit_list'),'M12E_unit_list')
            end
        
        end
        
        function obj = Single_Units(obj)
            global spike_times_all spike_clusters Cchannels Cids Cgroups SU CchanMap
            
            
            %[data1,timestamps1,info] = load_open_ephys_data_faster([fpath filesep file_type '_CH' num2str(ch) '.continuous' ]);
            %
            [data, timestamps, info] = load_open_ephys_data(fullfile(obj.params.fpath, 'all_channels.events'));
            
            start_stim_times = timestamps(find(info.eventId ==1))-obj.params.rec_start_time;
            end_stim_times = timestamps(find(info.eventId ==0))-obj.params.rec_start_time;
            
            if obj.manual == 0 
                SU = find(Cgroups == 2);
                SU = intersect(SU,find(obj.waveforms.SNR>10));
            end
            
            MU = find(Cgroups == 1);
            
            x = eval(obj.params.xbz_file_name);
            obj.analysis_code = x.analysis_code;
            obj.analysis_type = x.analysis_type;
            %             switch obj.analysis_code
            if obj.analysis_code < 99
                if length(unique(x.stimulus_ch1(:,3)))~= 1 && length(unique(x.stimulus_ch1(:,8)))~= 1
                    obj.analysis_type = 'FRA';
                else
                    obj.plot_type  = 'Tuning';
                end
            elseif obj.analysis_code >= 2300
                obj.analysis_type = 'VT';
                obj.plot_type = 'VT';
                
            elseif obj.analysis_code >= 100 && obj.analysis_code < 200 
                obj.analysis_type = 'User';
                obj.plot_type = 'User';
            end
            
            
            
            PreStim = x.pre_stimulus_record_time*1e-3; %s
            PostStim = x.post_stimulus_record_time*1e-3; %s
            StimDur = x.stimulus_ch1(:,5)*1e-3;
            
            stim_info = x.data(find(x.data(:,3) == 1 & x.data(:,4) == -1),:);
            
            
            data_new = stim_info;
            
            nreps = x.stimulus_ch1(1,4);
            nStim = max(x.stimulus_ch1(:,1));
            %             nreps = 10;
            
            %
            % nStim = nStim + max(x2.stimulus_ch1(:,1));
            %
            TotalReps = nStim*nreps;
            false_start = length(start_stim_times)-TotalReps;
            start_stim_times = start_stim_times(false_start+1:end);
            end_stim_times = end_stim_times(false_start+1:end);
            
            
            raster.stim = {};
            raster.rep = {};
            raster.spikes = {};
            spikes_pooled = {};
            spike_timesSU = {};
            rate.stim = {};
            rate.pre = {};
            
            for id = 1:length(SU)
                raster.stim{id} = [];
                raster.rep{id} = [];
                raster.spikes{id} = [];
                spikes_pooled{id} = [];
%                 rate_stim{id} = [];
                if obj.manual == 0         
                    spike_timesSU{id} = spike_times_all(find(spike_clusters == Cids(SU(id))))/obj.params.fs;
                else
                    spike_timesSU{id} = obj.spike_times{id};
                end
                for rep = 1:TotalReps
                    %for raster
%                     if id ==4
%                         spike_timesSU{id} = a.spike_times{4};
%                     end
                    spikes1 = spike_timesSU{id}(find(spike_timesSU{id}>=start_stim_times(rep)-PreStim & ...
                        spike_timesSU{id}<=end_stim_times(rep)+ PostStim)).';
                    spikes1 = spikes1 - start_stim_times(rep);
                    spikes_pooled{id} = [spikes_pooled{id} spikes1];
                    raster.stim{id} = [raster.stim{id} data_new(rep,1)*ones(size(spikes1))];
                    raster.rep{id} = [raster.rep{id} data_new(rep,2)*ones(size(spikes1))];
                    raster.spikes{id} = [raster.spikes{id} spikes1];
                    
                    %for rate
                    spikes2 = spike_timesSU{id}(find(spike_timesSU{id}>=start_stim_times(rep) & ...
                        spike_timesSU{id}<=end_stim_times(rep))).';
                    rate.stim{id}(data_new(rep,1),data_new(rep,2)) = length(spikes2)/StimDur(data_new(rep,1));
                    spikes3 = spike_timesSU{id}(find(spike_timesSU{id}<=start_stim_times(rep) & ...
                        spike_timesSU{id}>=start_stim_times(rep)-PreStim)).' ;
                    rate.pre{id}(data_new(rep,1),data_new(rep,2)) = length(spikes3)/PreStim;
                end
            end
            if obj.manual == 0 
                obj.spike_times = spike_timesSU;
            end
            stim_label  = x.stimulus_ch1(:,8);
            
            SUrate = {};
            % figure
            
            
            for id = 1:length(SU)
                figure
                suptitle(['cluster ' num2str(Cids(SU(id))) ' channel ' num2str(CchanMap(Cchannels(SU(id)))) ' SNR: ' num2str(obj.waveforms.SNR(SU(id)))])
                set(gcf, 'Position', [400 300 1800 1000]);
                
                %// this should change depending on stim_set. make it modular?
                SUrate{id}.mean = mean(rate.stim{id},2); %-mean(rate_pre{id},2);
                SUrate{id}.error = std(rate.stim{id},1,2)/sqrt(nreps);
                SUrate{id}.spont = mean(mean(rate.pre{id}));
                
                
                
                %waveform
                subplot(2,3,4)
                time_step = [1:60]/obj.params.fs*1e3;
                Nb_wave = min(size(obj.waveforms.raw{SU(id)},1),100);
                idx = 1:size(obj.waveforms.raw{SU(id)},1);
                rand_select = randsample(idx,Nb_wave);
                for n = 1:length(rand_select)
                    
                    plot(time_step,obj.waveforms.raw{SU(id)}(rand_select(n),:))
                    hold on
                end
                plot(time_step,obj.waveforms.mean{SU(id)},'Linewidth',2,'Color','k')
                title('waveform')
                xlabel('time (ms)')
                ylabel('uV');
                hold off
                
                %// PCA
                [Wi,score,latent] = pca(obj.waveforms.raw{SU(id)});
                
                X = [score(:,1) score(:,2)];
                idx = kmeans(X,2);
                subplot(2,3,2)
                gscatter(X(:,1),X(:,2),idx,'rb','+o',5)
                
                title('PCA')
                
                
                %// ISI plot
                subplot(2,3,5)
                ISI = obj.spike_times{id}(2:end)-obj.spike_times{id}(1:end-1);
                
                edges = -3.5:0.25:2;
                [N,~] = histcounts(ISI,10.^edges);
                refract = sum(N(1:2))*100/sum(N);
                histogram(ISI,10.^edges)
                set(gca,'xscale','log')
                xticks([1E-4 1E-3 1E-2 1E-1 1 10]);
                xticklabels({0.0001, 0.001, 0.01, 0.1, 1, 10})
                xlabel('time (ms)')
                title([num2str(refract) '%']);
                
                %                 histogram
                if strcmp(obj.plot_type,'Tuning') == 1
                    
                    %// Tuning curve
                    subplot(2,3,1)
                    

                    errorbar(stim_label,SUrate{id}.mean,SUrate{id}.error,'LineWidth',2);
                    hold on
                    
                    plot(stim_label,ones(1,length(stim_label))*SUrate{id}.spont,'--k');
                    set(gca,'xscale','log')
                    
                    xticks(round((stim_label(1:4:end).')*1e2)*1e-2);
                    xtickangle(45)
                    xlabel('Hz')
                    ylabel('firing rate (spikes/s)')
                    title('tuning curve')
                    
                    
                    %Raster plot
                    subplot(2,3,[3 6])
                    for st = 1:length(StimDur)
                        rectangle('Position',[0 nreps*(st-1),StimDur(st) nreps],'FaceColor',[0.9,0.9,0.9],'EdgeColor','none')
                    end
                    hold on
                    plot(raster.spikes{id},nreps*(raster.stim{id}-1)+raster.rep{id},'k.','MarkerSize',15);
                    %     pause
                    xlabel('time (s)')
                    ylabel('reps')
                    yticks([1:nreps*2:TotalReps]+10)
                    yticks([1:nreps*2:TotalReps]+10)
                    stim_ticks = {};
                    if length(stim_label)>2
                        for stim = 1:length(stim_label)/2
                            stim_ticks{stim}=num2str(round(stim_label(stim*2)*10)/10);
                        end
                    else
                        for stim = 1:length(stim_label)
                            stim_ticks{stim}=num2str(round(stim_label(stim)*10)/10);
                        end
                    end
                    yticklabels(stim_ticks)
                    axis([-PreStim max(StimDur) + PostStim 0 TotalReps+1])
                    hold off
                    title('rasterplot')
                    
                    drawnow()
                   
                elseif strcmp(obj.plot_type,'VT') == 1
                    stim_number = x.stimulus_ch1(:,1);
                    stim_length = x.stimulus_ch1(:,5);
                    for i = 9:22
                        VT_para(i-8) = length(unique(x.stimulus_ch1(:,i)));
                    end
                    stim_label = x.stimulus_ch1(:,find(VT_para ~= 1)+8);
                    if obj.analysis_code == 2320 || obj.analysis_code == 2335 
                         stim_label = x.stimulus_ch1(:,10);
                    end
                    
                    
                    %// Tuning Curve
                    SUrate{id}.mean_fr = SUrate{id}.mean./stim_length;            
                    subplot(2,3,1)
                    errorbar(stim_label,SUrate{id}.mean,SUrate{id}.error,'LineWidth',2);
                    hold on
                    plot(stim_label,ones(1,length(stim_label))*SUrate{id}.spont,'--k');
%                     xticks(stim_number);
%                     xticklabels(stim_label)
                    xtickangle(45)
                    %                 xlabel('Hz')
                    ylabel('firing rate (spikes/s)')
                    title('tuning curve')
                    
                    %// Rasterplot
                    subplot(2,3,[3 6])
                    for st = 1:length(StimDur)
                        rectangle('Position',[0 nreps*(st-1),StimDur(st) nreps],'FaceColor',[0.9,0.9,0.9],'EdgeColor','none')
                    end
                    hold on
                    plot(raster.spikes{id},nreps*(raster.stim{id}-1)+raster.rep{id},'k.','MarkerSize',15);
                    %     pause
                    xlabel('time (s)')
                    ylabel('reps')
                    yticks([1:nreps*2:TotalReps]+10)
                    yticks([1:nreps*2:TotalReps]+10)
                    stim_ticks = {};
                    if length(stim_label)>2
                        for stim = 1:length(stim_label)/2
                            stim_ticks{stim}=num2str(round(stim_label(stim*2)*10)/10);
                        end
                    else
                        for stim = 1:length(stim_label)
                            stim_ticks{stim}=num2str(round(stim_label(stim)*10)/10);
                        end
                    end
                    yticklabels(stim_ticks)
                    axis([-PreStim max(StimDur) + PostStim 0 TotalReps+1])
                    hold off
                    title('rasterplot')
                    
                    drawnow()
                    
                elseif strcmp(obj.plot_type,'User') == 1
                    %// Tuning Curve
                    stim_label = {};
                    stim_number = x.stimulus_ch1(:,1);
                    stim_length = x.stimulus_ch1(:,5);
                    for st = 1: length(stim_number)
                        idx1 = strfind(x.user_stimulus_desc_ch1{st},': ')+2;
                        idx2 = strfind(x.user_stimulus_desc_ch1{st},'_m')-1;
                        idx3 = strfind(x.user_stimulus_desc_ch1{st},'rev');
                        stim_label{st} = x.user_stimulus_desc_ch1{st}(idx1:idx2);
                        if isempty(idx3) ~=1
                            stim_label{st} = [stim_label{st} '_rev'];
                        end
                        
                    end
                    SUrate{id}.mean_fr = SUrate{id}.mean./stim_length;
                    subplot(2,3,1)
                    errorbar(stim_number,SUrate{id}.mean,SUrate{id}.error,'LineWidth',2);
                    hold on
                    plot(stim_number,ones(1,length(stim_label))*SUrate{id}.spont,'--k');
                    xticks(stim_number);
                    xticklabels(stim_label)
                    xtickangle(45)
                    %                 xlabel('Hz')
                    ylabel('firing rate (spikes/s)')
                    title('tuning curve')
                    
                    
                    
                    %// Rasterplot
                    
                    subplot(2,3,[3 6])
                    %                     area([0 StimDur StimDur 0],[0 TotalReps+5 0 TotalReps+5],'LineStyle','none','FaceColor',[.85 .85 1]);
                    for st = 1:length(StimDur)
                        rectangle('Position',[0 nreps*(st-1),StimDur(st) nreps],'FaceColor',[0.9,0.9,0.9],'EdgeColor','none')
                    end
                    hold on
                    plot(raster.spikes{id},nreps*(raster.stim{id}-1)+raster.rep{id},'k.','MarkerSize',15);
                    %     pause
                    xlabel('time (s)')
                    ylabel('reps')
                    yticks([1:nreps:TotalReps]+5)

                    yticklabels(stim_label)
                    axis([-PreStim max(StimDur) + PostStim 0 TotalReps+5])
                    hold off
                    title('rasterplot')
                    
                    drawnow()
                    
                end
                %             elseif strcmp(obj.plot_type,'FRA') == 1
                %             elseif strcmp(obj.plot_type,'User') == 1
                
            end
            
        end
        
        function obj = merge_split(obj,in,cid,k,varargin)
            
            
            global Cgroups Cids SU Cchannels
            

            
            switch in
                case 'split'
                    id = find(Cids == cid);
                    [~,score,~] = pca(obj.waveforms.raw{id} );
                    scatter(score(:,1),score(:,2));
                    X = [score(:,1) score(:,2)];
                    figure
                    subplot(2,3,1)
                    if ~isempty(varargin)
                        if strcmp(varargin{1},'gm')
                            gm = fitgmdist(X,k);
                            idx = cluster(gm,X);
                        end
                    else
                        idx = kmeans(X,k);
                        
                    end
                    gscatter(X(:,1),X(:,2),idx,'rb','+o',5)
                    drawnow
                    for s = 1:2
                        test{s} = randsample(find(idx ==s),min(length(find(idx ==s)),200));
                        subplot(2,3,s+1)
                         for t = 1:length(test{s})
                             plot(obj.waveforms.raw{id}(test{s}(t),:))
                             hold on
                         end
                         plot(mean(obj.waveforms.raw{id}(test{s}.',:)),'Linewidth',2,'Color','k')
%                          hold off
                         subplot(2,3,s+3)
                         spktimes{s} = obj.spike_times{find(SU == id)}(find(idx ==s));
                         ISI = spktimes{s}(2:end)-spktimes{s}(1:end-1);
                         
                         edges = -3.5:0.25:2;
                         [N,~] = histcounts(ISI,10.^edges);
                         refract = sum(N(1:2))*100/sum(N);
                         histogram(ISI,10.^edges)
                         set(gca,'xscale','log')
                         xticks([1E-4 1E-3 1E-2 1E-1 1 10]);
                         xticklabels({0.0001, 0.001, 0.01, 0.1, 1, 10})
                         xlabel('time (ms)')
                         title([num2str(refract) '%']);
                    end
                     
                    prompt = 'choose between "1","2","both","cancel"';
                    x = input(prompt);
                    switch x
                        case '1'
                            obj.spike_times{find(SU == id)} = spktimes{1};
                            obj.waveforms.raw{id} = obj.waveforms.raw{id}(find(idx ==1),:);
                            obj.manual = 1;
                        case '2'
                            obj.spike_times{find(SU == id)} = spktimes{2};
                            obj.waveforms.raw{id} = obj.waveforms.raw{id}(find(idx ==2),:);
                            obj.manual = 1;
                        case 'both'
                            
                                
                            obj.waveforms.raw{length(Cids)+1} = [];
                            obj.waveforms.raw{length(Cids)+1} = obj.waveforms.raw{id}(find(idx ==2),:);
                            obj.spike_times{length(SU)+1} = [];
                            obj.spike_times{length(SU)+1} = spktimes{2};
                            
                            obj.waveforms.mean{length(Cids)+1} = mean(obj.waveforms.raw{id}(find(idx ==2),:),1);
                            obj.waveforms.std{length(Cids)+1} = std(obj.waveforms.raw{id}(find(idx ==2),:),[],1);
                            
                            peak_to_peak = max(obj.waveforms.mean{length(Cids)+1})-min(obj.waveforms.mean{length(Cids)+1});
                            noise = mean(obj.waveforms.std{length(Cids)+1}(1:10));
                            
                            obj.waveforms.SNR(length(Cids)+1) = 20.*log10(peak_to_peak./noise);
                            
                            obj.spike_times{find(SU == id)} = spktimes{1};
                            obj.waveforms.raw{id} = obj.waveforms.raw{id}(find(idx ==1),:);
                            SU = [SU length(Cids)+1];
                            
%                             for t = 1:length(cid)
%                                 id = find(Cids == cid(t));
%                                 obj.spike_times(find(SU == id)) = [];
%                                 SU(find(SU == id)) =[];
%                             end
                            obj.manual = 1;
                            Cids = [Cids max(Cids)+1];
                            Cchannels = [Cchannels; Cchannels(find(Cids == cid(1)))];
                            
                            
                            
                            
                        case 'cancel'
                    end
                    
                case 'merge'
                    
                    %                                     SU = find(Cgroups == 2);
                    
                    temp_waveforms = [];
                    temp_spike_times = [];
                    for t = 1:length(cid)
                        id = find(Cids == cid(t));
                        temp_waveforms = [temp_waveforms; obj.waveforms.raw{id}];
                        temp_spike_times = [temp_spike_times; obj.spike_times{find(SU == id)}];
                    end
                    
                    [~,score,~] = pca(temp_waveforms );
                    scatter(score(:,1),score(:,2));
                    X = [score(:,1) score(:,2)];
                    figure
                    subplot(1,3,1)
                    if ~isempty(varargin)
                        if strcmp(varargin{1},'gm')
                            gm = fitgmdist(X,k);
                            idx = cluster(gm,X);
                        end
                    else
                        idx = kmeans(X,k);
                    end
                    gscatter(X(:,1),X(:,2),idx,'rb','+o',5)
                    drawnow
                    idxs = 1:length(temp_spike_times);
                    test = randsample(idxs,min(length(idxs),200));
                    subplot(1,3,2)
                    for t = 1:length(test)
                        plot(temp_waveforms(test(t),:))
                        hold on
                    end
                    temp_waveforms_mean = mean(temp_waveforms,1);
                    temp_waveforms_std = std(temp_waveforms,[],1);
                    peak_to_peak = max(temp_waveforms_mean)-min(temp_waveforms_mean);
                    noise = mean(temp_waveforms_std(1:10));
                    temp_waveforms_SNR = 20.*log10(peak_to_peak./noise);
                    plot(temp_waveforms_mean ,'Linewidth',2,'Color','k')
                    title(['SNR = ' num2str(temp_waveforms_SNR)])
                    %                          hold off
                    subplot(1,3,3)
                    spktimes = temp_spike_times(idxs);
                    ISI = spktimes(2:end)-spktimes(1:end-1);
                    
                    edges = -3.5:0.25:2;
                    [N,~] = histcounts(ISI,10.^edges);
                    refract = sum(N(1:2))*100/sum(N);
                    histogram(ISI,10.^edges)
                    set(gca,'xscale','log')
                    xticks([1E-4 1E-3 1E-2 1E-1 1 10]);
                    xticklabels({0.0001, 0.001, 0.01, 0.1, 1, 10})
                    xlabel('time (ms)')
                    title([num2str(refract) '%']);
                    
                    prompt = 'choose between "merge",cancel"';
                    x = input(prompt);
                    
                    switch x
                        case 'merge'
                           
                            obj.waveforms.raw{length(Cids)+1} = [];
                            obj.spike_times{length(SU)+1} = [];
                            for t = 1:length(cid)
                                id = find(Cids == cid(t));
                                obj.waveforms.raw{length(Cids)+1} = [obj.waveforms.raw{length(Cids)+1}; obj.waveforms.raw{id}];
                                obj.spike_times{length(SU)+1} = [obj.spike_times{length(SU)+1}; obj.spike_times{find(SU == id)}];
                            end
                            obj.waveforms.mean{length(Cids)+1} = temp_waveforms_mean;
                           obj.waveforms.std{length(Cids)+1} = temp_waveforms_std;
                            obj.waveforms.SNR(length(Cids)+1) = temp_waveforms_SNR;
                            SU = [SU length(Cids)+1];
                            
                            for t = 1:length(cid)
                                id = find(Cids == cid(t));
                                obj.spike_times(find(SU == id)) = [];
                                SU(find(SU == id)) =[];
                            end
                            obj.manual = 1;
                            Cids = [Cids max(Cids)+1];
                            Cchannels = [Cchannels; Cchannels(find(Cids == cid(1)))];
                        case 'cancel'
                    end
                        
            end
            %% PCA and clustering
%             [Wi,score,latent] = pca(obj.waveforms.raw{cluster} );
%             scatter(score(:,1),score(:,2));
%             
%             X = [score(:,1) score(:,2)];
%             idx = kmeans(X,3);
%             
%             
%             figure
%             gscatter(score(:,1),score(:,2),idx)

        end
        
        function obj = merge_more(obj)
            global Cgroups Cids SU Cchannels CchanMap
            for id = 1:length(SU)
                [Y,Sig,X] = svd(obj.waveforms.raw{SU(id)},'econ');
%                 sig = diag(Sig);%figure; semilogy(sig(sig>1),'kx-')
                k = 1:3;
                P = Y(:,k)*Sig(k,k)*X(:,k)';
                obj.template(id,:) = mean(P,1).';               
            end

            %// calculate distance between templates and channels
            
            chanmap = reshape(1:64,[4,16]);
            d_chan = zeros(length(SU));
            d_template = zeros(length(SU));
            
            for i = 1:length(SU)
                for j = 1:length(SU)
                    
                    [y_A,z_A] = find(chanmap == CchanMap(Cchannels(SU(i))));
                    [y_B,z_B] = find(chanmap == CchanMap(Cchannels(SU(j))));
                    d_chan(j,i) = sum([y_A-y_B,z_A-z_B].^2)^0.5;
                    d_template(j,i) = sqrt(sum((obj.template(i,:)-obj.template(j,:)).^2));
                end
            end
            d_chan = d_chan/max(max(d_chan));
            d_template = d_template/max(max(d_template));
            d_final = d_chan.*d_template;
            figure
%             imagesc(-d_final)
            
            
            % plot nearest neighbor and choose whether to merge or not
            j_list = [];
            merge_list =[];
            for i = 1:length(SU)
                d_final(i,i) = 1;
                
                if min(d_final(i,:))<0.05 && isempty(intersect(i,j_list))
                    j = find(d_final(i,:) == min(d_final(i,:)));
                    if length(j)>1
                        j = j(1);
                    end
                    temp_spike_times = [obj.spike_times{i}; obj.spike_times{j}];
                    tmp_wave = [obj.waveforms.raw{SU(i)};obj.waveforms.raw{SU(j)}];
                    [~,score,~] = pca(tmp_wave);
%                     [~,score1,~] = pca(obj.waveforms.raw{SU(i)});
%                     [~,score2,~] = pca(obj.waveforms.raw{SU(j)});
                    figure
                    subplot(2,2,1)
                    set(gcf, 'Position', [200 500 800 600]);
                    for t = 1:min(size((obj.waveforms.raw{SU(i)}),1),200)
                        plot(obj.waveforms.raw{SU(i)}(t,:))
                        hold on
                    end
                    plot(mean(obj.waveforms.raw{SU(i)},1),'LineWidth',2,'Color','k')
                    axis([0 60 -300 200])
                    title(['cid = ' num2str(Cids(SU(i))) 'ch = ' num2str( CchanMap(Cchannels(SU(i))))])
                    subplot(2,2,2)
                    for t = 1:min(size((obj.waveforms.raw{SU(j)}),1),200)
                        plot(obj.waveforms.raw{SU(j)}(t,:))
                        hold on
                    end
                    plot(mean(obj.waveforms.raw{SU(j)},1),'LineWidth',2,'Color','k')
                    axis([0 60 -300 200])
                    title(['cid = ' num2str(Cids(SU(j))) 'ch = ' num2str( CchanMap(Cchannels(SU(j))))])
                    
                    subplot(2,2,3)
                    scatter(score(1:length(obj.waveforms.raw{SU(i)}),1),score(1:length(obj.waveforms.raw{SU(i)}),2))
                    hold on
                    scatter(score(length(obj.waveforms.raw{SU(i)})+1:end,1),score(length(obj.waveforms.raw{SU(i)})+1:end,2));
                    
                    subplot(2,2,4)
                    idxs = 1:length(temp_spike_times);
                    spktimes = temp_spike_times(idxs);
                    ISI = spktimes(2:end)-spktimes(1:end-1);
                    
                    edges = -3.5:0.25:2;
                    [N,~] = histcounts(ISI,10.^edges);
                    refract = sum(N(1:2))*100/sum(N);
                    histogram(ISI,10.^edges)
                    set(gca,'xscale','log')
                    xticks([1E-4 1E-3 1E-2 1E-1 1 10]);
                    xticklabels({0.0001, 0.001, 0.01, 0.1, 1, 10})
                    xlabel('time (ms)')
                    title([num2str(refract) '%']);
                    
                    drawnow
                    
                    prompt = 'merge (m) or cancel (c)';
                    xx = input(prompt);
                    j_list = [j_list; j];
                    
                    switch xx
                        case 'm'
                            merge_list = [merge_list;i,j];
                            
                        case 'c'
                    end
                end
                
            end
            
            prompt = 'merge all? "y", "n"';
            
            xxx = input(prompt);
            if strcmp(xxx,'y')
                for m = 1:size((merge_list),1)
                    cid = merge_list(m,:);
                    obj.waveforms.raw{length(Cids)+1} = [];
                    obj.spike_times{length(SU)+1} = [];
                    for t = 1:length(cid)
                        %                     i = find(SU ==cid(t));
                        i = cid(t);
                        obj.waveforms.raw{length(Cids)+1} = [obj.waveforms.raw{length(Cids)+1}; obj.waveforms.raw{SU(i)}];
                        obj.spike_times{length(SU)+1} = [obj.spike_times{length(SU)+1}; obj.spike_times{i}];
                    end
                    
                    obj.waveforms.mean{length(Cids)+1} = mean(obj.waveforms.raw{length(Cids)+1},1);
                    obj.waveforms.std{length(Cids)+1} = std(obj.waveforms.raw{length(Cids)+1},[],1);
                    
                    peak_to_peak = max(obj.waveforms.mean{length(Cids)+1})-min(obj.waveforms.mean{length(Cids)+1});
                    noise = mean(obj.waveforms.std{length(Cids)+1}(1:10));
                    obj.waveforms.SNR(length(Cids)+1) = 20.*log10(peak_to_peak./noise);
                    
                    SU = [SU length(Cids)+1];
                    Cids = [Cids max(Cids)+1];
                    Cchannels = [Cchannels; Cchannels(SU(i))];
                    
                    
                    
                end
                
                merge_list = reshape(merge_list,[],1);
                merge_list = sort(unique(merge_list),'descend');
                for m = 1: length(merge_list)
                    obj.spike_times(merge_list(m)) = [];
                    SU(merge_list(m)) =[];
                    
                    obj.manual = 1;
                end
                
            end
            
            
        end
    end
end




%% Code for PCA analysis

% figure
% 
% for ls = 1:length(list)
%     [Wi,score,latent] = pca(a.waveforms.raw{list(ls)});
%     X = [score(:,1),score(:,2)];
%     t = 0;
%     Crit_values = [];
%     Opt_k = [];
%     while t <10
%     eva = evalclusters(X,'gmdistribution','Silhouette','KList',[1:3]);
%     Crit_values = [Crit_values; eva.CriterionValues];
%     Opt_k = [Opt_k eva.OptimalK];
%      t = t+1;
%     end
%     mean_Crit_values = mean(Crit_values,1);
%    [h,p] = ttest2(Crit_values(:,2),Crit_values(:,3));
%     subplot(2,2,ls)
%     gm = fitgmdist(X,mode(Opt_k));
%     idx = cluster(gm,X);
%     gscatter(X(:,1),X(:,2),idx,'rb','+o',5)
%     xlabel(['Critvalue = ' num2str(mean_Crit_values(2)) ', ' num2str(mean_Crit_values(3))]);
%     ylabel(['p =' num2str(p)]);
%     drawnow 
% end
%     
% 
% test{1} = randsample(find(idx ==1),min(length(find(idx ==1)),200));
% test{2} = randsample(find(idx ==2),min(length(find(idx ==2)),200));
% test{3} = randsample(find(idx ==3),min(length(find(idx ==3)),200));
% figure
% for s = 1:2
%     subplot(2,3,s)
%     for t = 1:length(test{s})
%         plot(a.waveforms.raw{34}(test{s}(t),:))
%         hold on
%     end
%     subplot(2,3,s+2)
%     spktimes{s} = a.spike_times{18}(find(idx ==s));
%     ISI = spktimes{s}(2:end)-spktimes{s}(1:end-1);
%     
%     edges = -3.5:0.25:2;
%     [N,~] = histcounts(ISI,10.^edges);
%     refract = sum(N(1:2))*100/sum(N);
%     h = histogram(ISI,10.^edges)
%     set(gca,'xscale','log')
%     xticks([1E-4 1E-3 1E-2 1E-1 1 10]);
%     xticklabels({0.0001, 0.001, 0.01, 0.1, 1, 10})
%     xlabel('time (ms)')
%     title([num2str(refract) '%']);
% end
    