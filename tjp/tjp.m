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
    end
    
    
    methods(Static)
        
        
        
    end
    
    methods
        function obj = initialize(obj)
            global spike_times_all spike_clusters Cchannels Cids Cgroups
            spike_times_all = readNPY(fullfile(obj.params.fpath, 'spike_times.npy'));
            spike_times_all = double(spike_times_all);
            spike_clusters = readNPY(fullfile(obj.params.fpath, 'spike_clusters.npy'));
            [Cids, Cgroups] = readClusterGroupsCSV(fullfile(obj.params.fpath, 'cluster_group.tsv'));
            spike_templates = readNPY(fullfile(obj.params.fpath, 'spike_templates.npy'));
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
            global spike_times_all spike_clusters Cchannels Cids Cgroups
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
            
            channels  = obj.params.chanMap(Cchannels);
            
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
        
        function obj = Single_Units(obj)
            global spike_times_all spike_clusters Cchannels Cids Cgroups
            
            
            %[data1,timestamps1,info] = load_open_ephys_data_faster([fpath filesep file_type '_CH' num2str(ch) '.continuous' ]);
            %
            [data, timestamps, info] = load_open_ephys_data(fullfile(obj.params.fpath, 'all_channels.events'));
            
            start_stim_times = timestamps(find(info.eventId ==1))-obj.params.rec_start_time;
            end_stim_times = timestamps(find(info.eventId ==0))-obj.params.rec_start_time;
            
            
            SU = find(Cgroups == 2);
            SU = intersect(SU,find(obj.waveforms.SNR>10));
            
            MU = find(Cgroups == 1);
            
            x = eval(obj.params.xbz_file_name);
            obj.analysis_code = x.analysis_code;
            obj.analysis_type = x.analysis_type;
            if obj.analysis_code ~= 100
                if length(unique(x.stimulus_ch1(:,3)))~= 1 && length(unique(x.stimulus_ch1(:,8)))~= 1
                    obj.analysis_type = 'FRA';
                else
                    obj.plot_type  = 'Tuning';
                end
            else
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
                
                spike_timesSU{id} = spike_times_all(find(spike_clusters == Cids(SU(id))))/obj.params.fs;
                for rep = 1:TotalReps
                    %for raster
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
            
            
            stim_label  = x.stimulus_ch1(:,8);
            
            SUrate = {};
            % figure
            
            
            for id = 1:length(SU)
                figure
                suptitle(['cluster ' num2str(Cids(SU(id))) ' channel ' num2str(Cchannels(id)) ' SNR: ' num2str(obj.waveforms.SNR(SU(id)))])
                set(gcf, 'Position', get(gcf,'Position').*[1 1 0 0] + [0 -600 1000 800]);
                
                %// this should change depending on stim_set. make it modular?
                SUrate{id}.mean = mean(rate.stim{id},2); %-mean(rate_pre{id},2);
                SUrate{id}.error = std(rate.stim{id},1,2)/sqrt(nreps);
                SUrate{id}.spont = mean(mean(rate.pre{id}));
                
                
                
                %waveform
                subplot(2,2,3)
                time_step = [1:60]/obj.params.fs*1e3;
                %     Nb_wave = min(size(waveforms.raw{find(Cchannels ==Cchannels(id))},1),100);
                Nb_wave = min(size(obj.waveforms.raw{SU(id)},1),100);
                
                for n = 1:Nb_wave
                    plot(time_step,obj.waveforms.raw{SU(id)}(n,:))
                    hold on
                end
                plot(time_step,obj.waveforms.mean{SU(id)},'Linewidth',2,'Color','k')
                title('waveform')
                xlabel('time (ms)')
                ylabel('mV');
                hold off
                
                if strcmp(obj.plot_type,'Tuning') == 1
                    
                    %// Tuning curve
                    subplot(2,2,1)
                    
%                     xs = 1:length(stim_label);
%                     h =1.0;
%                     for i = 1:length(stim_label)
%                         ys1(i) = gaussian_kern_reg(xs(i),xs,SUrate{id}.mean.',h);
%                     end
                    errorbar(stim_label,SUrate{id}.mean,SUrate{id}.error,'LineWidth',2);
                    hold on
                    
                    plot(stim_label,ones(1,length(stim_label))*SUrate{id}.spont,'--k');
                    xticks(stim_label(1:2:end).');
                    xtickangle(45) 
%                     labels  = {};
%                     for st = 1:length(stim_label)
%                         labels = {labels, num2str(stim_label(st))};
%                     end
%                     xticklabels(num2str(xticks))
                    xlabel('Hz')
                    ylabel('firing rate (spikes/s)')
                    title('tuning curve')
                    
                    
                    %Raster plot
                    subplot(2,2,[2 4])
%                     area([0 StimDur StimDur 0],[0 TotalReps+5 0 TotalReps+5],'LineStyle','none','FaceColor',[.85 .85 1]);
                    for st = 1:length(StimDur)
                        rectangle('Position',[0 nreps*(st-1),StimDur(st) nreps],'FaceColor',[0.9,0.9,0.9],'EdgeColor','none')
                    end
                    hold on
                    plot(raster.spikes{id},nreps*(raster.stim{id}-1)+raster.rep{id},'k.','MarkerSize',9);
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
                        SUrate{id}.mean_fr = SUrate{id}.mean./stim_length;
                    end
                    
                    subplot(2,2,1)
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
                    
                    subplot(2,2,[2 4])
                    %                     area([0 StimDur StimDur 0],[0 TotalReps+5 0 TotalReps+5],'LineStyle','none','FaceColor',[.85 .85 1]);
                    for st = 1:length(StimDur)
                        rectangle('Position',[0 nreps*(st-1),StimDur(st) nreps],'FaceColor',[0.9,0.9,0.9],'EdgeColor','none')
                    end
                    hold on
                    plot(raster.spikes{id},nreps*(raster.stim{id}-1)+raster.rep{id},'k.','MarkerSize',9);
                    %     pause
                    xlabel('time (s)')
                    ylabel('reps')
                    yticks([1:nreps:TotalReps]+5)
%                     yticks([1:nreps*2:TotalReps]+10)
%                     for stim = 1:length(stim_label)/2
%                         stim_ticks{stim}=num2str(round(stim_label(stim*2)*10)/10);
%                     end
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
        
        
        
        
        
    end
end
