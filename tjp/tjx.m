classdef tjx < handle
    
    properties(Constant)
    end
    
    properties
        params = eval('parameters_x');
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
    end
    
    
    methods(Static)
        
    end
    methods
        %%
        function obj = initialize(obj)
            global spike_times_all spike_clusters Cchannels Cids Cgroups new_channels
            
            % load sorted spikes for entire session
            spike_times_all = readNPY(fullfile(obj.params.fpath2, 'spike_times.npy'));
            spike_times_all = double(spike_times_all);
            spike_clusters = readNPY(fullfile(obj.params.fpath2, 'spike_clusters.npy'));
            [Cids, Cgroups] = readClusterGroupsCSV(fullfile(obj.params.fpath2, 'cluster_group.tsv'));
            %             spike_templates = readNPY(fullfile(obj.params.fpath2, 'spike_templates.npy'));
            new_channels = readNPY(fullfile(obj.params.fpath2, 'channel_map.npy'));
            new_channels = new_channels+1;
            Cchannels = [];
            temp = zeros(1,length(obj.params.list));
            parfor ls = 1:length(obj.params.list)
                
                [data1,~,~] = load_open_ephys_data_faster([obj.params.fpath obj.params.list{ls} filesep obj.params.file_type '_CH' num2str(1) '.continuous' ]);
                temp(ls) = length(data1);
            end
            obj.seg_length = cumsum(temp);
            
        end
        %%
        function obj = load_data_OE(obj)
            
            data_all = [];
            timestamps1 = {};
            info1 = {};
            file_type = '100';
            
            tic
            fprintf('Time %3.0fs. Loading data... \n', toc);
            parfor ch = 1:obj.params.Nb_ch
                %     disp(['loading data from channel ' num2str(ch)])
                [data1,timestamps1{ch},info1{ch}] = load_open_ephys_data_faster([obj.params.fpath obj.params.session_name filesep obj.params.file_type '_CH' num2str(ch) '.continuous' ]);
                data_all(ch,:) = data1;
            end
            obj.params.rec_start_time = timestamps1{1}(1);
            obj.raw_data = data_all;
            fprintf('Time %3.0fs. Loading data... Done \n', toc);
        end
        %%
        function obj = extract_units(obj)
            tic
            global spike_times_all spike_clusters Cchannels Cids Cgroups new_channels SU 
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
            filtData = filtData(obj.params.chanMap,:);
            
            fprintf('Time %3.0fs. Extracting waveforms and templates... \n', toc);
            %// extracting waveforms from processed data
            obj.waveforms.raw = {};
            obj.waveforms.mean = {};
            obj.waveforms.std = {};
            obj.waveforms.SNR = [];
            obj.templates = {};
            %             spike_times_all = {};
            spike_times_SU = {};
            SU = find(Cgroups == 2);
            obj.SU_good = ones(length(SU),1);
            %             obj.seg_length(1) =  12737536;
            %             obj.seg_length(2) = obj.seg_length(1) + 3003392;
            for id = 1:length(SU)
                disp(id)
                spike_times_SU{id} = spike_times_all(find(spike_clusters == Cids(SU(id)))); %/obj.params.fs;
                
                if obj.params.session_id > 1
                    spike_times_SU{id} = spike_times_SU{id}(find(spike_times_SU{id} < (obj.seg_length(obj.params.session_id)-40) & (obj.seg_length(obj.params.session_id-1)+20)< spike_times_SU{id}));
                    
                    spike_times_SU{id} = spike_times_SU{id} - obj.seg_length(obj.params.session_id-1);
                else
                    spike_times_SU{id} = spike_times_SU{id}(find(spike_times_SU{id} < obj.seg_length(obj.params.session_id)));
                end
                
                obj.waveforms.raw{id} = zeros(64,60,length(spike_times_SU{id}));
                
                for t = 1:length(spike_times_SU{id})
                    obj.waveforms.raw{id}(:,:,t) = filtData(:, ...
                        round(spike_times_SU{id}(t,1)-19):round(spike_times_SU{id}(t,1)+40)); %
                end
                if size(obj.waveforms.raw{id},3) < 10 % if number of spikes is less than 10, remove cluster
                    obj.SU_good(id) = 0;
                else
                    % extracting templates via SVD
                    obj.templates{id}.data = zeros(64,60);
                    %                  obj.templates{id}.best_ch = [];
                    for ch = 1:64
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
            
            %             for id = 1:length(SU)
            %                 plot(obj.templates{id}(2,:))
            %                 hold on
            %             end
            fprintf('Time %3.0fs. Extracting waveforms... Done \n', toc)
            
            
            %             for id = 1:length(SU)
            %                 figure
            %                 for ch = 1:64
            %                 subplot(16,4,ch)
            %                 plot(obj.templates{id}.data(ch,:))
            %                 axis([0 60 -10 5]);
            %                 end
            %                 drawnow
            % % %             end
            %             figure
            %             for n = 1:100
            %             plot(obj.waveforms.raw{id}(64,:,n))
            %             hold on
            %             end
            
            
            
        end
        %%
        
        function obj = save_Units(obj)
            global   Cchannels Cids  SU 
            
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
                M12E_unit_list.tags{5} = 'Cid';
                M12E_unit_list.data = {};
                save(fullfile(save_dir,'M12E_unit_list'),'M12E_unit_list');
            end
            
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
                    s_unit.xbz_file_name = obj.params.xbz_file_name;
                    s_unit.amplitude = max(abs(obj.waveforms.mean{id}));
                    s_unit.start_times = obj.start_times;
                    s_unit.end_times =obj.end_times;
                    unitname= 'M12Eu001';
                    if isempty(M12E_unit_list.data)
                        save(fullfile(save_dir,unitname),'s_unit')
                        M12E_unit_list.data{1,1} = 1;
                        M12E_unit_list.data{1,2} = s_unit.SNR;
                        M12E_unit_list.data{1,3} = s_unit.amplitude;
                        M12E_unit_list.data{1,4} = s_unit.xbz_file_name;
                        M12E_unit_list.data{1,5} = s_unit.cid;
                    else
                        unitname = [unitname(1:end-length(num2str(size((M12E_unit_list.data),1)+1))) num2str(size((M12E_unit_list.data),1)+1)];
                        save(fullfile(save_dir,unitname),'s_unit')
                        M12E_unit_list.data{size((M12E_unit_list.data),1)+1,1} = size((M12E_unit_list.data),1)+1;
                        M12E_unit_list.data{size((M12E_unit_list.data),1),2} = s_unit.SNR;
                        M12E_unit_list.data{size((M12E_unit_list.data),1),3} = s_unit.amplitude;
                        M12E_unit_list.data{size((M12E_unit_list.data),1),4} = s_unit.xbz_file_name;
                        M12E_unit_list.data{size((M12E_unit_list.data),1),5} = s_unit.cid;
                        
                    end
                    save(fullfile(save_dir,'M12E_unit_list'),'M12E_unit_list')
                end
            end
            
        end
        
        %%
        
        function obj = Single_units(obj)
            
            global SU  Cids  
            
            
            [~, timestamps, info] = load_open_ephys_data(fullfile(obj.params.fpath, obj.params.session_name, 'all_channels.events'));
            
            start_stim_times = timestamps(find(info.eventId ==1))-obj.params.rec_start_time;
            end_stim_times = timestamps(find(info.eventId ==0))-obj.params.rec_start_time;
            obj.start_times = start_stim_times;
            obj.end_times = end_stim_times; 
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
            
            %             if obj.analysis_code == 2367
            if strcmp(x.analysis_type,'User: Trill_list.txt')
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
            %             for id = 1:length(SU)
            %             obj.spike_times{id} = obj.spike_times{id}.*30000;
            %             end
            for id = 1:length(SU)
                if obj.SU_good(id) ==1
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
                        rate.pre{id}(data_new(rep,1),data_new(rep,2)) = length(spikes3)/PreStim;
                    end
                end
            end
            
            
            stim_label  = x.stimulus_ch1(:,8);
            
            for id = 1:length(SU)
                if obj.SU_good(id) == 1
                    figure
                    suptitle(['cluster ' num2str(Cids(SU(id))) ' channel ' num2str(obj.templates{id}.best_ch) ' SNR: ' num2str(obj.waveforms.SNR(id))])
                    set(gcf, 'Position', [400 300 1800 1000]);
                    
                    %// this should change depending on stim_set. make it modular?
                    SUrate{id}.mean = mean(rate.stim{id},2); %-mean(rate_pre{id},2);
                    SUrate{id}.error = std(rate.stim{id},1,2)/sqrt(nreps);
                    SUrate{id}.spont = mean(mean(rate.pre{id}));
                    
                    
                    
                    %waveform
                    subplot(2,3,4)
                    time_step = [1:60]/obj.params.fs*1e3;
                    Nb_wave = min(size(obj.spike_times{id},1),100);
                    idx = 1:size(obj.spike_times{id},1);
                    rand_select = randsample(idx,Nb_wave);
                    for n = 1:length(rand_select)
                        
                        plot(time_step,obj.waveforms.raw{id}(obj.templates{id}.best_ch,:,rand_select(n)))
                        hold on
                    end
                    plot(time_step,obj.waveforms.mean{id},'Linewidth',2,'Color','k')
                    title('waveform')
                    xlabel('time (ms)')
                    ylabel('uV');
                    hold off
                    
                    %// PCA
                    Y = [];
                    Y(:,:) = obj.waveforms.raw{id}(obj.templates{id}.best_ch,:,:);
                    Y = Y.';
                    [~,score,~] = pca(Y);
                    
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
            
            %          figure
            %          for id = 1:length(SU)
            %              plot(obj.templates{id}.data(obj.templates{id}.best_ch,:))
            %              hold on
            %              pause(1)
            %
            %          end
            %
            
            
            
            
            
        end
        
        
        
        
    end
end




















