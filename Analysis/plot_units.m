save_dir = 'D:\Data\M12E\Neurons';
load(fullfile(save_dir,'M12E_neuron_list.mat'));
u = 11;


ch = 1;
obj.params.xbz_file_name = M12E_neuron_list.list{u};
obj.params.fpath = ['C:\Data\OpenEphys\M12E\' '2019-07-25_14-25-19'];
obj.params.file_type = '100';

[data1,timestamps1{ch},info1{ch}] = load_open_ephys_data_faster([obj.params.fpath filesep obj.params.file_type '_CH' num2str(ch) '.continuous' ]);
obj.params.rec_start_time = timestamps1{1}(1);
[data, timestamps, info] = load_open_ephys_data(fullfile(obj.params.fpath, 'all_channels.events'));
start_stim_times = timestamps(find(info.eventId ==1))-obj.params.rec_start_time;
end_stim_times = timestamps(find(info.eventId ==0))-obj.params.rec_start_time;
addpath('C:\Data\OpenEphys\M12E');



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




A = unique(M12E_neuron_list.data(:,u));
for i = 2:length(A)
    unit_file_name = 'M12Eu000';
    unit_file_name = [unit_file_name(1:end-size(num2str(A(i)),2)) num2str(A(i)) '.mat'];
    y = load(unit_file_name);
    id= 1;
    raster.stim{id} = [];
    raster.rep{id} = [];
    raster.spikes{id} = [];
    spikes_pooled{id} = [];
    spike_timesSU{id} = y.s_unit.spiketimes;
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
    
    figure
    set(gcf, 'Position', [400 300 1800 1000]);
    suptitle(['cluster_id = ' num2str(A(i))])
    subplot(2,2,3)
    time_step = [1:60]/30000*1e3;
    Nb_wave = min(size(y.s_unit.waveforms{1},1),100);
    idx = 1:size(y.s_unit.waveforms{1},1);
    rand_select = randsample(idx,Nb_wave);
    for n = 1:length(rand_select)
        
        plot(time_step,y.s_unit.waveforms{1}(rand_select(n),:))
        hold on
    end
    plot(time_step,mean(y.s_unit.waveforms{1},1),'Linewidth',2,'Color','k')
    title('waveform')
    xlabel('time (ms)')
    ylabel('uV');
    hold off
    
    SUrate{id}.mean = mean(rate.stim{id},2); %-mean(rate_pre{id},2);
    SUrate{id}.error = std(rate.stim{id},1,2)/sqrt(nreps);
    SUrate{id}.spont = mean(mean(rate.pre{id}));
    stim_label  = x.stimulus_ch1(:,8);
    if strcmp(obj.plot_type,'Tuning') == 1
        
        %// Tuning curve
        subplot(2,2,1)
        
        
        errorbar(stim_label,SUrate{id}.mean,SUrate{id}.error,'LineWidth',2);
        hold on
        
        plot(stim_label,ones(1,length(stim_label))*SUrate{id}.spont,'--k');
        set(gca,'xscale','log')
        
        xticks(round((stim_label(1:4:end).')*1e2)*1e-2);
        xtickangle(45)
        xlabel('kHz')
        ylabel('firing rate (spikes/s)')
        title('tuning curve')
        
        
        %Raster plot
        subplot(2,2,[2 4])
        for st = 1:length(StimDur)
            rectangle('Position',[0 nreps*(st-1),StimDur(st) nreps],'FaceColor',[0.9,0.9,0.9],'EdgeColor','none')
        end
        hold on
        plot(raster.spikes{id},nreps*(raster.stim{id}-1)+raster.rep{id},'k.','MarkerSize',15);
        %     pause
        xlabel('time (s)')
        ylabel('kHz')
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
        subplot(2,2,1)
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
        subplot(2,2,[2 4])
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
    %
    
    
    
end


    
    
    
    
    
    
    
    
    
    
    
    
    
    