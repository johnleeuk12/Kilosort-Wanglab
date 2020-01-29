params = eval('parameters_distance');

addpath('D:\DATA\M12E\Experiments')
addpath('D:\DATA\M12E\Units')

load('M12E_unit_list.mat');
for i  = 1:length(M12E_unit_list.data)
    x = eval(M12E_unit_list.data{i,4});
    M12E_unit_list.data{i,6} = x.analysis_code;
    M12E_unit_list.data{i,7} = x.analysis_type;
end




%% The following code can be used universally, to extract spike and stim times from unit data.
dB = 40;
a_code = 2320;
u_list = find([M12E_unit_list.data{:,6}] == a_code);
p = 1;
Pool = {};
tic
for i = 1:length(u_list)
    y = eval(M12E_unit_list.data{u_list(i),4});
%     if y.stimulus_ch1(1,3) == dB
        unit_file_name = 'M12Eu000';
        unit_file_name = [unit_file_name(1:end-size(num2str(u_list(i)),2)) num2str(u_list(i)) '.mat'];
        x = load(unit_file_name);
        if length(x.s_unit.spiketimes)>10
            if ~strcmp(y.analysis_type,'User: Trill_list.txt')  %temp solution
                
                
                if ~isempty(x.s_unit.templates)
                    Pool(p).best_ch = x.s_unit.templates.best_ch;
                else
                    data = zeros(64,60);
                    for ch = 1:64
                        temp = [];
                        temp(:,:) = x.s_unit.waveforms{1}(ch,:,:);
                        [Y,Sig,X] = svd(temp,'econ');
                        %                 sig = diag(Sig);%figure; semilogy(sig(sig>1),'kx-')
                        k = 1:3;
                        P = Y(:,k)*Sig(k,k)*X(:,k)';
                        data(ch,:) = mean(P,2).';
                    end
                    data = abs(data);
                    [~,Pool(p).best_ch] = max(mean(data(:,10:40),2));
                end
                
                Pool(p).waveforms(:,:) = x.s_unit.waveforms{1}(Pool(p).best_ch,:,:);
                Pool(p).spiketimes = x.s_unit.spiketimes;
                Pool(p).xb = y;
                
                prev_unit_fn = [unit_file_name(1:end-4) '_prev.mat'];
                if isfield(x.s_unit,'start_times')
                    Pool(p).stim_times = [x.s_unit.start_times x.s_unit.end_times];
                else
                    z = load(prev_unit_fn);
                    Pool(p).stim_times = [z.s_unit.start_times z.s_unit.end_times];
                end
                Pool(p).unit_nb = u_list(i);
                
                for ll = 1:length(M12E_neurons_list)
                    if sum(M12E_neurons_list(ll,:) == u_list(i)) ==1
                        Pool(p).neuron_nb = ll;
                    end
                end
                
                p = p+1;
                
            end
        end
%     end
    toc
end
toc

%% Raster and Tuning curves

% p = 10;

raster.stim = {};
raster.rep = {};
raster.spikes = {};
spikes_pooled = {};
spike_timesSU = {};
rate.stim = {};
rate.pre = {};

unique_neuron_list = unique([Pool.neuron_nb]);
neuron_check_list = [];
figure
for p = 1:length(Pool)
    PreStim = Pool(p).xb.pre_stimulus_record_time*1e-3; %s
    PostStim = Pool(p).xb.post_stimulus_record_time*1e-3; %s
    
    stim_info = Pool(p).xb.data(find(Pool(p).xb.data(:,3) == 1 & Pool(p).xb.data(:,4) == -1),:);
    
    
    data_new = stim_info;
    StimDur = Pool(p).xb.stimulus_ch1(:,5)*1e-3;
    
    nreps = Pool(p).xb.stimulus_ch1(1,4);
    nStim = max(Pool(p).xb.stimulus_ch1(:,1));
    %             nreps = 10;
    
    %
    % nStim = nStim + max(x2.stimulus_ch1(:,1));
    %
    TotalReps = nStim*nreps;
    
    false_start = length(Pool(p).stim_times)-TotalReps;
    start_stim_times = Pool(p).stim_times(false_start+1:end,1);
    end_stim_times = Pool(p).stim_times(false_start+1:end,2);
    
    
    
    %             for id = 1:length(SU)
    
    
    
    raster.stim{p} = [];
    raster.rep{p} = [];
    raster.spikes{p} = [];
    spikes_pooled{p} = [];
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
        spikes1 = Pool(p).spiketimes(find(Pool(p).spiketimes>=start_stim_times(rep)-PreStim & ...
            Pool(p).spiketimes<=end_stim_times(rep)+ PostStim)).';
        spikes1 = spikes1 - start_stim_times(rep);
        spikes_pooled{p} = [spikes_pooled{p} spikes1];
        raster.stim{p} = [raster.stim{p} data_new(rep,1)*ones(size(spikes1))];
        raster.rep{p} = [raster.rep{p} data_new(rep,2)*ones(size(spikes1))];
        raster.spikes{p} = [raster.spikes{p} spikes1];
        
        %for rate
        spikes2 = Pool(p).spiketimes(find(Pool(p).spiketimes>=start_stim_times(rep) & ...
            Pool(p).spiketimes<=end_stim_times(rep))).';
        rate.stim{p}(data_new(rep,1),data_new(rep,2)) = length(spikes2)/StimDur(data_new(rep,1));
        spikes3 = Pool(p).spiketimes(find(Pool(p).spiketimes<=start_stim_times(rep) & ...
            Pool(p).spiketimes>=start_stim_times(rep)-PreStim)).' ;
        rate.pre{p}(data_new(rep,1),data_new(rep,2)) = length(spikes3)/PreStim;
    end
    
    if a_code == 1
        stim_label = Pool(p).xb.stimulus_ch1(:,8);                   
    elseif a_code >2300
        stim_number = Pool(p).xb.stimulus_ch1(:,1);
        stim_length = Pool(p).xb.stimulus_ch1(:,5);
        for i = 9:22
            VT_para(i-8) = length(unique(Pool(p).xb.stimulus_ch1(:,i)));
        end
        stim_label = Pool(p).xb.stimulus_ch1(:,find(VT_para ~= 1)+8);
        if a_code == 2320 || a_code == 2335
            stim_label = Pool(p).xb.stimulus_ch1(:,10);
        end
        
        
        
%         if length(unique(Pool(p).xb.stimulus_ch1(:,9)))>length(unique(Pool(p).xb.stimulus_ch1(:,10)))
%             stim_label = Pool(p).xb.stimulus_ch1(:,9);
%         else
%             stim_label = Pool(p).xb.stimulus_ch1(:,10);
%         end
    end

    
    SUrate{p}.mean = mean(rate.stim{p},2) -mean(rate.pre{p},2); %(Spikes/second)
    SUrate{p}.error = std(rate.stim{p},1,2)/sqrt(nreps);
    SUrate{p}.spont = mean(mean(rate.pre{p}));
    
    
%     if SUrate{p}.spont> 1 && isempty(intersect(Pool(p).neuron_nb,neuron_check_list))
%         if ~isempty(Pool(p).neuron_nb)
%             errorbar(stim_label,SUrate{p}.mean,SUrate{p}.error,'LineWidth',2);
%             hold on
%         end
%     end
    % plot(stim_label,ones(1,length(stim_label))*SUrate{p}.spont,'--k');
    neuron_check_list = [neuron_check_list, Pool(p).neuron_nb];
    drawnow
    
end


set(gca,'xscale','log')

xticks(round((stim_label(1:4:end).')*1e2)*1e-2);
xtickangle(45)
xlabel('Hz')
ylabel('firing rate (spikes/s)')
title('tuning curve')
best_fr = zeros(1,40);
tuning = [];
pp = 1;
for p = 1:length(Pool)
    if SUrate{p}.spont >1
        for st = 2:length(SUrate{p}.mean)-1
            best_fr(st) = SUrate{p}.mean(st-1) + SUrate{p}.mean(st) + SUrate{p}.mean(st+1);
        end
        temp = find(abs(best_fr) == max(abs(best_fr)));
        tuning(pp) = temp(1);
        pp = pp+1;
    end
end

figure
histogram(tuning)
xticklabels(round((stim_label(1:6:end).')*1e2)*1e-2)
xtickangle(45)
xlabel('Hz')


figure
pp = 1;
for p = 1:length(Pool)
    if SUrate{p}.spont >1
    subplot(3,7,pp)
    for st = 1:length(StimDur)
        rectangle('Position',[0 nreps*(st-1),StimDur(st) nreps],'FaceColor',[0.9,0.9,0.9],'EdgeColor','none')
    end
    hold on
    plot(raster.spikes{p},nreps*(raster.stim{p}-1)+raster.rep{p},'k.','MarkerSize',15);
    %     pause
    xlabel('time (s)')
    ylabel('IPI')
    yticks([1:nreps*2:TotalReps]+10)
    yticks([1:nreps*2:TotalReps]+10)
    stim_ticks = {};
    if length(stim_label)>2
        for stim = 1:length(stim_label)/2
%             stim_ticks{stim}=num2str(round(stim_label(stim*2)*10)/10);
            stim_ticks{stim}=num2str(round(stim_label(stim*2)*1000));
        end
    else
        for stim = 1:length(stim_label)
            %             stim_ticks{stim}=num2str(round(stim_label(stim)*10)/10);
            stim_ticks{stim}=num2str(round(stim_label(stim)*100));
        end
    end
    yticklabels(stim_ticks)
    axis([-PreStim max(StimDur) + PostStim 0 TotalReps+1])
    hold off
    title('rasterplot')
    
    drawnow()
    pp = pp+1;
    end
end

%% User 


raster.stim = {};
raster.rep = {};
raster.spikes = {};
spikes_pooled = {};
spike_timesSU = {};
rate.stim = {};
rate.pre = {};

unique_neuron_list = unique([Pool.neuron_nb]);
neuron_check_list = [];
figure
for p = 1:length(Pool)
    PreStim = Pool(p).xb.pre_stimulus_record_time*1e-3; %s
    PostStim = Pool(p).xb.post_stimulus_record_time*1e-3; %s
    StimDur = Pool(p).xb.stimulus_ch1(:,5)*1e-3;
    
    stim_info = Pool(p).xb.data(find(Pool(p).xb.data(:,3) == 1 & Pool(p).xb.data(:,4) == -1),:);
    
    
    data_new = stim_info;
    
    nreps = Pool(p).xb.stimulus_ch1(1,4);
    nStim = max(Pool(p).xb.stimulus_ch1(:,1));
    %             nreps = 10;
    
    %
    % nStim = nStim + max(x2.stimulus_ch1(:,1));
    %
    TotalReps = nStim*nreps;
    
    false_start = length(Pool(p).stim_times)-TotalReps;
    start_stim_times = Pool(p).stim_times(false_start+1:end,1);
    end_stim_times = Pool(p).stim_times(false_start+1:end,2);
    
    
    
    %             for id = 1:length(SU)
    
    
    
    raster.stim{p} = [];
    raster.rep{p} = [];
    raster.spikes{p} = [];
    spikes_pooled{p} = [];
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
        spikes1 = Pool(p).spiketimes(find(Pool(p).spiketimes>=start_stim_times(rep)-PreStim & ...
            Pool(p).spiketimes<=end_stim_times(rep)+ PostStim)).';
        spikes1 = spikes1 - start_stim_times(rep);
        spikes_pooled{p} = [spikes_pooled{p} spikes1];
        raster.stim{p} = [raster.stim{p} data_new(rep,1)*ones(size(spikes1))];
        raster.rep{p} = [raster.rep{p} data_new(rep,2)*ones(size(spikes1))];
        raster.spikes{p} = [raster.spikes{p} spikes1];
        
        %for rate
        spikes2 = Pool(p).spiketimes(find(Pool(p).spiketimes>=start_stim_times(rep) & ...
            Pool(p).spiketimes<=end_stim_times(rep))).';
        rate.stim{p}(data_new(rep,1),data_new(rep,2)) = length(spikes2)/StimDur(data_new(rep,1));
        spikes3 = Pool(p).spiketimes(find(Pool(p).spiketimes<=start_stim_times(rep) & ...
            Pool(p).spiketimes>=start_stim_times(rep)-PreStim)).' ;
        rate.pre{p}(data_new(rep,1),data_new(rep,2)) = length(spikes3)/PreStim;
    end
    

    stim_label = {};
    stim_number = Pool(p).xb.stimulus_ch1(:,1);
    stim_length = Pool(p).xb.stimulus_ch1(:,5);
    for st = 1: length(stim_number)
        idx1 = strfind(Pool(p).xb.user_stimulus_desc_ch1{st},': ')+2;
        idx2 = strfind(Pool(p).xb.user_stimulus_desc_ch1{st},'_m')-1;
        idx3 = strfind(Pool(p).xb.user_stimulus_desc_ch1{st},'rev');
        stim_label{st} = Pool(p).xb.user_stimulus_desc_ch1{st}(idx1:idx2);
        if isempty(idx3) ~=1
            stim_label{st} = [stim_label{st} '_rev'];
        end
        
    end
    
    
    SUrate{p}.mean = mean(rate.stim{p},2) -mean(rate.pre{p},2); %(Spikes/second)
    SUrate{p}.error = std(rate.stim{p},1,2)/sqrt(nreps);
    SUrate{p}.spont = mean(mean(rate.pre{p}));
    
    
    SUrate{p}.mean_fr = SUrate{p}.mean./stim_length;
    if SUrate{p}.spont> 1 && isempty(intersect(Pool(p).neuron_nb,neuron_check_list))
        if ~isempty(Pool(p).neuron_nb)
            errorbar(stim_number,SUrate{p}.mean,SUrate{p}.error,'LineWidth',2);
            hold on
            
        end
    end
    xticks(stim_number);
    xticklabels(stim_label)
    xtickangle(45)
    %                 xlabel('Hz')
    ylabel('firing rate (spikes/s)')
    title('tuning curve')
    
        neuron_check_list = [neuron_check_list, Pool(p).neuron_nb];
    drawnow
    
   
end


set(gca,'xscale','log')

xticks(round((stim_label(1:4:end).')*1e2)*1e-2);
xtickangle(45)
xlabel('Hz')
ylabel('firing rate (spikes/s)')
title('tuning curve')
best_fr = zeros(1,40);
tuning = [];
pp = 1;
for p = 1:length(Pool)
    if SUrate{p}.spont >1
        for st = 2:length(SUrate{p}.mean)-1
            best_fr(st) = SUrate{p}.mean(st-1) + SUrate{p}.mean(st) + SUrate{p}.mean(st+1);
        end
        temp = find(abs(best_fr) == max(abs(best_fr)));
        tuning(pp) = temp(1);
        pp = pp+1;
    end
end

figure
histogram(tuning)
xticklabels(round((stim_label(1:6:end).')*1e2)*1e-2)
xtickangle(45)
xlabel('Hz')


figure
pp = 1;
neuron_check_list = [];

for p = 1:length(Pool)
    if SUrate{p}.spont> 1 && isempty(intersect(Pool(p).neuron_nb,neuron_check_list))
        if ~isempty(Pool(p).neuron_nb)
            subplot(5,7,pp)
            
            for st = 1:length(StimDur)
                rectangle('Position',[0 nreps*(st-1),StimDur(st) nreps],'FaceColor',[0.9,0.9,0.9],'EdgeColor','none')
            end
            hold on
            plot(raster.spikes{p},nreps*(raster.stim{p}-1)+raster.rep{p},'k.','MarkerSize',15);
            %     pause
            xlabel('time (s)')
            ylabel('reps')
            yticks([1:nreps:TotalReps]+5)
            
            yticklabels(stim_label)
            axis([-PreStim max(StimDur) + PostStim 0 TotalReps+5])
            hold off
            title('rasterplot')
            drawnow()
            pp = pp+1;
        end
    end
end













%% Vocalization responses 

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

