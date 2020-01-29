



figure_on = 0;
N_list_1 = unique([Pool_1.neuron_nb]);
% N_list_2 = unique([Pool_2.neuron_nb]);

SUrate_1 = {};
BF_pool = [];

for n = 1:length(N_list_1)
    rec_list_1 = find([Pool_1.neuron_nb] ==N_list_1(n));
    p = rec_list_1(1); %try end as well
    nreps = size(rate_1.stim{1},2);
    SUrate_1{n}.mean = mean(rate_1.stim{p},2); % -mean(rate_1.pre{p},2); %(Spikes/second)
    SUrate_1{n}.error = std(rate_1.stim{p},1,2)/sqrt(nreps);
    SUrate_1{n}.spont = mean(mean(rate_1.pre{p}));
    StimDur = Pool_1(p).xb.stimulus_ch1(:,5)*1e-3;
    
    % Change stim label
    if Pool_1(p).xb.analysis_code >100 && Pool_1(p).xb.analysis_code< 1000
        plot_type = 'User';
        stim_label_1 = Pool_1(p).xb.stimulus_ch1(:,1);
        for st = 1: length(stim_label_1)
            idx1 = strfind(Pool_1(p).xb.user_stimulus_desc_ch1{st},': ')+2;
            idx2 = strfind(Pool_1(p).xb.user_stimulus_desc_ch1{st},'_m')-1;
            idx3 = strfind(Pool_1(p).xb.user_stimulus_desc_ch1{st},'rev');
            idx4 = strfind(Pool_1(p).xb.user_stimulus_desc_ch1{st},'.txt');
            stim_name{st} = Pool_1(p).xb.user_stimulus_desc_ch1{st}(idx1:idx4);
%             if isempty(idx3) ~=1
%                 stim_name{st} = [stim_name{st} '_rev'];
%             end
            
        end
    elseif  strcmp(Pool_1(p).xb.analysis_type,'User: Trill_list.txt') 
        plot_type = 'User';
        stim_label_1 = Pool_1(p).xb.stimulus_ch1(:,1);
        for st = 1: length(stim_label_1)
            idx1 = strfind(Pool_1(p).xb.user_stimulus_desc_ch1{st},': ')+2;
            idx2 = strfind(Pool_1(p).xb.user_stimulus_desc_ch1{st},'_m')-1;
            idx3 = strfind(Pool_1(p).xb.user_stimulus_desc_ch1{st},'rev');
            idx4 = strfind(Pool_1(p).xb.user_stimulus_desc_ch1{st},'.txt');
            stim_name{st} = Pool_1(p).xb.user_stimulus_desc_ch1{st}(idx1:idx4);
%             if isempty(idx3) ~=1
%                 stim_name{st} = [stim_name{st} '_rev'];
%             end
            
        end
        
        %     elseif Pool_1(p).xb.analysis_code >2000
        
    elseif      Pool_1(p).xb.analysis_code == 1
        plot_type = 'Tones';
        stim_label_1 = Pool_1(p).xb.stimulus_ch1(:,8);
        max_id = find(SUrate_1{n}.mean == max(SUrate_1{n}.mean));
        
        if SUrate_1{n}.mean(max_id(1)) > SUrate_1{n}.spont +2*SUrate_1{n}.error(max_id(1))
            BF_pool = [BF_pool stim_label_1(max_id(1))];
        end
    else
        plot_type = 'VT';
        stim_label_1 = Pool_1(p).xb.stimulus_ch1(:,10);
        
    end
    %
    nreps = Pool_1(p).xb.stimulus_ch1(1,4);
    nStim = max(Pool_1(p).xb.stimulus_ch1(:,1));
    %             nreps = 10;
    
    %
    % nStim = nStim + max(x2.stimulus_ch1(:,1));
    %
    TotalReps = nStim*nreps;
    
    % for real twitters
    nat_rev(n,1) = mean(SUrate_1{n}.mean(1:2:20));
    nat_rev(n,2) = mean(SUrate_1{n}.mean(2:2:20));
    
    
    if figure_on == 1
        if SUrate_1{n}.spont >0% && SUrate_2{n}.spont >1
            
            
            
            figure
            set(gcf, 'Position', [400 300 1800 1000]);
            
            subplot(1,2,1)
            errorbar(stim_label_1,SUrate_1{n}.mean,SUrate_1{n}.error,'LineWidth',2);
            hold on
            plot(stim_label_1,ones(1,length(stim_label_1))*SUrate_1{n}.spont,'--k');
            
            if strcmp(plot_type,'User')
                xticks(stim_label_1);
                xticklabels(stim_name)
                xtickangle(45)
            elseif strcmp(plot_type,'Tones')
                set(gca,'xscale','log')
                xticks(round((stim_label_1(1:4:end).')*1e-1)*1e1);
                xticklabels(round(stim_label_1(1:4:end)*1e-1)*1e-2)
                
                xtickangle(45)
                xlabel('kHz')
                ylabel('firing rate (spikes/s)')
                title('tuning curve')
                
            end
            
            
            
            subplot(1,2,2)
            PreStim = Pool_1(p).xb.pre_stimulus_record_time*1e-3; %s
            PostStim = Pool_1(p).xb.post_stimulus_record_time*1e-3; %s
            StimDur = Pool_1(p).xb.stimulus_ch1(:,5)*1e-3;
            
            nreps = Pool_1(p).xb.stimulus_ch1(1,4);
            nStim = max(Pool_1(p).xb.stimulus_ch1(:,1));
            %             nreps = 10;
            
            %
            % nStim = nStim + max(x2.stimulus_ch1(:,1));
            %
            TotalReps = nStim*nreps;
            
            
            
            for st = 1:length(StimDur)
                rectangle('Position',[0 nreps*(st-1),StimDur(st) nreps],'FaceColor',[0.9,0.9,0.9],'EdgeColor','none')
            end
            hold on
            plot(raster_1.spikes{p},nreps*(raster_1.stim{p}-1)+raster_1.rep{p},'k.','MarkerSize',15);
            %     pause
            xlabel('time (s)')
            ylabel('Hz')
            yticks([1:nreps*2:TotalReps]+10)
            %             yticks([1:nreps*2:TotalReps]+10)
            stim_ticks = {};
            
            for stim = 1:length(stim_label_1)
                %             stim_ticks{stim}=num2str(round(stim_label(stim)*10)/10);
                stim_ticks{stim}=num2str(round(stim_label_1(stim)*10)/10);
            end
            if strcmp(plot_type,'User')
                
                yticks([1:nreps:TotalReps]+5)
                
                yticklabels(stim_name)
            else
                yticklabels(stim_ticks(1:2:end))
            end
            axis([-PreStim max(StimDur) + PostStim 0 TotalReps+1])
            hold off
            title('rasterplot')
            
            drawnow
        end
    end
    
end

figure

scatter(nat_rev(:,1),nat_rev(:,2))
hold on
plot([1:20],[1:20])
% edges = stim_label_1;
% % histogram(stim_label(BF_pool),edges);
% histogram(BF_pool,stim_label_1)
%                     set(gca,'xscale','log')



