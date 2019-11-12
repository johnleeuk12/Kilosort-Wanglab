%% Code for comparing two stim set responses

N_list_1 = unique([Pool_1.neuron_nb]);
N_list_2 = unique([Pool_2.neuron_nb]);

N_list_all = intersect(N_list_1,N_list_2);
N_list_all = N_list_all(2:end); %remove -1
% stim_label_1 = Pool_1(1).xb.stimulus_ch1(:,8)/1000;
% stim_label_2 = Pool_2(1).xb.stimulus_ch1(:,9);



for n = 1:length(N_list_all)
%     stim_label_1 = Pool_1(p).xb.stimulus_ch1(:,8)/1000;

    rec_list_1 = find([Pool_1.neuron_nb] ==N_list_all(n));
    rec_list_2 = find([Pool_2.neuron_nb] ==N_list_all(n));
    
    if ~isempty(rec_list_1) && ~isempty(rec_list_2)
        
        p = rec_list_1(1); %try end as well
        nreps = size(rate_1.stim{1},2);
        SUrate_1{n}.mean = mean(rate_1.stim{p},2) -mean(rate_1.pre{p},2); %(Spikes/second)
        SUrate_1{n}.error = std(rate_1.stim{p},1,2)/sqrt(nreps);
        SUrate_1{n}.spont = mean(mean(rate_1.pre{p}));
        StimDur = Pool_1(p).xb.stimulus_ch1(:,5)*1e-3;
        
        stim_label_1 = Pool_1(p).xb.stimulus_ch1(:,8)/1000;
        
        nreps = Pool_1(p).xb.stimulus_ch1(1,4);
        nStim = max(Pool_1(p).xb.stimulus_ch1(:,1));
        %             nreps = 10;
        
        %
        % nStim = nStim + max(x2.stimulus_ch1(:,1));
        %
        TotalReps = nStim*nreps;
        
        
        
        
        
        p2 = rec_list_2(1); %try end as well
        nreps2 = size(rate_2.stim{1},2);
        SUrate_2{n}.mean = mean(rate_2.stim{p2},2) -mean(rate_2.pre{p2},2); %(Spikes/second)
        SUrate_2{n}.error = std(rate_2.stim{p2},1,2)/sqrt(nreps2);
        SUrate_2{n}.spont = mean(mean(rate_2.pre{p2}));
        
        stim_label_2 = Pool_2(p2).xb.stimulus_ch1(:,9);
        
        
        if SUrate_1{n}.spont >1% && SUrate_2{n}.spont >1
            figure
            set(gcf, 'Position', [400 300 1800 1000]);
            
            subplot(2,2,[1 2])
            errorbar(stim_label_1,SUrate_1{n}.mean,SUrate_1{n}.error,'LineWidth',2);
            hold on
            errorbar(stim_label_2,SUrate_2{n}.mean,SUrate_2{n}.error,'LineWidth',2);
            drawnow
            
            
                 
           % %  first rasterplot
            subplot(2,2,3)
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
            ylabel('kHz')
            yticks([1:nreps*2:TotalReps]+10)
%             yticks([1:nreps*2:TotalReps]+10)
            stim_ticks = {};
            
            for stim = 1:length(stim_label_1)
                %             stim_ticks{stim}=num2str(round(stim_label(stim)*10)/10);
                stim_ticks{stim}=num2str(round(stim_label_1(stim)*10)/10);
            end

            yticklabels(stim_ticks(1:2:end))
            axis([-PreStim max(StimDur) + PostStim 0 TotalReps+1])
            hold off
            title('rasterplot')
            % % second rasterplot
            subplot(2,2,4)
            
            
            PreStim = Pool_2(p2).xb.pre_stimulus_record_time*1e-3; %s
            PostStim = Pool_2(p2).xb.post_stimulus_record_time*1e-3; %s
            StimDur = Pool_2(p2).xb.stimulus_ch1(:,5)*1e-3;
            
            nreps = Pool_2(p2).xb.stimulus_ch1(1,4);
            nStim = max(Pool_2(p2).xb.stimulus_ch1(:,1));
            %             nreps = 10;
            
            %
            % nStim = nStim + max(x2.stimulus_ch1(:,1));
            %
            TotalReps = nStim*nreps;
            
            
            
            for st = 1:length(StimDur)
                rectangle('Position',[0 nreps*(st-1),StimDur(st) nreps],'FaceColor',[0.9,0.9,0.9],'EdgeColor','none')
            end
            hold on
            plot(raster_2.spikes{p2},nreps*(raster_2.stim{p2}-1)+raster_2.rep{p2},'k.','MarkerSize',15);
            %     pause
            xlabel('time (s)')
            ylabel('kHz')
%             yticks([1:nreps:TotalReps]+10)
            yticks([1:nreps*2:TotalReps]+10)
            stim_ticks = {};
            
            for stim = 1:length(stim_label_2)
                %             stim_ticks{stim}=num2str(round(stim_label(stim)*10)/10);
                stim_ticks{stim}=num2str(round(stim_label_2(stim)));
            end
            
            yticklabels(stim_ticks(1:2:end))
            axis([-PreStim max(StimDur) + PostStim 0 TotalReps+1])
            hold off
            title('rasterplot')
            
        end
    end
end


diff = [];

for n = 1:length(N_list_all)
    rec_list_1 = find([Pool_1.neuron_nb] ==N_list_all(n));
    p  = rec_list_1(1);
    stim_label_1 = Pool_1(p).xb.stimulus_ch1(:,8)/1000;
        
    if SUrate_1{n}.spont >1 % && SUrate_2{n}.spont >1
        [y1,Ty1] = resample(SUrate_1{n}.mean,stim_label_1,2.5);
        t1 = find(Ty1 == 4);
        t2 = find(Ty1 == 21.60);
        resamp = y1(t1:2:t2);
%         resamp = resamp/max(abs(resamp));
%         diff(n,:) = resamp- (SUrate_2{n}.mean)/max(abs(SUrate_2{n}.mean));
        diff(n,:) = resamp-SUrate_2{n}.mean;
%         plot(diff(n,:))
%         hold on
%         drawnow 
    end
end


p = 1;
diff2 =[];
for n = 1:length(diff)
    if sum(diff(n,:)) ~= 0
        diff2(p,:) = diff(n,:);
        p = p+1;
    end
end



figure
for n = 1:length(diff2)
%     if sum(diff(n,:)) ~= 0
        diff2(n,:) = diff2(n,:)/max(abs(diff2(n,:)));
        plot(diff2(n,:))
        hold on
%     end
end

av_diff = mean(diff2);
plot(av_diff,'LineWidth',2)

          
       