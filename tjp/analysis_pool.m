        params = eval('parameters_distance');

        addpath('D:\Data\M12E\Units')
        
        for i  = 1:length(M12E_unit_list.data)
            x = eval(M12E_unit_list.data{i,4});
            M12E_unit_list.data{i,6} = x.analysis_code;
            M12E_unit_list.data{i,7} = x.analysis_type;
        end
        
        
        
        
        
        dB = 40;
        
        u_list = find([M12E_unit_list.data{:,6}] == 1);
        p = 1;
        Pool = {};
        for i = 1:length(u_list)
            y = eval(M12E_unit_list.data{u_list(i),4});
            if y.stimulus_ch1(1,3) == dB
                unit_file_name = 'M12Eu000';
                unit_file_name = [unit_file_name(1:end-size(num2str(u_list(i)),2)) num2str(u_list(i)) '.mat'];
                x = load(unit_file_name);
                if length(x.s_unit.spiketimes)>10
                    if ~isempty(x.s_unit.templates)
                        Pool{p}.best_ch = x.s_unit.templates.best_ch;
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
                        [~,Pool{p}.best_ch] = max(mean(data(:,10:40),2));
                    end
                    
                    Pool{p}.waveforms(:,:) = x.s_unit.waveforms{1}(Pool{p}.best_ch,:,:);
                    Pool{p}.spiketimes = x.s_unit.spiketimes;
                    Pool{p}.xb = y;
                    Pool{p}.event_times = [x.s_unit.start_times x.s_unit.end_times];
                    p = p+1;
                end
            end
        end
        
        
        for p = 2:15
                        stim_label  =Pool{p}.xb.stimulus_ch1(:,8);

            
            PreStim = Pool{p}.xb.pre_stimulus_record_time*1e-3; %s
            PostStim = Pool{p}.xb.post_stimulus_record_time*1e-3; %s
            StimDur = Pool{p}.xb.stimulus_ch1(:,5)*1e-3;
            
            stim_info = Pool{p}.xb.data(find(Pool{p}.xb.data(:,3) == 1 & Pool{p}.xb.data(:,4) == -1),:);
            
            data_new = stim_info;
            nreps = Pool{p}.xb.stimulus_ch1(1,4);
            nStim = max(Pool{p}.xb.stimulus_ch1(:,1));
            %             nreps = 10;
            
            %
            % nStim = nStim + max(x2.stimulus_ch1(:,1));
            %
            TotalReps = nStim*nreps;
            false_start = length(Pool{p}.event_times)-TotalReps;
            start_stim_times = Pool{p}.event_times(false_start+1:end,1);
            end_stim_times = Pool{p}.event_times(false_start+1:end,2);
            spikes_pooled = [];
            raster.stim = [];
            raster.rep = [];
            raster.spikes = [];
            
            for rep = 1:TotalReps
                %for raster
                %                     if id ==4
                %                         spike_timesSU{id} = a.spike_times{4};
                %                     end
                spikes1 = Pool{p}.spiketimes(find(Pool{p}.spiketimes>=start_stim_times(rep)-PreStim & ...
                    Pool{p}.spiketimes<=end_stim_times(rep)+ PostStim)).';
                spikes1 = spikes1 - start_stim_times(rep);
                spikes_pooled = [spikes_pooled spikes1];
                raster.stim = [raster.stim data_new(rep,1)*ones(size(spikes1))];
                raster.rep = [raster.rep data_new(rep,2)*ones(size(spikes1))];
                raster.spikes = [raster.spikes spikes1];
                
                %for rate
                %             spikes2 = obj.spike_times{id}(find(obj.spike_times{id}>=start_stim_times(rep) & ...
                %                 obj.spike_times{id}<=end_stim_times(rep))).';
                %             rate.stim{id}(data_new(rep,1),data_new(rep,2)) = length(spikes2)/StimDur(data_new(rep,1));
                %             spikes3 = obj.spike_times{id}(find(obj.spike_times{id}<=start_stim_times(rep) & ...
                %                 obj.spike_times{id}>=start_stim_times(rep)-PreStim)).' ;
                %             rate.pre{id}(data_new(rep,1),data_new(rep,2)) = length(spikes3)/PreStim;
            end
            
            
            %            subplot(2,3,[3 6])
            figure
            for st = 1:length(StimDur)
                rectangle('Position',[0 nreps*(st-1),StimDur(st) nreps],'FaceColor',[0.9,0.9,0.9],'EdgeColor','none')
            end
            hold on
            plot(raster.spikes,nreps*(raster.stim-1)+raster.rep,'k.','MarkerSize',15);
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
        end
            
            
            