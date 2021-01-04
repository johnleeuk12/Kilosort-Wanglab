classdef tjxS_plots < handle
%{
Class for TJX GUI plots
    
This code is intended to serve as a external plotting code that can be called from the TJX GUI
    
%Principle author: Saul Meza
    
%%%%%%%%%%%%%% EDIT LOG %%%%%%%%%%%%%%%%%%%%%%%%
Last modification: Dec. 1, 2020
    %Created different functions for each type of plot
    
12/3/2020 JHL
    added functions raster2 and tuning2 for 2D Stim. function select stim
    is a utility function for these 2 functions
    %}
    
    
    properties
        
        out; %Output of tjxS
        out_su; %Output from Single_units only
        val_1; %Cluster selection 1
        val_2; %Cluster selection 2
        ax1;
        ax2;
        ax3;
        ax4;
        ax5;
        ax6;
        stim_info;
        %idx;
        %X;
    end
    
    methods
        %%
        function obj = tjxS_plots(stim_info,out,out_su,val1,val2,ax1,ax2,ax3,ax4,ax5,ax6)
            
                obj.out = out;
                obj.out_su = out_su;
                obj.val_1 = val1;
                obj.val_2 = val2;
                obj.ax1 = ax1;
                obj.ax2 = ax2;
                obj.ax3 = ax3;
                obj.ax4 = ax4;
                obj.ax5 = ax5;
                obj.ax6 = ax6;
                obj.stim_info = stim_info;
        end
        %%
        function obj = clear_all(obj) 
            % Function to clean axes titles and data
            
                delete(get(obj.ax1, 'children'));
                delete(get(obj.ax2, 'children'));
                delete(get(obj.ax3, 'children'));
                delete(get(obj.ax4, 'children'));
                delete(get(obj.ax5, 'children'));
                delete(get(obj.ax6, 'children'));
                
                delete(legend(obj.ax1));
                delete(legend(obj.ax2));
                delete(legend(obj.ax3));
                delete(legend(obj.ax4));
                
                delete(title(obj.ax1, ['Waveform Cluster ' obj.val_1]));
                delete(title(obj.ax2,['PCA Cluster' obj.val_1]));
                delete(title(obj.ax3,['Tuning Curve Cluster ' obj.val_1]));
                delete(title(obj.ax4,'ISI'));
                delete(title(obj.ax5,['Rasterplot Cluster ' obj.val_1]));
                delete(title(obj.ax6,['Rasterplot Cluster ' obj.val_2]));
            
            
        end
        %%
        function obj = waveform(obj,plt,ids)
            %Plotting of waveforms
            
            if length(ids) == 1
           
                if plt == 1
                    value = obj.val_1; %Cluster value
                elseif plt == 2
                    value = obj.val_2; %Cluster value
                end

                % WAVEFORM PLOT
                time_step = [1:60]/obj.out.params.fs*1e3;
                Nb_wave = min(size(obj.out.spike_times{ids},1),100);
                idx = 1:size(obj.out.spike_times{ids},1);
                rand_select = randsample(idx,Nb_wave);

                for n = 1:length(rand_select)

                    plot(obj.ax1,time_step,obj.out.waveforms.raw{ids}(obj.out.templates{ids}.best_ch,:,rand_select(n)));
                    hold on

                end
                plot(obj.ax1,time_step,obj.out.waveforms.mean{ids},'Linewidth',2,'Color','k');
                title(obj.ax1, ['Waveform Cluster: ', value, ', Best Ch: ', num2str(obj.out.templates{ids}.best_ch)]);
                xlabel(obj.ax1, 'Time (ms)');
                ylabel(obj.ax1, 'uV');
                hold off
            
            elseif length(ids) == 2
                
                time_step = [1:60]/obj.out.params.fs*1e3;
                plot(obj.ax1,time_step,obj.out.waveforms.mean{ids(1)},'Linewidth',2, 'Color', 'r');
                hold on
                plot(obj.ax1,time_step,obj.out.waveforms.mean{ids(2)},'Linewidth',2, 'Color', 'b');
                title(obj.ax1,'Waveform');
                xlabel(obj.ax1, 'Time (ms)');
                ylabel(obj.ax1, 'uV');
                legend(obj.ax1, obj.val_1, obj.val_2);
                hold off
                
            end
            
            
        end
        %%
        function obj = isi_plot(obj,plt, ids)
            %ISI plots
            if length(ids) == 1
            
                if plt == 1
                    value = obj.val_1; %Cluster value
                    color = 'r';
                elseif plt == 2
                    value = obj.val_2; %Cluster value
                    color = 'b';
                end

                ISI = obj.out.spike_times{ids}(2:end)-obj.out.spike_times{ids}(1:end-1);

                edges = -3.5:0.25:2;
                [N,~] = histcounts(ISI,10.^edges);
                refract = sum(N(1:2))*100/sum(N);
                histogram(obj.ax4,ISI,10.^edges, 'FaceColor', color)
                set(obj.ax4,'xscale','log')
                xticks(obj.ax4,[1E-4 1E-3 1E-2 1E-1 1 10]);
                xticklabels(obj.ax4,{0.0001, 0.001, 0.01, 0.1, 1, 10})
                xlabel(obj.ax4,'Time (ms)')
                title(obj.ax4,['ISI: ' num2str(refract) '%, Cluster ' value]);
                
            elseif length(ids) == 2
                
                ISI = obj.out.spike_times{ids(1)}(2:end)-obj.out.spike_times{ids(1)}(1:end-1);
                ISI2 = obj.out.spike_times{ids(2)}(2:end)-obj.out.spike_times{ids(2)}(1:end-1);
                
                edges = -3.5:0.25:2;
                [N,~] = histcounts(ISI,10.^edges);
                [N2,~] = histcounts(ISI2,10.^edges);
                refract = sum(N(1:2))*100/sum(N);
                refract2 = sum(N2(1:2))*100/sum(N2);
                histogram(obj.ax4,ISI,10.^edges, 'FaceColor', 'r');
                hold on
                histogram(obj.ax4,ISI2,10.^edges, 'FaceColor', 'b');
                set(obj.ax4,'xscale','log')
                xticks(obj.ax4,[1E-4 1E-3 1E-2 1E-1 1 10]);
                xticklabels(obj.ax4,{0.0001, 0.001, 0.01, 0.1, 1, 10})
                xlabel(obj.ax4,'Time (ms)')
                title(obj.ax4,'ISI');
                legend(obj.ax4,[num2str(refract) '%' ' (cluster ' obj.val_1 ')'],[num2str(refract2) '%' ' (cluster ' obj.val_2 ')']);
                hold off
                
            end
            
        end
        %%
        function obj = pca_plot(obj, plt, ids)
            %PCA plotting
            
            if length(ids) == 1
                if plt == 1
                    value = obj.val_1; %Cluster value
                    color = 'rs';
                elseif plt == 2
                    value = obj.val_2; %Cluster value
                    color = 'bs';
                end

                Y = [];
                Y(:,:) = obj.out.waveforms.raw{ids}(obj.out.templates{ids}.best_ch,:,:);
                Y = Y.';
                [~,score,~] = pca(Y);

                X = [score(:,1) score(:,2)];
                idx = kmeans(X,2);
                %obj.idx=idx;
                %obj.X=X;
                %disp(obj.X(obj.idx==2,1));
                %a=X(:,1);
                %b=X(:,2);
                %gscatter(obj.ax2,a,b,idx,'rb','+o',5);
                plot(obj.ax2,X(idx==1,1),X(idx==1,2),color,X(idx==2,1),X(idx==2,2),color);
                %plot(obj.ax2,X(idx==1,1),X(idx==1,2),'r+');
                title(obj.ax2, ['PCA Cluster ' value]);
                xlabel(obj.ax2,'a2');
                ylabel(obj.ax2,'b2');
                
            elseif length(ids) == 2
               
                Y1 = [];
                Y1(:,:) = obj.out.waveforms.raw{ids(1)}(obj.out.templates{ids(1)}.best_ch,:,:);
                Y1 = Y1.';
                [~,score1,~] = pca(Y1);
                
                X1 = [score1(:,1) score1(:,2)];
                idx1 = kmeans(X1,2);
                
                Y2 = [];
                Y2(:,:) = obj.out.waveforms.raw{ids(2)}(obj.out.templates{ids(2)}.best_ch,:,:);
                Y2 = Y2.';
                [~,score2,~] = pca(Y2);
                
                X2 = [score2(:,1) score2(:,2)];
                idx2 = kmeans(X2,2);

                plot(obj.ax2,X1(idx1==1,1),X1(idx1==1,2),'r*',X1(idx1==2,1),X1(idx1==2,2),'r*');
                %a1=X1(:,1);
                %b1=X1(:,2);
                %gscatter(obj.ax2,a1,b1,idx1,'r','+',5);
                hold on;
                plot(obj.ax2,X2(idx2==1,1),X2(idx2==1,2),'bo',X2(idx2==2,1),X2(idx2==2,2),'bo');
                %a2=X2(:,1);
                %b2=X2(:,2);
                %gscatter(obj.ax2,a2,b2,idx2,'b','o',5);
                
                title(obj.ax2, ['PCA']);
                xlabel(obj.ax2,'a2');
                ylabel(obj.ax2,'b2');
                
            end
            
            
        end
        %%       
        function obj = tuning_curve(obj, plt, ids)
            stim_label = obj.stim_info.tuning.x_data;
            if length(ids) == 1
                
                if plt == 1
                    value = obj.val_1; %Cluster value
                    color = 'r';
                elseif plt == 2
                    value = obj.val_2; %Cluster value
                    color = 'b';
                end
                
                errorbar(obj.ax3,stim_label,obj.out_su.SUrate{ids}.mean,obj.out_su.SUrate{ids}.error,'LineWidth',2, 'Color', color);
                hold on
                plot(obj.ax3,stim_label,ones(1,length(stim_label))*obj.out_su.SUrate{ids}.spont,'--k');
                
                legend(obj.ax3, value);
                
            elseif length(ids) == 2
                
                
                errorbar(obj.ax3,stim_label,obj.out_su.SUrate{ids(1)}.mean,obj.out_su.SUrate{ids(1)}.error,'LineWidth',2, 'Color', 'r');
                hold on
                errorbar(obj.ax3,stim_label,obj.out_su.SUrate{ids(2)}.mean,obj.out_su.SUrate{ids(2)}.error,'LineWidth',2, 'Color', 'b');
                plot(obj.ax3,stim_label,ones(1,length(stim_label))*obj.out_su.SUrate{ids(1)}.spont,'--k');
                plot(obj.ax3,stim_label,ones(1,length(stim_label))*obj.out_su.SUrate{ids(2)}.spont,'--k');
                legend(obj.ax3, obj.val_1, obj.val_2);
                title(obj.ax3,['Tuning Curve Cluster ' obj.val_1 ' and ' obj.val_2]);
                
            end
            xlabel(obj.ax3,obj.stim_info.tuning.x_label{1});
            xticks(obj.ax3,obj.stim_info.tuning.x_ticks(1:end));
            ylabel(obj.ax3,'Firing Rate (spikes/s)');
            xtickangle(obj.ax3,45);
            switch  obj.out.x.analysis_code
                case {1,10,62}
                    set(obj.ax3,'xscale','log')
%                     xtickangle(obj.ax3,45);
            end
            

%             if strcmp(plot_type,'Tuning') == 1
%                 
%                 %Tuning curve
%                 errorbar(obj.ax3,stim_label,SUrate{ids}.mean,SUrate{ids}.error,'LineWidth',2, 'Color', 'r');
%                 hold on
%                 
%                 plot(obj.ax3,stim_label,ones(1,length(stim_label))*SUrate{ids}.spont,'--k');
%                 set(obj.ax3,'xscale','log')
%                 %                             set(gca,'xscale','log')
%                 xticks(obj.ax3,round((stim_label(1:4:end).')*1e2)*1e-2);
%                 xtickangle(obj.ax3,45);
%                 xlabel(obj.ax3,'Hz');
%                 ylabel(obj.ax3,'Firing Rate (spikes/s)');
%                 title(obj.ax3,['Tuning Curve Cluster ' value]);
%                 %legend(obj.ax3, value);
%                 
%             elseif strcmp(plot_type,'VT') == 1
%                 stim_number = x.stimulus_ch1(:,1);
%                 stim_length = x.stimulus_ch1(:,5);
%                 for i = 9:22
%                     VT_para(i-8) = length(unique(x.stimulus_ch1(:,i)));
%                 end
%                 stim_label = x.stimulus_ch1(:,find(VT_para ~= 1)+8);
%                 if analysis_code == 2320 || analysis_code == 2335
%                     stim_label = x.stimulus_ch1(:,10);
%                 end
%                 
%                 %Tuning Curve
%                 SUrate{ids}.mean_fr = SUrate{ids}.mean./stim_length;
%                 errorbar(obj.ax3,stim_label,SUrate{ids}.mean,SUrate{ids}.error,'LineWidth',2, 'Color', 'r');
%                 hold on
%                 plot(obj.ax3,stim_label,ones(1,length(stim_label))*SUrate{ids}.spont,'--k');
%                 %                     xticks(stim_number);
%                 %                     xticklabels(stim_label)
%                 xtickangle(obj.ax3,45);
%                 %                 xlabel('Hz')
%                 ylabel(obj.ax3,'Firing Rate (spikes/s)');
%                 title(obj.ax3, ['Tuning Curve Cluster ' value]);
%                 legend(obj.ax3, value);
%                 
%             elseif strcmp(plot_type,'User') == 1
%                 
%                 %Tuning Curve
%                 stim_label = {};
%                 stim_number = x.stimulus_ch1(:,1);
%                 stim_length = x.stimulus_ch1(:,5);
%                 for st = 1: length(stim_number)
%                     idx1 = strfind(x.user_stimulus_desc_ch1{st},': ')+2;
%                     idx2 = strfind(x.user_stimulus_desc_ch1{st},'_m')-1;
%                     idx3 = strfind(x.user_stimulus_desc_ch1{st},'rev');
%                     stim_label{st} = x.user_stimulus_desc_ch1{st}(idx1:idx2);
%                     if isempty(idx3) ~=1
%                         stim_label{st} = [stim_label{st} '_rev'];
%                     end
%                     
%                 end
%                 SUrate{ids}.mean_fr = SUrate{ids}.mean./stim_length;
%                 errorbar(obj.ax3,stim_number,SUrate{ids}.mean,SUrate{ids}.error,'LineWidth',2, 'Color', 'r');
%                 hold on
%                 plot(obj.ax3,stim_number,ones(1,length(stim_label))*SUrate{ids}.spont,'--k');
%                 xticks(obj.ax3,stim_number);
%                 xticklabels(obj.ax3,stim_label);
%                 xtickangle(obj.ax3,45);
%                 %                 xlabel('Hz')
%                 ylabel(obj.ax3,'Firing Rate (spikes/s)');
%                 title(obj.ax3, ['Tuning Curve Cluster ' value]);
%                 %legend(obj.ax3, value);
%                     
%             end        
                    
            
        end
        %%
        function obj = raster(obj, plt, ids)
            
            nreps = obj.out_su.nreps;
            StimDur = obj.out_su.StimDur;
            raster = obj.out_su.raster;
            TotalReps = obj.out_su.TotalReps;
            PreStim = obj.out_su.PreStim;
            PostStim = obj.out_su.PostStim;
            
            if length(ids) == 1
                
                if plt == 1
                    rast_ax = obj.ax5; %Rasterplot axes
                    value = obj.val_1; %Cluster value
                elseif plt == 2
                    rast_ax = obj.ax6; %Rasterplot axes
                    value = obj.val_2; %Cluster value
                end
                
                
                for st = 1:length(StimDur)
                    rectangle(rast_ax,'Position',[0 nreps*(st-1),StimDur(st) nreps],'FaceColor',[0.9,0.9,0.9],'EdgeColor','none')
                end
                hold(rast_ax,'on')
                plot(rast_ax,raster.spikes{ids},nreps*(raster.stim{ids}-1)+raster.rep{ids},'k.','MarkerSize',15);
                xlabel(rast_ax,'Time (s)')
                axis(rast_ax,[-PreStim max(StimDur) + PostStim 0 TotalReps+1])
                hold(rast_ax,'off')
                title(rast_ax,['Rasterplot Cluster ' value]);
                ylabel(rast_ax,obj.stim_info.raster.y_label{1});
                yticks(rast_ax,obj.stim_info.raster.y_ticks);
                yticklabels(rast_ax,obj.stim_info.raster.y_ticklabels);
                
                drawnow()
                
            elseif length(ids) == 2
                
                % first plot
                for st = 1:length(StimDur)
                    rectangle(obj.ax5,'Position',[0 nreps*(st-1),StimDur(st) nreps],'FaceColor',[0.9,0.9,0.9],'EdgeColor','none')
                end
                hold(obj.ax5,'on')
                plot(obj.ax5,raster.spikes{ids(1)},nreps*(raster.stim{ids(1)}-1)+raster.rep{ids(1)},'k.','MarkerSize',15);
                xlabel(obj.ax5,'Time (s)')
                ylabel(obj.ax5,obj.stim_info.raster.y_label{1});
                yticks(obj.ax5,obj.stim_info.raster.y_ticks);
                yticklabels(obj.ax5,obj.stim_info.raster.y_ticklabels);
                axis(obj.ax5,[-PreStim max(StimDur) + PostStim 0 TotalReps+1])
                hold(obj.ax5,'off')
                title(obj.ax5,['Rasterplot Cluster ' obj.val_1]);
                drawnow()
                
                % second plot
                for st = 1:length(StimDur)
                    rectangle(obj.ax6,'Position',[0 nreps*(st-1),StimDur(st) nreps],'FaceColor',[0.9,0.9,0.9],'EdgeColor','none')
                end
                hold(obj.ax6,'on')
                plot(obj.ax6,raster.spikes{ids(2)},nreps*(raster.stim{ids(2)}-1)+raster.rep{ids(2)},'k.','MarkerSize',15);
                xlabel(obj.ax6,'Time (s)')
                ylabel(obj.ax6,obj.stim_info.raster.y_label{1});
                yticks(obj.ax6,obj.stim_info.raster.y_ticks);
                yticklabels(obj.ax6,obj.stim_info.raster.y_ticklabels);
                axis(obj.ax6,[-PreStim max(StimDur) + PostStim 0 TotalReps+1])
                hold(obj.ax6,'off')
                title(obj.ax6,['Rasterplot Cluster ' obj.val_2]);
                drawnow()
                
                
            end
            
            
        end
        %%
        function obj = raster2(obj, plt, ids,Stim1Para, Stim2Para) %raster for 2D stim
            nreps = obj.out_su.nreps;
            StimDur = obj.out_su.StimDur;
            raster = obj.out_su.raster;
            TotalReps = obj.out_su.TotalReps;
            PreStim = obj.out_su.PreStim;
            PostStim = obj.out_su.PostStim;
            
            %             if length(ids) == 1
            
            if plt == 1
                rast_ax = obj.ax5; %Rasterplot axes
                value = obj.val_1; %Cluster value
            elseif plt == 2
                rast_ax = obj.ax6; %Rasterplot axes
                value = obj.val_2; %Cluster value
            end
            
            [new_vec, new_vec2] = obj.select_stim(raster, ids, Stim1Para, Stim2Para);
            
            new_raster.spikes = raster.spikes{ids}(new_vec2);
            new_raster.stim1 = raster.stim{ids}(new_vec2);
            new_raster.rep = raster.rep{ids}(new_vec2);
            new_raster.stim2 = zeros(1,length(new_raster.stim1));
            for stt = 1:length(new_raster.stim1)
                new_raster.stim2(1,stt) = find(unique(new_raster.stim1) == new_raster.stim1(1,stt));
            end
            
            
            for st = 1: length(new_vec)
                rectangle(rast_ax,'Position',[0 nreps*(st-1),StimDur(new_vec(st)) nreps],'FaceColor',[0.9,0.9,0.9],'EdgeColor','none')
            end
            hold(rast_ax,'on')
            plot(rast_ax,new_raster.spikes,nreps*(new_raster.stim2-1)+new_raster.rep,'k.','MarkerSize',15);
            xlabel(rast_ax,'Time (s)')
            ylabel(rast_ax,'kHz')
            axis(rast_ax,[-PreStim max(StimDur) + PostStim 0 nreps*length(new_vec)+1])
            hold(rast_ax,'off')
            title(rast_ax,['Rasterplot Cluster ' value]);
            
            
            if Stim1Para ~= 0 % if the 1st parameter was set the 2nd parameter varies
                yticks(rast_ax,[1:nreps:nreps*length(new_vec)]+floor(nreps/2)) ;
                yticklabels(rast_ax, num2str(obj.stim_info.y_list))
                ylabel(rast_ax,obj.out.x.stimulus_tags_ch1{obj.stim_info.idx(2)});
                
                
            elseif Stim2Para ~= 0
                yticks(rast_ax,[1:nreps:nreps*length(new_vec)]+floor(nreps/2)) ;
                yticklabels(rast_ax, num2str(obj.stim_info.x_list))
                ylabel(rast_ax,obj.out.x.stimulus_tags_ch1{obj.stim_info.idx(1)});
            end
            
            drawnow()
                
%             elseif length(ids) == 2
%                 
%                 %first plot
%                 [new_vec, new_vec2] = obj.select_stim(raster, ids(1),Stim1Para, Stim2Para);
%                 
%                 new_raster.spikes = raster.spikes{ids(1)}(new_vec2);
%                 new_raster.stim1 = raster.stim{ids(1)}(new_vec2);
%                 new_raster.rep = raster.rep{ids(1)}(new_vec2);
%                 new_raster.stim2 = zeros(1,length(new_raster.stim1));
%                 for stt = 1:length(new_raster.stim1)
%                     new_raster.stim2(1,stt) = find(unique(new_raster.stim1) == new_raster.stim1(1,stt));
%                 end
%                 
%                 
%                 for st = 1: length(new_vec)
%                     rectangle(obj.ax5,'Position',[0 nreps*(st-1),StimDur(new_vec(st)) nreps],'FaceColor',[0.9,0.9,0.9],'EdgeColor','none')
%                 end
%                 hold(obj.ax5,'on')
%                 plot(obj.ax5,new_raster.spikes,nreps*(new_raster.stim2-1)+new_raster.rep,'k.','MarkerSize',15);
%                 xlabel(obj.ax5,'Time (s)')
%                 ylabel(obj.ax5,'kHz')
%                 axis(obj.ax5,[-PreStim max(StimDur) + PostStim 0 nreps*length(new_vec)+1])
%                 hold(obj.ax5,'off')
%                 title(obj.ax5,['Rasterplot Cluster ' obj.val_1]);
%                 
%                 
%                 if obj.Stim1Para.Value ~= 0 % if the 1st parameter was set the 2nd parameter varies
%                     yticks(obj.ax5,[1:nreps:nreps*length(new_vec)]+floor(nreps/2)) ;
%                     yticklabels(obj.ax5, num2str(obj.stim_info.y_list))
%                     ylabel(obj.ax5,obj.out.x.stimulus_tags_ch1{obj.stim_info.idx(2)});
%                     
%                     
%                 elseif app.Stim2Para.Value ~= 0
%                     yticks(obj.ax5,[1:nreps:nreps*length(new_vec)]+floor(nreps/2)) ;
%                     yticklabels(obj.ax5, num2str(obj.stim_info.x_list))
%                     ylabel(obj.ax5,obj.out.x.stimulus_tags_ch1{obj.stim_info.idx(1)});
%                 end
%                 
%                 drawnow()
%                 
%                 
%                 
%                 % second plot
%                 [new_vec, new_vec2] = obj.select_stim(raster, ids(2),Stim1Para, Stim2Para);
%                 
%                 new_raster.spikes = raster.spikes{ids(2)}(new_vec2);
%                 new_raster.stim1 = raster.stim{ids(2)}(new_vec2);
%                 new_raster.rep = raster.rep{ids(2)}(new_vec2);
%                 new_raster.stim2 = zeros(1,length(new_raster.stim1));
%                 for stt = 1:length(new_raster.stim1)
%                     new_raster.stim2(1,stt) = find(unique(new_raster.stim1) == new_raster.stim1(1,stt));
%                 end
%                 
%                 
%                 for st = 1: length(new_vec)
%                     rectangle(obj.ax6,'Position',[0 nreps*(st-1),StimDur(new_vec(st)) nreps],'FaceColor',[0.9,0.9,0.9],'EdgeColor','none')
%                 end
%                 hold(obj.ax6,'on')
%                 plot(obj.ax6,new_raster.spikes,nreps*(new_raster.stim2-1)+new_raster.rep,'k.','MarkerSize',15);
%                 xlabel(obj.ax6,'Time (s)')
%                 ylabel(obj.ax6,'kHz')
%                 axis(obj.ax6,[PreStim max(StimDur) + PostStim 0 nreps*length(new_vec)+1])
%                 hold(obj.ax6,'off')
%                 title(obj.ax6,['Rasterplot Cluster ' obj.val_1]);
%                 
%                 
%                 if obj.Stim1Para.Value ~= 0 % if the 1st parameter was set the 2nd parameter varies
%                     yticks(obj.ax6,[1:nreps:nreps*length(new_vec)]+floor(nreps/2)) ;
%                     yticklabels(obj.ax6, num2str(obj.stim_info.y_list))
%                     ylabel(obj.ax6,obj.out.x.stimulus_tags_ch1{obj.stim_info.idx(2)});
%                     
%                     
%                 elseif app.Stim2Para.Value ~= 0
%                     yticks(obj.ax6,[1:nreps:nreps*length(new_vec)]+floor(nreps/2)) ;
%                     yticklabels(obj.ax6, num2str(obj.stim_info.x_list))
%                     ylabel(obj.ax6,obj.out.x.stimulus_tags_ch1{obj.stim_info.idx(1)});
%                 end
%                 
%                 drawnow()
%                 
%                 
%                 
%             end
            
            
        end
        %%
        function obj = tuning2(obj, plt, ids,a,b) %tuning plot for 2D stim
            if length(ids) == 1
%                 
%                 if plt == 1
%                     value = obj.val_1; %Cluster value
%                     color = 'r';
%                 elseif plt == 2
%                     value = obj.val_2; %Cluster value
%                     color = 'b';
%                 end
%                 
                
                M  = obj.out_su.SUrate{ids}.mean;
                M = reshape(M,[length(obj.stim_info.x_list),length(obj.stim_info.y_list)]);
                if length(ids) == 1
                    %              imagesc(obj.ax3,obj.stim_info.x_list,obj.stim_info.y_list,M);
                    imagesc(obj.ax3,obj.stim_info.y_list,obj.stim_info.x_list,M.');
                    colorbar(obj.ax3);
                    
                end
                
            else
            end
            
        end
        
        function obj = RHSplot(obj, plt,ids)
            freq_set_ind = find(strcmp(obj.out.x.stimulus_tags_ch1,' freq sets')==1);
            freq_set = obj.out.x.stimulus_ch1(1,freq_set_ind);
            atten_ind = find(strcmp(obj.out.x.stimulus_tags_ch1,' levels')==1);
            D = obj.out.x.stimulus_ch1(1:end,atten_ind); %atten
            R = obj.out_su.SUrate{ids(1)}.mean; % stim rate
            R = R.';
            W = (D.'*D)\D.'*R.';
            CF = obj.out.x.stimulus_ch1(1,8);
            BF = obj.out.x.stimulus_ch1(1,end);
            if ~isempty(W)
                plot(obj.ax3,freq_set,W,'-b','LineWidth',2)
                hold(obj.ax3,'on')
                plot(obj.ax3,freq_set, zeros(1,length(freq_set)),'--k','LineWidth',1);
                drawnow
                plot(obj.ax3,CF,W(find(freq_set == CF)),'or','MarkerSize',10);
                hold(obj.ax3,'off')
            end
            
        
        end
        
        
        function [new_vec, new_vec2] = select_stim(obj,raster, ids,Stim1Para, Stim2Para)
%             
%             if isempty(obj.Stim1Para.Value) && isempty(obj.Stim2Para)
%                 new_vec2 = [1:length(raster.stim{ids})].';
%                 new_vec = obj.out.x.stimulus_ch1(:,1);
            if Stim1Para == 0 && Stim2Para == 0
                new_vec2 = [1:length(raster.stim{ids})].';
                new_vec = obj.out.x.stimulus_ch1(:,1);
            elseif Stim1Para ~= 0
                new_vec = find(obj.out.x.stimulus_ch1(:,obj.stim_info.idx(1)) == Stim1Para).';
                new_vec2 = [];
                for v = new_vec
                    new_vec2 = [new_vec2, find(raster.stim{ids} == v)];
                end
            elseif Stim2Para ~= 0
                new_vec = find(obj.out.x.stimulus_ch1(:,obj.stim_info.idx(2)) == Stim2Para).';
                new_vec2 = [];
                for v = new_vec
                    new_vec2 = [new_vec2, find(raster.stim{ids} == v)];
                end
            end
            
            
        end
    end
end