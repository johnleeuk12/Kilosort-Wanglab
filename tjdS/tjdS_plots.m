classdef tjdS_plots < handle
%{
Class for TJD GUI plots
    
This code is intended to serve as a external plotting code that can be called from the TJD GUI
    
%Principle author: Saul Meza
    
%%%%%%%%%%%%%% EDIT LOG %%%%%%%%%%%%%%%%%%%%%%%%
Last modification: May 4th 2021
Saul: Added functions to plot single units or two units in the GUI.
    
    %}
    
    
    properties
        channels = [];
        unit = [];
        waveforms_all = [];
        spike_times = [];
        ax_pca = [];
        ax_waveform = [];
        ax_waveform_2 = [];
        ax_isi = [];
        score = []; 
        idx_cls = []; 
        
    end
    
    methods
        %%
        function obj = tjdS_plots(channels,unit, score, idx_cls, waveforms_all, spike_times, ax_pca, ax_waveform, ax_waveform_2, ax_isi)
            
                obj.channels = channels;
                obj.unit = unit;
                obj.score = score;
                obj.idx_cls = idx_cls;
                obj.waveforms_all = waveforms_all;
                obj.spike_times = sort(spike_times);
                obj.ax_pca = ax_pca;
                obj.ax_waveform = ax_waveform;
                obj.ax_waveform_2 = ax_waveform_2;
                obj.ax_isi = ax_isi;
                
        end
        %%
        function obj = clear_all(obj) 
            % Function to clean axes titles and data
            
            cla(obj.ax_pca);
            cla(obj.ax_waveform);
            cla(obj.ax_isi);
            cla(obj.ax_waveform_2);
            
            delete(get(obj.ax_pca, 'children'));
            delete(get(obj.ax_waveform, 'children'));
            delete(get(obj.ax_isi, 'children'));
            delete(get(obj.ax_waveform_2, 'children'));
            
            delete(legend(obj.ax_pca));
            delete(legend(obj.ax_isi));
            
        end
        %%
        function obj = waveforms(obj)
            %Waveforms for splitting units
            id1=find(obj.idx_cls==0);
            id2=find(obj.idx_cls==1);
            
            if length(id1)<120
                idx1 = id1;
            else
                idx1 = randsample(id1,120);
            end

            for a = idx1
                plot(obj.ax_waveform,obj.waveforms_all(:,a));
            end
            plot(obj.ax_waveform,mean(obj.waveforms_all(:,id1),2), 'LineWidth',2,'Color','k');
            title(obj.ax_waveform, ['Unit ' num2str(obj.unit(1)) ', Best Channel: ' num2str(obj.channels(1))]);
            
            if length(id2)<120
                idx2 = id2;
            else
                idx2 = randsample(id2,120);
            end

            for b = idx2
                plot(obj.ax_waveform_2,obj.waveforms_all(:,b));
            end
            plot(obj.ax_waveform_2,mean(obj.waveforms_all(:,id2),2), 'LineWidth',2,'Color','k');
            title(obj.ax_waveform_2, ['Unit ' num2str(obj.unit(2)) ', Best Channel: ' num2str(obj.channels(2))]);

            disp('Waveform plot done');
        end
        %%
        function obj = waveform(obj)
            id1=find(obj.idx_cls==1);
            index = id1;
            
            if length(index)<120
                idx1 = index;
            else
                idx1 = randsample(index,120);
            end

            for a = idx1
                plot(obj.ax_waveform,obj.waveforms_all(:,a));
            end
            disp(obj.unit);
            plot(obj.ax_waveform,mean(obj.waveforms_all(:,index),2), 'LineWidth',2,'Color','k');
            title(obj.ax_waveform, ['Unit ' num2str(obj.unit(1)) ', Best Channel ' num2str(obj.channels)]);
        end
        %%
        function obj = isis(obj)
            %ISI plots for splitting units
            id1 = find(obj.idx_cls==0);
            id2 = find(obj.idx_cls==1);
            
            spike_times_1 = obj.spike_times(id1);
            spike_times_2 = obj.spike_times(id2);
            
            ISI_1 = spike_times_1(2:end)-spike_times_1(1:end-1);
            ISI_2 = spike_times_2(2:end)-spike_times_2(1:end-1);
            
            edges = -3.5:0.25:2;
            [N1,~] = histcounts(ISI_1,10.^edges);
            [N2,~] = histcounts(ISI_2,10.^edges);
            
            refract1 = sum(N1(1:2))*100/sum(N1);
            refract2 = sum(N2(1:2))*100/sum(N2);
            
            histogram(obj.ax_isi,ISI_1,10.^edges);
            hold on;
            histogram(obj.ax_isi,ISI_2,10.^edges);
            
            set(obj.ax_isi,'xscale','log')
            xticks(obj.ax_isi,[1E-4 1E-3 1E-2 1E-1 1 10]);
            xticklabels(obj.ax_isi,{0.0001, 0.001, 0.01, 0.1, 1, 10})
            xlabel(obj.ax_isi,'time (ms)')
            title(obj.ax_isi,'ISI');
            %legend(obj.ax_isi, ['Cluster 1: ' num2str(refract1) '% Cluster 2: ' num2str(refract2)]);
            legend(obj.ax_isi,[num2str(refract1) '%' ' Unit ' num2str(obj.unit(1))],[num2str(refract2) '%' ' Unit ' num2str(obj.unit(2))], 'Location', 'best');
            disp('ISI plot done');
            
        end
        %%
        function obj = isi(obj)
            ISI = obj.spike_times(2:end)-obj.spike_times(1:end-1);
            
            edges = -3.5:0.25:2;
            [NN,~] = histcounts(ISI,10.^edges);
            refract = sum(NN(1:2))*100/sum(NN);
            histogram(obj.ax_isi,ISI,10.^edges);
            set(obj.ax_isi,'xscale','log');
            xticks(obj.ax_isi,[1E-4 1E-3 1E-2 1E-1 1 10]);
            xticklabels(obj.ax_isi,{0.0001, 0.001, 0.01, 0.1, 1, 10});
            xlabel(obj.ax_isi,'time (ms)');
            title(obj.ax_isi,['ISI ' num2str(refract) '%']);
        end
        %%
        function obj = pcas(obj)
            %PCA plotting for splitting units
            
            gscatter(obj.ax_pca,obj.score(:,1),obj.score(:,2),obj.idx_cls, 'rb');

            title(obj.ax_pca, 'PCA');
            legend(obj.ax_pca, ['Unit ' num2str(obj.unit(1))], ['Unit ' num2str(obj.unit(2))], 'Location', 'best');
            xlabel(obj.ax_pca,'a2');
            ylabel(obj.ax_pca,'b2');
            
            disp('PCA plot done');
            
        end   
        %%
        function obj = pca(obj)
            gscatter(obj.ax_pca,obj.score(:,1),obj.score(:,2),obj.idx_cls, 'r');
            
            title(obj.ax_pca, 'PCA');
            legend(obj.ax_pca, ['Unit ' num2str(obj.unit(1))], 'Location', 'best');
            xlabel(obj.ax_pca,'a2');
            ylabel(obj.ax_pca,'b2');
            
            disp('PCA plot done');
            
        end
    end
end