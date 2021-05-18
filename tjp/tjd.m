classdef tjd < handle
    
    properties(Constant)
    end
    
    properties
        params = eval('parameters_distance');
        list_data = {};
        %         d_template = {};
        merge_pair = [];
        SU_keep = [];
        SU_good = [];
        L = [];
        N = [];
        SU_clear = [];
        Na = [];
        Nb = [];
        sequence = [];
        cid = [];
    end
    methods(Static)
    end
    
    methods
        function obj = Cdistance(obj)
            gpuDevice()
            tic
            obj.list_data = load([obj.params.animal_name '_unit_list.mat']);
            list_session = unique(obj.list_data.unit_list.data(:,4));
            d_template = {};
            for seg_ls = 1:length(obj.params.segment_list)-1
                fprintf(['Time %3.0fs. calculating distance between ' obj.params.segment_list{seg_ls} ' and ' obj.params.segment_list{seg_ls+1}  ' ... \n'], toc);
                seg_1 = find(strcmp(list_session,obj.params.segment_list{seg_ls}));
                seg_2 = find(strcmp(list_session,obj.params.segment_list{seg_ls+1}));
                
                A = find(strcmp(obj.list_data.unit_list.data(:,4),list_session{seg_1}));
                B = find(strcmp(obj.list_data.unit_list.data(:,4),list_session{seg_2}));
                %             centroids = {};
                pool.waveforms = [];
                pool.index = [];
                pool.channels = [];
                pool.templates = {};
                pool.templates.A = {};
                pool.templates.B = {};
                %                 pool.id = [];
                
                
                for i = 1:length(A)
                    unit_file_name = [obj.params.animal_name 'u00000'];
                    unit_file_name = [unit_file_name(1:end-size(num2str(A(i)),2)) num2str(A(i)) '.mat'];
                    x = load(unit_file_name);
                    try
                        pool.waveforms = cat(3,pool.waveforms,x.s_unit.waveforms{1});
                    catch                  
                        tempmat = x.s_unit.waveforms;
                        x.s_unit.waveforms = {};
                        x.s_unit.waveforms{1} = [];
                        x.s_unit.waveforms{1} = tempmat;
                        pool.waveforms = cat(3,pool.waveforms,x.s_unit.waveforms{1}); % weirdness where some waveforms are saved in cell and others in double...
                    end
                    pool.index = [pool.index;A(i)*ones(size(x.s_unit.waveforms{1},1),1)];
                    pool.id(i) = x.s_unit.cid;
                    
                    if ~isempty(x.s_unit.templates)
                        pool.channels = [pool.channels; [A(i) x.s_unit.templates.best_ch]];
                        if exist('x.s_unit.templates.data')
                            pool.templates.A{i} = x.s_unit.templates.data;
                        else
                            data = zeros(64,60);
                            %                  obj.templates{id}.best_ch = [];
                            for ch = 1:64
                                temp = [];
                                temp(:,:) = x.s_unit.waveforms{1}(ch,:,:);
                                [Y,Sig,X] = svd(gpuArray(temp),'econ');
                                %                 sig = diag(Sig);%figure; semilogy(sig(sig>1),'kx-')
                                k = 1:3;
                                P = Y(:,k)*Sig(k,k)*X(:,k)';
                                data(ch,:) = mean(gather(P),2).';
                            end
                            pool.templates.A{i}  = data;
                        end
                    end
                    %                 SpikeData = x.s_unit.waveforms{1};
                    %                 SpikeData= SpikeData';
                    %                 SpikeNum = length(SpikeData);
                    %                 % SpikeData = SpikeData - repmat(mean(SpikeData,2),1,SpikeNum);
                    %                 [Y,Sig,X] = svd(SpikeData,'econ');
                    %                 sig = diag(Sig);%figure; semilogy(sig(sig>1),'kx-')
                    %                 k = 1:3;
                    %                 P = Y(:,k)*Sig(k,k)*X(:,k)';
                    %                 P = P';
                    %                 WTempA(i,:) = mean(P,1);
                end
                
                obj.cid = pool.id; 
                x = load(unit_file_name);
                obj.SU_good.A = x.s_unit.SU_good;
                
                %                 fprintf('Time %3.0fs. decompositing segment 2... \n', toc);
                for i = 1:length(B)
                    unit_file_name = [obj.params.animal_name 'u00000'];
                    unit_file_name = [unit_file_name(1:end-size(num2str(B(i)),2)) num2str(B(i)) '.mat'];
                    x = load(unit_file_name);
                    
                    try
                        pool.waveforms = cat(3,pool.waveforms,x.s_unit.waveforms{1});
                    catch
                        tempmat = x.s_unit.waveforms;
                        x.s_unit.waveforms = {};
                        x.s_unit.waveforms{1} = [];
                        x.s_unit.waveforms{1} = tempmat;
                        pool.waveforms = cat(3,pool.waveforms,x.s_unit.waveforms{1}); % weirdness where some waveforms are saved in cell and others in double...
                    end
                    
                    pool.index = [pool.index;B(i)*ones(size(x.s_unit.waveforms{1},1),1)];
                    %                     pool.id = [pool.id ; x.s_unit.cid];
                    
                    if ~isempty(x.s_unit.templates)
                        pool.channels = [pool.channels; [B(i) x.s_unit.templates.best_ch]];
                        if exist('x.s_unit.templates.data')
                            pool.templates.B{i} = x.s_unit.templates.data;
                        else
                            data = zeros(64,60);
                            %                  obj.templates{id}.best_ch = [];
                            for ch = 1:64
                                temp = [];
                                temp(:,:) = x.s_unit.waveforms{1}(ch,:,:);
                                [Y,Sig,X] = svd(gpuArray(temp),'econ');
                                %                 sig = diag(Sig);%figure; semilogy(sig(sig>1),'kx-')
                                k = 1:3;
                                P = Y(:,k)*Sig(k,k)*X(:,k)';
                                data(ch,:) = mean(gather(P),2).';
                            end
                            pool.templates.B{i}  = data;
                        end
                        
                    end
                    %                 SpikeData = x.s_unit.waveforms{1};
                    %                 SpikeData= SpikeData';
                    %                 SpikeNum = length(SpikeData);
                    %                 % SpikeData = SpikeData - repmat(mean(SpikeData,2),1,SpikeNum);
                    %                 [Y,Sig,X] = svd(SpikeData,'econ');
                    %                 sig = diag(Sig);%figure; semilogy(sig(sig>1),'kx-')
                    %                 k = 1:3;
                    %                 P = Y(:,k)*Sig(k,k)*X(:,k)';
                    %                 P = P';
                    %                 WTempB(i,:) = mean(P,1);
                    
                end
                
                x = load(unit_file_name);
                obj.SU_good.B = x.s_unit.SU_good;
                
                
                pool.peak = obj.list_data.unit_list.data(:,3);
                
                %             B = B+39;
                %             pool.channels(40:end,1) = B;
                %
                
                
                %% Calculate template distance
                
                obj.N = length(obj.SU_good.A);
                
                d_template{seg_ls} = zeros(obj.N);
                
                obj.Na = find(obj.SU_good.A == 1);
                obj.Nb = find(obj.SU_good.B == 1);
                obj.L = length(obj.params.segment_list)-1;
                for i = 1:length(obj.Na)
                    for j =  1:length(obj.Nb)
                        if obj.SU_good.A(obj.Na(i))*obj.SU_good.B(obj.Nb(j)) ~= 0 %if SU on both segments, calculate distance
                            for ch = 1:64
                                temp = sqrt(sum((pool.templates.A{i}(ch,:)- pool.templates.B{j}(ch,:)).^2)); %euclidean distance between waveform template
                                d_template{seg_ls}(obj.Na(i),obj.Nb(j)) = d_template{seg_ls}(obj.Na(i),obj.Nb(j)) + temp;
                            end
                        else
                            
                        end
                    end
                end
                d_template{seg_ls}(find(d_template{seg_ls} == 0)) = NaN;
                fprintf('Time %3.0fs. calculating distance... done! \n', toc);
                
                seq = obj.SU_good.A;
                for s = 1:length(obj.Na)
                    seq(obj.Na(s)) = A(s);
                end
                obj.sequence = [obj.sequence seq];
                
                if seg_ls == length(obj.params.segment_list)-1 % adding last segment of the neurons_list
                    seq2 = obj.SU_good.B;
                    for s = 1:length(obj.Nb)
                        seq2(obj.Nb(s)) = B(s);
                    end
                    obj.sequence = [obj.sequence seq2];
                end
                
%                 obj.sequence = reshape(obj.sequence,7,[] );
            end
            %
            for seg_ls = 1:obj.L
                figure
                imagesc(d_template{seg_ls})
                
            end
            
            %Combine template distance with physical distance. 
            
            
            
            
            pool1 = [];
            pool2 = [];
            for seg_ls = 1:obj.L
                pool1 = [pool1 d_template{seg_ls}];
                pool2 = [pool2; d_template{seg_ls}];
            end
            pool2 = pool2.';
            pool1 = pool1/max(max(pool1));
            pool2 = pool2/max(max(pool2));
            idx =[];
            for i = 1:obj.L
                idx = [idx 1:obj.N];
            end
            cmap = colormap('jet');
            cmap = cmap(1:2:end,:);
            Neuron_good = zeros(1,obj.N);
            Neuron_merge = {};
            for u = 1:obj.N
                figure
                %                 subplot(6,5,u)
                set(gcf, 'Position', [1000 300 1000 800]);
                gscatter(pool1(u,:),pool2(u,:),idx,cmap,'o+*^x',10) %,'doleg','off')
                hold on
                test = sqrt(pool1(u,:).^2 + pool2(u,:).^2);
                grp1 = test(idx == u);
                grp2 = test(idx ~= u);
                %                 grp3 = test(idx == 7);
                idx2 = idx(idx ~=u);
                Y1 = prctile(grp2(find(~isnan(grp2))),5); %5 percentile of group 2
                Y2 = prctile(grp2,10);
                plot([sqrt(Y1^2+Y1^2),0],[0,sqrt(Y1^2+Y1^2)],'--')
                drawnow
                if grp1(find(~isnan(grp1))) < Y1
                    Neuron_good(u) = 1;
                elseif sum(grp1(find(~isnan(grp1))) < Y1)/length(grp1(find(~isnan(grp1)))) >=0.9
                    Neuron_good(u) = 2;
                elseif length(grp1(find(~isnan(grp1)))) <10 %Tenth percentile because sample size is too small.
                    if length(grp1(find(~isnan(grp1)))) - sum(grp1(find(~isnan(grp1))) < Y1)<2
                        Neuron_good(u) = 2;
                    end
                end
                C1 = unique(idx2(find(grp2<Y2)));
                C2= unique(idx2(find(grp2 > Y2)));
                
                %
                Neuron_merge{u} = setdiff(C1,C2);
                title(['unit number ' num2str(u)])
                axis([0 1 0 1])
            end
            
            %             test = sqrt(pool1(u,:).^2 + pool2(u,:).^2);
            %             le = ones(1,length(test));
            %             gscatter(test,le,idx,cmap,'o+*^x',10,'doleg','off');
            
            obj.merge_pair = [];
            for n = 1:obj.N
                if ~isempty(Neuron_merge{n})
                    temp_list = Neuron_merge{n};
                    for t = 1:length(temp_list)
                        if Neuron_merge{temp_list(t)} == n
                            obj.merge_pair = [obj.merge_pair;n,temp_list(t)];
                        end
                    end
                end
            end
            
            %             obj.merge_pair = merge_pair;
            i = 1;
            while i < size(obj.merge_pair,1)
                if ~isempty(intersect(obj.merge_pair(:,2),obj.merge_pair(i,1)))
                    obj.merge_pair(find(obj.merge_pair(:,2) == obj.merge_pair(i,1)),:) = [];
                end
                i = i+1;
            end
            
            obj.SU_keep = ones(obj.N,1);
            obj.SU_clear = Neuron_good;
            
        end
        %%
        
        function obj = clean_cluster(obj,SU_nb)
            list_session = unique(obj.list_data.unit_list.data(:,4));
            for seg_ls = 1:obj.L
                seg_1 = find(strcmp(list_session,obj.params.segment_list{seg_ls}));
                %                 seg_2 = find(strcmp(list_session,obj.params.segment_list{seg_ls+1}));
                
                A = obj.sequence(SU_nb,:);
                %                 B = find(strcmp(obj.list_data.M12E_unit_list.data(:,4),list_session{seg_2}));
                pool.waveforms = [];
                pool.index = [];
                pool.channels = [];
                pool.spike_times = [];
                pool.waveforms_all = [];
                if A(seg_ls) ~=0
                    unit_file_name = [obj.params.animal_name 'u00000'];
                    unit_file_name = [unit_file_name(1:end-size(num2str(A(seg_ls)),2)) num2str(A(seg_ls)) '.mat'];
                    x = load(unit_file_name);
                    if ~isempty(x.s_unit.templates)
%                         x.s_unit.templates.best_ch = 60;
                        pool.channels = [pool.channels; [0 x.s_unit.templates.best_ch]];
                        pool.waveforms = cat(3,pool.waveforms, x.s_unit.waveforms{1}(x.s_unit.templates.best_ch,:,:));
                        pool.spike_times = [pool.spike_times; x.s_unit.spiketimes];
                    end
                    pool.waveforms_all(:,:) = pool.waveforms(1,:,:);
                    [~,score,~] = pca(pool.waveforms_all.');
                    
                    
                    good_split = 0;
                    
                    while good_split == 0
                        
                        if size(pool.waveforms_all(1,:)) > 40
                            gm = fitgmdist(score(:,1:2),2);
                            idx_cls = cluster(gm,score(:,1:2));
                            
                        else
                            idx_cls = kmeans(score(:,1:2),2);
                        end
                        
                        
                        figure
                        subplot(2,2,1)
                        A_id = find(idx_cls ==1);
                        if length(A_id)<100
                            idx = A_id;
                        else
                            idx = randsample(A_id,100);
                        end
                        for i = idx
                            plot(pool.waveforms_all(:,i))
                            hold on
                        end
                        plot(mean(pool.waveforms_all(:,A_id),2),'LineWidth',2,'Color','k');
                        axis([0 60 -300 200])
                        
                        hold off
                        subplot(2,2,3)
%                         B_id = 1:size(find(idx_cls ==2),1);
                        B_id = find(idx_cls ==2);
                        if length(B_id)<100
                            idx = B_id;
                        else
                            idx = randsample(B_id,100);
                        end
                        for i = idx
                            plot(pool.waveforms_all(:,i))
                            hold on
                        end
                        plot(mean(pool.waveforms_all(:,B_id),2),'LineWidth',2,'Color','k');
                        axis([0 60 -300 200])
                        hold off
                        subplot(2,2,2)
                        gscatter(score(:,1),score(:,2),idx_cls)
                        
                        subplot(2,2,4)
                        ISI = pool.spike_times(2:end)-pool.spike_times(1:end-1);
                        
                        edges = -3.5:0.25:2;
                        [NN,~] = histcounts(ISI,10.^edges);
                        refract = sum(NN(1:2))*100/sum(NN);
                        histogram(ISI,10.^edges)
                        set(gca,'xscale','log')
                        xticks([1E-4 1E-3 1E-2 1E-1 1 10]);
                        xticklabels({0.0001, 0.001, 0.01, 0.1, 1, 10})
                        xlabel('time (ms)')
                        title([num2str(refract) '%']);
                        drawnow
                        
                        prompt = 'Is split good? yes(1) no(0) int';
                        good_split = input(prompt);
%                         repmat(sprintf('\b'),1,length(prompt));
                        
                    end
                    
                    
                    prompt = 'split? yes(1 2) cancel(c) int';
                    p = input(prompt);
%                     repmat(sprintf('\b'),1,length(prompt));
                    
                    switch p
                        case 1
                            save_dir = fullfile('D:\Data\Units', filesep, obj.params.animal_name);
                            
                            unitname = [unit_file_name(1:end-4) '_prev.mat'];
                            if exist(unitname)==0
                                
                                s_unit = {};
                                s_unit = x.s_unit;
                                save(fullfile(save_dir,unitname),'s_unit')
                                
                            end
                            s_unit = x.s_unit;
                            
                            %                         s_unit.waveforms = [];
                            s_unit.waveforms{1} = x.s_unit.waveforms{1}(:,:,A_id);
                            s_unit.spiketimes = pool.spike_times(A_id,1);
                            
                            data = zeros(64,60);
                            %                  obj.templates{id}.best_ch = [];
                            for ch = 1:64
                                temp = [];
                                temp(:,:) = s_unit.waveforms{1}(ch,:,:);
                                [Y,Sig,X] = svd(temp,'econ');
                                %                 sig = diag(Sig);%figure; semilogy(sig(sig>1),'kx-')
                                k = 1:3;
                                P = Y(:,k)*Sig(k,k)*X(:,k)';
                                data(ch,:) = mean(P,2).';
                            end
                            
                            s_unit.templates.data = data;
                            
                            save(fullfile(save_dir,unit_file_name),'s_unit')
                            
                            
                        case 2
                            save_dir = fullfile('D:\Data\Units', filesep, obj.params.animal_name);
                            
                            unitname = [unit_file_name(1:end-4) '_prev.mat'];
                            if exist(unitname)==0
                                
                                s_unit = {};
                                s_unit = x.s_unit;
                                save(fullfile(save_dir,unitname),'s_unit')
                                
                            end
                            s_unit = x.s_unit;
                            
                            %                         s_unit.waveforms = [];
                            s_unit.waveforms{1} = x.s_unit.waveforms{1}(:,:,B_id);
                            s_unit.spiketimes = pool.spike_times(B_id,1);
                            
                            data = zeros(64,60);
                            %                  obj.templates{id}.best_ch = [];
                            for ch = 1:64
                                temp = [];
                                temp(:,:) = s_unit.waveforms{1}(ch,:,:);
                                [Y,Sig,X] = svd(temp,'econ');
                                %                 sig = diag(Sig);%figure; semilogy(sig(sig>1),'kx-')
                                k = 1:3;
                                P = Y(:,k)*Sig(k,k)*X(:,k)';
                                data(ch,:) = mean(P,2).';
                            end
                            
                            s_unit.templates.data = data;
                            
                            save(fullfile(save_dir,unit_file_name),'s_unit')
                            
                            
                            
                            
                        case 'c'
                            
                    end
                    

                    
                    
                    
                end
            end
            
            prompt = 'accept cluster?';
            pp = input(prompt);
            switch pp
                case 'y'
                    obj.SU_clear(SU_nb) = 1;
                case 'n'
            end
        end
        
        
        
        %%
        function obj = merge_cluster(obj,pair_nb)
            %             pair_nb = 1;
            list_session = unique(obj.list_data.unit_list.data(:,4));
            
            cluster1 = obj.merge_pair(pair_nb,1);
            cluster2 = obj.merge_pair(pair_nb,2);
            
            
            to_split = zeros(obj.L,1);
            
            for seg_ls = 1:obj.L
                %                 seg_1 = find(strcmp(list_session,obj.params.segment_list{seg_ls}));
                %                 seg_2 = find(strcmp(list_session,obj.params.segment_list{seg_ls+1}));
                
                %                 A = find(strcmp(obj.list_data.M12E_unit_list.data(:,4),list_session{seg_1}));
                A = obj.sequence(:,seg_ls);
                
                %                 B = find(strcmp(obj.list_data.M12E_unit_list.data(:,4),list_session{seg_2}));
                pool.waveforms = [];
                pool.index = [];
                pool.channels = [];
                pool.spike_times = [];
                pool.waveforms_all = [];
                
                if A(cluster1)~= 0 && A(cluster2)~= 0
                    
                    unit_file_name = [obj.params.animal_name 'u00000'];
                    unit_file_name = [unit_file_name(1:end-size(num2str(A(cluster1)),2)) num2str(A(cluster1)) '.mat'];
                    x = load(unit_file_name);
                    if ~isempty(x.s_unit.templates)
                        pool.channels = [pool.channels; [0 x.s_unit.templates.best_ch]];
                        pool.waveforms = cat(3,pool.waveforms, x.s_unit.waveforms{1}(x.s_unit.templates.best_ch,:,:));
                        
                    end
                    
                    pool.index = [pool.index;0*ones(size(x.s_unit.waveforms{1},3),1)];
                    pool.spike_times = [pool.spike_times; x.s_unit.spiketimes];
                    
                    
                    
                    %                 fprintf('Time %3.0fs. decompositing segment 2... \n', toc);
                    unit_file_name = [obj.params.animal_name 'u00000'];
                    unit_file_name = [unit_file_name(1:end-size(num2str(A(cluster2)),2)) num2str(A(cluster2)) '.mat'];
                    y = load(unit_file_name);
                    

                     
                     if ~isempty(y.s_unit.templates)
                         pool.channels = [pool.channels; [1 y.s_unit.templates.best_ch]];
                         
                         
                         try
                             pool.waveforms = cat(3,pool.waveforms,y.s_unit.waveforms{1}(y.s_unit.templates.best_ch,:,:));
                             
                         catch
                             tempmat = y.s_unit.waveforms;
                             y.s_unit.waveforms = {};
                             y.s_unit.waveforms{1} = [];
                            y.s_unit.waveforms{1} = tempmat;
                             pool.waveforms = cat(3,pool.waveforms,y.s_unit.waveforms{1}(y.s_unit.templates.best_ch,:,:));
                        end
                        
%                         pool.waveforms = cat(3,pool.waveforms,y.s_unit.waveforms{1}(y.s_unit.templates.best_ch,:,:));
                        
                    end
                    pool.index = [pool.index;ones(size(y.s_unit.waveforms{1},3),1)];
                    pool.spike_times = [pool.spike_times; y.s_unit.spiketimes];
                    
                    
                    pool.spike_times_all = sort(pool.spike_times);
                    pool.waveforms_all(:,:) = pool.waveforms(1,:,:);
                    
                    
                    figure
                    subplot(2,2,1)
                    
                    A_id = find(pool.index == 0);
                    if length(A_id)<100
                        idx = A_id;
                    else
                        idx = randsample(A_id,100);
                    end
                    for i = idx
                        plot(pool.waveforms_all(:,i))
                        hold on
                    end
                    plot(mean(pool.waveforms_all(:,A_id),2),'LineWidth',2,'Color','k');
                    axis([0 60 -300 200])
                    if ~isempty(x.s_unit.templates)
                        
                        title(['channel  ' num2str(x.s_unit.templates.best_ch)]);
                    end
                    subplot(2,2,3)
                    B_id = find(pool.index == 1);
                    if length(B_id)<100
                        idx = B_id;
                    else
                        idx = randsample(B_id,100);
                    end
                    for i = idx
                        plot(pool.waveforms_all(:,i))
                        hold on
                    end
                    plot(mean(pool.waveforms_all(:,B_id),2),'LineWidth',2,'Color','k');
                    axis([0 60 -300 200])
                    if ~isempty(y.s_unit.templates)
                        
                        title(['channel  ' num2str(y.s_unit.templates.best_ch)]);
                    end
                    
                    % PCA
                    subplot(2,2,2)
                    [~,score,~] = pca(pool.waveforms_all.');
                    gscatter(score(:,1),score(:,2),pool.index)
                    
                    
                    
                    %ISI
                    subplot(2,2,4)
                    
                    ISI = pool.spike_times_all(2:end)-pool.spike_times_all(1:end-1);
                    
                    edges = -3.5:0.25:2;
                    [NN,~] = histcounts(ISI,10.^edges);
                    refract = sum(NN(1:2))*100/sum(NN);
                    histogram(ISI,10.^edges)
                    set(gca,'xscale','log')
                    xticks([1E-4 1E-3 1E-2 1E-1 1 10]);
                    xticklabels({0.0001, 0.001, 0.01, 0.1, 1, 10})
                    xlabel('time (ms)')
                    title([num2str(refract) '%']);
                    drawnow
                    
                    prompt = 'split pair? yes (1) no (0)';
                    
                    to_split(seg_ls) = input(prompt);
                end
                
            end
            
            prompt = 'merge (m),split (s) or cancel (c)';
            p = input(prompt);
            close all
            switch p
                case 'm'
                    save_dir = fullfile('D:\Data\Units', filesep, obj.params.animal_name);
                    
                    
                    for seg_ls = 1:obj.L
                        
                        
                        seg_1 = find(strcmp(list_session,obj.params.segment_list{seg_ls}));
                        %                         seg_2 = find(strcmp(list_session,obj.params.segment_list{seg_ls+1}));
                        
                        A = obj.sequence(:,seg_ls);
                        if A(cluster1)~= 0 && A(cluster2)~= 0
                            
                            %                         B = find(strcmp(obj.list_data.M12E_unit_list.data(:,4),list_session{seg_2}));
                            pool.waveforms = [];
                            pool.index = [];
                            pool.channels = [];
                            pool.spike_times = [];
                            pool.waveforms_all = [];
                            pool.waveforms_all2 = [];
                            
                            %                         pool.waveforms_all = [];
                            
                            
                            unit_file_name = [obj.params.animal_name 'u00000'];
                            unit_file_name = [unit_file_name(1:end-size(num2str(A(cluster1)),2)) num2str(A(cluster1)) '.mat'];
                            x = load(unit_file_name);
                            if ~isempty(x.s_unit.templates)
                                pool.channels = [pool.channels; [0 x.s_unit.templates.best_ch]];
                            end
                            
                            pool.waveforms = cat(3,pool.waveforms, x.s_unit.waveforms{1}(:,:,:));
                            pool.index = [pool.index;0*ones(size(x.s_unit.waveforms{1},3),1)];
                            pool.spike_times = [pool.spike_times; x.s_unit.spiketimes];
                            pool.waveforms_all(:,:) =  pool.waveforms(x.s_unit.templates.best_ch,:,:);
                            
                            unitname = [unit_file_name(1:end-4) '_prev.mat'];
                            if exist(unitname)==0
                                
                                s_unit = {};
                                s_unit = x.s_unit;
                                save(fullfile(save_dir,unitname),'s_unit')
                            end
                            %                 fprintf('Time %3.0fs. decompositing segment 2... \n', toc);
                            unit_file_name2 = [obj.params.animal_name 'u00000'];
                            unit_file_name2 = [unit_file_name2(1:end-size(num2str(A(cluster2)),2)) num2str(A(cluster2)) '.mat'];
                            y = load(unit_file_name2);
                            if ~isempty(y.s_unit.templates)
                                pool.channels = [pool.channels; [1 y.s_unit.templates.best_ch]];
                            end
                            
                            try
                                pool.waveforms = cat(3,pool.waveforms,y.s_unit.waveforms{1}(:,:,:));
                                
                            catch
                                tempmat = y.s_unit.waveforms;
                                y.s_unit.waveforms = {};
                                y.s_unit.waveforms{1} = [];
                                y.s_unit.waveforms{1} = tempmat;
                                pool.waveforms = cat(3,pool.waveforms,y.s_unit.waveforms{1}(:,:,:));
                            end
                            
                            
                            pool.index = [pool.index;ones(size(y.s_unit.waveforms{1},3),1)];
                            pool.spike_times = [pool.spike_times; y.s_unit.spiketimes];
                            pool.waveforms_all2(:,:) = y.s_unit.waveforms{1}(x.s_unit.templates.best_ch,:,:); % changed to look at same channel
                            
                            pool.waveforms_all = [pool.waveforms_all pool.waveforms_all2];
                            
                            pool.spike_times_all = sort(pool.spike_times);
                            %                         pool.waveforms_all(:,:) = pool.waveforms(1,:,:);
                            
                            unitname = [unit_file_name2(1:end-4) '_prev.mat'];
                            if exist(unitname)==0
                                
                                s_unit = {};
                                s_unit = y.s_unit;
                                save(fullfile(save_dir,unitname),'s_unit')
                            end
                            
                            [~,score,~] = pca(pool.waveforms_all.');
                            
                            gm = fitgmdist(score(:,1:2),2);
                            idx_cls = cluster(gm,score(:,1:2));
                            
                            
                            data1 = zeros(64,60);
                            for ch = 1:64
                                temp = [];
                                temp(:,:) = y.s_unit.waveforms{1}(ch,:,:);
                                [Y,Sig,X] = svd(temp,'econ');
                                %                 sig = diag(Sig);%figure; semilogy(sig(sig>1),'kx-')
                                k = 1:3;
                                P = Y(:,k)*Sig(k,k)*X(:,k)';
                                data1(ch,:) = mean(P,2).';
                            end
                            
                            
                            
                            figure
                            subplot(2,2,1)
                            A_id = find(idx_cls == 1);
                            if length(A_id)<100
                                idx = A_id;
                            else
                                idx = randsample(A_id,100);
                            end
                            for i = idx
                                plot(pool.waveforms_all(:,i))
                                hold on
                            end
                            plot(mean(pool.waveforms_all(:,A_id),2),'LineWidth',2,'Color','k');
                            axis([0 60 -300 200])
                            title('cluster 1');
                            
                            subplot(2,2,3)
                            B_id = find(idx_cls == 2);
                            if length(B_id)<100
                                idx = B_id;
                            else
                                idx = randsample(B_id,100);
                            end
                            for i = idx
                                plot(pool.waveforms_all(:,i))
                                hold on
                            end
                            plot(mean(pool.waveforms_all(:,B_id),2),'LineWidth',2,'Color','k');
                            axis([0 60 -300 200])
                            title('cluster 2');
                            
                            
                            % PCA
                            subplot(2,2,2)
                            
                            gscatter(score(:,1),score(:,2),idx_cls)
                            
                            
                            
                            %ISI
                            subplot(2,2,4)
                            
                            ISI = pool.spike_times_all(2:end)-pool.spike_times_all(1:end-1);
                            
                            edges = -3.5:0.25:2;
                            [NN,~] = histcounts(ISI,10.^edges);
                            refract = sum(NN(1:2))*100/sum(NN);
                            histogram(ISI,10.^edges)
                            set(gca,'xscale','log')
                            xticks([1E-4 1E-3 1E-2 1E-1 1 10]);
                            xticklabels({0.0001, 0.001, 0.01, 0.1, 1, 10})
                            xlabel('time (ms)')
                            title([num2str(refract) '%']);
                            drawnow
                            
                            
                            % merging
                            prompt = 'keep cluster1 (c1) cluster2(c2) both(m) or cancel (c)';
                            
                            xx = input(prompt);
                            
                            switch xx
                                case 'c1'
                                    s_unit.waveforms = pool.waveforms(:,:,A_id);
                                    s_unit.spiketimes = pool.spike_times(A_id);
                                    s_unit.xbz_file_name = x.s_unit.xbz_file_name;
                                    s_unit.SU_good = x.s_unit.SU_good;
                                    s_unit.templates = x.s_unit.templates;
                                    save(fullfile(save_dir,unit_file_name),'s_unit')
                                    obj.SU_keep(cluster2) = 0;
                                case 'c2'
                                    s_unit.waveforms = pool.waveforms(:,:,B_id);
                                    s_unit.spiketimes = pool.spike_times(B_id);
                                    s_unit.xbz_file_name = x.s_unit.xbz_file_name;
                                    s_unit.SU_good = x.s_unit.SU_good;
                                    s_unit.templates = x.s_unit.templates;
                                    
                                    save(fullfile(save_dir,unit_file_name),'s_unit')
                                    obj.SU_keep(cluster2) = 0;
                                    
                                case 'm'
                                    s_unit = {};
                                    s_unit.waveforms = pool.waveforms;
                                    s_unit.spiketimes = pool.spike_times_all;
                                    s_unit.xbz_file_name = x.s_unit.xbz_file_name;
                                    s_unit.SU_good = x.s_unit.SU_good;
                                    s_unit.templates = x.s_unit.templates;
                                    
                                    save(fullfile(save_dir,unit_file_name),'s_unit')
                                    
                                    obj.SU_keep(cluster2) = 0;
                                case 'c'
                            end
                            
                        end
                        
                    end
                case 's'
                    
                    save_dir = fullfile('D:\Data\Units', filesep, obj.params.animal_name);
                    
                    
                    for seg_ls = 1:obj.L
                        if to_split(seg_ls) ==1
                            seg_1 = find(strcmp(list_session,obj.params.segment_list{seg_ls}));
                            %                         seg_2 = find(strcmp(list_session,obj.params.segment_list{seg_ls+1}));
                            
                            A = obj.sequence(:,seg_ls);
                            %                         B = find(strcmp(obj.list_data.M12E_unit_list.data(:,4),list_session{seg_2}));
                            pool.waveforms = [];
                            pool.waveforms_all = [];
                            pool.waveforms_all2 = [];
                            pool.index = [];
                            pool.channels = [];
                            pool.spike_times = [];
                            pool.spike_times_all = [];
                            %                         pool.waveforms_all = [];
                            
                            
                            unit_file_name = 'M12Eu00000';
                            unit_file_name = [unit_file_name(1:end-size(num2str(A(cluster1)),2)) num2str(A(cluster1)) '.mat'];
                            x = load(unit_file_name);
                            if ~isempty(x.s_unit.templates)
                                pool.channels = [pool.channels; [0 x.s_unit.templates.best_ch]];
                            end
                            
                            pool.waveforms = cat(3,pool.waveforms, x.s_unit.waveforms{1}(:,:,:));
                            pool.index = [pool.index;0*ones(size(x.s_unit.waveforms{1},3),1)];
                            pool.spike_times = [pool.spike_times; x.s_unit.spiketimes];
                            pool.waveforms_all(:,:) =  pool.waveforms(x.s_unit.templates.best_ch,:,:);
                            
                            
                            unitname = [unit_file_name(1:end-4) '_prev.mat'];
                            if exist(unitname)==0
                                s_unit = {};
                                s_unit = x.s_unit;
                                save(fullfile(save_dir,unitname),'s_unit')
                            end
                            %                 fprintf('Time %3.0fs. decompositing segment 2... \n', toc);
                            unit_file_name2 = 'M12Eu00000';
                            unit_file_name2 = [unit_file_name2(1:end-size(num2str(A(cluster2)),2)) num2str(A(cluster2)) '.mat'];
                            y = load(unit_file_name2);
                            if ~isempty(y.s_unit.templates)
                                pool.channels = [pool.channels; [1 y.s_unit.templates.best_ch]];
                            end
                            pool.waveforms = cat(3,pool.waveforms,y.s_unit.waveforms{1}(:,:,:));
                            pool.index = [pool.index;ones(size(y.s_unit.waveforms{1},3),1)];
                            pool.spike_times = [pool.spike_times; y.s_unit.spiketimes];
                            
                            pool.waveforms_all2(:,:) = y.s_unit.waveforms{1}(x.s_unit.templates.best_ch,:,:);
                            
                            
                            pool.waveforms_all = [pool.waveforms_all pool.waveforms_all2];
                            pool.spike_times_all = sort(pool.spike_times);
                            
                            unitname = [unit_file_name2(1:end-4) '_prev.mat'];
                            if exist(unitname)==0
                                s_unit = {};
                                s_unit = y.s_unit;
                                save(fullfile(save_dir,unitname),'s_unit')
                            end
                            
                            [~,score,~] = pca(pool.waveforms_all.');
                            
                            gm = fitgmdist(score(:,1:2),2);
                            idx_cls = cluster(gm,score(:,1:2));
                            
                            
                            good_fit = 0;
                            
                            figure
                            
                            while good_fit == 0
                                subplot(2,2,1)
                                A_id = find(idx_cls == 1);
                                if length(A_id)<100
                                    idx = A_id;
                                else
                                    idx = randsample(A_id,100);
                                end
                                for i = idx
                                    plot(pool.waveforms_all(:,i))
                                    hold on
                                end
                                plot(mean(pool.waveforms_all(:,A_id),2),'LineWidth',2,'Color','k');
                                axis([0 60 -300 200])
                                hold off
                                title('cluster id = 1')
                                
                                
                                
                                subplot(2,2,3)
                                B_id = find(idx_cls == 2);
                                if length(B_id)<100
                                    idx = B_id;
                                else
                                    idx = randsample(B_id,100);
                                end
                                for i = idx
                                    plot(pool.waveforms_all(:,i))
                                    hold on
                                end
                                plot(mean(pool.waveforms_all(:,B_id),2),'LineWidth',2,'Color','k');
                                axis([0 60 -300 200])
                                hold off
                                title('cluster id = 2')
                                
                                % PCA
                                subplot(2,2,2)
                                
                                gscatter(score(:,1),score(:,2),idx_cls)
                                
                                
                                
                                
                                %ISI
                                subplot(2,2,4)
                                
                                ISI = pool.spike_times_all(2:end)-pool.spike_times_all(1:end-1);
                                
                                edges = -3.5:0.25:2;
                                [NN,~] = histcounts(ISI,10.^edges);
                                refract = sum(NN(1:2))*100/sum(NN);
                                histogram(ISI,10.^edges)
                                set(gca,'xscale','log')
                                xticks([1E-4 1E-3 1E-2 1E-1 1 10]);
                                xticklabels({0.0001, 0.001, 0.01, 0.1, 1, 10})
                                xlabel('time (ms)')
                                title([num2str(refract) '%']);
                                drawnow
                                
                                prompt = 'run fit again? (y) (n) (k)';
                                xx = input(prompt);
                                
                                switch xx
                                    case 'y'
                                        gm = fitgmdist(score(:,1:2),2);
                                        idx_cls = cluster(gm,score(:,1:2));
                                    case 'n'
                                        good_fit = 1;
                                        
                                    case 'k'
                                        idx_cls = kmeans(score(:,1:2),2);
                                end
                            end
                            
                            
                            
                            
                            prompt ='Assign gm_index 1 to cluster1(c1) or cluster2(c2)';
                            xx = input(prompt);
                            switch xx
                                case 'c1'
                                    s_unit = {};
                                    s_unit.waveforms{1} = pool.waveforms(:,:,A_id);
                                    s_unit.spiketimes = pool.spike_times(A_id);
                                    s_unit.xbz_file_name = x.s_unit.xbz_file_name;
                                    s_unit.SU_good = x.s_unit.SU_good;
                                    s_unit.templates.best_ch = x.s_unit.templates.best_ch;
                                    save(fullfile(save_dir,unit_file_name),'s_unit')
                                    
                                    s_unit = {};
                                    s_unit.waveforms{1} = pool.waveforms(:,:,B_id);
                                    s_unit.spiketimes = pool.spike_times(B_id);
                                    s_unit.xbz_file_name = y.s_unit.xbz_file_name;
                                    s_unit.SU_good = y.s_unit.SU_good;
                                    s_unit.templates.best_ch = y.s_unit.templates.best_ch;
                                    save(fullfile(save_dir,unit_file_name2),'s_unit')
                                    
                                    
                                case 'c2'
                                    s_unit = {};
                                    s_unit.waveforms{1} = pool.waveforms(:,:,A_id);
                                    s_unit.spiketimes = pool.spike_times(A_id);
                                    s_unit.xbz_file_name = y.s_unit.xbz_file_name;
                                    s_unit.SU_good = y.s_unit.SU_good;
                                    s_unit.templates.best_ch = y.s_unit.templates.best_ch;
                                    save(fullfile(save_dir,unit_file_name2),'s_unit')
                                    
                                    s_unit = {};
                                    s_unit.waveforms{1} = pool.waveforms(:,:,B_id);
                                    s_unit.spiketimes = pool.spike_times(B_id);
                                    s_unit.xbz_file_name = x.s_unit.xbz_file_name;
                                    s_unit.SU_good = x.s_unit.SU_good;
                                    s_unit.templates.best_ch = x.s_unit.templates.best_ch;
                                    save(fullfile(save_dir,unit_file_name),'s_unit')
                                case 'c'
                                    
                                    
                            end
                        end
                        
                    end
                    
                case 'c'
            end
            
            
           close all 
        end
        
        
        
    end
    
end

















