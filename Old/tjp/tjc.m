classdef tjc < handle
    
    properties(Constant)
    end
    
    properties
        params = eval('parameters_connect');
        list_data = load('M12E_unit_list.mat');
        raw_data = {};
        filt_data = {};
        waveforms = {};
        analysis_code = [];
        analysis_type = {};
        plot_type = {};
        spike_times = {};
        connect = [];
        manual = 0; % variable to indicate whether the clusters were sorted manually
        best_cluster = [];
        connect_add = [];
        connect_master = [];
    end
    
    
    methods(Static)
    end
    
    methods
        
        
        
        function obj = Cconnect4(obj)
            
            
            list_session = unique(obj.list_data.M12E_unit_list.data(:,4));
            
            seg_1 = find(strcmp(list_session,obj.params.segment1));
            seg_2 = find(strcmp(list_session,obj.params.segment2));
            
            A = find(strcmp(obj.list_data.M12E_unit_list.data(:,4),list_session{seg_1}));
            B = find(strcmp(obj.list_data.M12E_unit_list.data(:,4),list_session{seg_2}));
%             centroids = {};
            pool.waveforms = [];
            pool.index = [];
            pool.channels = [];
            pool.templates = {};
            pool.templates.A = {};
            pool.templates.B = {};
            pool.id = [];
            WTempA = zeros(length(A),60);
            tic
            fprintf('Time %3.0fs. decompositing segment 1... \n', toc);
            for i = 1:length(A)
                unit_file_name = 'M12Eu000';
                unit_file_name = [unit_file_name(1:end-size(num2str(A(i)),2)) num2str(A(i)) '.mat'];
                x = load(unit_file_name);
                pool.waveforms = cat(3,pool.waveforms,x.s_unit.waveforms{1});
                pool.index = [pool.index;A(i)*ones(size(x.s_unit.waveforms{1},1),1)];
                pool.id = [pool.id ; x.s_unit.cid];
                
                if ~isempty(x.s_unit.templates)
                    pool.channels = [pool.channels; [A(i) x.s_unit.templates.best_ch]];
                    pool.templates.A{i} = x.s_unit.templates.data;
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
            
            x = load(unit_file_name);
            SU_good.A = x.s_unit.SU_good;
            
            WTempB = zeros(length(B),60);
            
            fprintf('Time %3.0fs. decompositing segment 2... \n', toc);
            for i = 1:length(B)
                unit_file_name = 'M12Eu000';
                unit_file_name = [unit_file_name(1:end-size(num2str(B(i)),2)) num2str(B(i)) '.mat'];
                x = load(unit_file_name);
                pool.waveforms = cat(3,pool.waveforms,x.s_unit.waveforms{1});
                pool.index = [pool.index;B(i)*ones(size(x.s_unit.waveforms{1},1),1)];
                pool.id = [pool.id ; x.s_unit.cid];
                
                if ~isempty(x.s_unit.templates)
                    pool.channels = [pool.channels; [A(i) x.s_unit.templates.best_ch]];
                    pool.templates.B{i} = x.s_unit.templates.data;
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
            SU_good.B = x.s_unit.SU_good;
            
            
            fprintf('Time %3.0fs. decompositing segments ... done! \n', toc);
            pool.peak = obj.list_data.M12E_unit_list.data(:,3);
            
            %             B = B+39;
            %             pool.channels(40:end,1) = B;
            %
            
            
                
                
                
                %% calculate distance
            % channel
%             chanmap = reshape(1:64,[4,16]);
            chanmap = [1 0 3 0; 0 2 0 4];
            for r = 1:15
                chanmap = [chanmap; [1 0 3 0; 0 2 0 4]+[4 0 4 0; 0 4 0 4]*r];
            end

                    
            d_chan = zeros(length(A),length(B));
            d_pca = zeros(length(A),length(B));
            d_peak = zeros(length(A),length(B));
            d_template = zeros(length(A),length(B));
            % PCA
            % PCA
            [~,score,~] = pca(pool.waveforms );
            %
            % mdl = rica(pool.waveforms,3);
            % score = transform(mdl,pool.waveforms);
            
            centroids = {};
            for i = 1:length(A)
                centroids{1}(i,:) =  mean(score(find(pool.index==A(i)),1:3));
            end
            for j = 1:length(B)
                centroids{2}(j,:) =  mean(score(find(pool.index==B(j)),1:3));
            end
            
            
            for i = 1:length(B)
                for j = 1:length(A)
                    [y_A,z_A] = find(chanmap == pool.channels(find(pool.channels(:,1) == A(j)),2));
                    [y_B,z_B] = find(chanmap == pool.channels(find(pool.channels(:,1) == B(i)),2));
                    d_chan(j,i) = sum([y_A-y_B,z_A-z_B].^2)^0.5;
%                     d_pca(j,i) = (sum(centroids{1}(j,:)-centroids{2}(i,:)).^2)^1/3;
                    d_peak(j,i) = abs(pool.peak{find(pool.channels(:,1) == A(j))}- pool.peak{find(pool.channels(:,1) == B(i))});
                    d_template(j,i) = sqrt(sum((WTempB(i,:)-WTempA(j,:)).^2));
                end
            end
            
            % d_chan1 = exp(d_chan/4);
            d_chan1 = 1./(1/3 + exp(-d_chan+3));
            % d_chan1 = d_chan;
            imagesc(-d_chan1)
            d_chan1 = d_chan1/max(max(d_chan1));
            d_pca = d_pca/max(max(d_pca));
            d_peak = d_peak/max(max(d_peak));
            d_template = d_template/max(max(d_template));
            % d_peak = 1./1+exp(-d_chan +3);
            % d_final = (d_chan1).^10.*d_pca.*d_peak;
%             d_final = d_chan1.*d_pca.*d_peak;
            d_final = d_template.*d_chan1;
            figure
            imagesc(-d_final)
            best_clusterB = zeros(length(B),2);
            best_clusterB(:,2) = 1:length(B);
            % colormap(parula)
            best_clusterA = zeros(length(A),2);
            best_clusterA(:,1) = 1:length(A);
            
            for i = 1:length(B)
                best_clusterB(i,1) = find(d_final(:,i) == min(d_final(:,i)));
            end
            
            
            for i = 1:length(A)
                best_clusterA(i,2) = find(d_final(i,:) == min(d_final(i,:)));
                nearest_list = sort(d_final(i,:));
                nearest_A{i} = [];
                for p = 1:4
                nearest_A{i} = [nearest_A{i} find(d_final(i,:) == nearest_list(p))];
                 
                end
            end
            %find best mutual cluster
            
            mutual_cluster = [];
            for i = 1:length(A)
                r = find(best_clusterB(:,1) == best_clusterA(i,1) & best_clusterB(:,2) == best_clusterA(i,2));
                if ~isempty(r)
                    mutual_cluster = [mutual_cluster; best_clusterA(i,1:2)];
                end
            end
            
            
            
            for i = 1:length(mutual_cluster)
                figure
                set(gcf, 'Position', [700 700 1200 500]);
                subplot(1,3,1)
                B_id = find(pool.index== B(mutual_cluster(i,2)));
                A_id = find(pool.index == A(mutual_cluster(i,1)));
                scatter(score(find(pool.index==B(mutual_cluster(i,2))),1),score(find(pool.index==B(mutual_cluster(i,2))),2))
                hold on
                scatter(score(find(pool.index == A(mutual_cluster(i,1))),1),score(find(pool.index == A(mutual_cluster(i,1))),2))
                subplot(1,3,2)
                for p = 1:min(100,length(B_id))
                    plot(pool.waveforms(B_id(p),:))
                    hold on
                end
                axis([0 60 -300 200])
                title(num2str(pool.channels(find(pool.channels(:,1) == B(mutual_cluster(i,2))),2)))
                
                subplot(1,3,3)
                for p = 1:min(100,length(A_id))
                    plot(pool.waveforms(A_id(p),:))
                    hold on
                end
                axis([0 60 -300 200])
                title(num2str(pool.channels(find(pool.channels(:,1) == A(mutual_cluster(i,1))),2)))
                
                drawnow
            end
            
            
            nearest_cluster = [];
            for i = 1:length(A)
                if isempty(find(mutual_cluster(:,1)==i))
                    figure
                    set(gcf, 'Position', [200 500 800 600]);
                    
                    subplot(2,3,1)
                    
                    A_id = find(pool.index == A(i));
                    for p = 1:min(100,length(A_id))
                        plot(pool.waveforms(A_id(p),:))
                        hold on
                    end
                    plot(mean(pool.waveforms(A_id,:),1),'LineWidth',2,'Color','k')
                    title(num2str(pool.channels(find(pool.channels(:,1) == A(i)),2)))
                    axis([0 60 -300 200]) 

                    for pp = 1:4
                        B_id = find(pool.index == B(nearest_A{i}(pp)));
                        subplot(2,3,2+pp)
                        
                        for p = 1:min(100,length(B_id))
                            plot(pool.waveforms(B_id(p),:))
                            hold on
                        end
                        plot(mean(pool.waveforms(B_id,:),1),'LineWidth',2,'Color','k')
                        
                        title([num2str(pp) ' - ' num2str(pool.channels(find(pool.channels(:,1) == B(nearest_A{i}(pp))),2))])
                        axis([0 60 -300 200])
                        hold off
                    end
                    subplot(2,3,2)
                    hold off
                    for pp = 1:4
                        scatter3(score(find(pool.index==B(nearest_A{i}(pp))),1),score(find(pool.index==B(nearest_A{i}(pp))),2),score(find(pool.index==B(nearest_A{i}(pp))),3))
                        hold on
                    end
                    scatter3(score(find(pool.index == A(i)),1),score(find(pool.index == A(i)),2),score(find(pool.index == A(i)),3),'+')
                    
                    
                    prompt = 'choose or cancel';
                    xx = input(prompt);
                    if strcmp(class(xx),'double')
%                         if isempty(intersect(mutual_cluster(:,2),nearest_A{i}(xx)))
                            nearest_cluster = [nearest_cluster; i,nearest_A{i}(xx)];
%                         else
%                             fprintf('double exists');
%                         end
                    
                    end
                    
                end

                    
                % self merge first)
                
            end

            obj.connect = [A(mutual_cluster(:,1)) B(mutual_cluster(:,2))];
            obj.connect = [obj.connect ones(length(obj.connect),1)];
            if ~isempty(nearest_cluster)
                obj.connect_add = [A(nearest_cluster(:,1)) B(nearest_cluster(:,2)),2*ones(size(nearest_cluster(:,1),1),1)];
            end
            %             obj.best_cluster = [A(best_clusterA(:,1)) B(best_clusterA(:,2))];
            obj.connect_master = [obj.connect;obj.connect_add];
            
            
            
            
            
            
        end
        function obj = save_connect(obj)
            idx_mut = find(obj.connect_master(:,3) ==1);
            idx_near = find(obj.connect_master(:,3) ==2);
            save_dir = 'D:\Data\M12E\Neurons';
            if exist(fullfile(save_dir,filesep,'M12E_neuron_list.mat'))
                load(fullfile(save_dir,'M12E_neuron_list.mat'));
                
                
                %                 if  ~isempty(find(strcmp(M12E_neuron_list.data(:,4),obj.params.xbz_file_name)))
                %                     error('Units already saved')
                %                 end
                %                 if ~isempty(M12E_neuron_list)
                %                 n = intersect(intersect(find([M12E_neuron_list(:).Hole] == obj.params.hole_number), ...
                %                     find([M12E_neuron_list(:).Track] == obj.params.track_number)), ...
                %                     find([M12E_neuron_list(:).Depth] == obj.params.depth));
                %                 else
                %                     error('error')
                %                 end
                
                
                % look for existing data
                n = intersect(intersect(find([M12E_neuron_list(:).Hole] == obj.params.hole_number), ...
                    find([M12E_neuron_list(:).Track] == obj.params.track_number)), ...
                    find([M12E_neuron_list(:).Depth] == obj.params.depth));
                
                if isempty(n)
                    
                    n = length(M12E_neuron_list)+1;
                    M12E_neuron_list(n).Hole = obj.params.hole_number;
                    M12E_neuron_list(n).Track = obj.params.track_number;
                    M12E_neuron_list(n).Hemi = obj.params.hemi;
                    M12E_neuron_list(n).Depth = obj.params.depth;
                    M12E_neuron_list(n).list = {};
                    M12E_neuron_list(n).data = [];
                    M12E_neuron_list(n).list{1} = obj.params.segment1;
                    M12E_neuron_list(n).list{2} = obj.params.segment2;
                    M12E_neuron_list(n).data(:,1) = obj.connect(:,1);
                    M12E_neuron_list(n).data(:,2) = obj.connect(:,2);
                else
                    if ~isempty(find(strcmp(M12E_neuron_list.list(:),obj.params.segment1)))
                        m = length(M12E_neuron_list.list);
                        if length(M12E_neuron_list.list) == find(strcmp(M12E_neuron_list.list(:),obj.params.segment1))
                            M12E_neuron_list.list{m+1} = obj.params.segment2;
                            % connecting mutual nearest neighbours
                            
                            for u = 1:length(obj.connect_master(:,1))
                                if ~isempty(find(M12E_neuron_list.data(:,m) == obj.connect_master(u,1)))
                                    u_0 = find(M12E_neuron_list.data(:,m) == obj.connect_master(u,1));
                                    M12E_neuron_list.data(u_0,m+1) = obj.connect_master(u,2);
                                    M12E_neuron_list.data(u_0,m) = obj.connect_master(u,1);
                                else
                                    if  obj.connect_master(u,3) ==1 % if unit doesnt exist from previous segment
                                        u_0 = length(M12E_neuron_list.data)+1;
                                        M12E_neuron_list.data(u_0,m+1) = obj.connect_master(u,2);
                                        M12E_neuron_list.data(u_0,m) = obj.connect_master(u,1);
                                    else
                                        if isempty(intersect(obj.connect_master(idx_mut,2),obj.connect_master(u,2)))
                                            % if unit from segment 2 is not
                                            % already connected to segment
                                            % 1 by mutual nearest neighbor
                                            u_0 = length(M12E_neuron_list.data)+1;
                                            M12E_neuron_list.data(u_0,m+1) = obj.connect_master(u,2);
                                            M12E_neuron_list.data(u_0,m) = obj.connect_master(u,1);
                                        end
                                    end
                                end
                            end
                            
                            
                            
                            
                            
                            
                            
%                             
%                             
%                             
%                             for u = 1:length(obj.connect(:,1))
%                                 if ~isempty(find(M12E_neuron_list.data(:,m) == obj.connect(u,1)))
%                                     u_0 = find(M12E_neuron_list.data(:,m) == obj.connect(u,1));
%                                     M12E_neuron_list.data(u_0,m+1) = obj.connect(u,2);
%                                     M12E_neuron_list.data(u_0,m) = obj.connect(u,1);
%                                 else
%                                     if isempty(intersect(M12E_neuron_list.data(:,2),obj.connect(find(obj.connect(:,2) == obj.connect(u,2)),1)))
%                                         u_0 = length(M12E_neuron_list.data)+1;
%                                         M12E_neuron_list.data(u_0,m+1) = obj.connect(u,2);
%                                         M12E_neuron_list.data(u_0,m) = obj.connect(u,1);
%                                     end
%                                 end
%                             end
%                             for u = 1:length(obj.connect_add)
%                                 %                                 if isempty(intersect(M12E_neuron_list.data(:,2),obj.connect(find(obj.connect(:,2) == obj.connect_add(u,2)),1)))
%                                 if isempty(find(obj.connect(:,2) == obj.connect_add(u,2)))
%                                     if ~isempty(find(M12E_neuron_list.data(:,m) == obj.connect_add(u,1)))
%                                         u_0 = find(M12E_neuron_list.data(:,m) == obj.connect_add(u,1));
%                                         M12E_neuron_list.data(u_0,m+1) = obj.connect_add(u,2);
%                                         M12E_neuron_list.data(u_0,m) = obj.connect_add(u,1);
%                                     else
%                                         
%                                         u_0 = length(M12E_neuron_list.data)+1;
%                                         M12E_neuron_list.data(u_0,m+1) = obj.connect(u,2);
%                                         M12E_neuron_list.data(u_0,m) = obj.connect(u,1);
%                                     end
%                                 end
%                             end
%                             % connecting onesided neighbors
%                             not_connected = intersect(find(M12E_neuron_list.data(:,m) ~= 0), find(M12E_neuron_list.data(:,m+1) == 0));
%                             
%                             for u_0 = 1:length(not_connected)
%                                 M12E_neuron_list.data(not_connected(u_0),m+1) = obj.best_cluster(find(obj.best_cluster(:,1)== M12E_neuron_list.data(not_connected(u_0),m)),2);
%                             end
                        else
                            
                            error('connection already added')
                        end
                    else
                        error('first segment not found')
                    end
                end
                
                
                
                
            else
                M12E_neuron_list = {};
                M12E_neuron_list(1).Hole = obj.params.hole_number;
                M12E_neuron_list(1).Track = obj.params.track_number;
                M12E_neuron_list(1).Hemi = obj.params.hemi;
                M12E_neuron_list(1).Depth = obj.params.depth;
                M12E_neuron_list(1).list = {};
                M12E_neuron_list(1).data = [];
                M12E_neuron_list(1).list{1} = obj.params.segment1;
                M12E_neuron_list(1).list{2} = obj.params.segment2;
                M12E_neuron_list(1).data(:,1) = obj.connect(:,1);
                M12E_neuron_list(1).data(:,2) = obj.connect(:,2);
                
            end
%             pause
            save(fullfile(save_dir,'M12E_neuron_list'),'M12E_neuron_list');
            
            
            
            
            % write data
            
            
        end
        
    end
    
    
end
