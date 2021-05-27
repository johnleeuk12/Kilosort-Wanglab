classdef tjdS < handle
    %{
 
%%%%%%%%%%%%%% EDIT LOG %%%%%%%%%%%%%%%%%%%%%%%%

May 21st: 
    - Modified code to work with inputs from parameters_distanceS file
    - Animal name, PC name and save directory must be specified in the
    parameters_distanceS file before running the code
    -Modified code to display error message when unit has 0 waveforms 
    %}
    
    properties(Constant)
    end
    
    properties
        params = eval('parameters_distanceS');
        list_data = [];
        merge_pair = [];
        SU_keep = [];
        SU_good = [];
        L = [];
        N = [];
        SU_clear = [];
        Na = [];
        Nb = [];
        sequence = [];
        xbz_L = []; %Number of xbz files on parameters file
        Naa=[]; %New variable, same as Na but inside the new for loop of 6 iterations
        waveforms_all = {}; %For storing all waveforms from 6 xbz files
        spike_times_all = {}; %For storing all spike times
        channels = []; %Channels
        unit=[]; %Unit numbers
        index_all={};
        d_temp;
    end
    methods(Static)
    end
    
    methods

        function obj = Cdistance(obj)
            tic
            obj.list_data = load([obj.params.animal_name '_unit_list.mat']);
            list_session = unique([obj.list_data.unit_list.data(:,4)]);
            
            d_template = {};
            for seg_ls = 1:length(obj.params.segment_list)-1
                
                fprintf(['Time %3.0fs. calculating distance between ' obj.params.segment_list{seg_ls} ' and ' obj.params.segment_list{seg_ls+1}  ' ... \n'], toc);
                seg_1 = find(strcmp(list_session,obj.params.segment_list{seg_ls}));
                seg_2 = find(strcmp(list_session,obj.params.segment_list{seg_ls+1}));

                
                %A = find(strcmp(obj.list_data.M12E_unit_list.data(:,4),list_session{seg_1}));
                A = find(strcmp(obj.list_data.unit_list.data(:,4),list_session{seg_1}));
                %B = find(strcmp(obj.list_data.M12E_unit_list.data(:,4),list_session{seg_2}));
                B = find(strcmp(obj.list_data.unit_list.data(:,4),list_session{seg_2}));
               
                
                pool.waveforms = [];
                pool.index = [];
                pool.channels = [];
                pool.templates = {};
                pool.templates.A = {};
                pool.templates.B = {};


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
                    %                     pool.id = [pool.id ; x.s_unit.cid];
                    
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
                                [Y,Sig,X] = svd(temp,'econ');
                                %                 sig = diag(Sig);%figure; semilogy(sig(sig>1),'kx-')
                                k = 1:3;
                                P = Y(:,k)*Sig(k,k)*X(:,k)';
                                data(ch,:) = mean(P,2).';
                            end
                            pool.templates.A{i}  = data;
                        end
                    end

                end
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
                                [Y,Sig,X] = svd(temp,'econ');
                                %                 sig = diag(Sig);%figure; semilogy(sig(sig>1),'kx-')
                                k = 1:3;
                                P = Y(:,k)*Sig(k,k)*X(:,k)';
                                data(ch,:) = mean(P,2).';
                            end
                            pool.templates.B{i}  = data;
                        end
                        
                    end
                    
                end
                
                x = load(unit_file_name);
                obj.SU_good.B = x.s_unit.SU_good;
                
                
                pool.peak = obj.list_data.unit_list.data(:,3);

                
                %% Calculate template distance
                
                obj.N = length(obj.SU_good.A);
                
                d_template{seg_ls} = zeros(obj.N);

                obj.Na = find(obj.SU_good.A == 1);
                obj.Nb = find(obj.SU_good.B == 1);
                obj.L = length(obj.params.segment_list)-1;
                obj.xbz_L = length(obj.params.segment_list);
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
               

            end
            %% NEW CODE SECTION TO OBTAIN SEQUENCE INFO FROM ALL XBZ FILES
            for seg_ls_2 = 1:length(obj.params.segment_list)
                seg = find(strcmp(list_session,obj.params.segment_list{seg_ls_2}));
                
                Aa = find(strcmp(obj.list_data.unit_list.data(:,4),list_session{seg}));
                
                for i = 1:length(Aa)
                    a_unit_file_name = [obj.params.animal_name 'u00000'];
                    a_unit_file_name = [a_unit_file_name(1:end-size(num2str(Aa(i)),2)) num2str(Aa(i)) '.mat'];
                end
                xa = load(a_unit_file_name);
                obj.SU_good.Aa = xa.s_unit.SU_good;
                seq = obj.SU_good.Aa;
                
                obj.Naa = find(obj.SU_good.Aa == 1);

                for s = 1:length(obj.Naa)
                    seq(obj.Naa(s)) = Aa(s);
                end
                obj.sequence = [obj.sequence seq];
            end
            
            %%

            
            obj.d_temp = d_template;
            
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
            Neuron_good = zeros(1,obj.N);
            Neuron_merge = {};
            for u = 1:obj.N

                test = sqrt(pool1(u,:).^2 + pool2(u,:).^2);
                grp1 = test(idx == u);
                grp2 = test(idx ~= u);

                idx2 = idx(idx ~=u);
                Y1 = prctile(grp2(find(~isnan(grp2))),5); %5 percentile of group 2
                Y2 = prctile(grp2,10);

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
                
                Neuron_merge{u} = setdiff(C1,C2);
            end
            
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

            i = 1;
            while i < size(obj.merge_pair,1)
                if ~isempty(intersect(obj.merge_pair(:,2),obj.merge_pair(i,1)))
                    obj.merge_pair(find(obj.merge_pair(:,2) == obj.merge_pair(i,1)),:) = [];
                end
                i = i+1;
            end
            
            obj.SU_keep = ones(obj.N,1);
            obj.SU_clear = Neuron_good;
            
            disp('Opening GUI');
            GUI_TJD(obj.params.segment_list, obj.d_temp, obj.N, obj.sequence, obj.params.animal_name, obj.params.save_dir);
            
            
        end
        %%
        
        function [obj, out] = clean_cluster(obj,best_ch,SU_nb,sequence,current)
            
            seg_ls = current;
            
            A = sequence(SU_nb,:);
            
            pool.waveforms = [];
            pool.index = [];
            pool.channels = [];
            pool.spike_times = [];
            pool.waveforms_all = [];
            pool.waveforms_chan_all = [];
            
            if A(seg_ls) ~=0
                unit_file_name = [obj.params.animal_name 'u00000'];
                unit_file_name = [unit_file_name(1:end-size(num2str(A(seg_ls)),2)) num2str(A(seg_ls)) '.mat'];
                x = load(unit_file_name);

                if x.s_unit.waveforms{1} == 0
                    disp(['ERROR: File ' unit_file_name ' waveforms are equal to 0']);
                else
                    
                    if ~isempty(x.s_unit.templates)
                        pool.waveforms = cat(3,pool.waveforms, x.s_unit.waveforms{1}(best_ch,:,:));
                        pool.spike_times = [pool.spike_times; x.s_unit.spiketimes];
                    end
                    pool.waveforms_chan_all = cat(3,pool.waveforms_chan_all, x.s_unit.waveforms{1}(:,:,:));
                    pool.waveforms_all(:,:) = pool.waveforms(1,:,:);
                    pool.index = [pool.index;ones(size(x.s_unit.waveforms{1},3),1)];
                    
                    
                    %OUTPUTS
                    
                    out.waveforms_all = pool.waveforms_all(:,:);
                    
                    out.spike_times = pool.spike_times;
                    
                    out.index_all = pool.index;
                    
                    out.waveforms_chan_all = pool.waveforms_chan_all;
                    
                    disp('Out');
                end
            end
            
        end
    


        %%
        %The output is waveform information from 2 clusters
        function [obj, out] = merge_cluster(obj,best_ch,pair_nb,sequence,current)
            disp(best_ch);
            seg_ls=current;
            cluster1 = pair_nb(1);
            cluster2 = pair_nb(2);
            
            A=sequence(:,seg_ls);
            
            pool.waveforms = [];
            pool.index = [];
            pool.channels = [];
            pool.spike_times = [];
            pool.waveforms_all = [];
            pool.waveforms_all2 = [];
            
            
            if A(cluster1)~= 0 && A(cluster2)~= 0
                unit_file_name1 = [obj.params.animal_name 'u00000'];
                unit_file_name1 = [unit_file_name1(1:end-size(num2str(A(cluster1)),2)) num2str(A(cluster1)) '.mat'];
                x = load(unit_file_name1);
                if x.s_unit.waveforms{1} == 0
                    disp(['ERROR: File ' unit_file_name1 ' waveforms are equal to 0']);
                else
                    
                    if ~isempty(x.s_unit.templates)
                        
                        pool.channels = [pool.channels; [0 x.s_unit.templates.best_ch]];
                        
                    end
                    pool.waveforms = cat(3,pool.waveforms, x.s_unit.waveforms{1}(:,:,:));
                    pool.index = [pool.index;0*ones(size(x.s_unit.waveforms{1},3),1)];
                    pool.spike_times = [pool.spike_times; x.s_unit.spiketimes];
                    pool.waveforms_all(:,:) =  pool.waveforms(best_ch,:,:);
                end
                
                unit_file_name2 = [obj.params.animal_name 'u00000'];
                unit_file_name2 = [unit_file_name2(1:end-size(num2str(A(cluster2)),2)) num2str(A(cluster2)) '.mat'];
                y = load(unit_file_name2);
                if y.s_unit.waveforms{1} == 0
                    disp(['ERROR: File ' unit_file_name2 ' waveforms are equal to 0']);
                else

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
                    pool.waveforms_all2(:,:) = y.s_unit.waveforms{1}(best_ch,:,:);
                end
                
                pool.waveforms_all = [pool.waveforms_all pool.waveforms_all2];

                pool.spike_times_all = sort(pool.spike_times);
                
                
                %OUTPUTS
                
                out.waveforms_all = pool.waveforms_all;
               
                out.spike_times = pool.spike_times;
                
                out.index_all = pool.index;
                
                out.waveforms_chan_all = pool.waveforms;
                
                disp('Out');
                
                
            end

        end
        
        
        
    end
    
end

















