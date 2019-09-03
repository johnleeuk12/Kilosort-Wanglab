save_dir = 'D:\Data\M12E\Neurons';
load(fullfile(save_dir,'M12E_neuron_list.mat'));

%% re matching duplicates

% ChanMap = reshape(1:64,[4,16]).';
Data_mat = M12E_neuron_list.data;
ChanMap = [1 0 3 0; 0 2 0 4];
for r = 1:15
    ChanMap = [ChanMap; [1 0 3 0; 0 2 0 4]+[4 0 4 0; 0 4 0 4]*r];
end


for s = 1:11
    fprintf(['segment ' num2str(s) '\n'])
    
    for i  =1:size(M12E_neuron_list.data,1)
        A = [];
        B = [];
        %find duplications in segment s+1;
        dupli = find(Data_mat(:,s+1) == Data_mat(i,s+1));
        if length(dupli) >1
            B = Data_mat(i,s+1);
            A = Data_mat(dupli,s);
            if isempty(intersect([A; B],0))
                %duplicate on segment s+1 exists
                pool.waveforms = [];
                pool.index = [];
                pool.channels = [];
                for f = 1:length(A)
                    
                    unit_file_name = 'M12Eu000';
                    unit_file_name = [unit_file_name(1:end-size(num2str(A(f)),2)) num2str(A(f)) '.mat'];
                    x = load(unit_file_name);
                    pool.waveforms = [pool.waveforms ;x.s_unit.waveforms{1}];
                    pool.index = [pool.index;f*ones(size(x.s_unit.waveforms{1},1),1)];
                    pool.channels = [pool.channels; [f x.s_unit.ch]];
                    
                end
                unit_file_name = 'M12Eu000';
                unit_file_name = [unit_file_name(1:end-size(num2str(B),2)) num2str(B) '.mat'];
                x = load(unit_file_name);
                pool.waveforms = [pool.waveforms ;x.s_unit.waveforms{1}];
                pool.index = [pool.index;f+1*ones(size(x.s_unit.waveforms{1},1),1)];
                pool.channels = [pool.channels; [f+1 x.s_unit.ch]];
                
                
                [~,score,~] = pca(pool.waveforms );
                figure
                set(gcf, 'Position', [200 500 800 600]);
                
                for f = 1:length(dupli)+1
                    subplot(2,3,f)
                    idx = find(pool.index == f);
                    nb_waves = min(length(idx),100);
                    idx_plot = randsample(idx,nb_waves);
                    for p = 1:length(idx_plot)
                        plot(pool.waveforms(idx_plot(p),:))
                        hold on
                    end
                    plot(mean(pool.waveforms(idx,:),1),'LineWidth',2,'Color','k')
                    hold off
                    axis([0 60 -300 200])
                    title(['ch = ' num2str(pool.channels(f,2)) ' id = ' num2str(f)])
                end
                subplot(2,3,6)
                gscatter(score(:,1),score(:,2),pool.index)
                prompt = 'remove 1, 2, both';
                xx = input(prompt);
                
                switch xx
                    case 1
                        Data_mat(dupli(xx),s+1:end) = 0;
                        
                    case 2
                        Data_mat(dupli(xx),s+1:end) = 0;
                        
                    case both
                end
            end
        end
        
    end
end

Data_mat_copy = Data_mat;
Data_mat = Data_mat_copy;
%% connecting disconnected segments

for s = 1:10
    fprintf(['segment ' num2str(s) '\n'])
    
    for i  =1:size(M12E_neuron_list.data,1)
        if Data_mat(i,s) ~= 0
            A = [];
            B = [];
            if Data_mat(i,s+1) == 0 && Data_mat(i,s) ~= 0
                seg2_idx = find(Data_mat(:,s) == 0 & Data_mat(:,s+1) ~= 0);
                
                B = Data_mat(i,s);
                A = Data_mat(seg2_idx,s+1);
                pool.waveforms = [];
                pool.index = [];
                pool.channels = [];
                unit_file_name = 'M12Eu000';
                unit_file_name = [unit_file_name(1:end-size(num2str(B),2)) num2str(B) '.mat'];
                x = load(unit_file_name);
                pool.waveforms = [pool.waveforms ;x.s_unit.waveforms{1}];
                pool.index = [pool.index;0*ones(size(x.s_unit.waveforms{1},1),1)];
                pool.channels = [pool.channels; [0 x.s_unit.ch]];
                
                SpikeData = x.s_unit.waveforms{1};
                SpikeData= SpikeData';
                SpikeNum = length(SpikeData);
                % SpikeData = SpikeData - repmat(mean(SpikeData,2),1,SpikeNum);
                [Y,Sig,X] = svd(SpikeData,'econ');
                sig = diag(Sig);%figure; semilogy(sig(sig>1),'kx-')
                k = 1:3;
                P = Y(:,k)*Sig(k,k)*X(:,k)';
                P = P';
                WTempB = mean(P,1);
                for f = 1:length(A)
                    
                    unit_file_name = 'M12Eu000';
                    unit_file_name = [unit_file_name(1:end-size(num2str(A(f)),2)) num2str(A(f)) '.mat'];
                    x = load(unit_file_name);
                    pool.waveforms = [pool.waveforms ;x.s_unit.waveforms{1}];
                    pool.index = [pool.index;f*ones(size(x.s_unit.waveforms{1},1),1)];
                    pool.channels = [pool.channels; [f x.s_unit.ch]];
                    
                    SpikeData = x.s_unit.waveforms{1};
                    SpikeData= SpikeData';
                    SpikeNum = length(SpikeData);
                    % SpikeData = SpikeData - repmat(mean(SpikeData,2),1,SpikeNum);
                    [Y,Sig,X] = svd(SpikeData,'econ');
                    sig = diag(Sig);%figure; semilogy(sig(sig>1),'kx-')
                    k = 1:3;
                    P = Y(:,k)*Sig(k,k)*X(:,k)';
                    P = P';
                    WTempA(f,:) = mean(P,1);
                end
                
                [~,score,~] = pca(pool.waveforms );
                
                
                % calculate distance metric
                d_chan = [];
                d_template = [];
                for f = 1:length(A)
                    [y_A,z_A] = find(ChanMap == pool.channels(find(pool.channels(:,1) == f),2));
                    [y_B,z_B] = find(ChanMap == pool.channels(find(pool.channels(:,1) == 0),2));
                    d_chan(f) = sum([y_A-y_B,z_A-z_B].^2)^0.5;
                    d_template(f) = sqrt(sum((WTempB(1,:)-WTempA(f,:)).^2));
                end
                d_chan1 = 1./(1/3 + exp(-d_chan+3));
                d_chan1 = d_chan1/max(max(d_chan1));
                d_template = d_template/max(max(d_template));
                d_final = d_template.*d_chan1;
                
                min_d_final = sort(d_final);
                % plot 3 with smallest distance
                near_list = [];
                for pp = 1:min(3,length(min_d_final))
                    near_list = [near_list find(d_final == min_d_final(pp))]; 
                end
                figure
                set(gcf, 'Position', [200 500 800 600]);
                
                subplot(2,3,1)
                idx =  find(pool.index == 0);
                nb_waves = min(length(idx),100);
                idx_plot = randsample(idx,nb_waves);
                for p = 1:length(idx_plot)
                    plot(pool.waveforms(idx_plot(p),:))
                    hold on
                end
                plot(mean(pool.waveforms(idx,:),1),'LineWidth',2,'Color','k')
                hold off
                axis([0 60 -300 200])
                title(['ch = ' num2str(pool.channels(1,2))])
                    subplot(2,3,2)
                    scatter(score(idx,1),score(idx,2))
                hold on
                
                for f = 1:min(3,length(min_d_final))
                    subplot(2,3,f+3)
                    idx = find(pool.index ==near_list(f));
                    
                    nb_waves = min(length(idx),100);
                    idx_plot = randsample(idx,nb_waves);
                    for p = 1:length(idx_plot)
                        plot(pool.waveforms(idx_plot(p),:))
                        hold on
                    end
                    plot(mean(pool.waveforms(idx,:),1),'LineWidth',2,'Color','k')
                    hold off
                    axis([0 60 -300 200])
                    title(['ch = ' num2str(pool.channels(find(pool.channels(:,1) == near_list(f)),2)) ' d = ' num2str(d_final(near_list(f)))])
                    
                    subplot(2,3,2)
                    scatter(score(idx,1),score(idx,2))
                    hold on
                end
                
                
                
                prompt = 'choose or cancel';
                xx = input(prompt);
                if strcmp(class(xx),'double')
                    r = find(Data_mat(:,s+1) == A(xx));
                    Data_mat(i,s+1:end) = Data_mat(r,s+1:end);
                    Data_mat(r,s+1:end) = 0;
                end
                
            end
        end
        
    end
end
    
        
        
    
    
    
 



% 
% 

%% Plotting Single neuron groups

% for u = 1:71

u = 39;
pool.waveforms = [];
pool.index = [];
pool.channels = [];
% A = Data_mat_copy(u,:);
% A = M12E_neuron_list.data(u,:);
A = [0,20,39,69,87,98,122,157,169,205,217,244];
%
list = A;
for i = 1:length(A)
    
    unit_file_name = 'M12Eu000';
    if A(i) ==0
        pool.waveforms = [pool.waveforms ;zeros(1,60)];
        pool.index = [pool.index;i*ones(1,1)];
        pool.channels = [pool.channels; [i 0]];
    else
        unit_file_name = [unit_file_name(1:end-size(num2str(A(i)),2)) num2str(A(i)) '.mat'];
        x = load(unit_file_name);
        pool.waveforms = [pool.waveforms ;x.s_unit.waveforms{1}];
        pool.index = [pool.index;i*ones(size(x.s_unit.waveforms{1},1),1)];
        pool.channels = [pool.channels; [i x.s_unit.ch]];
    end
end

% % WTempA = zeros(length(A),60);
% [Y,Sig,X] = svd(pool.waveforms,'econ');
% sig = diag(Sig);%figure; semilogy(sig(sig>1),'kx-')
% k = 1:3;
% P = Y(:,k)*Sig(k,k)*X(:,k)';
% P = P';
% % WTempA(i,:) = mean(P,1);

[~,score,~] = pca(pool.waveforms );

figure
c = colormap('jet');
color_list = 1:floor(64/length(list)):64; 
color_list = color_list(1:length(list));
set(gcf, 'Position', [100 700 2200 500]);

for i = 1:length(list)
    subplot(2,length(list),i)
    idx = find(pool.index == i);
    nb_waves = min(length(idx),100);
    idx_plot = randsample(idx,nb_waves);
    for p = 1:length(idx_plot)
        plot(pool.waveforms(idx_plot(p),:),'Color',c(color_list(i),:))
        hold on
    end
    plot(mean(pool.waveforms(idx,:),1),'LineWidth',2,'Color','k')
    axis([0 60 -300 200])
    title(['ch = ' num2str(pool.channels(i,2))])
    
    subplot(2,length(list),[length(list)+1 length(list)*2])
    scatter(idx,score(idx,1),[],c(color_list(i),:))
    hold on
    drawnow
end
% pause
% end


% 
% 
% 













 
            
            
            
            
            
            
            