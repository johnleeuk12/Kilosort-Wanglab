



%% load data
addpath(genpath(fullfile('D:\Data\M12E\Units')));

list_data = load('M12E_unit_list.mat');
% x = load('M12Eu001.mat');
% [~,score,~] = pca(x.s_unit.waveforms{1} );
% scatter3(x.s_unit.spiketimes,score(:,1),score(:,2));


list_session = unique(list_data.M12E_unit_list.data(:,4));
A = find(strcmp(list_data.M12E_unit_list.data(:,4),list_session{2}));
B = find(strcmp(list_data.M12E_unit_list.data(:,4),list_session{3}));

centroids = {};
pool.waveforms = []; 
pool.index = [];
pool.channels = [];
WTempA = zeros(length(A),60);
for i = 1:length(A)
    unit_file_name = 'M12Eu000';
    unit_file_name = [unit_file_name(1:end-size(num2str(A(i)),2)) num2str(A(i)) '.mat'];
    x = load(unit_file_name);
    pool.waveforms = [pool.waveforms ;x.s_unit.waveforms{1}];
    pool.index = [pool.index;A(i)*ones(size(x.s_unit.waveforms{1},1),1)];
    pool.channels = [pool.channels; [A(i) x.s_unit.ch]];
    
    SpikeData = x.s_unit.waveforms{1};
    SpikeData= SpikeData';
    SpikeNum = length(SpikeData);
    % SpikeData = SpikeData - repmat(mean(SpikeData,2),1,SpikeNum);
    [Y,Sig,X] = svd(SpikeData);
    sig = diag(Sig);%figure; semilogy(sig(sig>1),'kx-')
    k = 1:3;
    P = Y(:,k)*Sig(k,k)*X(:,k)';
    P = P';
    WTempA(i,:) = mean(P,1);
end

WTempB = zeros(length(B),60);

for i = 1:length(B)
    unit_file_name = 'M12Eu000';
    unit_file_name = [unit_file_name(1:end-size(num2str(B(i)),2)) num2str(B(i)) '.mat'];
    x = load(unit_file_name);
    pool.waveforms = [pool.waveforms ;x.s_unit.waveforms{1}];
    pool.index = [pool.index;B(i)*ones(size(x.s_unit.waveforms{1},1),1)];
    pool.channels = [pool.channels; [B(i) x.s_unit.ch]];
    
        SpikeData = x.s_unit.waveforms{1};
    SpikeData= SpikeData';
    SpikeNum = length(SpikeData);
    % SpikeData = SpikeData - repmat(mean(SpikeData,2),1,SpikeNum);
    [Y,Sig,X] = svd(SpikeData);
    sig = diag(Sig);%figure; semilogy(sig(sig>1),'kx-')
    k = 1:3;
    P = Y(:,k)*Sig(k,k)*X(:,k)';
    P = P';
    WTempB(i,:) = mean(P,1);
    
end

pool.peak = list_data.M12E_unit_list.data(:,3);



%% calculate distance
% channel
chanmap = reshape(1:64,[4,16]);
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
        d_pca(j,i) = (sum(centroids{1}(j,:)-centroids{2}(i,:))^2)^1/3;
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
d_final = d_chan1.*d_pca.*d_peak;
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



%% calculate d_PCA





imagesc(d_pca)



channel = {};
for i = 1:length(A)
    channel{1}(i,:) =  mean(score(find(pool.index==A(i)),1:5));
end
for j = 1:length(B)
    channel{2}(j,:) =  mean(score(find(pool.index==B(j)),1:5));
end




for i = 1:length(B)
    unit_file_name = 'M12Eu000';
    unit_file_name = [unit_file_name(1:end-size(num2str(B(i)),2)) num2str(B(i)) '.mat'];
    x = load(unit_file_name);
    [~,score,~] = pca(x.s_unit.waveforms{1} );
    centroids{2}(i,:) =  mean(score(:,1:5));
end











%% testing svd deocm and reconstruction

waves = x.s_unit.waveforms{1};
waves = waves';
SpikeData = [SpikeData waves];

SpikeNum = length(SpikeData);

SpikeData = SpikeData - repmat(mean(SpikeData,2),1,SpikeNum); 
[Y,Sig,X] = svd(SpikeData); 
sig = diag(Sig); figure; semilogy(sig(sig>1),'kx-')      % plot the significant singular values xlabel('index','fontsize',14); ylabel('singular value','fontsize',14)l 





