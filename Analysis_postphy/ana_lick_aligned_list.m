function [list_pre2,list_post] = ana_lick_aligned_list(C2, spont, Pool,ops)

%%

% Calculate number of neurons with less than 1spk per sec

% bin into 20ms
bin = ops.bin;
edges = ops.edges;

% Normalizing C
C2norm = {};
C2norm.hit = zeros(length(edges)-1,length(Pool)); 
C2norm.FA = zeros(length(edges)-1,length(Pool)); 
% C2norm = zeros(length(edges),length(Pool));

w = gausswin(5); %gaussian smoothing
C3 = {};
% C3.hit = C2.hit;
% C3.FA = C2.FA;
C3.hit = imfilter(C2.hit,w,'replicate');
% C3.hit = smoothdata(C2.hit,1,'gaussian',10);

C3.FA = imfilter(C2.FA,w,'replicate');
% C3.FA = smoothdata(C2.FA,1,'gaussian',10);
% C3.hit = imgaussfilt(C2.hit,[2,1]);
% C3.FA = imgaussfilt(C2.FA,[2,1]);
% C3.hit = C2.hit;
% C3.FA = C2.FA;


max_ind = {};
max_ind.hit = [];
max_ind.FA = [];
thresh = zeros(1,length(Pool));
% m_ind = {};
for n =1:length(Pool)
    [~,m1] = max(C3.hit(:,n));
    [~,m2] = max(C3.FA(:,n));
%     try
%         maxC2 = max([mean(C3.hit(m1-1:m1+1,n)),mean(C3.hit(m1-1:m1+1,n))]);
%     catch
        maxC2 = max([mean(C3.hit(m1,n)),mean(C3.hit(m1,n))]);
%     end
    %     maxC2 = max([mean(C3.hit(m1-1:m1+1,n)),mean(C3.hit(m1-1:m1+1,n))]);
    %     C2norm(:,n) = (C3(:,n)-prctile(C3(:,n),25))/(maxC2-prctile(C3(:,n),25)+2);
    C2norm.hit(:,n) = (C3.hit(:,n)-spont.mean(n,1))/(maxC2-spont.mean(n,1)+5);
    %     C2norm.hit(:,n) = (C3.hit(:,n)-spont.hit(n,1))/(maxC2+5);
    
    C2norm.FA(:,n) = (C3.FA(:,n)-spont.mean(n,1))/(maxC2-spont.mean(n,1)+5);
    %     C2norm.FA(:,n) = (C3.FA(:,n)-spont.FA(n,1))/(maxC2+5);
    thresh(1,n) = spont.std(n,1)/(maxC2-spont.mean(n,1)+5);

    try
        if mean(C3.hit(m1-1:m1+1,n))> 2*spont.std(n,1)+spont.mean(n,1)
            max_ind.hit = [max_ind.hit, n];
        end
        if max(C3.FA(m2-1:m2+1,n))> 2*spont.std(n,1)+spont.mean(n,1)
            max_ind.FA = [max_ind.FA, n];
        end
    catch
    end
end




max_ind.all = intersect(max_ind.hit,max_ind.FA);


Av = {};
max_ind.all2 = unique([max_ind.hit, max_ind.FA]);
Av.pre= zeros(length(Pool),2);
Av.post = zeros(length(Pool),2);
% win = floor(200/ops.bin); %(100ms before)
% onset = floor((length(ops.edges)-1)/2);

onset = abs(ops.edges(1)/ops.bin);
win = floor(200/ops.bin);


for n = 1:length(max_ind.all2)
    Av.pre(max_ind.all2(n),1) = mean(C2norm.hit(onset-win:onset,max_ind.all2(n)));
    Av.pre(max_ind.all2(n),2) = mean(C2norm.FA(onset-win:onset,max_ind.all2(n)));
    Av.post(max_ind.all2(n),1) = mean(C2norm.hit(onset+1:onset+win,max_ind.all2(n)));
    Av.post(max_ind.all2(n),2) = mean(C2norm.FA(onset+1:onset+win,max_ind.all2(n)));    
end

list_pre2 = zeros(4,length(Pool));
for n= 1:length(Pool)
    if Av.pre(n,1) > 2*thresh(1,n) % Hit
        list_pre2(1,n) = 1;
    end
    if Av.pre(n,2) > 2*thresh(1,n) % FA
        list_pre2(2,n) = 1;
    end
end

list_pre2(3,:) = list_pre2(1,:).*list_pre2(2,:);
list_pre2(4,:) =  (list_pre2(1,:)>0) + (list_pre2(2,:)>0);
list_pre2 = list_pre2.' ;
list_pre2 = logical(list_pre2);

list_post = zeros(4,length(Pool));
for n= 1:length(Pool)
    if Av.post(n,1) > 2*thresh(1,n) % Hit
        list_post(1,n) = 1;
    end
    if Av.post(n,2) > 2*thresh(1,n) % FA
        list_post(2,n) = 1;
    end
end

list_post(3,:) = list_post(1,:).*list_post(2,:);
list_post(4,:) =  (list_post(1,:)>0) + (list_post(2,:)>0);
list_post = list_post.' ;
list_post = logical(list_post);




