function [list_pre2,list_post] = ana_lick_aligned_list2(C2, spont, Pool,ops)
% 
% bin = ops.bin;
% edges = ops.edges;


C3 = {};
C3.hit = smoothdata(C2.hit,1,'gaussian',10);
C3.FA = smoothdata(C2.FA,1,'gaussian',10);
C3.hit = C2.hit;
C3.FA;


Av = {};
Av.pre= zeros(length(Pool),2);
Av.post = zeros(length(Pool),2);
onset = abs(ops.edges(1)/ops.bin);
win = floor(100/ops.bin);



for n = 1:length(Pool)
    Av.pre(n,1) = mean(C3.hit(onset-win:onset,n));
    Av.pre(n,2) = mean(C3.FA(onset-win:onset,n));
    Av.post(n,1) = mean(C3.hit(onset+1:onset+win,n));
    Av.post(n,2) = mean(C3.FA(onset+1:onset+win,n));    
end

list_pre2 = zeros(4,length(Pool));
for n= 1:length(Pool)
    if Av.pre(n,1) > 2*spont.std(n,1)+spont.mean(n,1) % Hit
        list_pre2(1,n) = 1;
    end
    if Av.pre(n,2) > 2*spont.std(n,1)+spont.mean(n,1)
        list_pre2(2,n) = 1;
    end
end

list_pre2(3,:) = list_pre2(1,:).*list_pre2(2,:);
list_pre2(4,:) =  (list_pre2(1,:)>0) + (list_pre2(2,:)>0);
list_pre2 = list_pre2.' ;
list_pre2 = logical(list_pre2);


list_post = zeros(4,length(Pool));
for n= 1:length(Pool)
    if Av.post(n,1) > 1*spont.std(n,1)+spont.mean(n,1) % Hit
        list_post(1,n) = 1;
    end
    if Av.post(n,2) > 1*spont.std(n,1)+spont.mean(n,1)
        list_post(2,n) = 1;
    end
end


list_post(3,:) = list_post(1,:).*list_post(2,:);
list_post(4,:) =  (list_post(1,:)>0) + (list_post(2,:)>0);
list_post = list_post.' ;
list_post = logical(list_post);



