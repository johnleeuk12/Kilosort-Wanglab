
ops = {};
ops.bin = 5;
ops.edges = -2e3:ops.bin:4e3;

[raster2, rate, Lick2, Lick_rate] = gather_raster_ephys(2, 6, Pool);
[C2_1, spont_1, L_1] = ana_lick_aligned(Pool,1,ops);
[C2_2, spont_2, L_2] = ana_lick_aligned(Pool,2,ops);

% figure
% plot(mean(C2_1.hit,2))
% hold on 
% plot(mean(C2_2.hit,2))
% 
% figure
% plot(mean(C2_1.FA,2))
% hold on 
% plot(mean(C2_2.FA,2))

list_pre2_r1 =ana_lick_aligned_list(C2_1,spont_1,Pool,ops);
list_pre2_r2 =ana_lick_aligned_list(C2_2,spont_2,Pool,ops);

list3 = (list_pre2_r1(:,4)) | (list_pre2_r2(:,4));
x_axis = ops.edges(1:end-1)*1e-3;
figure
plot(x_axis,mean(C2_1.hit(:,list_pre2_r1(:,1)),2))
hold on 
plot(x_axis,mean(C2_2.hit(:,list_pre2_r2(:,1)),2))


figure
plot(mean(C2_1.FA(:,list_pre2_r1(:,2)),2))
hold on 
plot(mean(C2_2.FA(:,list_pre2_r2(:,2)),2))

%% normalize R1 and R2

C3_1 = normC2(C2_1,Pool,spont_1);
C3_2 = normC2(C2_2,Pool,spont_2);

figure
shadedErrorBar(x_axis,smoothdata(mean(C3_1.hit(:,list_pre2_r1(:,1)),2),'gaussian',2), ...
    std(C3_1.hit(:,list_pre2_r1(:,1)),[],2)/sqrt(size(C3_1.hit(:,list_pre2_r1(:,1)),2)),'lineProps','b')
% plot(x_axis,mean(C3_1.hit(:,list_pre2_r1(:,1)),2))
hold on 
shadedErrorBar(x_axis,smoothdata(mean(C3_2.hit(:,list_pre2_r2(:,1)),2),'gaussian',2), ...
    std(C3_2.hit(:,list_pre2_r2(:,1)),[],2)/sqrt(size(C3_2.hit(:,list_pre2_r2(:,1)),2)),'lineProps','r')

% figure
% L_1.hit = smoothdata(L_1.hit,"gaussian",10); 
% L_2.hit = smoothdata(L_2.hit,"gaussian",10); 

figure
shadedErrorBar(x_axis,smoothdata(mean(L_1.hit(:,list_pre2_r1(:,1)),2),'gaussian',2), ...
    std(L_1.hit(:,list_pre2_r1(:,1)),[],2)/sqrt(size(L_1.hit(:,list_pre2_r1(:,1)),2)),'lineProps','b')
% plot(x_axis,mean(C3_1.hit(:,list_pre2_r1(:,1)),2))
hold on 
shadedErrorBar(x_axis,smoothdata(mean(L_2.hit(:,list_pre2_r2(:,1)),2),'gaussian',2), ...
    std(L_2.hit(:,list_pre2_r2(:,1)),[],2)/sqrt(size(L_2.hit(:,list_pre2_r2(:,1)),2)),'lineProps','r')


% plot(x_axis,mean(L_1.hit(:,list_pre2_r1(:,1)),2))
% hold on
% plot(x_axis,mean(L_2.hit(:,list_pre2_r2(:,1)),2))
ylim([0,20])
%%
ops = {};
ops.bin = 5;
ops.edges = -2e3:ops.bin:4e3;

[C2, spont, L] = ana_lick_aligned(Pool,0, ops);
% list_pre2 =ana_lick_aligned_list(C2,spont,Pool,ops);
[list_pre2, list_post] =ana_lick_aligned_list(C2,spont,Pool,ops);
figure
plot(mean(C2.hit(:,list_pre2(:,3)),2))
hold on
plot(mean(C2.FA(:,list_pre2(:,3)),2))


% ops.edges = -2e3:ops.bin:4e3;
% [C2, spont, L] = ana_lick_aligned(Pool,0, ops);
%% 

% list_pre2 = list_post;
% list_pre2(:,4) = (list_pre2(:,4)+ list_post(:,4))>0;
% list_pre2(:,1) = list_pre2(:,4);
% list_pre2(:,2) = list_pre2(:,4);
% list_pre2 = list_post;
% Change bin size
C3 = {};
L3 = {};
sz2=  length(C2.hit)/10;


C3.hit = zeros(sz2, sum(list_pre2(:,1)));
L3.hit = zeros(sz2, sum(list_pre2(:,1)));
C3.FA = zeros(sz2, sum(list_pre2(:,2)));
L3.FA = zeros(sz2, sum(list_pre2(:,2)));



for t = 1:sz2
    C3.hit(t,:) = mean(C2.hit((t-1)*10+1:t*10,list_pre2(:,1)),1);
    C3.FA(t,:) = mean(C2.FA((t-1)*10+1:t*10,list_pre2(:,2)),1);
    L3.hit(t,:) = mean(L.hit((t-1)*10+1:t*10,list_pre2(:,1)),1);
    L3.FA(t,:) = mean(L.FA((t-1)*10+1:t*10,list_pre2(:,2)),1);
end
    
coef = {}; 
coef.hit= zeros(1,sum(list_pre2(:,1)));
coef.FA = zeros(1,sum(list_pre2(:,2)));
for n = 1:sum(list_pre2(:,1))
    ans = corrcoef(C3.hit(:,n),L3.hit(:,n));
    coef.hit(n) = ans(2,1);
end
for n = 1:sum(list_pre2(:,2))
    ans = corrcoef(C3.FA(:,n),L3.FA(:,n));
    coef.FA(n) = ans(2,1);
end


% figure
% scatter(zeros(1,length(coef.hit)),coef.hit)
% hold on
% scatter(ones(1,length(coef.FA)),coef.FA)




% figure
% plot(mean(C3.hit,2));
% hold on
% plot(mean(C3.FA,2));

% figure
% imagesc(C3.hit.')
% figure
% imagesc(C2.hit(:,list_pre2(:,1)))


C4 = {};
C4.hit = zeros(size(C3.hit));
C4.FA = zeros((size(C3.FA)));
for n = 1:size(C3.hit,2)
    C4.hit(:,n) = (C3.hit(:,n)-spont.mean(n,1))/(max(abs(C3.hit(:,n)-spont.mean(n,1)))+1);
end

for n = 1:size(C3.FA,2)
    C4.FA(:,n) = (C3.FA(:,n)-spont.mean(n,1))/(max(abs(C3.FA(:,n)-spont.mean(n,1)))+1);
end

list_hit = (mean(C4.hit,1)>0);
list_FA = (mean(C4.FA,1)>0);


figure
imagesc(C4.hit(:,list_hit).')
caxis([-0.4,0.9])

figure
imagesc(C4.FA(:,list_FA).')
caxis([-0.4,0.9])


x_axis = -2:5*1e-2:4;
x_axis = x_axis(1:end-1);
figure
% plot(mean(C4.hit,2))
shadedErrorBar(x_axis,smoothdata(mean(C4.hit(:,list_hit),2),'gaussian',5),std(C4.hit(:,list_hit),[],2)/sqrt(size(C4.hit(:,list_hit),2)),'lineProps','b')

hold on
shadedErrorBar(x_axis,smoothdata(mean(C4.FA(:,list_FA),2),'gaussian',5),std(C4.FA(:,list_FA),[],2)/sqrt(size(C4.FA(:,list_FA),2)),'lineProps','r')

%%

onset = {};
onset.hit = mean(C4.hit(35:39,list_hit),1);
onset.FA = mean(C4.FA(35:39,list_FA),1);
sus = {};
sus.hit = mean(C4.hit(50:60,list_hit),1);
sus.FA = mean(C4.FA(50:60,list_FA),1);

figure(5)
x = [onset.hit,onset.FA];
x = x.';
g = [ones(length(onset.hit),1);2*ones(length(onset.FA),1)];
boxplot(x,g)
hold on 
scatter(g,x,20,'k','filled');
xlim([0,3]);
ylim([0,1]);
hold off


figure(6)
x = [sus.hit,sus.FA];
x = x.';
g = [ones(length(sus.hit),1);2*ones(length(sus.FA),1)];
boxplot(x,g)
hold on 
scatter(g,x,20,'k','filled');
xlim([0,3]);
ylim([0,1]);
hold off


%% 

list_H = {};
list_H.IC = list_hit;
list_F = {};
list_F.IC = list_FA;
C4_IC = C4;

list_H.AC = list_hit;
list_F.AC = list_FA;
C4_AC = C4;


%%

figure
shadedErrorBar(x_axis,smoothdata(mean(C4_AC.hit(:,list_H.AC),2),'gaussian',2),std(C4_AC.hit(:,list_H.AC),[],2)/sqrt(size(C4_AC.hit(:,list_H.AC),2)),'lineProps','b')
hold on
shadedErrorBar(x_axis,smoothdata(mean(C4_IC.hit(:,list_H.IC),2),'gaussian',2),std(C4_IC.hit(:,list_H.IC),[],2)/sqrt(size(C4_IC.hit(:,list_H.IC),2)),'lineProps','r')
xlim([-0.5,0.5])



%%
[C2, spont, L] = ana_lick_aligned(Pool,0);
list_pre2 =ana_lick_aligned_list(C2,spont,Pool);
C3 = {};

edges = -2e3:5:2e3;
edges = edges(1:end-1);

C2.hit = smoothdata(C2.hit,"gaussian",20); 
C2.FA = smoothdata(C2.FA,"gaussian",20);

C3.hit = zeros(size(C2.hit,1),size(C2.hit,2));
C3.FA = zeros(size(C2.hit,1),size(C2.hit,2));


for n = 1:length(Pool)
    C3.hit(:,n) = (C2.hit(:,n)-spont.mean(n,1))/(max(abs(C2.hit(:,n)-spont.mean(n,1)))+0.5);
    C3.FA(:,n) = (C2.FA(:,n)-spont.mean(n,1))/(max(abs(C2.FA(:,n)-spont.mean(n,1)))+0.5);
end

figure
plot(edges,mean(C3.hit,2))
hold on
plot(edges,mean(C3.FA,2))

% coef = {}; 
% coef.hit= zeros(1,length(Pool));
% coef.FA = zeros(1,length(Pool));
% for n = 1:length(Pool)
%     ans = corrcoef(C3.hit(:,n),smoothdata(L.hit(:,n),"gaussian",20));
%     coef.hit(n) = ans(2,1);
% end
% for n = 1:length(Pool)
%     ans = corrcoef(C3.FA(:,n),smoothdata(L.FA(:,n),"gaussian",20));
%     coef.FA(n) = ans(2,1);
% end
% 
% figure
% scatter(zeros(1,sum(list_pre2(:,4))),abs(coef.hit(list_pre2(:,4))))
% hold on
% scatter(ones(1,sum(list_pre2(:,4))),abs(coef.FA(list_pre2(:,4))))





imagesc(C2.all.')
% % caxis([0,1])

% %%

% [coeff,score,lat] = pca(((C3.hit(:,list_pre2(:,4))+C3.FA(:,list_pre2(:,4)))/2).');
% figure
% plot(coeff(1,:))
% 
% 
% traj_hit = C3.hit(:,list_pre2(:,4)) * score;
% traj_FA = C3.FA(:,list_pre2(:,4)) * score;
% 
% plot3(traj_hit(:,1),traj_hit(:,2),traj_hit(:,3));
% hold on
% plot3(traj_FA(:,1),traj_FA(:,2),traj_FA(:,3));
% 
% figure
% t1 = tiledlayout(1,5);
% for tr = 1:5
%     nexttile
%     plot(traj_hit(:,tr))
%     hold on
%     plot(traj_FA(:,tr))
% end
% 

%% figure 
figure
t = tiledlayout(2,2);
nexttile
imagesc(edges,1:sum(list_pre2(:,4)),C3.hit(:,list_pre2(:,4)).')
caxis([0,1])
xlim([-1000,1500])
nexttile
imagesc(edges,1:sum(list_pre2(:,4)),C3.FA(:,list_pre2(:,4)).')
caxis([0,1])
xlim([-1000,1500])

nexttile
imagesc(edges,1:sum(list_pre2(:,4)),L.hit(:,list_pre2(:,4)).')
caxis([0,50])
xlim([-1000,1500])

nexttile
imagesc(edges,1:sum(list_pre2(:,4)),L.FA(:,list_pre2(:,4)).')
caxis([0,50])
xlim([-1000,1500])

t.TileSpacing = 'compact';
t.Padding = 'compact';


figure
plot(edges,mean(C3.hit(:,list_pre2(:,4)),2))
hold on
plot(edges,mean(C3.FA(:,list_pre2(:,4)),2))





