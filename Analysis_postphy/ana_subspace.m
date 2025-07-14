SS = {};
SS.mean = zeros(1,4);
SS.std = zeros(1,4);
%% 

O = load('D:\DATA\subspace_OFC.mat');
ind = 4;



S = {};
S.mean = zeros(1,6);
S.mean(1,4) = mean(O.hit1);
S.mean(1,5) = mean(O.hit2);
S.mean(1,6) = mean(O.hit3);
S.mean(1,1) = mean(O.lick1);
S.mean(1,2) = mean(O.lick2);
S.mean(1,3) = mean(O.lick3);

S.std = zeros(1,6);
S.std(1,4) = std(O.hit1);
S.std(1,5) = std(O.hit2);
S.std(1,6) = std(O.hit3);
S.std(1,1) = std(O.lick1);
S.std(1,2) = std(O.lick2);
S.std(1,3) = std(O.lick3);

% figure
% bar([1:6],S.mean,"facecolor","none")
% hold on
% errorbar(1:6,S.mean,S.std/sqrt(20),"LineStyle","none","LineWidth",2,"Color",'k')
% ylim([0.3,1])



SS.mean(1,ind) = mean(S.mean(1,1:3));
SS.std(1,ind) = mean(S.std(1,1:3));


%% 

figure
bar([1:4],SS.mean,"facecolor","none")
hold on
errorbar(1:4,SS.mean,SS.std,"LineStyle","none","LineWidth",2,"Color",'k')
ylim([0.3,1])
