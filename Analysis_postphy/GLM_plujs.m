

% Get data from 

load('D:\Python\LeeLab\saved_struct.mat');
W = {};
list1 = list0 +1;
W{1} = CD0(list1,:); % R1
W{2} = CD1(list1,:); % RT
W{3} = CD2(list1,:); % R2
wbin = 21; % from stim onset + 10, then average response of 4 seconds (over the response period + 1sec;

A = mean(W{1}(:,wbin:wbin+40),2);
B = mean(W{2}(:,wbin:wbin+40),2);


[c,S] = polyfit(A,B,1);
[y_fit,delta] = polyval(c,A,S);
R2 = 1-(S.normr/norm(A)-mean(A))^2;

figure
hold on
scatter(A,B)
xline([0,0.1,-0.1])
yline([0,0.1,-0.1])
plot(A,y_fit,'r--','LineWidth',1)
text(0.5,-0.5,['R^2 = ',num2str(R2,'%.2f')])
axis([-1,1,-1,1])

hold off
