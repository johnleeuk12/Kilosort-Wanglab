%{
EDIT LOG
2024/01/26 JHL

Within the Experiments/animal_name folder, create a folder named Lick
this folder should originally contain all the lick.mat files (raw analog
signal for lick detection) from the post-phy analysis code.
change files to YYYY-MM-DD_animal_name_region.mat

%} 

animal_name= 'RLn609';

expfpath = fullfile('D:\DATA', filesep,'Experiments',filesep,'RL',animal_name,filesep,'Lick');
addpath(expfpath);
list = dir(expfpath);


fs = 3e4;
for n = 3:length(list)
    raw_L = load(list(n).name);
    L = downsample(raw_L.lick,fs/1e3);
    gL = gradient(double(L));
    pks = (gL > .5*1e4);
    locs = 1:length(L);
    locs = locs(pks);
    ISI = locs(2:end)-locs(1:end-1);
    ind = (ISI>100);
    locs2 = locs(2:end);
    locs = [locs(1), locs2(ISI>100)];
    
    save(fullfile(expfpath,filesep,[list(n).name(1:10),'_licks.mat']),'locs')
end
