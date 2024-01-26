

animal_name= 'RLn602';

expfpath = fullfile('D:\DATA', filesep,'Experiments',filesep,'RL',animal_name,filesep,'Lick');
addpath(expfpath);
list = dir(expfpath);


fs = 3e4;
for n = 3:length(list)
    n = 11
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
