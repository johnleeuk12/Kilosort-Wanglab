
function Pool = pool_data_ephys(animal_name,varargin)



% animal_name = 'RLn602';
addpath(fullfile('D:\DATA', filesep,'Units',filesep,'RL',animal_name)); % path to Units


expfpath = fullfile('D:\DATA', filesep,'Experiments',filesep,'RL',animal_name);
addpath(fullfile('D:\DATA', filesep,'Experiments',filesep,'RL',animal_name)); % path to Units

load([animal_name, '_unit_list.mat']);




% X = event_time_presentation(1,expfpath, fname);


%% Pool data 


% varargin{1} is area name : 'PPC', 'A1', 'IC'
% varargin{1} = 'PPC';
u_list = unit_list.data(find(strcmp(unit_list.data(:,7),varargin{1})),1);

Pool = {};
tic
for i = 1:length(u_list)
    if mod(i,10) ==1
        fprintf(['%4d /' num2str(length(u_list)) ' time : %6.2f sec \n'],i,toc')
    end
    unit_file_name = [animal_name 'u00000'];
    unit_file_name = [unit_file_name(1:end-size(num2str(u_list{i}),2)) num2str(u_list{i}) '.mat'];
    x = load(unit_file_name);
    fname = unit_list.data(u_list{i},5);
    y = event_time_presentation(1,expfpath, fname);
%     l = 
    for st =  1: length(y.event_code)
        ind = find(y.trial_type(:,2) == y.event_code(st));
        y.trial_type(ind,3) = 1:length(ind);
        y.trial_type(ind,4) = st;
    end
    y.trial_tag{3} = 'rep';
    

    Pool(i).neuron_nb = u_list(i);
    Pool(i).best_ch = x.s_unit.best_ch;
    Pool(i).waveforms(:,:) = x.s_unit.waveforms.waveForms(1,:,Pool(i).best_ch,:);
    %extract spiketimes and other information
    Pool(i).spiketimes = x.s_unit.spiketimes;
    Pool(i).licktimes = y.licks.locs*1e-3;
    Pool(i).licktimes = Pool(i).licktimes.';
    Pool(i).xb = y;
    Pool(i).eventtimes = x.s_unit.event_time;
    




end











































