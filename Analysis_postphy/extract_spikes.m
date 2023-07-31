function extract_spikes()


fpath = fullfile('D:\DATA\Reversal Learning\RLn602\2023-07-12_15-27-21');



addpath(genpath('D:\GitHub\Kilosort2')) % path to kilosort folder

addpath(genpath('D:\\GitHub\npy-matlab')) % path to npy-matlab scripts
addpath(genpath('D:\\GitHub\spikes')) % path to npy-matlab scripts


sp = loadKSdir(fpath);

sp.spikeTemplates(1)