clear ops
ops.GPU                 = 1; % whether to run this code on an Nvidia GPU (much faster, mexGPUall first)		
ops.parfor              = 1; % whether to use parfor to accelerate some parts of the algorithm		
ops.verbose             = 1; % whether to print command line progress		
ops.showfigures         = 0; % whether to plot figures during optimization		
		
ops.datatype            = 'tev';  % binary ('dat', 'bin') or 'openEphys'		
% ops.fbinary             = 'test_binary.dat'; % will be created for 'openEphys'	
% ops.ephys_type          = file_type; % OpenEphys file type. check header number of one "Continuous" file
ops.fproc               = fullfile(fpath, 'temp_wh.dat'); % residual from RAM of preprocessed data		
ops.root                = fullfile(fpath); % 'openEphys' only: where raw files are		
		
% ops.chanMap             = fullfile('D:\DATA\Spikes\HarrisProbe\', 'chanMap.mat'); % make this file using createChannelMapFile.m		

ops.fs                  = 24414;        % sampling rate		(omit if already in chanMap file)
ops.NchanTOT            = 32;           % total number of channels (omit if already in chanMap file)
ops.Nchan               = 32;           % number of active channels (omit if already in chanMap file)
ops.Nfilt               = 128;           % number of clusters to use (2-4 times more than Nchan, should be a multiple of 32)     		
ops.nNeighPC            = 12; % visualization only (Phy): number of channnels to mask the PCs, leave empty to skip (12)		
ops.nNeigh              = 16; % visualization only (Phy): number of neighboring templates to retain projections of (16)		
		
% options for channel whitening		
ops.whitening           = 'full'; % type of whitening (default 'full', for 'noSpikes' set options for spike detection below)		
ops.nSkipCov            = 1; % compute whitening matrix from every N-th batch (1)		
ops.whiteningRange      = 32; % how many channels to whiten together (Inf for whole probe whitening, should be fine if Nchan<=32)		
		
% define the channel map as a filename (string) or simply an array		
ops.chanMap             = fullfile(fpath, 'chanMap.mat'); % make this file using createChannelMapFile.m		
ops.criterionNoiseChannels = 0.2; % fraction of "noise" templates allowed to span all channel groups (see createChannelMapFile for more info). 		
% ops.chanMap = 1:ops.Nchan; % treated as linear probe if a chanMap file		
		
% other options for controlling the model and optimization		
ops.Nrank               = 3;    % matrix rank of spike template model (3)		
ops.nfullpasses         = 6;    % number of complete passes through data during optimization (6)		
ops.maxFR               = 10000;  % maximum number of spikes to extract per batch (20000)		
ops.fshigh              = 300;   % frequency for high pass filtering		
ops.fslow               = 3000;   % frequency for low pass filtering (optional)
ops.ntbuff              = 64;    % samples of symmetrical buffer for whitening and spike detection		
ops.scaleproc           = 200;   % int16 scaling of whitened data		
ops.NT                  = 32*1024+ ops.ntbuff;% this is the batch size (try decreasing if out of memory) 		
% for GPU should be multiple of 32 + ntbuff		
		
% the following options can improve/deteriorate results. 		
% when multiple values are provided for an option, the first two are beginning and ending anneal values, 		
% the third is the value used in the final pass. 		
% ops.Th               = [4 10 10];    % threshold for detecting spikes on template-filtered data ([6 12 12])		
% ops.lam              = [5 10 10];   % large means amplitudes are forced around the mean ([10 30 30])		
ops.nannealpasses    = 4;            % should be less than nfullpasses (4)		
% ops.momentum         = 1./[20 1000];  % start with high momentum and anneal (1./[20 1000])		
ops.shuffle_clusters = 1;            % allow merges and splits during optimization (1)		
ops.mergeT           = .1;           % upper threshold for merging (.1)		
ops.splitT           = .1;           % lower threshold for splitting (.1)		
		
% options for initializing spikes from data		
ops.initialize      = 'no'; %'fromData' or 'no'		
ops.spkTh           = - 6;      % spike threshold in standard deviations (4)		
ops.loc_range       = [3  1];  % ranges to detect peaks; plus/minus in time and channel ([3 1])		
ops.long_range      = [30  6]; % ranges to detect isolated peaks ([30 6])		
ops.maskMaxChannels = 5;       % how many channels to mask up/down ([5])		
ops.crit            = .65;     % upper criterion for discarding spike repeates (0.65)		
ops.nFiltMax        = 10000;   % maximum "unique" spikes to consider (10000)		
		
% load predefined principal components (visualization only (Phy): used for features)		
dd                  = [];%load(fullfile(fpath,'PCspikes_3.mat')); % you might want to recompute this from your own data		
ops.wPCA            = [];%dd.Wi(:,1:7);   % PCs 		
		
% options for posthoc merges (under construction)		
ops.fracse  = 0.1; % binning step along discriminant axis for posthoc merges (in units of sd)		
ops.epu     = Inf;		
		
ops.ForceMaxRAMforDat   = 20e9; % maximum RAM the algorithm will try to use; on Windows it will autodetect.


%% New config options for Kilosort2

% minimum firing rate on a "good" channel (0 to skip)
ops.minfr_goodchannels = 0.1; 

ops.Th =  [15 10]; %[25 10];
% how important is the amplitude penalty (like in Kilosort1, 0 means not used, 10 is average, 50 is a lot) 
ops.lam = 15;  

% splitting a cluster at the end requires at least this much isolation for each sub-cluster (max = 1)
ops.AUCsplit = 0.9; 

% minimum spike rate (Hz), if a cluster falls below this for too long it gets removed
ops.minFR = 1/10; 

% number of samples to average over (annealed from first to second value) 
ops.momentum = [20 400]; 

% spatial constant in um for computing residual variance of spike
ops.sigmaMask = 10; 

% threshold crossings for pre-clustering (in PCA projection space)
ops.ThPre = 7; 

%% danger, changing these settings can lead to fatal errors
% options for determining PCs
ops.spkTh           = -6;      % spike threshold in standard deviations (-6)
ops.reorder         = 1;       % whether to reorder batches for drift correction. 
ops.nskip           = 25;  % how many batches to skip for determining spike PCs

ops.GPU                 = 1; % has to be 1, no CPU version yet, sorry
% ops.Nfilt               = 1024; % max number of clusters
ops.nfilt_factor        = 4; % max number of clusters per good channel (even temporary ones)
ops.ntbuff              = 64;    % samples of symmetrical buffer for whitening and spike detection
ops.NT                  = 32*1024+ ops.ntbuff; % must be multiple of 32 + ntbuff. This is the batch size (try decreasing if out of memory). 
ops.whiteningRange      = 32; % number of channels to use for whitening each channel
ops.nSkipCov            = 25; % compute whitening matrix from every N-th batch
ops.scaleproc           = 200;   % int16 scaling of whitened data
ops.nPCs                = 3; % how many PCs to project the spikes into
ops.useRAM              = 0; % not yet available

