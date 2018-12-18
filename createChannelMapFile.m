%  create a channel map file

function make_HarrisChannelMap(fpath)
% function make_HarrisChannelMap()
% create a channel Map file for simulated data (eMouse)

Nchannels = 64;





% here I know a priori what order my channels are in.  So I just manually 
% make a list of channel indices (and give
% an index to dead channels too). chanMap(1) is the row in the raw binary
% file for the first channel. chanMap(1:2) = [33 34] in my case, which happen to
% be dead channels. 


% chanMap = [26 36 41 43 35 38 39 30 33 25 27 20 24 34 40 45 ...
%     37 21 23 29 28 22 17 50 48 47 18 32 46 19 42 44 54 56 ...
%     13 52 2 16 49 31 1 15 3 51 14 4 53 12 11 64 7 60 62 9 58 ...
%     5 57 55 6 59 10 63 61 8];

% chanMap = [23 7 14 5 13 9 12 19 15 22 24 29 25 10 16 3 11 26 ...
%     28 20 21 32 27 64 2 31 1 17 4 8 30 6 60 36 58 62 47 63 33 18 ...
%     48 46 34 61 35 59 45 37 38 42 50 54 52 56 40 44 55 43 57 53 ...
%     39 51 39 31];

chanMap = [23 8 14 6 14 10 11 19 16 22 24 29 25 9 15 4 12 26 ...
    28 20 21 32 27 63 1 31 2 17 3 7 30 5 59 36 57 61 ...
    47 64 33 18 48 46 34 62 35 60 45 37 38 42 49 53 ...
    51 55 40 44 56 43 58 54 39 52 50 41]; 

% chanMap = [21 18 25 3 23 20 26 5 24 27 32 7 30 19 28 9 ...
%     29 17 22 11 16 13 1 15 31 12 14 10 8 4 6 2 ...
%     63 59 61 57 55 51 53 34 50 54 52 56 58 62 60 64 ...
%     43 33 37 39 40 46 48 38 45 49 47 36 35 42 41 44];   
%     
% chanMap = [28 31 24 46 26 29 23 44 25 22 17 42 19 20 21 40 ...
%     20 32 27 28 33 36 48 34 18 37 35 39 41 45 43 47 ...
%     49 53 51 55 57 61 59 16 64 60 62 58 56 52 54 50 ...
%     5 15 11 9 10 4 2 12 3 63 1 14 13 8 7 6];

% the first thing Kilosort does is reorder the data with data = data(chanMap, :).
% Now we declare which channels are "connected" in this normal ordering, 
% meaning not dead or used for non-ephys data

connected = true(Nchannels, 1);

% now we define the horizontal (x) and vertical (y) coordinates of these
% 34 channels. For dead or nonephys channels the values won't matter. Again
% I will take this information from the specifications of the probe. These
% are in um here, but the absolute scaling doesn't really matter in the
% algorithm. 

xcoords   = repmat([1 2 3 4]', 1, Nchannels/4);
xcoords   = 20*xcoords(:);
ycoords   = repmat(1:Nchannels/4, 4, 1);
% ycoords(2,:) = ycoords(2,:)+1;
% ycoords(4,:) = ycoords(4,:)+1;
ycoords   = 20*ycoords;
ycoords(2,:) = ycoords(2,:)+10;
ycoords(4,:) = ycoords(4,:)+10;
ycoords = ycoords(:);
ycoords = abs(ycoords-330)+10;

% Often, multi-shank probes or tetrodes will be organized into groups of
% channels that cannot possibly share spikes with the rest of the probe. This helps
% the algorithm discard noisy templates shared across groups. In
% this case, we set kcoords to indicate which group the channel belongs to.
% In our case all channels are on the same shank in a single group so we
% assign them all to group 1. 

kcoords   = ones(Nchannels,1);

% grouping of channels (i.e. tetrode groups)
% at this point in Kilosort we do data = data(connected, :), ycoords =
% ycoords(connected), xcoords = xcoords(connected) and kcoords =
% kcoords(connected) and no more channel map information is needed (in particular
% no "adjacency graphs" like in KlustaKwik). 
% Now we can save our channel map for the eMouse. 

% would be good to also save the sampling frequency here
%channelmap, from right to left, top to bottom
fs = 30000; % sampling frequency


save(fullfile(fpath, 'chanMap.mat'), ... 
    'chanMap','connected', 'xcoords', 'ycoords', 'kcoords', 'fs')

% save('D:\DATA\Spikes\HarrisProbe\chanMap.mat', ...
%     'chanMap','connected', 'xcoords', 'ycoords', 'kcoords', 'fs')

% kcoords is used to forcefully restrict templates to channels in the same
% channel group. An option can be set in the master_file to allow a fraction 
% of all templates to span more channel groups, so that they can capture shared 
% noise across all channels. This option is

% ops.criterionNoiseChannels = 0.2; 

% if this number is less than 1, it will be treated as a fraction of the total number of clusters

% if this number is larger than 1, it will be treated as the "effective
% number" of channel groups at which to set the threshold. So if a template
% occupies more than this many channel groups, it will not be restricted to
% a single channel group. 