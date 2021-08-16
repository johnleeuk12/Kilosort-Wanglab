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



% % old batch of 64ch
% chanMap = [27 32 21 3 25 30 19 5 23 17 24 7 20 29 26 9 ...
%     22 31 28 11 16 13 1 15 18 12 14 10 8 4 6 2 ...
%     64 60 62 58 56 52 54 48 49 53 51 55 57 61 59 63 ...
%     38 42 40 45 43 35 33 47 36 50 34 44 46 39 41 37];


% new batch 2021
chanMap = [27 32 21 3 25 30 19 5 23 17 24 7 20 29 26 9 ...
    22 31 28 11 16 13 1 15 18 12 14 10 8 4 6 2 ...
    64 59 61 57 55 51 53 49 48 54 52 56 58 62 60 64 ...
    37 41 39 46 44 34 50 36 47 33 35 43 45 40 42 38];

% the first thing Kilosort does is reorder the data with data = data(chanMap, :).
% Now we declare which channels are "connected" in this normal ordering, 
% meaning not dead or used for non-ephys data

connected = true(Nchannels, 1);
% connected(find(chanMap == 21)) = false;
% connected(find(chanMap == 63)) = false;
% connected(find(chanMap == 15)) = false;

% connected(64) = false;

% connected(1:2) = false;
% connected(17,7,11,15,38,55) = false;
% connected(17,1) = false;
% connected([1:4:end],1) = false;
% now we define the horizontal (x) and vertical (y) coordinates of these
% 34 channels. For dead or nonephys channels the values won't matter. Again
% I will take this information from the specifications of the probe. These
% are in um here, but the absolute scaling doesn't really matter in the
% algorithm. 

xcoords   = repmat([1 2 3 4]', 1, Nchannels/4);
xcoords   = 20*xcoords(:);
ycoords   = repmat(1:Nchannels/4, 4, 1);
%  ycoords(2,:) = ycoords(2,:)+1;
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
% kcoords(1:4:end,1) =NaN;



% comment: what you see in openEphys is already channel arranged. ch17 is
% kcoords 17.
% dead = find(chanMap == 40); %change here for bad channels

% kcoords(51) =0;

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