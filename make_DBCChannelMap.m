%  create a channel map file

function make_DBCChannelMap(fpath)
% function make_HarrisChannelMap()
% create a channel Map file for simulated data (eMouse)

Nchannels = 64;





% here I know a priori what order my channels are in.  So I just manually 
% make a list of channel indices (and give
% an index to dead channels too). chanMap(1) is the row in the raw binary
% file for the first channel. chanMap(1:2) = [33 34] in my case, which happen to
% be dead channels. 
% 
% chanMap = [ 3 17 13 31 5 1 11 19 15 29 7 6 9 21 12 27 ...
%     2 8 16 23 10 25 4 22 14 20 28 30 18 24 32 26];

% chanMap = [1:64];
% -- Samtec Adaptor --
% chanMap = [38 34 44 61 40 36 46 60 42 48 41 58 45 35 39 56 ...
%            43 33 37 54 49 52 63 50 47 53 51 55 57 62 59 64 ...
%            1 6 3 8 10 14 12 16 17 11 13 9 7 4 5 2 ... 
%            28 24 26 19 21 32 15 30 18 31 29 22 20 25 23 27];
% -- DBC adaptor correct -- (MPM probeB)      
chanMap = [59 64 53 35 57 62 51 37 55 49 56 39 52 61 58 41 ...
           54 63 60 43 48 45 33 47 50 44 46 42 40 36 38 34 ...
           31 27 29 25 23 19 21 17 16 22 20 24 26 30 28 32 ...
           5 9 7 14 12 2 18 4 155 1 3 11 13 8 10 6];       
% %        ...
% %        65,66,67,68,69,70,71,72];
% -- DBC adaptor mirrored --
% chanMap  = [37 34 43 61 39 36 45 59 41 47 42 57 46 35 40 55 ...
%             44 33 38 53 50 51 63 49 48 54 52 56 58 62 60 64 ...
%             1 5 3 7 9 13 11 15 18 12 14 10 8 4 6 2 27 23 25 20 ...
%             22 32 16 30 17 31 29 21 19 26 24 28];
% chanMap = [27 32 21 3 25 30 19 5 23 17 24 7 20 29 26 9 ...
%            22 31 28 11 16 13 1 15 18 12 14 10 8 4 6 2 ...
%            63 59 61 57 55 51 53 49 48 54 52 56 58 62 60 64 ...
%            37 41 39 46 44 34 50 36 47 33 35 43 45 40 42 38];

       
% xcoords = xcoords/50;
% ycoords = ycoords/50;
% 
% ycoords = ycoords +1;
% xcoords = xcoords +2;

% T = zeros(12,3);
% 
% for c = 1:32
%     T(ycoords(c),xcoords(c)) = c;
% end


% the first thing Kilosort does is reorder the data with data = data(chanMap, :).
% Now we declare which channels are "connected" in this normal ordering, 
% meaning not dead or used for non-ephys data

connected = true(Nchannels, 1);
% connected  = [connected; false; false; false];
% connected(find(chanMap == 2)) = false;


% connected(64) = false;

% connected(1:2) = false;


% now we define the horizontal (x) and vertical (y) coordinates of these
% 34 channels. For dead or nonephys channels the values won't matter. Again
% I will take this information from the specifications of the probe. These
% are in um here, but the absolute scaling doesn't really matter in the
% algorithm. 
xcoords   = repmat([1 2 3 4]', 1, Nchannels/4);
xcoords   = 20*xcoords(:);

% xcoords = [xcoords; 1;1;1];
ycoords   = repmat(1:Nchannels/4, 4, 1);
%  ycoords(2,:) = ycoords(2,:)+1;
% ycoords(4,:) = ycoords(4,:)+1;
ycoords   = 20*ycoords;
ycoords(2,:) = ycoords(2,:)+10;
ycoords(4,:) = ycoords(4,:)+10;
ycoords = ycoords(:);
ycoords = abs(ycoords-330)+10;
% ycoords = flipud(ycoords);
% ycoords = [ycoords;1;1;1];
% Often, multi-shank probes or tetrodes will be organized into groups of
% channels that cannot possibly share spikes with the rest of the probe. This helps
% the algorithm discard noisy templates shared across groups. In
% this case, we set kcoords to indicate which group the channel belongs to.
% In our case all channels are on the same shank in a single group so we
% assign them all to group 1. 

kcoords   = ones(Nchannels,1);
if Nchannels> 64
    kcoords(65:end) = 0;
    connected(65:end) = false;

end
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