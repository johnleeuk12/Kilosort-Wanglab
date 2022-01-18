# User Instructions
1. Edit "make_HarrisChannelMap" to modify channel mapping of used probe (64 or 128ch), and define which channels to remove from spike-sorting. The code is well commented for this.
2. Edit "concat_OE" to define the list of OE_folders (eg. '2021-04-20_17-12-18') to concatanate for spike-sorting. For each hole, track and session, there should be a list of these folders, where each folder corresponds to a recording segment. Also remember to edit the path to Kilosort repositories if needed.
3. Edit "Session_name" variable in "master_Harris3" eg(H4T1S1_concat) to define the name of the concatanated folder. Also remember to edit the path to Kilosort repositories if needed.
4. Run master_Harris3

