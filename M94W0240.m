function x = M94W0240()
% Data File Format Version 4.10

x = []; % initialize structure

% BEGIN HEADER
x.chamber = 2;
x.data_format = 4.10;
x.crown_amplifier_setting = 24;
x.ampl_dBSPL_at_1000Hz_0dB_att = 90;
x.analysis_type = 'WB-noise';
x.analysis_code = [61];
x.animal = 'M94W';
x.datetime = '05-Apr-2019 14:15:51';
x.hemisphere = 'Left';
x.hole_number = 3;
x.track_number = 4;
x.starting_depth = -1;
x.first_spike = -1;
x.unit_depth = 339;
x.unit_number = 2;
x.cf = -1;
x.threshold = -1;
x.latency = -1;
x.bandwidth = -1;
x.spikesort_ch1 = '';
x.spikesort_ch2 = '';
x.spikesort_ch3 = '';
x.spikesort_ch4 = '';
x.spikesort_ch5 = '';
x.spikesort_ch6 = '';
x.spikesort_ch7 = '';
x.spikesort_ch8 = '';
x.independent_var = -1;
x.nstim = 1;
  x.stimulus_tags_ch1 = { 'stim_num' 'spkr' 'attn(dB)' 'rep' 'len(ms)' 'event' 'delay(ms)' ' NOISE Center Frequency' ' Bandwidth' ' randn_seed' ' randn_seed' };
x.stimulus_ch1 = [
 % BEGIN STIMULUS CH 1
	1.0000	1.0000	50.0000	10.0000	100.0000	-1.0000	-1.0000	0.0080	99.0000	300101.0000	300101.0000
 % END STIMULUS CH 1
 ];
x.user_stimulus_desc_ch1 = {
	'Stimulus 1 : NOISE: Center Frequency  0.008       Bandwidth  99, randn_seed 300101 300101'
 };
x.user_parms_formatted_ch1 = 'true';
x.attenuation_ch1 = [	50.00	];
x.spkr_number_ch1 = [	1	];
 x.stimulus_tags_ch2 = { };
x.stimulus_ch2 = [
 % BEGIN STIMULUS CH 2
 % END STIMULUS CH 2
 ];
x.user_stimulus_desc_ch2 = {
 };
x.user_parms_formatted_ch2 = 'true';
x.attenuation_ch2 = [	-1.00	];
x.spkr_number_ch2 = [	-1	];
 x.presentation_mode = 'Random';
x.pre_stimulus_record_time = 200;
x.post_stimulus_record_time = 300;
x.iti_min = 200;
x.iti_max = 500;
x.iti = [ ];
x.spkr_tags = { 'speaker' 'azimuth' 'elevation' 'type' };
x.spkr = {
 % BEGIN SPEAKER
	'    Speaker 	ID 		  Azimuth 	 Elevation'
	'    1			B&W 		 +0			 0'
	'    2			EStim		-30			45'
 % END SPEAKER
 };
x.comments = '%';
x.data_tags = { 'stim_num' 'stim_rep' 'ch_num' 'event_time_microsecs' };
 x.data = [
 % BEGIN DATA
1 1 1 -1
1 1 1 3.072000e+01
1 1 2 -1
1 1 2 5.120000e+01
1 1 3 -1
1 1 3 4.096000e+01
1 1 4 -1
1 1 5 -1
1 1 6 -1
1 1 7 -1
1 2 1 -1
1 2 1 3.072000e+01
1 2 2 -1
1 2 2 5.120000e+01
1 2 3 -1
1 2 3 4.096000e+01
1 2 3 4.096000e+01
1 2 4 -1
1 2 5 -1
1 2 6 -1
1 2 7 -1
1 3 1 -1
1 3 1 3.072000e+01
1 3 2 -1
1 3 2 5.120000e+01
1 3 3 -1
1 3 3 4.096000e+01
1 3 3 4.096000e+01
1 3 3 4.096000e+01
1 3 4 -1
1 3 5 -1
1 3 6 -1
1 3 7 -1
1 4 1 -1
1 4 1 3.072000e+01
1 4 2 -1
1 4 2 5.120000e+01
1 4 3 -1
1 4 3 4.096000e+01
1 4 3 4.096000e+01
1 4 3 4.096000e+01
1 4 3 4.096000e+01
1 4 4 -1
1 4 5 -1
1 4 6 -1
1 4 7 -1
1 5 1 -1
1 5 1 3.072000e+01
1 5 2 -1
1 5 2 5.120000e+01
1 5 3 -1
1 5 3 4.096000e+01
1 5 3 4.096000e+01
1 5 3 4.096000e+01
1 5 3 4.096000e+01
1 5 3 4.096000e+01
1 5 4 -1
1 5 5 -1
1 5 6 -1
1 5 7 -1
1 6 1 -1
1 6 1 3.072000e+01
1 6 2 -1
1 6 2 5.120000e+01
1 6 3 -1
1 6 3 4.096000e+01
1 6 3 4.096000e+01
1 6 3 4.096000e+01
1 6 3 4.096000e+01
1 6 3 4.096000e+01
1 6 3 4.096000e+01
1 6 4 -1
1 6 5 -1
1 6 6 -1
1 6 7 -1
1 7 1 -1
1 7 1 3.072000e+01
1 7 2 -1
1 7 2 5.120000e+01
1 7 3 -1
1 7 3 4.096000e+01
1 7 3 4.096000e+01
1 7 3 4.096000e+01
1 7 3 4.096000e+01
1 7 3 4.096000e+01
1 7 3 4.096000e+01
1 7 3 4.096000e+01
1 7 4 -1
1 7 5 -1
1 7 6 -1
1 7 7 -1
1 8 1 -1
1 8 1 3.072000e+01
1 8 2 -1
1 8 2 5.120000e+01
1 8 3 -1
1 8 3 4.096000e+01
1 8 3 4.096000e+01
1 8 3 4.096000e+01
1 8 3 4.096000e+01
1 8 3 4.096000e+01
1 8 3 4.096000e+01
1 8 3 4.096000e+01
1 8 3 4.096000e+01
1 8 4 -1
1 8 5 -1
1 8 6 -1
1 8 7 -1
1 9 1 -1
1 9 1 3.072000e+01
1 9 2 -1
1 9 2 5.120000e+01
1 9 3 -1
1 9 3 4.096000e+01
1 9 3 4.096000e+01
1 9 3 4.096000e+01
1 9 3 4.096000e+01
1 9 3 4.096000e+01
1 9 3 4.096000e+01
1 9 3 4.096000e+01
1 9 3 4.096000e+01
1 9 3 4.096000e+01
1 9 4 -1
1 9 5 -1
1 9 6 -1
1 9 7 -1
1 10 1 -1
1 10 1 3.072000e+01
1 10 2 -1
1 10 2 5.120000e+01
1 10 3 -1
1 10 3 4.096000e+01
1 10 3 4.096000e+01
1 10 3 4.096000e+01
1 10 3 4.096000e+01
1 10 3 4.096000e+01
1 10 3 4.096000e+01
1 10 3 4.096000e+01
1 10 3 4.096000e+01
1 10 3 4.096000e+01
1 10 3 4.096000e+01
1 10 4 -1
1 10 5 -1
1 10 6 -1
1 10 7 -1
 % END DATA
 ];

x.trial_complete = 'true';

% END DATA FILE