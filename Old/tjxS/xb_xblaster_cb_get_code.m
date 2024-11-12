function analysis_code = xb_xblaster_cb_get_analysis_code(analysis_type)
%
% Last edited: 01/30/2007 (SS)
%
% *** YOU MUST DOCUMENT ANY CHANGES MADE BELOW & MODIFY "Last edited" ABOVE ***
% This file is shared between M and C and should remain identical in two computers
% (72) 5/24/2019 added VirtualCall -JHL
% (71) 1/30/2018 added EStimPulseTrains - SDK
% (70) 03/9/2015 added DualChannelStimuli - sdk
% (69) 10/6/2013 added TORC - dwg
% (69) 5/6/2013 added 'User: pump.txt' - DG
% (68) 12/5/2011 added RHS_control -LF
% (67) 3/10/2010 added CIstimpanel_plusAcoustic -LAJ
% (66) 10/21/2009 Added CIstimpanel -LAJ
% (65) 07/28/2009 Added general clicks and mixed from XB2 - MJ
% (64) 07/02/2008 Added ABR stimulation panel - MJ
% (63) 03/18/2008 Added visual flash stimulation - MJ
% (62) 02/2008    Added combo_splinter GUI - MJ
% (61) 2007       Added spontaneous recording - MJ
% (60) 2007       Added acoustic/electric stimuluation interaction - MJ
% (59) 2007       Added electrical stimuluation - MJ
% (58) 09/06/2007 Added vocalization segments in conditioning stimulus - ER
% (57) 07/22/2006 Added 1201, 1210 (covered in STRFer block) - SS
% (56) 01/13/2006 Added 1203, 1204 (covered in STRFer block) - SS
% (55) 06/23/2005 Added 1205 (covered in STRFer block): twopip - SS
% (54) 06/21/2005 Added 1701-1800: General Clicks, EI
% (53) 05/12/2005 Added 1601 - 1700: FM_Interrupted, Poppy (SS)
% (52) 03/16/2005 Added 1208: Ripples, SS
% (51) 02/25/2004 Added 1207: Sparse Pips, SS
% (50) 10/10/2004 Add 1400-1600 (XB_Voc_Decom), S.Xu
% (49) 5/03/2004 Added TonePip codes (covered in STRFer block 1208-1209) - S. Sadagopan
% (48) 4/26/2004 Updated VV codes (version 1.2) C. DiMattina
% (47) 3/7/2004 Updated VV codes (version 1.1) C. DiMattina
% (46) 1/14/2004  Updated VV codes C. DiMattina
% (45) 1/06/2004  Updated VV codes C. DiMattina
% (44) 12/12/2003 Added more codes (1270-1279) for STRFer E. Issa
% (43) 11/30/2003 Recombined analysis codes from chambers 1 & 2 so that
% this file is now consistent between the two computers, E. Issa
% (42) 11/30/2003 Add 1301 ('linear_fm'), CD
% (41) 9/30/2003 Add 196 ('User: colonynoise&control_list.txt'), S.Xu
% (40) 9/30/2003 Add 191 ('User: babycries&control_list.txt'), S.Xu
% (39) 9/30/2003 Add 113 ('User: call_types&control_list.txt'), S.Xu
% (38) 9/11/2003 Add 1191 - Equal Spectrum noise for Pblaster, (D. Bendor)
% (37) 7/25/2003 Add 195-199 ('User: colonynoise_list.txt'), S. Xu
% (36) 7/25/2003 Add 190-194 ('User: babycries_list.txt'), S. Xu
% (35) 7/16/2003 Add 1201-1300 - STRFer (reconstruction stimuli), S. Sadagopan
% (34) 6/22/2003 Changed name of analyis types for virtual phee call manipulations (C. DiMattina)
% (33) 6/18/2003 Add 45 - linfm (C. DiMattina)
% (32) 5/25/2003 Add 1100-1200 - Pblaster (pitch stimuli), (D. Bendor)
% (31) 4/16/2003 Add 906-999 - Virtual Vocalizations (C. DiMattina)
% (30) 4/16/2003 Add 900-905 - STAR (program removed)
% (29) 4/22/2002 Add 891 gamma envelope modulation, E. Bartlett
% (28) 4/14/2002 Add 890 beta envelope modulation, E. Bartlett
% (27) 4/14/2002 Add 890-899 envelope modulation, E. Bartlett
% (26) 10/20/2001 871-877: Updated conditioning codes to include vocalizations, E. Bartlett
% (25) 8/22/2001 Upgraded comodulation to version 1.1, D. Barbour
% (24) 4/22/2001 Add new cases for FM,AMFM conditioning stimuli, E. Bartlett
% (23) 1/15/2001 Changed names of cases for conditioning stimuli, E. Bartlett
% (22) 1/14/2001 Add specific cases for conditioning stimuli, E. Bartlett
% (21) 1/7/2001 Add: 800-888 (Conditioning stimuli), E. Bartlett
% (20) 1/1/2001 Add: 771-2 (MWBEL v2.0, v2.1), D. Barbour
% (19) 1/1/2001 Add: 753-4 (MWBE v2.0, v2.1), D. Barbour
% (18) 12/12/2000 Add: 650 (Comodulation-Combination), D. Barbour
% (17) 12/12/2000 Add: 620 (Comodulation-Noise), D. Barbour
% (16) 12/11/2000 Add: 752 (MWBE v1.3), D. Barbour
% (15) 12/11/2000 Add: 751 (MWBE v1.2), D. Barbour
% (14) 11/20/2000 Add: 102 ('User: Twitter1_list.dat'), X. Wang
% (13) 11/19/2000 Add: 44 (fmtrain), 112 (TrilTwit_list.txt), X. Wang
% (12) 10/3/2000 Add: 750 (MWBE v1.1), D. Barbour
% (11) 4/24/2000 Add: 610 (Comodulation-Trill), D. Barbour
% (10) 4/6/2000 Add: 705 (wideband estimator v2.2), D. Barbour
% (9)  3/23/2000 Add: 704 (wideband estimator v2.10), D. Barbour
% (8)  3/14/2000 Add: 702-3 (wideband estimator v1.2, 2.0), D. Barbour
% (7)  3/13/2000 Add: 701 (wideband estimator v1.1), D. Barbour
% (6)  2/2/2000 Add: 700-799 (wideband estimator), D. Barbour
% (5)  12/14/1999 Add: 400-420 (ppclick, et. al), T. Lu
% (4)  12/1/1999 Add: 320-339, D. Barbour
% (3)  9/13/1999 Add: 600-699, D. Barbour
% (2)  9/3/1999  Add: 151-154, 161-164, 171-174, 180
% (1)  9/1/1999  Add: 300-319, 500-599
%switch lower(analysis_type)
switch (analysis_type)
    
    % unspecified data types
    case 'unspecified'
        analysis_code = -1;
    case 'User: temp_list.txt'  % e.g. M1k1562.dat
        analysis_code = -2;
    case 'User: out_list.txt'  % e.g. M1k1562.dat
        analysis_code = -3;
        
        % characterizations (1-9):
    case 'bandwidth'  % Freq-seq data
        analysis_code = 1;
    case 'ratelevel'  % Rate-level data
        analysis_code = 2;
    case 'latency'  % latency data
        analysis_code = 3;
        
        % Tones (10-19):
    case 'tones'
        analysis_code = 10;
        
        % Two stimuli groups (20-29):
    case 'Two_Tone'
        analysis_code = 20;
        
        % click stimuli (30-39):
    case 'click_vnum'
        analysis_code = 31;
    case 'clicks'
        analysis_code = 32;
    case 'gaussian noise click'
        analysis_code = 33;
    case 'click_test1'
        analysis_code = 34;
    case 'click_test2'
        analysis_code = 35;
    case 'rclick'
        analysis_code = 36;
    case 'gclick'
        analysis_code = 37;
    case 'gnclick'
        analysis_code = 38;
        
        % linear FM stimuli (40-49):
    case 'fmsweep'
        analysis_code = 41;
    case 'fmrdtrain'
        analysis_code = 42;
    case 'lFM'
        analysis_code = 43;
    case 'fmtrain'
        analysis_code = 44;
    case 'linfm'
        analysis_code = 45;
        
        % sAM/sFM stimuli (50-59):
    case 'amfm'
        analysis_code = 50;
        % later further classified as:
        %  sAM=51, sFM=52, sMM=53, BP-AM=54, WB-AM=55
    case 'sAM'
        analysis_code = 51;
    case 'sFM'
        analysis_code = 52;
    case 'sMM'
        analysis_code = 53;
    case 'BP-AM'
        analysis_code = 54;
    case 'WB-AM'
        analysis_code = 55;
        
        % noise stimuli (60-69):
    case 'noise'
        analysis_code = 60;
        % later further classified as:
        %  WB=61, BP=62, BR=63, LP=64, HP=65
    case 'WB-noise'
        analysis_code = 61;
    case 'BP-noise'
        analysis_code = 62;
    case 'BR-noise'
        analysis_code = 63;
    case 'LP-noise'
        analysis_code = 64;
    case 'HP-noise'
        analysis_code = 65;
        
        % envelope[ramped/damped], rise/fall-time stimuli (70-79):
    case 'ramped tone'
        analysis_code = 71;
    case 'rdamped'
        analysis_code = 72;
    case 'risetime'
        analysis_code = 73;
    case 'nrdamped'
        analysis_code = 74;
        
        % multiple sine components stimuli (80-89):
    case 'fpitch'
        analysis_code = 81;
        
        % *** User-panel stimuli ***
        
        % vocalizations (100-129):
    case 'User: Call_Types_list.txt'
        analysis_code = 100;
    case 'User: marmocalls_list.dat'   % = "Vocal 1"
        analysis_code = 101;
    case 'User: Twitter1_list.dat'
        analysis_code = 102;
    case 'User: Twitter2_list.dat'
        analysis_code = 103;
    case 'User: phee_list.txt'
        analysis_code = 104;
    case 'User: pheepeep_list.txt'
        analysis_code = 105;
    case 'User: trilpeep_list.txt'
        analysis_code = 106;
    case 'User: Trill_list.txt'
        analysis_code = 107;
    case 'User: twitphee_list.txt'
        analysis_code = 108;
    case 'User: trilphee_list.txt'
        analysis_code = 109;
    case 'User: tsikbark_list.txt'
        analysis_code = 110;
    case 'User: pheestrg_list.txt'
        analysis_code = 111;
    case 'User: triltwit_list.txt'
        analysis_code = 112;
    case 'User: call_types&control_list.txt'
        analysis_code = 113;
    case 'User: Call_Types_Lei_list.txt'
        analysis_code = 114;
    case 'User: pump.txt'
        analysis_code = 115;
        
        % sectioned calls (130-139):
    case 'User: Sectioned_Call_list.txt'
        analysis_code = 130;
        
        % twitter call discrimination analyses (140-159):
    case {'User: twitter_discrim_list1.txt','User: call_discrim_list1.txt'}
        analysis_code = 141;
    case {'User: twitter_discrim_list2.txt','User: call_discrim_list2.txt'}
        analysis_code = 142;
    case {'User: twitter_discrim_list3.txt','User: call_discrim_list3.txt'}
        analysis_code = 143;
    case {'User: twitter_discrim_list4.txt','User: call_discrim_list4.txt'}
        analysis_code = 144;
    case {'User: twitter_discrim_list5.txt','User: call_discrim_list5.txt'}
        analysis_code = 145;
    case {'User: twitter_discrim_list6.txt','User: call_discrim_list6.txt'}
        analysis_code = 146;
    case {'User: twitter_discrim_list7.txt','User: call_discrim_list7.txt'}
        analysis_code = 147;
        
    case {'User: twitter_discrim2_list1.txt'}
        analysis_code = 151;
    case {'User: twitter_discrim2_list2.txt'}
        analysis_code = 152;
    case {'User: twitter_discrim2_list3.txt'}
        analysis_code = 153;
    case {'User: twitter_discrim2_list4.txt'}
        analysis_code = 154;
        
        % phee call discrimination analyses (160-169):
    case {'User: phee_discrim_list1.txt'}
        analysis_code = 161;
    case {'User: phee_discrim_list2.txt'}
        analysis_code = 162;
    case {'User: phee_discrim_list3.txt'}
        analysis_code = 163;
    case {'User: phee_discrim_list4.txt'}
        analysis_code = 164;
        
        % trill call discrimination analyses (170-179):
    case {'User: trill_discrim_list1.txt'}
        analysis_code = 171;
    case {'User: trill_discrim_list2.txt'}
        analysis_code = 172;
    case {'User: trill_discrim_list3.txt'}
        analysis_code = 173;
    case {'User: trill_discrim_list4.txt'}
        analysis_code = 174;
        
        % trillphee call discrimination analyses (180-189):
    case {'User: trilphee_discrim_list.txt'}
        analysis_code = 180;
    case {'User: trilphee_discrim_list1.txt'}
        analysis_code = 181;
    case {'User: trilphee_discrim_list2.txt'}
        analysis_code = 182;
    case {'User: trilphee_discrim_list3.txt'}
        analysis_code = 183;
    case {'User: trilphee_discrim_list4.txt'}
        analysis_code = 184;
        
        % deaf and normal babies' cries analyses (190-194):
    case {'User: babycries_list.txt'}
        analysis_code = 190;
    case {'User: babycries&control_list.txt'}
        analysis_code = 191;
        
        % colony noise analyses (195-199):
    case {'User: colonynoise_list.txt'}
        analysis_code = 195;
    case {'User: colonynoise&control_list.txt'}
        analysis_code = 196;
        
        
        % ***** Vocalization panel stimuli (200-): *****
    case 'vocal'
        analysis_code = 200;
        
        % modulated call discrimination analysis (300-319):
    case 'User: Call_Types_Mod_list.txt'
        analysis_code = 300;
        
        % temporal call discrimination analysis (320-339):
    case 'User: Call_Types_Temp_list.txt'
        analysis_code = 320;
        
        % random interval click train stimuli (400-420):
    case 'ppclick'
        analysis_code = 400;
        
        % multichannel synchrony experiment (500-599):
    case 'Synchrony'
        analysis_code = 500;
    case 'Synchrony-Grouping'
        analysis_code = 510;
    case 'Synchrony-Capture'
        analysis_code = 520;
    case 'Synchrony-Comod/Mod'
        analysis_code = 530;
        
        % multichannel comodulation experiment(600-699):
    case 'Comodulation'  % tone-based stimuli
        analysis_code = 600;
    case 'Comodulation-Trill'
        analysis_code = 610;
    case 'Comodulation-Noise'
        analysis_code = 620;
    case 'Comodulation-Combination'
        analysis_code = 650;
    case 'Comodulation v2.0'  % tone-based stimuli
        analysis_code = 601;
    case 'Comodulation-Trill v2.0'
        analysis_code = 611;
    case 'Comodulation-Noise v2.0'
        analysis_code = 621;
    case 'Comodulation-Combination v2.0'
        analysis_code = 651;
        
        % wideband estimator (700-799)
    case 'Wideband Estimator'
        analysis_code = 700;
    case 'Wideband Estimator v1.1'
        analysis_code = 701;
    case 'Wideband Estimator v1.2'
        analysis_code = 702;
    case 'Wideband Estimator v2.0'
        analysis_code = 703;
    case 'Wideband Estimator v2.1'
        analysis_code = 704;
    case 'Wideband Estimator v2.2'
        analysis_code = 705;
        
    case 'Modulated Wideband Estimator (MWBE) v1.1'
        analysis_code = 750;
    case 'Modulated Wideband Estimator (MWBE) v1.2'
        analysis_code = 751;
    case 'Modulated Wideband Estimator (MWBE) v1.3'
        analysis_code = 752;
    case 'Modulated Wideband Estimator (MWBE) v2.0'
        analysis_code = 753;
    case 'Modulated Wideband Estimator (MWBE) v2.1'
        analysis_code = 754;
        
    case 'Modulated Wideband Estimator Levels (MWBEL) v2.0'
        analysis_code = 771;
    case 'Modulated Wideband Estimator Levels (MWBEL) v2.1'
        analysis_code = 772;
        
        % Conditioning stimuli(800-888)
        % Numbers 801-810 and 820,830,etc are open analysis codes
    case 'Conditioning stimuli'
        analysis_code = 800;
    case 'Conditioning stimulus:tone-tone'
        analysis_code = 811;
    case 'Conditioning stimulus:tone-noise'
        analysis_code = 812;
    case 'Conditioning stimulus:tone-AM tone'
        analysis_code = 813;
    case 'Conditioning stimulus:tone-AM noise'
        analysis_code = 814;
    case 'Conditioning stimulus:tone-FM tone'
        analysis_code = 815;
    case 'Conditioning stimulus:tone-AMFM tone'
        analysis_code = 816;
    case 'Conditioning stimulus:tone-vocalization'
        analysis_code = 817;
    case 'Conditioning stimulus:tone-vocalization segment'
        analysis_code = 818;
    case 'Conditioning stimulus:noise-tone'
        analysis_code = 821;
    case 'Conditioning stimulus:noise-noise'
        analysis_code = 822;
    case 'Conditioning stimulus:noise-AM tone'
        analysis_code = 823;
    case 'Conditioning stimulus:noise-AM noise'
        analysis_code = 824;
    case 'Conditioning stimulus:noise-FM tone'
        analysis_code = 825;
    case 'Conditioning stimulus:noise-AMFM tone'
        analysis_code = 826;
    case 'Conditioning stimulus:noise-vocalization'
        analysis_code = 827;
    case 'Conditioning stimulus:noise-vocalization segment'
        analysis_code = 828;
    case 'Conditioning stimulus:AM tone-tone'
        analysis_code = 831;
    case 'Conditioning stimulus:AM tone-noise'
        analysis_code = 832;
    case 'Conditioning stimulus:AM tone-AM tone'
        analysis_code = 833;
    case 'Conditioning stimulus:AM tone-AM noise'
        analysis_code = 834;
    case 'Conditioning stimulus:AM tone-FM tone'
        analysis_code = 835;
    case 'Conditioning stimulus:AM tone-AMFM tone'
        analysis_code = 836;
    case 'Conditioning stimulus:AM tone-vocalization'
        analysis_code = 837;
    case 'Conditioning stimulus:AM tone-vocalization segment'
        analysis_code = 838;
    case 'Conditioning stimulus:AM noise-tone'
        analysis_code = 841;
    case 'Conditioning stimulus:AM noise-noise'
        analysis_code = 842;
    case 'Conditioning stimulus:AM noise-AM tone'
        analysis_code = 843;
    case 'Conditioning stimulus:AM noise-AM noise'
        analysis_code = 844;
    case 'Conditioning stimulus:AM noise-FM tone'
        analysis_code = 845;
    case 'Conditioning stimulus:AM noise-vocalization'
        analysis_code = 847;
    case 'Conditioning stimulus:AM noise-AMFM tone'
        analysis_code = 846;
    case 'Conditioning stimulus:AM noise-vocalization segment'
        analysis_code = 848;
    case 'Conditioning stimulus:FM tone-tone'
        analysis_code = 851;
    case 'Conditioning stimulus:FM tone-noise'
        analysis_code = 852;
    case 'Conditioning stimulus:FM tone-AM tone'
        analysis_code = 853;
    case 'Conditioning stimulus:FM tone-AM noise'
        analysis_code = 854;
    case 'Conditioning stimulus:FM tone-FM tone'
        analysis_code = 855;
    case 'Conditioning stimulus:FM tone-AMFM tone'
        analysis_code = 856;
    case 'Conditioning stimulus:FM tone-vocalization'
        analysis_code = 857;
    case 'Conditioning stimulus:FM tone-vocalization segment'
        analysis_code = 858;
    case 'Conditioning stimulus:FM tone-tone'
        analysis_code = 861;
    case 'Conditioning stimulus:AMFM tone-noise'
        analysis_code = 862;
    case 'Conditioning stimulus:AMFM tone-AM tone'
        analysis_code = 863;
    case 'Conditioning stimulus:AMFM tone-AM noise'
        analysis_code = 864;
    case 'Conditioning stimulus:AMFM tone-FM tone'
        analysis_code = 865;
    case 'Conditioning stimulus:AMFM tone-AMFM tone'
        analysis_code = 866;
    case 'Conditioning stimulus:AMFM tone-vocalization'
        analysis_code = 867;
    case 'Conditioning stimulus:AMFM tone-vocalization segment'
        analysis_code = 868;
    case 'Conditioning stimulus:vocalization-tone'
        analysis_code = 871;
    case 'Conditioning stimulus:vocalization-noise'
        analysis_code = 872;
    case 'Conditioning stimulus:vocalization-AM tone'
        analysis_code = 873;
    case 'Conditioning stimulus:vocalization-AM noise'
        analysis_code = 874;
    case 'Conditioning stimulus:vocalization-FM tone'
        analysis_code = 875;
    case 'Conditioning stimulus:vocalization-AMFM tone'
        analysis_code = 876;
    case 'Conditioning stimulus:vocalization-vocalization'
        analysis_code = 877;
    case 'Conditioning stimulus:vocalization-vocalization segment'
        analysis_code = 878;
    case 'Conditioning stimulus:vocalization segment-tone'
        analysis_code = 881;
    case 'Conditioning stimulus:vocalization segment-noise'
        analysis_code = 882;
    case 'Conditioning stimulus:vocalization segment-AM tone'
        analysis_code = 883;
    case 'Conditioning stimulus:vocalization segment-AM noise'
        analysis_code = 884;
    case 'Conditioning stimulus:vocalization segment-FM tone'
        analysis_code = 885;
    case 'Conditioning stimulus:vocalization segment-AMFM tone'
        analysis_code = 886;
    case 'Conditioning stimulus:vocalization segment-vocalization'
        analysis_code = 887;
    case 'Conditioning stimulus:vocalization segment-vocalization segment'
        analysis_code = 888;
        
        
        % Modulation envelope (890-899)
    case 'beta'
        analysis_code = 890;
    case 'gamma'
        analysis_code = 891;
        
        
        % Virtual Vocalizations (900-1100)
        
        % 901-909 : Token Sets
    case {'xb_vv_1.0_real_virtual','xb_vv_1.1_real_virtual','xb_vv_1.2_real_virtual'},
        analysis_code=901;
    case {'xb_vv_1.0_tsik_egg_token','xb_vv_1.1_tsik_egg_token','xb_vv_1.2_tsik_egg_token'},
        analysis_code=902;
        
        % 910-911 : Twit Space
    case {'xb_vv_1.0_twit_space','xb_vv_1.1_twit_space','xb_vv_1.2_twit_space'},
        analysis_code=910;
        
        % 912-920 : Probe Panels
    case {'vblaster_probe_all','xb_vv_1.0_probe','xb_vv_1.1_probe','xb_vv_1.2_probe'}
        analysis_code=912;
    case 'vblaster_probe_twitter'
        analysis_code=913;
    case 'vblaster_probe_trill'
        analysis_code=914;
    case 'vblaster_probe_trillphee'
        analysis_code=915;
    case 'vblaster_probe_phee'
        analysis_code=916;
    case {'vblaster_probe_hcs','xb_vv_1.0_harmonic_probe','xb_vv_1.1_harmonic_probe','xb_vv_1.2_harmonic_probe'}
        analysis_code=917;
        
        
        % 921-930 : Stimulus Decompositions
    case {'vblaster_amfm_rev','xb_vv_1.0_twitter_rev','xb_vv_1.1_twitter_rev','xb_vv_1.2_twitter_rev'}
        analysis_code=921;
    case {'vblaster_amfm_esn','xb_vv_1.0_phee_ssn', ...
            'xb_vv_1.0_trill_ssn','xb_vv_1.0_twitter_ssn'}
        analysis_code=922;
    case {'vblaster_amfm_twitter','xb_vv_1.0_trill_decompose', ...
            'xb_vv_1.0_twitter_decompose','xb_vv_1.0_phee_decompose', ...
            'xb_vv_1.1_trill_decompose', 'xb_vv_1.1_twitter_decompose','xb_vv_1.1_phee_decompose', ...
            'xb_vv_1.2_trill_decompose', 'xb_vv_1.2_twitter_decompose','xb_vv_1.2_phee_decompose'}
        analysis_code=923;
        % 931-939 : General Paradigms
    case {'vblaster_hcs','xb_vv_1.0_harmonic_series'}
        analysis_code=931;
    case {'xb_vv_1.0_harmonic_series_f2f1','xb_vv_1.1_harmonic_series_f2f1','xb_vv_1.2_harmonic_series_f2f1'}
        analysis_code=932;
    case {'xb_vv_1.0_harmonic_series_a2a1','xb_vv_1.1_harmonic_series_a2a1','xb_vv_1.2_harmonic_series_a2a1'}
        analysis_code=933;
    case 'xb_vv_1.0_snr_series'
        analysis_code=936;
    case {'xb_vv_1.0_conditioning','xb_vv_1.1_conditioning'}
        analysis_code=937;
        
        % 940-959 : Trill/Trillphee
    case {'vblaster_trill_selectivity','xb_vv_1.0_trill_select','xb_vv_1.1_trill_select','xb_vv_1.2_trill_select'}
        analysis_code=941;
    case {'vblaster_trill_multiple','xb_vv_1.0_trill_domain','xb_vv_1.1_trill_domain'}
        analysis_code=942;
    case 'vblaster_trill_fm_rate'
        analysis_code=943;
    case {'vblaster_trill_fm_depth','xb_vv_1.0_trill_single_fm'}
        analysis_code=944;
    case 'vblaster_trill_am_rate'
        analysis_code=945;
    case {'vblaster_trill_am_depth','xb_vv_1.0_trill_single_am'}
        analysis_code=946;
    case {'vblaster_trill_am1_fm1_phase','xb_vv_1.0_trill_single_th'}
        analysis_code=947;
    case 'vblaster_trill_am2_fm2_phase'
        analysis_code=948;
    case {'vblaster_trill_trans','xb_vv_1.0_trillphee_transition','xb_vv_1.1_trillphee_transition','xb_vv_1.2_trillphee_transition'}
        analysis_code=949;
    case {'vblaster_trill_trillrate','xb_vv_1.0_trill_single_rt'}
        analysis_code=950;
    case {'xb_vv_1.0_trill_single','xb_vv_1.1_trill_single','xb_vv_1.2_trill_single'}
        analysis_code=951;
    case {'xb_vv_1.0_trill_dust','xb_vv_1.1_trill_dust','xb_vv_1.2_trill_dust'}
        analysis_code=952;
        
        % 960-976 : Twitter
    case {'vblaster_twitter_ipi','xb_vv_1.0_twitter_single_ipi'}
        analysis_code=966;
    case {'vblaster_twitter_tsw','xb_vv_1.0_twitter_single_tsw'}
        analysis_code=967;
    case 'vblaster_twitter_fkn'
        analysis_code=968;
    case {'vblaster_twitter_selectivity','xb_vv_1.0_twitter_select','xb_vv_1.1_twitter_select','xb_vv_1.2_twitter_select'}
        analysis_code=969;
    case {'vblaster_twitter_multiple','xb_vv_1.0_twitter_domain','xb_vv_1.1_twitter_domain'}
        analysis_code=970;
    case {'xb_vv_1.0_twitter_single','xb_vv_1.1_twitter_single','xb_vv_1.2_twitter_single'}
        analysis_code=971;
    case {'xb_vv_1.0_twitter_dust','xb_vv_1.1_twitter_dust','xb_vv_1.2_twitter_dust'}
        analysis_code=972;
        
        % 977-979 : Phee Call
    case 'vblaster_phee_dur'
        analysis_code=977;
    case 'vblaster_phee_mod'
        analysis_code=978;
    case {'vblaster_phee_selectivity','xb_vv_1.0_phee_dust','xb_vv_1.1_phee_dust','xb_vv_1.2_phee_dust'}
        analysis_code=979;
        %case 'vblaster_phee_multiple'
        %   analysis_code=980;
        
        % 980-983 : Morphing (Vblaster)
    case 'vblaster_morph_twit_chimera'
        analysis_code=980;
    case 'vblaster_morph_twit_morph'
        analysis_code=981;
    case 'vblaster_morph_trill_chimera'
        analysis_code=982;
    case 'vblaster_morph_trill_morph'
        analysis_code=983;
        
        % 984-990 : Morphing (VV)
    case {'xb_vv_1.0_trill_phee_chimera','xb_vv_1.1_trill_phee_chimera','xb_vv_1.2_trill_phee_chimera'}
        analysis_code=984;
    case {'xb_vv_1.0_trill_phee_morph','xb_vv_1.1_trill_phee_morph','xb_vv_1.2_trill_phee_morph'}
        analysis_code=985;
        
        % 990-1000 : Twitter Caller Discrimination
    case 'xb_vv_1.2_twit_disc_test'
        analysis_code=990;
    case 'xb_vv_1.2_twit_disc_ab'
        analysis_code=991;
    case 'xb_vv_1.2_twit_disc_ac'
        analysis_code=992;
    case 'xb_vv_1.2_twit_disc_bc'
        analysis_code=993;
        
        % 1000-1010 : Trill/Phee Call Discrimination
    case 'xb_vv_1.2_trill_phee_disc'
        analysis_code=1000;
        
        % 1010 - 1100 : Free
        
        % 1099 : End VV Block
        
        %1100-1200 : Pblaster panel (pitch)
    case 'pblaster_IteratedRippleNoise'
        analysis_code=1101;
        
    case 'pblaster_GaussianClicktrain_ToneCarrier'
        analysis_code=1111;
    case 'pblaster_GaussianClicktrain_NoiseCarrier'
        analysis_code=1112;
        
        
    case 'pblaster_RampDampClicktrain_ToneCarrier'
        analysis_code=1121;
    case 'pblaster_RampDampClicktrain_NoiseCarrier'
        analysis_code=1122;
        
    case 'pblaster_RectangularClicktrain_SamePolarity'
        analysis_code=1131;
    case 'pblaster_RectangularClicktrain_AltPolarity'
        analysis_code=1132;
        
        
    case 'pblaster_FMsweepClicktrain_SamePolarity'
        analysis_code=1141;
    case 'pblaster_FMsweepClicktrain_AltPolarity'
        analysis_code=1142;
        
    case 'pblaster_TwoToneParadigm'
        analysis_code=1161;
        
    case 'pblaster_HarmonicComplexTone_SinePhase'
        analysis_code=1171;
    case 'pblaster_HarmonicComplexTone_CosinePhase'
        analysis_code=1172;
    case 'pblaster_HarmonicComplexTone_AltPhase'
        analysis_code=1173;
    case 'pblaster_HarmonicComplexTone_SchroederPhase'
        analysis_code=1174;
        
    case 'pblaster_NoiseMaskMissingFund_SinePhase'
        analysis_code=1181;
    case 'pblaster_NoiseMaskMissingFund_CosinePhase'
        analysis_code=1182;
    case 'pblaster_NoiseMaskMissingFund_AltPhase'
        analysis_code=1183;
    case 'pblaster_NoiseMaskMissingFund_SchroederPhase'
        analysis_code=1184;
        
    case 'pblaster_EqualSpectrumNoise'
        analysis_code=1191;
        
        % 1201 - 1300 (S. Sadagopan)
        
        % 1201: ga_optim
    case 'ga_optim'
        analysis_code = 1201;
        
        % 1202: fl_path
    case 'fl_path'
        analysis_code = 1202;
        
        % 1203: fl_ellipse
    case 'fl_ellipse'
        analysis_code = 1203;
        
        % 1204: fl_plane
    case 'fl_plane'
        analysis_code = 1204;
        
        % 1205: twopip
    case 'twopip'
        analysis_code = 1205;
        
        % 1206: Ripples
    case 'ripples'
        analysis_code = 1206;
        
        % 1207 : Sparse tonepips
    case 'sparse_pips'
        analysis_code = 1207;
        
        % 1208 - 1209 : Tone Pip Reverse Correlation
    case 'tonepip_strf'
        analysis_code = 1208;
    case 'tonepip_validate'
        analysis_code = 1209;
        
        % 1210 : ga2_optim
    case 'ga2_optim'
        analysis_code = 1210;
        
        % 1211 : xbts
    case 'xbts'
        analysis_code = 1211;
        % 1212 : linfm2
    case 'xblinfm2'
        analysis_code = 1212;
        % 1213: xstep (contrast adaptation)
    case 'xstep'
        analysis_code = 1213;
        % 1214 - 1219 : Clicks (reserved for future use)
        
        % 1221 - 1239 : Tone-based stimuli
        % 1220 - space
    case 'tone_moving_flat'
        analysis_code=1221;
    case 'tone_moving_LP'
        analysis_code=1222;
    case 'tone_moving_SH'
        analysis_code=1223;
    case 'tone_moving_comod'
        analysis_code=1224;
        % 1225 - space
    case 'harm_moving_flat'
        analysis_code=1226;
    case 'harm_moving_LP'
        analysis_code=1227;
    case 'harm_moving_SH'
        analysis_code=1228;
    case 'harm_moving_comod'
        analysis_code=1229;
        % 1230 - space
    case 'tone_null_flat'
        analysis_code=1231;
    case 'tone_null_LP'
        analysis_code=1232;
    case 'tone_null_SH'
        analysis_code=1233;
    case 'tone_null_comod'
        analysis_code=1234;
        % 1235 - space
    case 'harm_null_flat'
        analysis_code=1236;
    case 'harm_null_LP'
        analysis_code=1237;
    case 'harm_null_SH'
        
        analysis_code=1238;
    case 'harm_null_comod'
        analysis_code=1239;
        
        % 1240 - 1259 : Noise-based stimuli
        % 1240 - space
    case 'noise_moving_flat'
        analysis_code=1241;
    case 'noise_moving_LP'
        analysis_code=1242;
    case 'noise_moving_SH'
        analysis_code=1243;
    case 'noise_moving_comod'
        analysis_code=1244;
        % 1245 - space
    case 'noise_null_flat'
        analysis_code=1246;
    case 'noise_null_LP'
        analysis_code=1247;
    case 'noise_null_SH'
        analysis_code=1248;
    case 'noise_null_comod'
        analysis_code=1249;
        % 1250 - 1259: reserved for future use.
        
        % 1260 - 1279 : IRN-based stimuli
        % 1260 - space
    case 'irn_moving_flat'
        analysis_code=1261;
    case 'irn_moving_LP'
        analysis_code=1262;
    case 'irn_moving_SH'
        analysis_code=1263;
    case 'irn_moving_comod'
        analysis_code=1264;
        % 1265 - space
    case 'irn_null_flat'
        analysis_code=1266;
    case 'irn_null_LP'
        analysis_code=1267;
    case 'irn_null_SH'
        analysis_code=1268;
    case 'irn_null_comod'
        analysis_code=1269;
        
        % 1270 - 1279: FM-based stimuli
        % 1270 - space
    case 'fm_moving_flat'
        analysis_code=1271;
    case 'fm_moving_LP'
        analysis_code=1272;
    case 'fm_moving_SH'
        analysis_code=1273;
    case 'fm_moving_comod'
        analysis_code=1274;
        % 1275 - space
    case 'fm_null_flat'
        analysis_code=1276;
    case 'fm_null_LP'
        analysis_code=1277;
    case 'fm_null_SH'
        analysis_code=1278;
    case 'fm_null_comod'
        analysis_code=1279;
        
        % 1280 - 1299 : Reconstruction test stimuli
        
        % 1301 (Linear FM, CD)
    case 'linear_fm'
        analysis_code=1301;
        
        % 1400-1600 (XB_Voc_Decom, SX)
    case 'Voc_Decom: Original Envelope'
        analysis_code = 1401;
    case 'Voc_Decom: Original Fine Structure'
        analysis_code = 1451;
    case 'Voc_Decom: Both Envelope and Fine Structure Processing'
        analysis_code = 1501;
    case 'Voc_Decom: Band Pass Signal Processing'
        analysis_code=1551;
        
        % 1601 - 1700 (fm_interrupted, PC)
    case 'fm_int_Sweep_Noise_Tone_PF'
        analysis_code=1601;
    case 'fm_int_Absent_Noise_Tone_PF'
        analysis_code=1602;
    case 'fm_int_Sweep_Noise_Tone_SF'
        analysis_code=1603;
    case 'fm_int_Absent_Noise_Tone_SF'
        analysis_code=1604;
    case 'fm_int_Sweep_Noise_Tone_NseL'
        analysis_code=1605;
    case 'fm_int_Absent_Noise_Tone_NseL'
        analysis_code=1606;
    case 'fm_int_Sweep_Silence_Tone_PF'
        analysis_code=1611;
    case 'fm_int_Absent_Silence_Tone_PF'
        analysis_code=1612;
    case 'fm_int_Sweep_Silence_Tone_SF'
        analysis_code=1613;
    case 'fm_int_Absent_Silence_Tone_SF'
        analysis_code=1614;
    case 'fm_int_Sweep_Silence_Tone_NseL'
        analysis_code=1615;
    case 'fm_int_Absent_Silence_Tone_NseL'
        analysis_code=1616;
    case 'fm_int_Sweep_Noise_Sweep_PF'
        analysis_code=1621;
    case 'fm_int_Absent_Noise_Sweep_PF'
        analysis_code=1622;
    case 'fm_int_Sweep_Noise_Sweep_SF'
        analysis_code=1623;
    case 'fm_int_Absent_Noise_Sweep_SF'
        analysis_code=1624;
    case 'fm_int_Sweep_Noise_Sweep_NseL'
        analysis_code=1625;
    case 'fm_int_Absent_Noise_Sweep_NseL'
        analysis_code=1626;
    case 'fm_int_Sweep_Silence_Sweep_PF'
        analysis_code=1631;
    case 'fm_int_Absent_Silence_Sweep_PF'
        analysis_code=1632;
    case 'fm_int_Sweep_Silence_Sweep_SF'
        analysis_code=1633;
    case 'fm_int_Absent_Silence_Sweep_SF'
        analysis_code=1634;
    case 'fm_int_Sweep_Silence_Sweep_NseL'
        analysis_code=1635;
    case 'fm_int_Absent_Silence_Sweep_NseL'
        analysis_code=1636;
        
        
        % GENERAL CLICKS AND MIXED FROM XB2
    case 'general_clicks'
        analysis_code = 1701;
    case 'tones-masked'
        analysis_code = 1704;
    case 'sAM-masked'
        analysis_code = 1706;
    case 'BP-AM-mixed'
        analysis_code = 1703;
    case 'BP-noise-masked'
        analysis_code = 1705;
    case 'BP-AM-masked'
        analysis_code = 1707;
    case 'sFM-mixed'
        analysis_code = 1710;
    case 'sMM-mixed'
        analysis_code = 1711;
    case 'mixed'
        analysis_code = 1702;
    case 'BR-AM'
        analysis_code = 1712;
    case 'LP-AM'
        analysis_code = 1713;
    case 'HP-AM'
        analysis_code = 1714;
    case 'Two_Noise'
        analysis_code = 1708;
    case 'ABA'
        analysis_code = 1709;
        
        
        % Marcus: 1800 - 1899
    case 'xb_estim_01'
        analysis_code = 1801;
    case 'xb_estim_astim_01'
        analysis_code = 1821;
    case 'xb_spontaneousactivity_01'
        analysis_code = 1841;
    case 'xb_flash_01'
        analysis_code = 1851;
    case 'xb_Click_03'
        analysis_code = 1861;
        
        % Marcus: ABR panels   1901 - 1920
    case 'xb_ABR_tone'
        analysis_code = 1901;
    case 'xb_ABR_click'
        analysis_code = 1902;
        
        %% luke johnson 2101-2200
        %% Cochlear Implant Stimuli (CIstim)
    case 'CIstim_rectpulse'
        analysis_code = 2101;
    case 'CIstim_rectpulse_plusAcoustic'
        analysis_code = 2102;
    case 'CIstim_rectpulse_plusAcoustic2'
        analysis_code = 2103;
        
        %%EVAN taking 2201-2300
        
        %% John 2301-2400
        
    case 'puretone_cf'
        analysis_code = 2301;
    case 'puretone_attn'
        analysis_code = 2302;
        
    case 'puretone_dur'
        analysis_code = 2303;
    case 'noise_cf'
        analysis_code = 2304;
    case 'noise_attn'
        analysis_code = 2305;
    case 'noise_dur'
        analysis_code = 2306;
    case 'noise_bw'
        analysis_code = 2307;
    case 'noise_seed'
        analysis_code = 2308;
    case 'harmonic_cf'
        analysis_code = 2309;
    case 'harmonic_attn'
        analysis_code = 2310;
    case 'harmonic_dur'
        analysis_code = 2311;
    case 'harmonic_rattn'
        analysis_code = 2312;
    case 'harmonic_order'
        analysis_code = 2313;
    case 'sFM_rate'
        analysis_code = 2314;
    case 'sFM_depth'
        analysis_code = 2315;
    case 'sAM_rate'
        analysis_code = 2316;
    case 'sAM_depth'
        analysis_code = 2317;
    case 'sFMAM_rate'
        analysis_code = 2318;
    case 'complextone'
        analysis_code = 2319;
    case 'phee_cf'
        analysis_code = 2320;
    case 'phee_attn'
        analysis_code = 2321;
    case 'phee_dur'
        analysis_code = 2322;
    case 'phee_FMmod'
        analysis_code = 2323;
    case 'phee_FMrate'
        analysis_code = 2324;
    case 'phee_FMdepth'
        analysis_code = 2325;
    case 'phee_AMmod'
        analysis_code = 2326;
    case 'phee_AMrate'
        analysis_code = 2327;
    case 'phee_AMdepth'
        analysis_code = 2328;
    case 'phee_FMAMrate'
        analysis_code = 2329;
    case 'phee_Hrattn'
        analysis_code = 2330;
    case 'phee_Horder'
        analysis_code = 2331;
    case 'phee_Hbw'
        analysis_code = 2332;
    case 'phee_SNR'
        analysis_code = 2333;
    case 'phee'
        analysis_code = 2334;
    case 'trill_cf'
        analysis_code =2335;
    case 'trill_attn'
        analysis_code = 2336;
    case 'trill_dur'
        analysis_code = 2337;
    case 'trill_FMmod'
        analysis_code = 2338;
    case 'trill_FMrate'
        analysis_code =2339;
    case 'trill_FMdepth'
        analysis_code = 2340;
    case 'trill_AMmod'
        analysis_code = 2341;
    case 'trill_AMrate'
        analysis_code = 2342;
    case 'trill_AMdepth'
        analysis_code =2343;
    case 'trill_FMAMrate'
        analysis_code = 2344;
    case 'trill_Hrattn'
        analysis_code = 2345;
    case 'trill_Horder'
        analysis_code = 2346;
    case 'trill_Hbw'
        analysis_code =2347;
    case 'trill_SNR'
        analysis_code = 2348;
    case 'trill'
        analysis_code = 2349;
    case 'trillphee_cf'
        analysis_code = 2350;
    case 'trillphee_attn'
        analysis_code =2351;
    case 'trillphee_dur'
        analysis_code = 2352;
    case 'trillphee_FMmod'
        analysis_code = 2353;
    case 'trillphee_FMrate'
        analysis_code = 2354;
    case 'trillphee_FMdepth'
        analysis_code =2356;
    case 'trillphee_AMmod'
        analysis_code = 2357;
    case 'trillphee_AMrate'
        analysis_code = 2358;
    case 'trillphee_AMdepth'
        analysis_code =  2359;
    case 'trillphee_FMAMrate'
        analysis_code = 2360;
    case 'trillphee_tTrans'
        analysis_code = 2361;
    case 'trillphee_Hrattn'
        analysis_code = 2362;
    case 'trillphee_Horder'
        analysis_code = 2363;
    case 'trillphee_Hbw'
        analysis_code =2364;
    case 'trillphee_SNR'
        analysis_code = 2365;
    case 'trillphee'
        analysis_code = 2366;
    case 'twitter_cf'
        analysis_code = 2367;
    case 'twitter_attn'
        analysis_code = 2368;
    case 'twitter_nphr'
        analysis_code = 2369;
    case 'twitter_IPI'
        analysis_code = 2370;
    case 'twitter_tphr'
        analysis_code = 2371;
    case 'twitter_bwphr'
        analysis_code = 2372;
    case 'twitter_fknee'
        analysis_code = 2373;
    case 'twitter_tknee'
        analysis_code = 2374;
    case 'twitter_Hrattn'
        analysis_code = 2375;
    case 'twitter_SNR'
        analysis_code = 2376;
    case 'twitter'
        analysis_code = 2377;
        
        
        %% Darik 2501-2600
    case {'TORC', 'TORCs'}
        analysis_code = 2501;
    case 'ArbitraryAudio'
        analysis_code = 2502;
    case 'Ripples'
        analysis_code = 2503;
    case {'GeneticRipplesDummy'}
        analysis_code = 2504;
        
        %%%%%% Lei Feng 2601 - 2700
    case 'RSS'
        analysis_code = 2601;
    case 'RHS'
        analysis_code = 2602;
    case 'Pair_RHS'
        analysis_code = 2603;
    case 'Tone_Stacks'
        analysis_code = 2604;
    case 'Two_Tone_Stacks'
        analysis_code = 2605;
    case 'Two_BP'
        analysis_code = 2606;
    case 'jitter_stacks'
        analysis_code = 2607;
    case 'IHTC'
        analysis_code = 2608;
    case 'optimal_model1'
        analysis_code = 2609;
    case 'optimal_model2'
        analysis_code = 2610;
        
        %% Seth 2701 - 2800
    case 'DualChannelStimuli'
        analysis_code = 2701;
        
        %%%%% Seth Koehler (VNS) 3000 - 3100
    case 'EStimPulseTrains'
        analysis_code=3000;
        
    otherwise
        analysis_code = 0;
        disp('*** xb_xblaster_cb_get_code: analysis_code NOT defined <Pause> ***');
        pause;
end
