% 
% Eunji Jung: neuneunji@gmail.com / eunji_jung@kaist.ac.kr
% Next step of running [pool_data_ephys.m] 
% 1. jitter correction: intan input timestamp vs. presentation log port code timestamp
% 2. sorting eventtype
% 3. Block change time, TTR 
% 4. lickbout onset, offset and duration
% ver1. 2024-Aug-22
% ver2. 2025-May-20, jitter correction edited 



% clear; clc; close all; 
% roi =3; 
% if roi == 1
%     load PPCPool.mat
%     datasource = PPCPool;
%     clear PPCPool
%     roi2 = ['PPC'];
% elseif roi == 2
%     load ACPool.mat
%     datasource = ACPool; 
%     clear ACPool
%     roi2 = ['AC'];
% elseif roi == 3
%     load ICPool.mat
%     datasource = ICPool; 
%     clear ICPool
%     roi2 = ['IC'];
% end

datasource = Pool;
%%%

stim = 500;         delay = 1000;       rw = 2000;      
lickboutinterval = 2000;                itilick = 4000; 
lickpre = 2000; 
stimperiod = stim+delay+rw; 
antperiod = stim+delay; 
CS_code = [1211 1212 12122 12112]; % Go - Nogo
Tra_code = [211 212 2122 2112]; % Go - Nogo 
% Tra_R1_5kHz = 211;      Tra_R1_10kHz = 212;     Cond_R1_5kHz = 1211;    Cond_R1_10kHz = 1212;

for i =1:size(datasource,2)
    event_type = datasource(i).xb.trial_type(:,1:2);
    event = datasource(i).eventtimes;
    num = length(datasource(i).xb.trial_type);
    temp_block = event(find(abs(diff(event_type(:,2))) > 500)+1,1);
    if isempty(datasource(i).eventtimes(find(event_type(1:num,2) == CS_code(1,3)),1)) == 0
        CS_type(i,1) = 1;
        bb = 1;
        for b = 1:2:size(temp_block,1)
            rule_change{i,1}(bb,1) = temp_block(b,1);
            bb = bb + 1; 
        end
    else
         CS_type(i,1) = 0; 
        % sessions with Task only 
        rule_change{i,1} = temp_block;
    end

    rule_change{i,1}(end+1,1) = max(event(:,1)) + 1000;
    rev_log(i,1) = rule_change{i,1}(1,1);
end

log_session = vertcat(1,find(diff(rev_log) ~= 0)+1); 
%%
thr = 50; % jitter threshold, in ms. 
clear event_session
for i = 1:size(log_session,1)
    event_session{i,1} = datasource(log_session(i)).xb.licks.locs';

    n_eventtimes = length(datasource(log_session(i)).eventtimes(:,1));
    n_trialtype = length(datasource(log_session(i)).xb.trial_type(:,1));

    if n_eventtimes - n_trialtype ~= 0
        disp(['mismatched session, #' num2str(i)])
        clear tnum behavlog intanlog
        tnum = min(n_eventtimes, n_trialtype);
        behavevent = diff(datasource(log_session(i)).xb.trial_type(1:tnum,1));
        intanevent = diff(datasource(log_session(i)).eventtimes(1:tnum,1)*1000);

        % jitter
        while true
            jitter = diff(horzcat(behavevent, intanevent),1,2);
            if isempty(find(jitter > thr, 1))
                break
            end

            n = length(behavevent);
            intanevent2 = [];
            skip = false;
            for j = 1:n
                if skip
                    skip = false;
                    continue;
                end
                if j < n && abs(behavevent(j) - intanevent(j)) > thr
                    combined = intanevent(j) + intanevent(j+1);
                    if abs(behavevent(j) - combined) < thr
                        intanevent2 = [intanevent2, combined];
                        skip = true;
                        continue;
                    end
                end
                intanevent2 = [intanevent2, intanevent(j)];
            end
            intanevent = intanevent2';
            tnum = length(intanevent) + 1;
            behavevent = diff(datasource(log_session(i)).xb.trial_type(1:tnum,1));
        end

        event_session{i,2}(1:tnum,1) = datasource(log_session(i)).eventtimes(1:tnum,1)*1000;
        event_session{i,2}(1:tnum,3) = datasource(log_session(i)).xb.trial_type(1:tnum,2);
        event_session{i,2}(1:tnum,2) = datasource(log_session(i)).xb.trial_type(1:tnum,1);

    else
        % from intan log file
        event_session{i,2}(:,1) = datasource(log_session(i)).eventtimes(:,1)*1000;
        event_session{i,2}(:,3) = datasource(log_session(i)).xb.trial_type(:,2);
        % from presentation log
        event_session{i,2}(:,2) = datasource(log_session(i)).xb.trial_type(:,1);
    end

    % Block change point
    BlockChangeTime{i,1}(:,1) = event_session{i,2}(find(abs(diff(event_session{i,2}(:,3))) > 100)+1, 1);
    BlockChangeTime{i,1}(:,2) = event_session{i,2}(find(abs(diff(event_session{i,2}(:,3))) > 100)+1, 2);

end




%%
%%% lick : ms
lickboutinterval = 2000;                itilick = 4000; 
lickpre = 2000;  
%%% perform: sec 
stim = 500;         delay = 1000;       rw = 2000;       pre = 1000;       
stimperiod = stim+delay+rw; 
antperiod = stim+delay; 

CS_code = [1211 1212 12122 12112]; % Go - Nogo
Tra_code = [211 212 2122 2112]; % Go - Nogo 
All_stim = [Tra_code CS_code]; 

smth = 11; performcri = 0.7; consec = 2; consec_tnum = 3; grace = 0; 

% Tra_R1_5kHz = 211;      Tra_R1_10kHz = 212;     Cond_R1_5kHz = 1211;    Cond_R1_10kHz = 1212;
% Tra_R2_5kHz = 2112;     Tra_R2_10kHz = 2122;    Cond_R2_5kHz = 12112;   Cond_R2_10kHz =12122;

clear behav 
for i = 1:size(event_session,1)

    %%% lickbout and spont_lick trial sorting 
    clear lick_ms tempcri lickbout_on lickbout_off lickbout_duration event_ms_mat event_ms
    lick_ms = event_session{i,1}(:,1);
    tempcri = find(diff(lick_ms) > lickboutinterval);       % finding out lick bout onset
    for k = 1:length(tempcri)-1
        lickbout_on(k,1) = lick_ms(tempcri(k)+1);              % lick bout onset
        lickbout_off(k,1) = lick_ms(tempcri(k+1));             % lick bout offset
    end

    lickbout_on = vertcat(lick_ms(1), lickbout_on, lick_ms(tempcri(end)+1));
    lickbout_off = vertcat(lick_ms(tempcri(1)), lickbout_off, lick_ms(end));
    lickbout_duration(:,1) = lickbout_off - lickbout_on;
    behav.lick{i,1}(:,1) = lickbout_on;
    behav.lick{i,1}(:,2) = lickbout_off;
    behav.lick{i,1}(:,3) = lickbout_duration;
    behav.lickms{i,1} = lick_ms; 

    % for t = 1:stimperiod+3000
    %     event_ms_mat(:,t) = lick_ms(:,1)-(lickpre+1)+t;
    % end
    % 
    % %%% spont. lick trials w/o other lickbout contamination 
    % spont = 1;
    % for s = 1:size(lickbout_on,1)
    %     if length(find(lickbout_on(s) < event_ms_mat & lickbout_off(s) > event_ms_mat )) == 0
    %         if lickbout_on(s) ~= lickbout_off(s)
    %             behav.spont{i,1}(spont,1) = lickbout_on(s);
    %             behav.spont{i,2}(spont,2) = lickbout_off(s);
    %             spont = spont + 1; 
    %         end
    %     end
    % end
    % clear lick_ms tempcri lickbout_on lickbout_off lickbout_duration event_ms_mat event_ms

    %%% Go/Nogo Sorting: Tra{R1go R1ngo R2go R2ngo} CS{R1go R1ngo R2go R2ngo}
    for j = 1:size(Tra_code,2)
        behav.Tra{i,j} = event_session{i,2}(find(event_session{i,2}(:,3)==Tra_code(:,j)),:);
        behav.CS{i,j} = event_session{i,2}(find(event_session{i,2}(:,3)==CS_code(:,j)),:);        
    end

    for j = [1 3]
        lick = 1; nolick = 1;
        for t = 1:size(behav.Tra{i,j},1)
            onset = behav.Tra{i,j}(t,1);
            if sum(find(lick_ms > onset+antperiod & lick_ms < onset+stimperiod)) > 0
                behav.Perform.Tra{i,j*2-1}{lick,1}(1,1) = onset;
                behav.Perform.Tra{i,j*2-1}{lick,1}(1,2) = behav.Tra{i,j}(t,2);
                behav.Perform.Tra{i,j*2-1}{lick,2} = lick_ms(find(lick_ms > onset-pre & lick_ms < onset+stimperiod)) - onset;
                behav.Perform.Tra{i,j*2-1}{lick,3} = 1;
                behav.Perform.Tra{i,j*2-1}{lick,4} = 1;
                lick = lick + 1;
                behav.Tra{i,j}(t,4) = 1; % Contingency
                behav.Tra{i,j}(t,5) = 1; % Lick
                behav.Tra{i,j}(t,6) = 1; % correct
                behav.Tra{i,j}(t,7) = 1; % Stage
            elseif sum(find(lick_ms > onset+antperiod & lick_ms < onset+stimperiod)) == 0
                behav.Perform.Tra{i,j*2}{nolick,1}(1,1) = onset;
                behav.Perform.Tra{i,j*2}{nolick,1}(1,2) = behav.Tra{i,j}(t,2);
                behav.Perform.Tra{i,j*2}{nolick,2} = lick_ms(find(lick_ms > onset-pre & lick_ms < onset+stimperiod)) - onset;
                behav.Perform.Tra{i,j*2}{nolick,3} = 0;
                behav.Perform.Tra{i,j*2}{nolick,4} =0;
                nolick = nolick + 1;
                behav.Tra{i,j}(t,4) = 1; % Contingency
                behav.Tra{i,j}(t,5) = 0; % Lick
                behav.Tra{i,j}(t,6) = 0; % correct
                behav.Tra{i,j}(t,7) = 1; % Stage 
            end
        end

        Antlick = 1; stay = 1;
        if isempty(behav.CS{i,j}) == 0
            for t = 1:size(behav.CS{i,j},1)
                onset = behav.CS{i,j}(t,1);
                if sum(find(lick_ms > onset & lick_ms < onset+antperiod)) > 0
                    behav.Perform.CS{i,j*2-1}{Antlick,1}(1,1) = onset;
                    behav.Perform.CS{i,j*2-1}{Antlick,1}(1,2) = behav.CS{i,j}(t,2);
                    behav.Perform.CS{i,j*2-1}{Antlick,2} = lick_ms(find(lick_ms > onset-pre & lick_ms < onset+stimperiod)) - onset;
                    behav.Perform.CS{i,j*2-1}{Antlick,3} = 1;
                    behav.Perform.CS{i,j*2-1}{Antlick,4} = 1;
                    Antlick = Antlick + 1;
                    behav.CS{i,j}(t,4) = 1; % Contingency
                    behav.CS{i,j}(t,5) = 1; % Lick
                    behav.CS{i,j}(t,6) = 1; % correct
                    behav.CS{i,j}(t,7) = 0; % Stage
                elseif sum(find(lick_ms > onset & lick_ms < onset+antperiod)) == 0
                    behav.Perform.CS{i,j*2}{stay,1}(1,1) = onset;
                    behav.Perform.CS{i,j*2}{stay,1}(1,2) = behav.CS{i,j}(t,2);
                    behav.Perform.CS{i,j*2}{stay,2} = lick_ms(find(lick_ms > onset-pre & lick_ms < onset+stimperiod)) - onset;
                    behav.Perform.CS{i,j*2}{stay,3} = 0;
                    behav.Perform.CS{i,j*2}{stay,4} = 0;
                    stay = stay + 1;
                    behav.CS{i,j}(t,4) = 1; % Contingency
                    behav.CS{i,j}(t,5) = 0; % Lick
                    behav.CS{i,j}(t,6) = 0; % correct
                    behav.CS{i,j}(t,7) = 0; % Stage
                end
            end
        end
    end

    for j = [2 4]
        lick = 1; nolick = 1;
        for t = 1:size(behav.Tra{i,j},1)
            onset = behav.Tra{i,j}(t,1);
            if sum(find(lick_ms > onset+antperiod & lick_ms < onset+stimperiod)) > 0
                behav.Perform.Tra{i,j*2-1}{lick,1}(1,1) = onset;
                behav.Perform.Tra{i,j*2-1}{lick,1}(1,2) = behav.Tra{i,j}(t,2);
                behav.Perform.Tra{i,j*2-1}{lick,2} = lick_ms(find(lick_ms > onset-pre & lick_ms < onset+stimperiod)) - onset;
                behav.Perform.Tra{i,j*2-1}{lick,3} = 0;
                behav.Perform.Tra{i,j*2-1}{lick,4} = 1;
                lick = lick + 1;
                behav.Tra{i,j}(t,4) = 0; % Contingency
                behav.Tra{i,j}(t,5) = 1; % Lick
                behav.Tra{i,j}(t,6) = 0; % correct 
                behav.Tra{i,j}(t,7) = 1; % Stage 
            elseif sum(find(lick_ms > onset+antperiod & lick_ms < onset+stimperiod)) == 0
                behav.Perform.Tra{i,j*2}{nolick,1}(1,1) = onset;
                behav.Perform.Tra{i,j*2}{nolick,1}(1,2) = behav.Tra{i,j}(t,2);
                behav.Perform.Tra{i,j*2}{nolick,2} = lick_ms(find(lick_ms > onset-pre & lick_ms < onset+stimperiod)) - onset;
                behav.Perform.Tra{i,j*2}{nolick,3} = 1;
                behav.Perform.Tra{i,j*2}{nolick,4} = 0;
                nolick = nolick + 1;
                behav.Tra{i,j}(t,4) = 0; % Contingency
                behav.Tra{i,j}(t,5) = 0; % Lick
                behav.Tra{i,j}(t,6) = 1; % correct 
                behav.Tra{i,j}(t,7) = 1; % Stage 
            end
        end

        Antlick = 1; stay = 1;
        if isempty(behav.CS{i,j}) == 0
            for t = 1:size(behav.CS{i,j},1)
                onset = behav.CS{i,j}(t,1);
                if sum(find(lick_ms > onset & lick_ms < onset+antperiod)) > 0
                    behav.Perform.CS{i,j*2-1}{Antlick,1}(1,1) = onset;
                    behav.Perform.CS{i,j*2-1}{Antlick,1}(1,2) = behav.CS{i,j}(t,2);
                    behav.Perform.CS{i,j*2-1}{Antlick,2} = lick_ms(find(lick_ms > onset-pre & lick_ms < onset+stimperiod)) - onset;
                    behav.Perform.CS{i,j*2-1}{Antlick,3} = 0;
                    behav.Perform.CS{i,j*2-1}{Antlick,4} = 1;
                    Antlick = Antlick + 1;
                    behav.CS{i,j}(t,4) = 0; % Contingency
                    behav.CS{i,j}(t,5) = 1; % Lick
                    behav.CS{i,j}(t,6) = 0; % correct
                    behav.CS{i,j}(t,7) = 0; % Stage
                elseif sum(find(lick_ms > onset & lick_ms < onset+antperiod)) == 0
                    behav.Perform.CS{i,j*2}{stay,1}(1,1) = onset;
                    behav.Perform.CS{i,j*2}{stay,1}(1,2) = behav.CS{i,j}(t,2);
                    behav.Perform.CS{i,j*2}{stay,2} = lick_ms(find(lick_ms > onset-pre & lick_ms < onset+stimperiod)) - onset;
                    behav.Perform.CS{i,j*2}{stay,3} = 1;
                    behav.Perform.CS{i,j*2}{stay,4} = 0;
                    stay = stay + 1;
                    behav.CS{i,j}(t,4) = 0; % Contingency
                    behav.CS{i,j}(t,5) = 0; % Lick
                    behav.CS{i,j}(t,6) = 1; % correct
                    behav.CS{i,j}(t,7) = 0; % Stage
                end
            end
        end
    end


    clear temp temp1 
    tnum = 0; temp2  = [];
    temp = sortrows(vertcat(vertcat(behav.CS{i,3:4}),vertcat(behav.Tra{i,3:4})),2);
    temp1 = movmean(temp(:,6), floor(smth)); 
    for t = 1:length(temp1)-floor(smth/2)+1
        if sum(temp1(t:t+consec,1) >= performcri) == consec_tnum
            tnum = tnum + 1;
            temp2(tnum,:) = temp(t,:);
            temp2(tnum,end) = t+floor(smth/2)-1;
        end
    end
    if isempty(temp2) == 1
        temp2(1,:) = temp(end,:);
        temp2(1,end) = 1;
    end
    behav.TR_time(i,:) = temp(temp2(1,7)+grace,:);
    behav.TTR(i,:) = temp2(1,7);

    behav.All{i,1} = sortrows(vertcat(vertcat(behav.CS{i,:}), vertcat(behav.Tra{i,:})),2);
    behav.BlockChange{i,1} = BlockChangeTime{i,1};

end

%%%
log_session = vertcat(log_session, length(datasource)+1); 
clear dataset
for s = 1:size(log_session,1)-1
    clear target
    target = log_session(s,1):log_session(s+1,1)-1; 
    for i = 1:size(target,2)
        dataset{target(i),1} = datasource(target(i)).spiketimes;
        dataset{target(i),2} = behav.lick{s,1};
        dataset{target(i),3} = behav.All{s,1}; 
        dataset{target(i),4} = behav.TR_time(s,:);
        dataset{target(i),5} = behav.TTR(s,:); 
        dataset{target(i),6} = s;
        dataset{target(i),7} = behav.BlockChange{s,1}; 
        dataset{target(i),8} = behav.lickms{s,1}; 
    end
end
%%%
clearvars -except Pool dataset BlockChangeTime All_stim Antlick antperiod CS_code CS_type delay event_session log_session rev_log rule_change roi roi2 stim stimperiod Tra_code
%%%
% if roi == 1
%     save PPCPool_Event_TTR_250520.mat -v7.3
% elseif roi == 2
%     save ACPool_Event_TTR_250520.mat -v7.3
% elseif roi == 3
%     save ICPool_Event_TTR_250520.mat -v7.3
% end
% disp('save temp dataset')

%% Post code, add to SUdata for just the stuff I need 

for n = 1:length(Pool)
    GLM_dataset{n,3} = dataset{n,3}(:,1)*1e-3;
%     GLM_dataset{n,4} = dataset{n,3}(:,4:7);
    GLM_dataset{n,7} = dataset{n,5}; % TTR
%     test= dataset{n,3};
%     test2 = GLM_dataset{n,5};
    tempdata = [];
    tempdata2 = [];
    for t = 1:length(dataset{n,3})
        ind = find(GLM_dataset{n,5}(:,1) == dataset{n,3}(t,2)); 
        tempdata(t,:) = GLM_dataset{n,5}(ind,:);
        tempdata2(t,:) = GLM_dataset{n,4}(ind,:);
    end
    GLM_dataset{n,5} = tempdata;
    GLM_dataset{n,4} = tempdata2;

end
    












