clc;
clear all; close all;

rats = {'IM1413' 'IM1414' 'IM1415'}; 

for rat_idx = 1:3
    clearvars -except rat_idx rats

    rat_idx      
    dir1 = '../rat_data_Apr2022/';
    file1 = append(rats{rat_idx}, '.mat');
    load([dir1, file1]);

    Fs=250;
    tWin =[-3 9]; 
    t_axis = linspace(tWin(1),tWin(2),1+(tWin(2)-tWin(1))*Fs);
    t_axis = t_axis(1:end-1);
    responseFood=nan(numel(t_axis),9,15);
    response=nan(4,9,15,numel(t_axis));
    Fs_session  = nan(3,15);
    zscore_sigma=nan(3,15,2); %
    event = 1; % for cue time    
        
    for session=1:15 %length(Session)
                   
            boxts = Session(session).boxts;
            box =   Session(session).box;
            boxts.event = round(Fs*boxts.event / (box.SamplingRate));
            FoodLine = downsample(boxts.FoodLine,box.SamplingRate/Fs);
            ratConditioning=Session(session).ratConditioning;
            outcome=ratConditioning.Food;
            trial_num = length(outcome); 
    
            for trial=1:trial_num               
                if ~isnan(boxts.event(event,trial)) && (boxts.event(event,trial)+tWin(1)*Fs>0) && (boxts.event(event,trial)+tWin(2)*Fs<length(FoodLine))
                    foodEntry(trial,:) = FoodLine(boxts.event(event,trial)+tWin(1)*Fs:boxts.event(event,trial)+tWin(2)*Fs-1);
                else
                    foodEntry(trial,:) = nan(1,1+Fs*(tWin(2)-tWin(1)));
                end                
            end
                    
            ind = -1*ones(size(ratConditioning.Tone));
            ind(ratConditioning.Tone==0) = 1;
            ind(ratConditioning.Tone==ratConditioning.Cue_Contingency(ratConditioning.Cue_Contingency(:,2)==75)) = 2;
            ind(ratConditioning.Tone==ratConditioning.Cue_Contingency(ratConditioning.Cue_Contingency(:,2)==25)) = 3;
            ind(ratConditioning.Tone==ratConditioning.Cue_Contingency(ratConditioning.Cue_Contingency(:,2)==0))  = 4;
            % 1: click; 2: 75%; 3: 25%; 4: 0%
            
    
            responseFood(:,1,session) = squeeze(nanmean(foodEntry(ind==4,:))); % 0%
            responseFood(:,2,session) = squeeze(nanmean(foodEntry(ind==3,:))); % 25%
            responseFood(:,3,session) = squeeze(nanmean(foodEntry(ind==2,:))); % 75%
            responseFood(:,4,session) = squeeze(nanmean(foodEntry(ind==3 & outcome==1,:))); % 25% r=1
            responseFood(:,5,session) = squeeze(nanmean(foodEntry(ind==2 & outcome==1,:))); % 75% r=1
            responseFood(:,6,session) = squeeze(nanmean(foodEntry(ind==1 & outcome==1,:))); % click, r=1
            responseFood(:,7,session) = squeeze(nanmean(foodEntry(ind==3 & outcome==0,:))); % 25% r=0
            responseFood(:,8,session) = squeeze(nanmean(foodEntry(ind==2 & outcome==0,:))); % 75% r=0
            responseFood(:,9,session) = squeeze(nanmean(foodEntry(:,:))); % all
            
            for j2=1:9
               response(4,j2,session,:)=responseFood(:,j2,session);                    
            end
    
            cues= {Session(session).DLS.cue_ts, Session(session).DMS.cue_ts, Session(session).VS.cue_ts};
            tss = {Session(session).DLS.signal.ts, Session(session).DMS.signal.ts, Session(session).VS.signal.ts};
            sigs= {Session(session).DLS.signal.data,Session(session).DMS.signal.data,Session(session).VS.signal.data};
            refs= {Session(session).DLS.ref.data, Session(session).DMS.ref.data, Session(session).VS.ref.data};
            
            for i=1:3 % for DLS, DMS, VS
                clear cue ts sig ref dF1 dF1_z
                cue = cues{i};
                ts  = tss{i};
                %sig = sigs{i};
                %ref = refs{i};
                sig = medfilt1(sigs{i});
                ref = medfilt1(refs{i});
    
                if ((i==3) && (rat_idx == 1) && (session == 12) || ((i==3) && (rat_idx == 3) && (session == 9)))
                    continue
                else
                    dt = ts(2)-ts(1);
                    if abs(1/dt -25)<0.5
                        Fs_orig = 25;
                        length_orig = 300;
                        p = 10;
                        q = 1;
                    elseif abs(1/dt -33.3)<0.5
                        Fs_orig = 33.3;
                        length_orig = 400;
                        p = 15;
                        q = 2;
                    elseif abs(1/dt -250)<0.5
                        Fs_orig = 250;
                        length_orig = 3000;
                        p = 1;
                        q = 1;
                    else
                        print('something wrong');
                    end
                    Fs_session(i, session)=Fs_orig; 
    
                    if ((p==10) || (p==15))
                        % shift the data righhand by two
                        if size(sig,1)==1
                            sig=[NaN NaN sig];
                            ref=[NaN NaN ref];
                        else
                            sig=[NaN; NaN; sig];
                            ref=[NaN; NaN; ref];
                        end
                    end
                    
                    % seesion 5 sig and ref length is not the same
                    start1 = 10;
                    length_whosession = min(length(sig),length(ref))-start1;
                    sig = sig(start1:length_whosession);
                    ref = ref(start1:length_whosession);
                    ts  = ts(start1:length_whosession);
                    [blue,~,~] = isosbestic_correction_Fs(sig,ref,'fs',Fs_orig);
                          
                    idx_whole = find(ts>cue(1)+tWin(1) & ts<cue(trial_num)+tWin(2));
                    trialAvg_sig = mean(blue(idx_whole));
                    trialStd_sig = std(blue(idx_whole));
                    
                    for trial=1:trial_num
                        idx = find(ts>=cue(trial)+tWin(1) & ts<cue(trial)+tWin(1)+0.1);
                        idx0 = idx(1);
                        test1 = blue(idx0:idx0+length_orig-1);
                        dF1(trial,:) = resample(test1,p,q);
                    end
    
                    % z score
                    dF1_z = (dF1-trialAvg_sig)/trialStd_sig;
                    normFactor=max(mean(dF1_z(ind==1 & outcome==1,:)));
                            
                    % 0: 0%; 1: 25%; 2: 75%; 3: 25% r= 1; 4: 75% r = 1; 5, click; 
                    % 6: 25% r=0; 7: 75% r=0; 8: all
                    response(i,1,session,:) = mean(dF1_z(ind==4,:))/normFactor; % 0%
                    response(i,2,session,:) = mean(dF1_z(ind==3,:))/normFactor; % 25%
                    response(i,3,session,:) = mean(dF1_z(ind==2,:))/normFactor; % 75%
                    response(i,4,session,:) = mean(dF1_z(ind==3 & outcome==1,:))/normFactor; % 25% r=1
                    response(i,5,session,:) = mean(dF1_z(ind==2 & outcome==1,:))/normFactor; % 75% r=1
                    response(i,6,session,:) = mean(dF1_z(ind==1 & outcome==1,:))/normFactor; % click, r=1
                    response(i,7,session,:) = mean(dF1_z(ind==3 & outcome==0,:))/normFactor; % 25% r=0
                    response(i,8,session,:) = mean(dF1_z(ind==2 & outcome==0,:))/normFactor; % 75% r=0
                    response(i,9,session,:) = mean(dF1_z(:,:,1))/normFactor; % all  
    
                    zscore_sigma(i,session,1) = normFactor;
                    zscore_sigma(i,session,2) = trialStd_sig;
                  
                end     
                
            end
    end    
    data1 = append('response_',rats{rat_idx},'_1_15_','zc_per_session','.mat');
    save(data1, 'Fs_session','zscore_sigma','t_axis','response');
end   