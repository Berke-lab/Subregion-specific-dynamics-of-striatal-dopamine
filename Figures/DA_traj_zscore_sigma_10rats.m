clc;
clear all; close all;

rats = {'IM1277' 'IM1278' 'IM1299' 'IM1300' 'IM1301' 'IM1358' 'IM1359' ...
        'IM1366' 'IM1367' 'IM1381' 'IM1382' };

for rat_idx = 1:11
    clearvars -except rat_idx rats

    rat_idx   
    dir1 = '../rat_data_Apr2022/';
    file1 = append(rats{rat_idx}, '.mat');
    Subject(1) = load([dir1, file1]);

    Fs=250;
    tWin =[-3 9];
    t_axis = linspace(tWin(1),tWin(2),1+(tWin(2)-tWin(1))*Fs);
    t_axis = t_axis(1:end-1);

    responseFood=nan(numel(t_axis),9,15);
    response=nan(3,9,15,numel(t_axis));
    Fs_session  = nan(1,15);
    zscore_sigma=nan(2,15,2); % chanel 1, 3
    event = 1; % for cue time    
    
    subject = 1;
    Session = Subject(subject).Session;
    for session=1:15            
            
            boxts = Session(session).boxts;
            Fs = Session(session).Fs;
            ratConditioning = Session(session).behavData;
            FoodLine = downsample(boxts.FoodLine,Session(session).phot.SamplingRate/Fs);
            outcome=ratConditioning.Food;
            trial_num = length(outcome); 

            ind = -1*ones(size(ratConditioning.Tone));
            ind(ratConditioning.Tone==0) = 1;
            ind(ratConditioning.Tone==ratConditioning.Cue_Contingency(ratConditioning.Cue_Contingency(:,2)==75)) = 2;
            ind(ratConditioning.Tone==ratConditioning.Cue_Contingency(ratConditioning.Cue_Contingency(:,2)==25)) = 3;
            ind(ratConditioning.Tone==ratConditioning.Cue_Contingency(ratConditioning.Cue_Contingency(:,2)==0)) = 4;
            
            for trial=1:trial_num               
                if ~isnan(boxts.event(event,trial)) && (boxts.event(event,trial)+tWin(1)*Fs>0) && (boxts.event(event,trial)+tWin(2)*Fs<length(FoodLine))
                    foodEntry(trial,:) = FoodLine(boxts.event(event,trial)+tWin(1)*Fs:boxts.event(event,trial)+tWin(2)*Fs-1);
                else
                    foodEntry(trial,:) = nan(1,1+Fs*(tWin(2)-tWin(1)));
                end                
            end
            responseFood(:,1,session,subject) = squeeze(nanmean(foodEntry(ind==4,:,1))); % 0%
            responseFood(:,2,session,subject) = squeeze(nanmean(foodEntry(ind==3,:,1))); % 25%
            responseFood(:,3,session,subject) = squeeze(nanmean(foodEntry(ind==2,:,1))); % 75%
            responseFood(:,4,session,subject) = squeeze(nanmean(foodEntry(ind==3 & outcome==1,:,1))); % 25% r=1
            responseFood(:,5,session,subject) = squeeze(nanmean(foodEntry(ind==2 & outcome==1,:,1))); % 75% r=1
            responseFood(:,6,session,subject) = squeeze(nanmean(foodEntry(ind==1 & outcome==1,:,1))); % click, r=1
            responseFood(:,7,session,subject) = squeeze(nanmean(foodEntry(ind==3 & outcome==0,:,1))); % 25% r=0
            responseFood(:,8,session,subject) = squeeze(nanmean(foodEntry(ind==2 & outcome==0,:,1))); % 75% r=0
            responseFood(:,9,session,subject) = squeeze(nanmean(foodEntry(:,:,1))); % all

            for j2=1:9
               response(3,j2,session,:)=squeeze(responseFood(:,j2,session));                    
            end

            sigs= {Session(session).Channel(1).data, Session(session).Channel(3).data};
            refs= {Session(session).Channel(2).data, Session(session).Channel(4).data};
            
            Fs_s(session) = Fs;
            for i=1:2 % for recordings
                clear sig ref dF1 dF1_z                
                %sig = sigs{i};
                %ref = refs{i};
                sig = medfilt1(sigs{i});
                ref = medfilt1(refs{i});

                start1 = 10;
                length_whosession = min(length(sig),length(ref))-start1;
                sig = sig(start1:length_whosession);
                ref = ref(start1:length_whosession);

                [dF_F1,~,~] = isosbestic_correction_Fs(sig,ref,'fs',Fs);
                          
                idx_whole = boxts.event(event,1)+tWin(1)*Fs:boxts.event(event,size(boxts.event,2))+tWin(2)*Fs;
                trialAvg_sig = mean(dF_F1(idx_whole));
                trialStd_sig = std(dF_F1(idx_whole));                
                
                for trial=1:size(boxts.event,2)                    
                    idx0 = boxts.event(1,trial)+tWin(1)*Fs;
                    dF1(trial,:) = dF_F1(idx0:idx0+length(t_axis)-1);                        
                end

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
    data1 = append('response_',rats{rat_idx},'_1_15_','zc_per_session','.mat');
    save(data1, 'Fs_session','zscore_sigma','t_axis','response');
end
