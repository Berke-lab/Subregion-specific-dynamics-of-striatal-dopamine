clc;
clear all; close all;


load ('RL_dLight_pruned_triRegion_2022-04-25.mat')


k=0;
for session=1:length(Session)
    
    clear response baseline minresponse trace traceNext
    
    
    k=k+1;
    sig=Session(session).sig;
    
    Location = Session(session).Location';
    
    boxts = Session(session).boxts;
    Fs = Session(session).Fs;
    behavData = Session(session).behavData;
    validTrial = find(diff([0; behavData.trial]) == 1);
    ratexperience = behavData.ratexperience;
    
    dF_F1(1:500)=dF_F1(501:1000);
    [b,a] = butter(3,30/(Fs/2));
    dF_F1 = filter(b,a,dF_F1);
    dF_F1 = medfilt1(dF_F1,25);
       
    [RewRate Corr tau] = estimate_rewrate_AM(ratexperience,boxts,55);
    RewRate = RewRate(1:length(dF_F1));
    rewRateTrial = RewRate(boxts.event(2,validTrial));
    
    
    outcome = diff([0;behavData.food]);
    CNIlatency = behavData.CNIlatency;
    
    outcome = outcome(validTrial);
    CNIlatency = CNIlatency(validTrial);
    
    clear preSidein_Value DA dF
    event=5;
    tWin =[-5 5];
    
    for trial=1:length(outcome)
        if ~isnan(boxts.event(event,validTrial(trial)))
            dF(trial,:) = dF_F1(boxts.event(event,validTrial(trial))+tWin(1)*Fs:boxts.event(event,validTrial(trial))+tWin(2)*Fs);
        else
            dF(trial,:) = nan(1,1+Fs*(tWin(2)-tWin(1)));
        end
    end
    
    t_axis=linspace(tWin(1),tWin(2),(tWin(2)-tWin(1))*Fs+1);
    

    peak=max(dF(:,t_axis>0 & t_axis<1),[],2);
    valley=min(dF(:,t_axis>0 & t_axis<1),[],2);

    jj=0; alphaC=0:0.01:1;
    for aa=alphaC
        jj=jj+1;
        [RPE V]= AC_return(ratexperience, aa);
        CCp(jj,session)=corr(peak(outcome==1),RPE(outcome==1)');
        CCn(jj,session)=corr(valley(outcome==0),RPE(outcome==0)');
    end
    [M I]=max(CCp(:,session));
    alphaOptimal_p(session)=alphaC(I);
    [M I]=max(CCn(:,session));
    alphaOptimal_n(session)=alphaC(I);

    [RPE V]= AC_return(ratexperience, alphaOptimal_p(session));

    [N,edges,bin] = histcounts((V),quantile((V),[0 1/3 2/3 1]));

    dF1=dF;
    peak=max(nanmean(dF1(bin==1 & outcome'==1,t_axis>0 & t_axis<1)));
    dF1=dF1/peak;

    for i=1:3
        sideInResponse(:,i,1,k)= nanmean(dF1(bin==i & outcome'==1,:));
        sideInResponse(:,i,2,k)= nanmean(dF1(bin==i & outcome'==0,:));
        sideInResponse(:,i,3,k)= nanmean(dF1(bin==i,:));
    end
    
    
    for t=1:size(dF,2)
        [b1(:,t,k,event),bint1]= regress(dF(outcome==1,t),[ones(sum(outcome==1),1)';rewRateTrial(outcome==1)]',0.01);
        [b2(:,t,k,event),bint2]= regress(dF(outcome==0,t),[ones(sum(outcome==0),1)';rewRateTrial(outcome==0)]',0.01);
        [b3(:,t,k,event),bint3]= regress(dF(:,t),[ones(numel(rewRateTrial),1)';rewRateTrial;outcome']',0.05);
        p1(2,t,k,event,1) = double((bint1(2,1)<0 & bint1(2,2)<0));
        p2(2,t,k,event,1) = double((bint2(2,1)<0 & bint2(2,2)<0));
        p3(2,t,k,event,1) = double((bint3(2,1)>0 & bint3(2,2)>0));
    end
    
end

for i=1:length(Session)
    if strcmp(Session(i).Location,'DLS'),Region(i)=1;end
    if strcmp(Session(i).Location,'DMS'),Region(i)=2;end
    if strcmp(Session(i).Location,'VS'),Region(i)=3;end
end

DLS = find(Region==1);DMS = find(Region==2);NAc = find(Region==3);



region{1}=DLS;region{2}=DMS;region{3}=NAc;

h=figure('units','normalized','outerposition',[0.2 0.2 .5 .7]); %[left bottom width height]
t_axis=linspace(tWin(1),tWin(2),(tWin(2)-tWin(1))*Fs+1);
COLORS = [100 0 0;175 0 0;255 0 0;0 0 100;0 0 175;0 0 255;175 175 175;90 90 90;0 0 0]/255;
medFiltSize=3;
for i=1:3
    subplot(2,3,i)
    idx=find(t_axis>0);
    bar(t_axis(idx),medfilt1(100*sum(p3(2,idx,region{i},5),3)/numel(region{i}),medFiltSize),1,'FaceColor','k','EdgeColor','none');hold on
    bar(t_axis(idx),medfilt1(100*sum(p1(2,idx,region{i},5),3)/numel(region{i}),medFiltSize),1,'FaceColor','r','EdgeColor','none');hold on
    bar(t_axis(idx),medfilt1(100*sum(p2(2,idx,region{i},5),3)/numel(region{i}),medFiltSize),1,'FaceColor','b','EdgeColor','none')
    idx=find(t_axis<=0);
    bar(t_axis(idx),medfilt1(100*sum(p3(2,idx,region{i},5),3)/numel(region{i}),medFiltSize),1,'FaceColor','k','EdgeColor','none');hold on
    bar(t_axis(idx),medfilt1(100*sum(p1(2,idx,region{i},5),3)/numel(region{i}),medFiltSize),1,'FaceColor','r','EdgeColor','none');hold on
    bar(t_axis(idx),medfilt1(100*sum(p2(2,idx,region{i},5),3)/numel(region{i}),medFiltSize),1,'FaceColor','b','EdgeColor','none')
    
    for j=1:length(t_axis)
        if sum(p3(2,j,region{i},5),3)>2
            plot(t_axis(j),100,'k.')
        end
        if sum(p2(2,j,region{i},5),3)>2
            plot(t_axis(j),110,'b.')
        end
        if sum(p1(2,j,region{i},5),3)>2
            plot(t_axis(j),120,'r.')
        end
    end
    
    
    xlim([-1.5 1.5])
    
    COLORS_pre=[0 0 0; 70 70 70;150 150 150]/255;
    subplot(2,3,i+3)
    for j=1:3
        idx =find(t_axis>=0);
        XX = permute(sideInResponse,[4 1 2 3]);
        plotSEM(t_axis(idx),(XX(region{i},idx,j,1)),COLORS(j,:),[.8 .8 .8]);hold on
        plotSEM(t_axis(idx),(XX(region{i},idx,j,2)),COLORS(j+3,:),[.8 .8 .8]);hold on

        idx =find(t_axis<=0);
        plotSEM(t_axis(idx),(XX(region{i},idx,j,3)),COLORS_pre(j,:),[.8 .8 .8]);hold on

    end
    xlim([-1.5 1.5])
    line([0 0],[-0.2 1],'linestyle','--','color',[.1 .1 .1])
    xlabel('Side-in (s)');ylabel('dF/F')
    set(gca,'FontSize',16)
end
set(h,'color','w')
set(h,'Units','inches');
screenposition = get(h,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',[screenposition(3:4)]);
print('-dpdf','-vector',['sideIn_average2_' datestr(now, 'yyyy-mm-dd')])

