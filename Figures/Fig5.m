clc;
clear all; close all;

load('Pav_Delay_2023-05-15_Summary.mat')

session=1;
saveDirectory = '/Users/ali/Library/CloudStorage/OneDrive-UCSF/Berke Lab/workflow/2023/05/Timing'
figureExport=1;

k=0;
for subject=1:length(Subject)
    for session=1:length(Subject(subject).Session)
        k=k+1;
        dF1=Subject(subject).Session(session).dF1;
        dF2=Subject(subject).Session(session).dF2;
        foodEntry= Subject(subject).Session(session).foodEntry;
        ind=Subject(subject).Session(session).ind;
        t_axis=Subject(subject).Session(session).t_axis;
        rew = Subject(subject).Session(session).rew;
        tWin=[0 1];t_idx=find(t_axis>tWin(1) & t_axis<tWin(2));
        for cond=1:4
            for region=1:3
                M1=max(nanmean(dF1(ind==1,t_idx,region)));
                M2=max(nanmean(dF2(ind==1,t_idx,region)));
                cue_raw(k,:,cond,region) = nanmean(dF1(ind==cond & rew==1,:,region));
                Click_raw(k,:,cond,region) = nanmean(dF2(ind==cond & rew==1,:,region));
                if M1>1
                    cue(k,:,cond,region) = nanmean(dF1(ind==cond & rew==1,:,region))/M1;
                    Click(k,:,cond,region) = nanmean(dF2(ind==cond & rew==1,:,region))/M2;
                    valid{subject}(region,session)=1;
                    x=strfind(Subject(subject).Session(session).Address,filesep);
                else
                    cue(k,:,cond,region) = nan(size(t_axis));
                    Click(k,:,cond,region) = nan(size(t_axis));
                    valid{subject}(region,session)=0;
                end
            end
            food(k,:,cond) = nanmean(foodEntry(ind==cond & rew==1,:,1))*100;
        end
    end
end
%% Poor placement for IM-1556:: VS

cue(7:9,:,:,3)=nan(3,size(cue,2),size(cue,3));
Click(7:9,:,:,3)=nan(3,size(Click,2),size(Click,3));

%%
COLORS_region = [0 0 0;122 179 84;85 152 73;43 87 62]/255;
COLORS = [0 0 0;251 221 89;247 179 76;187 149 64]/255;
Region_Name ={'DLS','DMS','VS'};

M=nan(size(Click,3),size(Click,4));

for region=1:3
    for cond=1:4
        tWin=[0 .5];t_idx=find(t_axis>tWin(1) & t_axis<tWin(2));
        if cond==1
            dummy=squeeze(nanmean(Click));
            if isnan(max(dummy(t_idx,cond,region),[],1))
                M(cond,region)= nan;
            else
                M(cond,region)= 1;
            end
        else
            dummy=squeeze(nanmean(cue));
            M(cond,region)= max(dummy(t_idx,cond,region),[],1)/.75;
        end
    end
end




%% Exponential
h=figure('units','normalized','outerposition',[0 0 1 1]); %[left bottom width height]
subplot(121);hold on
for region=1:3
    f = @(a,b,c,x) b+a*exp(-x/c);
    [fitresult x y]= fit([0 1 3 12]',[squeeze((M(1:4,region)))],f,'StartPoint',[1 1 1]);
    xfit = linspace(0,50,500);
    yfit = feval(fitresult,xfit);
    regionTau(region)=fitresult.c;
    hh(region)=plot(xfit,yfit,'linewidth',2,'Color',COLORS_region(region+1,:))
    plot([1]',squeeze((M(2,region)))','o','MarkerEdgeColor',COLORS_region(region+1,:),'MarkerFaceColor',COLORS(2,:),'MarkerSize',10,'linewidth',3)
    plot([3]',squeeze((M(3,region)))','o','MarkerEdgeColor',COLORS_region(region+1,:),'MarkerFaceColor',COLORS(3,:),'MarkerSize',10,'linewidth',3)
    plot([12]',squeeze((M(4,region)))','o','MarkerEdgeColor',COLORS_region(region+1,:),'MarkerFaceColor',COLORS(4,:),'MarkerSize',10,'linewidth',3)
end
title(['f=b+a*exp(-x/\tau)'])
set(gca,'FontSize',16)
xlim([0 20]);ylim([0 1.1])
subplot(122);hold on;
for region=1:3
    f = @(a,b,c,x) a + b./(1+x/c);
    [fitresult x y]= fit([0 1 3 12]',[squeeze((M(1:4,region)))],f,'StartPoint',[1 1 1]);
    xfit = linspace(0,50,500);
    yfit = feval(fitresult,xfit);
    regionTau(region)=fitresult.c;
    hh(region)=plot(xfit,yfit,'linewidth',2,'Color',COLORS_region(region+1,:))
    plot([1]',squeeze((M(2,region)))','o','MarkerEdgeColor',COLORS_region(region+1,:),'MarkerFaceColor',COLORS(2,:),'MarkerSize',10,'linewidth',3)
    plot([3]',squeeze((M(3,region)))','o','MarkerEdgeColor',COLORS_region(region+1,:),'MarkerFaceColor',COLORS(3,:),'MarkerSize',10,'linewidth',3)
    plot([12]',squeeze((M(4,region)))','o','MarkerEdgeColor',COLORS_region(region+1,:),'MarkerFaceColor',COLORS(4,:),'MarkerSize',10,'linewidth',3)

end
title(['f=a + b./(1+x/c)'])
set(gca,'FontSize',16)
xlim([0 20]);ylim([0 1.1])

if figureExport
    set(h,'Units','inches');
    screenposition = get(h,'Position');
    set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',[screenposition(3:4)]);
    print('-dpdf','-vector',[saveDirectory filesep  'Average_Fit_75pCorrected_'  datestr(now, 'yyyy-mm-dd')])
end
print('-dpdf','-vector',['Average_Fit_75pCorrected_'  datestr(now, 'yyyy-mm-dd')])
