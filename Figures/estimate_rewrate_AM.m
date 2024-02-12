function [rewRate Corr tau] = estimate_rewrate_AM(ratexperience,boxts,tau)


% [tau,fval,exitflag,output,lambda,grad,hessian] = fmincon(@(x)dummy_estimate_rewards(x,ratexperience,boxts),50,[],[],[],[],0.5,1000);
%
% Fs = 250;
% endTime=max(max(boxts.event))+Fs*120;
% kernel = exp((-linspace(0,5*tau,Fs*tau)/tau));
% behavData = ratexperience.behavData;
% outcome = diff([0;behavData.food]);
%
% rewRate = zeros(1,endTime);
% rewIDX = find(outcome==1);
% for trial=1:length(rewIDX)
%     idx = boxts.event(5,rewIDX(trial));
%     rewRate(idx:idx+length(kernel)-1) = rewRate(idx:idx+length(kernel)-1)  + kernel;
% end
%
% validTrial = find(diff([0; behavData.trial]) == 1);
% CNIlatency = behavData.CNIlatency;
% CNIlatency = CNIlatency(validTrial);
% CNIlatency(CNIlatency==0)=nan;
%
% Corr = nancorr(log10(CNIlatency),rewRate(boxts.event(2,validTrial))');

pad=3500;

if tau==0
    
    for tau=1:200
        Fs = 250;
        endTime=max(max(boxts.event))+Fs*pad;
        kernel = exp((-linspace(0,5*tau,Fs*tau)/tau));
        behavData = ratexperience.behavData;
        outcome = diff([0;behavData.food]);
        
        rewRate = zeros(1,(endTime));
        rewIDX = find(outcome==1);
        for trial=1:length(rewIDX)
            idx = (boxts.event(7,rewIDX(trial)));
            if ~isnan(idx)
                rewRate(idx:idx+length(kernel)-1) = rewRate(idx:idx+length(kernel)-1)  + kernel;
            end
        end
        
        validTrial = find(diff([0; behavData.trial]) == 1);
        CNIlatency = behavData.CNIlatency;
        CNIlatency = CNIlatency(validTrial);
        CNIlatency(CNIlatency==0)=nan;
        
        Corr(tau) = nancorr(log10(CNIlatency),rewRate((boxts.event(2,validTrial)))');
    end
    
    [~,tau]=min(Corr);
    
else
    Fs = 250;
    endTime=max(max(boxts.event))+Fs*pad;
    kernel = exp((-linspace(0,5*tau,Fs*tau)/tau));
    behavData = ratexperience.behavData;
    outcome = diff([0;behavData.food]);
    
    t_axis = linspace(0,endTime/Fs,(endTime));
    
    rewRate = zeros(1,(endTime));
    rewIDX = find(outcome==1);
    for trial=1:length(rewIDX)
        idx = (boxts.event(7,rewIDX(trial)));
        if ~isnan(idx)
            rewRate(idx:idx+length(kernel)-1) = rewRate(idx:idx+length(kernel)-1)  + kernel;
        end
    end
    
    %figure,plot(t_axis/60,rewRate)
    
    validTrial = find(diff([0; behavData.trial]) == 1);
    CNIlatency = behavData.CNIlatency;
    CNIlatency = CNIlatency(validTrial);
    CNIlatency(CNIlatency==0)=nan;
    
    Corr = nancorr(log10(CNIlatency),rewRate((boxts.event(2,validTrial)))');
    
end
    
%     figure;
%     plot(rewRate)
%     for trial=1:length(rewIDX)
%         idx = (boxts.event(7,rewIDX(trial)));hold on;
%         line([idx idx],[0 max(rewRate)],'color','r');hold on;
%     end
%
%     function C=dummy_estimate_rewards(x,ratexperience,boxts)
%         Fs = 250;
%         endTime=max(max(boxts.event))+Fs*120;
%         tau=x(1);
%         kernel = exp((-linspace(0,5*tau,Fs*tau)/tau));
%         behavData = ratexperience.behavData;
%         outcome = diff([0;behavData.food]);
%
%         rewRate = zeros(1,endTime);
%         rewIDX = find(outcome==1);
%         for trial=1:length(rewIDX)
%             idx = boxts.event(5,rewIDX(trial));
%             rewRate(idx:idx+length(kernel)-1) = rewRate(idx:idx+length(kernel)-1)  + kernel;
%         end
%
%         validTrial = find(diff([0; behavData.trial]) == 1);
%         CNIlatency = behavData.CNIlatency;
%         CNIlatency = CNIlatency(validTrial);
%         CNIlatency(CNIlatency==0)=nan;
%
%         C = nancorr(log10(CNIlatency),rewRate(boxts.event(2,validTrial))');
%
%     end

end




