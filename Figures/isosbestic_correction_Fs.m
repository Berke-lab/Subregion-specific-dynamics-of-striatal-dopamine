function [dF_F,ref_fitted,slope] = isosbestic_correction(signal,ref,varargin)

type = 'linear';
Fs = 250;
win=2;
filterorder=5;

for i = 1:2:length(varargin)
    switch lower(varargin{i})
        case 'type'
            type = varargin{i+1};
        case 'fs'
            Fs= varargin{i+1};
        case 'win'
            win = varargin{i+1};
        otherwise
            disp('invalid optional argument passed to plots_unitStability');
    end
end

% iter=1;
%
% [x,fval,exitflag,output,lambda,grad,hessian] = fmincon(@(x)fit_isosbestic(x,signal,ref),rand(1,2),[],[],[],[],[-10 -10],[10 10]);

switch type
    
    case 'zscore'
        %         [b,a] = butter(filterorder,.1/(Fs/2),'high');
        %         z = filter(b,a,signal);
        
        zSignal = zscore(signal);
        zRef = zscore(ref);
        
        mdl = fitlm(zRef,zSignal,'linear');
        x = [mdl.Coefficients.Estimate];
        
        ref_fitted = x(2)*zRef + x(1);
        dF_F = 100*(zSignal - ref_fitted);
        
        slope= x;
        
    case 'ratio'
        
        dF_F = signal./ref;
        ref_fitted=[];
        slope=[];
        
    case 'smooth'
              
        ref_fitted = filter(ones(1,bw)./bw,1,signal);
        %S=gpuArray(signal);
        S=signal;
        ref_fitted = medfilt1(S,bw);
        dF_F = 100*(signal - ref_fitted)./ ref_fitted;
        slope = nan;
        
        %tic;A=smooth(S,bw);t1=toc;
        %tic;A=smooth(signal,bw);t2=toc;
        
    case 'linear'
        if abs(Fs-33.3)<0.5
            win=3;
            bw = 100;
        else
            bw = Fs*win;
        end
        start = round(Fs*250); % starting from 250 sec

        ref = medfilt1(ref,bw)';
        mdl = fitlm(ref(start:end),signal(start:end),'linear');
        x = [mdl.Coefficients.Estimate];
        
        ref_fitted = x(2)*ref + x(1); ref_fitted = medfilt1(ref_fitted,bw)';
        dF_F = 100*(signal - ref_fitted)./ ref_fitted;
        
        slope= x;
        
    case 'exponential'
        
        x_axis = linspace(0,length(signal)/Fs,length(signal))';
        
        est = fminsearch(@(x) fit_isosbestic(x,signal,x_axis,'type','exponential'),rand(3,1));
        Aest = est(1);Best = est(2);Cest = est(3);
        fitted_470= Aest + Best.*exp(-Cest*x_axis);
        dBleachSignal = signal - fitted_470;
        
        
        est = fminsearch(@(x) fit_isosbestic(x,ref,x_axis,'type','exponential'),rand(3,1));
        Aest = est(1);Best = est(2);Cest = est(3);
        fitted_405= Aest + Best.*exp(-Cest*x_axis);
        dBleachRef= ref - fitted_405;
        
        
        dBleachSignal = dBleachSignal + 2*(round(abs(min(dBleachSignal))));
        dBleachRef = dBleachRef + 2*(round(abs(min(dBleachRef))));
        
        
        dF_F = (dBleachSignal - dBleachRef )./dBleachRef;
        
        figure,
        subplot(211),plot(dBleachSignal);hold on;plot(signal)
        subplot(212),plot(dBleachRef);hold on;plot(ref)
        %%
        
        
        bleach_subtracted_470 = (signal-fitted_470);
        bleach_subtracted_405 = (ref-fitted_405);
        out_470 = bleach_subtracted_470-bleach_subtracted_405;
        out_405 = bleach_subtracted_405;
        scaling = nan;
        offset = nan;
        
end


end


