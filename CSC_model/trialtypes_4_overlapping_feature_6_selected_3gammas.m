% code modified from Elliot A. Ludvig's matlab code for the paper:
% Evaluating the TD model of classical conditioning 2012
%% Model Parameters
%numms = 6;          %Number of microstimuli per stimulus. 6, 20
numstimuli = 5;     %4 for trial type. 1 for overlapping
%alpha = 0.05;       %Step-size
alpha = 0.01;
%decay = 0.985;      %Memory Trace decay rate

lambda = 0.98; %0.95;      %Eligibility trace decay rate. lamda =0 for no egibility trace
stimulustime = 10;           %Conditioned Stimulus (CS) time
rewardt = stimulustime+31;                   %Reward Time
numtrials = 5000;                 %Last 3 trials are different probes    
numtimesteps = 60;
stimrep = 1;                     %Current code does not allow for mixed representations
%1 = CSC; 2 = Pres; 3 = MS
if stimrep == 1 %CSC             %Ensures that vectors are big enough to store data
    numms = numtimesteps-1;
end
% dt = 0.1;
gammas = [0.95 0.99 0.9999];
seeds  = 0:9;
for i = 1:length(gammas)
    for seed_i = 1:length(seeds)
        seed = seeds(seed_i);
        rng(seed);
        cue_types = randi(4,numtrials,1);
    
        
        gamma = gammas(i);
        % Initialize Data Vectors:
        x = zeros(1, (numms + 1)*numstimuli);       %Stimulus representation
        w = zeros(1, (numms + 1)*numstimuli);       %Weight vector
        e = zeros(1, (numms + 1)*numstimuli);       %Eligibility Traces
        delta = zeros (numtrials, numtimesteps);    %TD Errors
        value = zeros (numtrials, numtimesteps);    %Value functions
        rewards = zeros (numtrials,1);
        %
        % cue type: 1: 0%, 2: %25, 3: 75%, 4: click; 


        % Run Simuluation:
        for trial = 1:numtrials
            oldvalue = 0;               % Reset on every trial
            e = zeros(1, (numms + 1)*numstimuli);
            triallen = numtimesteps;

            rewardtime = rewardt;
            cue_idx = cue_types(trial);

            for timestep = 1:triallen
                reward = 0;
                if (((cue_idx==2) && (timestep == rewardtime) && (rand<0.25)) ...
                        || ((cue_idx==3) && (timestep == rewardtime) && (rand<0.75)) ...
                        || ((cue_idx==4) && (timestep == rewardtime)))
                    reward = 1;
                    rewards(trial) = 1;

                end
                % CSC Code:
                if stimrep == 1
                    x = zeros(1, (numms + 1)*numstimuli);

                    if timestep >= stimulustime 
                        if cue_idx<=3
                            x ((cue_idx-1)*numtimesteps + (timestep-stimulustime)+1) = 1; 
                            x (4*numtimesteps + (timestep-stimulustime)+1) = 1; % the last feature share by the 3 cues.            
                        end
                    end
                end
                % Value Calculation          
                value (trial, timestep) = dot(x,w);
                delta (trial, timestep) = reward + (gamma * value (trial, timestep)) - oldvalue; %TD Learning
                w = w + (alpha * delta (trial, timestep) * e);
                e = x + (gamma * lambda * e);
                oldvalue = dot(x,w); % Or oldvalue = value. Difference in which weights are used.
                %oldreward = reward;
            end

        end
        
        if ~exist('gamma3', 'dir')
           mkdir('gamma3')
        end
        data1 = './gamma3/seed_%d_gamma_%d_alpha_%0.2f_overlap_60.mat';
        file1 = sprintf(data1,seed,i,alpha);
        save(file1, 'cue_types','rewards', 'value','delta');
    end
end
