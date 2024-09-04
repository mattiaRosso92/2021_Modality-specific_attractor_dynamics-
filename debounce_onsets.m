%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mattia Rosso - IPEM Institute for Systematic Musicology - 06/05/2020 
% Project ID: 2019005
% 
% This script serves to import the raw onsets from the 'Drifting
% Metronomes' experiment, check that no mistakes occurred during data
% acquisition and correct eventual false positives by 'debouncing' the
% signal offline.
%
% The format of raw .csv files is the following: 
% timestamp(ms) , trialID , subjectID(1,2)/metronomeID(3,4) , amplitude
%
% The first timestamp of metronome1 represents t=0, and coincides with t=0 
% of the eeg dataset.
%
% The processed data are saved in a .m file, to be imported in integrated 
% workspace for the analysis pipeline.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear 
close all
clc

% Sample size
nsample = 14;

% Design settings
nconds      = 4;
ntrials     = 10;
cond_labels = {'Visual Coupling' , 'Visual Control' , 'Audio Coupling' , 'Audio Control'}; 

% Settings for matching eeg time series
srate = 1000; %as eeg


%% Initialize stuff

nclicks1   = 65; %metronomes clicks per trial
nclicks2   = 64;
bpm1       = 100; %metronomes speed
bpm2       = nclicks2/nclicks1 * 100;
hz1        = bpm1/60; %metronomes frequency
hz2        = bpm2/60;
ioi_ref1   = 1/hz1 * 1000; %reference for expected inter-onset intervals 
ioi_ref2   = 1/hz2 * 1000;
time_trial = ioi_ref1 * nclicks1;
time_block = time_trial * ntrials;
tolerance  = 5; %tolerance of jittering in ms

%initialize timestamps vectors for the whole sample
%pre-allocate cells because we don't know the size a priori
onset_metro1 = {[]};  
onset_metro2 = {[]};   
onset_sub1   = {[]};    
onset_sub2   = {[]}; 
false_positive = 400; %threshold in milliseconds to detect extra-logging

%initialize inter-onset intervals vectors
ioi_metro1 = {[]}; 
ioi_metro2 = {[]}; 
ioi_sub1   = {[]}; 
ioi_sub2   = {[]}; 




%% Import and process

for pairi = 1:nsample %loop over pairs, or select one pair to visualize

    %import raw data grouped for conditions
    raw = cell(1,nconds);
    
    for condi = 1:nconds % loop over conditions
            
        raw{condi} = csvread(['Pair' num2str(pairi) '_Condition' num2str(condi) '.csv']);
         

        %assign timestamps; codes: 1=sub1 2=sub2 3=metro1 4=metro2    
        onset_metro1{condi} = raw{condi}( find(raw{condi}(:,3) == 3) , 1 ); 
        %align to first click of metronome1: metronomes takes precedence
        time_zero = onset_metro1{condi}(1,1);
        onset_metro1{condi}(:,1) = onset_metro1{condi} - time_zero; 
        %transpose vector to match EEG data format
        onset_metro1{condi} = onset_metro1{condi}';

        %repeat for metronome2
        onset_metro2{condi} = raw{condi}( find(raw{condi}(:,3) == 4) , 1 );
        onset_metro2{condi} = onset_metro2{condi}(:,1) - time_zero;
        onset_metro2{condi} = onset_metro2{condi}';

        %now for the subjects
        onset_sub1{condi} = raw{condi}( find(raw{condi}(:,3) == 1) , 1 );
        onset_sub1{condi}(:,1) = onset_sub1{condi}(:,1) - time_zero;
        onset_sub1{condi} = onset_sub1{condi}';        

        onset_sub2{condi} = raw{condi}( find(raw{condi}(:,3) == 2) , 1 ); 
        onset_sub2{condi}(:,1) = onset_sub2{condi}(:,1) - time_zero;
        onset_sub2{condi} = onset_sub2{condi}';


        %Detect false positives (i.e. extra-inputs from sensors)      
        for i = 2:length(onset_sub1{condi})

            if onset_sub1{condi}(i)-onset_sub1{condi}(i-1) < false_positive         
                onset_sub1{condi}(i) = nan;       
            end

        end
        %and remove them
        onset_sub1{condi}(find(isnan(onset_sub1{condi}))) = [];

        %Repeat for subject2
        for i = 2:length(onset_sub2{condi})

            if onset_sub2{condi}(i)-onset_sub2{condi}(i-1) < false_positive         
                onset_sub2{condi}(i) = nan;       
            end

        end
        onset_sub2{condi}(find(isnan(onset_sub2{condi}))) = [];


        %Check there are not anomalies in metronomes data logging (e.g. very 
        %few times, a metronome might be logged twice       
        [~,anomaly1] = find(ioi_metro1{condi} > ioi_ref1 + tolerance | ioi_metro1{condi} < ioi_ref1 - tolerance);
        [~,anomaly2] = find(ioi_metro2{condi} > ioi_ref2 + tolerance | ioi_metro2{condi} < ioi_ref2 - tolerance);
        if ~isempty(anomaly1) 
            ioi_metro1{condi}(anomaly1) = [];
            onset_metro1{condi}(anomaly1+1) = [];
            if length(onset_metro1{condi}) ~= nclicks1*ntrials %double-check
                sprintf(['Anomaly not corrected! check metronome1 in condition ' num2str(condi) ])
                pause; 
            end
        end
        if ~isempty(anomaly2)
            ioi_metro2{condi}(anomaly2) = [];
            onset_metro2{condi}(anomaly2+1) = [];
            if length(onset_metro2{condi}) ~= nclicks2*ntrials %double-check
                sprintf(['Anomaly not corrected! check metronome2 in condition ' num2str(condi) ])
                pause;
            end
        end

        
    end  

end



%%














