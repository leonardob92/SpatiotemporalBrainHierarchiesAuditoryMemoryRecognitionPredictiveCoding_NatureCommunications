%%

%% SPATIOTEMPORAL BRAIN HIERARCHIES OF AUDITORY MEMORY RECOGNITION AND PREDICTIVE CODING (DATA COLLECTION: AUDITORY PATTERN RECOGNITION 2020) - LEONARDO BONETTI - NATURE COMMUNICATIONS

% Aarhus University, University of Oxford, University of Bologna, Massachusetts Institute of Technology (MIT)

%%

%%

% The code provided below shows the full pipeline connected to the paper.
% In some instances, when relevant, it shows both the code for the original 
% submission and the code developed during the revision process.
% When the revision made completely obsolete the original code, such code
% has been removed from this script (but it is still avaialable in another
% script within this repository for full disclosure).
% Please, note that in some instances the source data files for the figures
% was saved manually.

%%

%% *** START UP FUNCTIONS.. (LBPD_startup_D) ***

%starting up some functions for LBPD toolbox.

%starting up some of the functions that I wrote for the following analysis
pathl = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD'; %path to stored functions (THIS MUST BECOME A FOLDER ON YOUR OWN COMPUTER)
addpath(pathl);
LBPD_startup_D(pathl);

%%

%% *** PREPROCESSING ***

%%% OBS!! the preprocessing was computed for (nearly) all the data
%%% collected and not only for the data that was actually analyzed and
%%% reported in this paper.

%%

%% Maxfilter

%OBS! before running maxfilter you need to close matlab, open the terminal and write: 'use anaconda', then open matlab and run maxfilter script

maxfilter_path = '/neuro/bin/util/maxfilter';
project = 'MINDLAB2020_MEG-AuditoryPatternRecognition';
maxDir = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2'; %output path

path = '/raw/sorted/MINDLAB2020_MEG-AuditoryPatternRecognition'; %path with all the subjects folders
jj = dir([path '/0*']); %list all the folders starting with '0' in order to avoid hidden files
for ii = 90:length(jj) %over subjects (STARTING FROM 13, SINCE THE PREVIOUS ONES WERE THE PILOTS (IF THEY DID ALSO THE PROPER EXPERIMENT, THE MAXFILTER COMPUTATION WILL BE DONE IN THE NEXT  SECTION))
    cart = [jj(ii).folder '/' jj(ii).name]; %create a path that combines the folder name and the file name
    pnana = dir([cart '/2*']); %search for folders starting with '2'
    for pp = 1:length(pnana) %loop to explore ad analyze all the folders inside the path above
        cart2 = [pnana(1).folder '/' pnana(pp).name];
        pr = dir([cart2 '/ME*']); %looks for meg folder
        if ~isempty(pr) %if pr is not empty, proceed with subfolders inside the meg path
            pnunu = dir([pr(1).folder '/' pr(1).name '/00*']);
            %if length(pnunu) > 1
            %warning(['subj ' num2str(ii) ' has files nuber = ' num2str(length(pnunu))]) %show a warning message if any subj has more thatn 1 meg sub-folder
            %end
            for dd = 1:length(pnunu)
                if strcmp(pnunu(dd).name(5:6),'re') || strcmp(pnunu(dd).name(5:6),'sa') || strcmp(pnunu(dd).name(5:6),'vi') || strcmp(pnunu(dd).name(5:6),'pd') %checks whether characters 5 to 6 are equal to 're', 'vi, 'sa' or 'pd'; the loop continues if this is true (1) and it stops if this is false (0)
                    %idx2 = strfind(pnunu(1).name,'Mus'); % search for musmelo folder in order to avoid other projects
                    %if ~isempty(idx2)
                    fpath = dir([pnunu(1).folder '/' pnunu(dd).name '/files/*.fif']); % looks for .fif file
                    rawName = ([fpath.folder '/' fpath.name]); %assigns the final path of the .fif file to the rawName path used in the maxfilter command
                    maxfName = ['SUBJ' jj(ii).name '_' fpath.name(1:end-4)]; %define the output name of the maxfilter processing
                    %movement compensation
                    cmd = ['submit_to_cluster -q maxfilter.q -n 4 -p ' ,project, ' "',maxfilter_path,' -f ',[rawName],' -o ' [maxDir '/' maxfName '_tsssdsm.fif'] ' -st 4 -corr 0.98 -movecomp -ds 4 ',' -format float -v | tee ' [maxDir '/log_files/' maxfName '_tsssdsm.log"']];
                    %no movement compensation (to be used if HPI coils did not work properly)
%                     cmd = ['submit_to_cluster -q maxfilter.q -n 4 -p ' ,project, ' "',maxfilter_path,' -f ',[rawName],' -o ' [maxDir '/' maxfName '_tsssdsm.fif'] ' -st 4 -corr 0.98 -ds 4 ',' -format float -v | tee ' [maxDir '/log_files/' maxfName '_tsssdsm.log"']];
                    system(cmd);
                end
            end
        end
    end
end

%% Maxfilter - Subjects who did the experiment (version 2.0) after the first pilot (version 1.0)

%%% DONE SUBJECTS 2 10 11 7 %%%

%OBS! before running maxfilter you need to close matlab, open the terminal and write: 'use anaconda', then open matlab and run maxfilter script

%settings
maxfilter_path = '/neuro/bin/util/maxfilter';
project = 'MINDLAB2020_MEG-AuditoryPatternRecognition';
maxDir = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2'; %output path

clear rawName maxfName
%path to raw files for these particular subjects
pathraw{1} = '/raw/sorted/MINDLAB2020_MEG-AuditoryPatternRecognition/0002/20210413_000000/MEG'; %assigns the final path of the .fif file to the rawName path used in the maxfilter command
pathraw{2} = '/raw/sorted/MINDLAB2020_MEG-AuditoryPatternRecognition/0010/20210416_000000/MEG';
pathraw{3} = '/raw/sorted/MINDLAB2020_MEG-AuditoryPatternRecognition/0011/20210423_000000/MEG';
pathraw{4} = '/raw/sorted/MINDLAB2020_MEG-AuditoryPatternRecognition/0007/20210511_000000/MEG';
pathraw{5} = '/raw/sorted/MINDLAB2020_MEG-AuditoryPatternRecognition/0012/20210622_000000/MEG';

for ii = 5:length(pathraw) %over subjects
    list = dir([pathraw{ii} '/0*']);
    for jj = 1:length(list) %over experimental blocks
        if strcmp(list(jj).name(5:6),'re') || strcmp(list(jj).name(5:6),'sa') || strcmp(list(jj).name(5:6),'vi') || strcmp(list(jj).name(5:6),'pd') %checks whether characters 5 to 6 are equal to 're', 'vi, 'sa' or 'pd'; the loop continues if this is true (1) and it stops if this is false (0)
            
            fpath = dir([list(jj).folder '/' list(jj).name '/files/*.fif']); % looks for .fif file
            rawName = [fpath(1).folder '/' fpath(1).name];
            maxfName = ['SUBJ' pathraw{ii}(56:59) '_' fpath.name(1:end-4) '_bis']; %define the output name of the maxfilter processing
            %movement compensation
            cmd = ['submit_to_cluster -q maxfilter.q -n 4 -p ' ,project, ' "',maxfilter_path,' -f ',[rawName],' -o ' [maxDir '/' maxfName '_tsssdsm.fif'] ' -st 4 -corr 0.98 -movecomp -ds 4 ',' -format float -v | tee ' [maxDir '/log_files/' maxfName '_tsssdsm.log"']];
            %no movement compensation
%             cmd = ['submit_to_cluster -q maxfilter.q -n 4 -p ' ,project, ' "',maxfilter_path,' -f ',[rawName],' -o ' [maxDir '/' maxfName '_tsssdsm.fif'] ' -st 4 -corr 0.98 -ds 4 ',' -format float -v | tee ' [maxDir '/log_files/' maxfName '_tsssdsm.log"']];
            system(cmd);
        end
    end
end

%% Converting the .fif files into SPM objects

%OBS! remember to run 'starting up OSL' first

%setting up the cluster
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/Cluster_ParallelComputing') %add the path to the function that submits the jobs to the cluster
clusterconfig('scheduler', 'cluster');
clusterconfig('long_running', 1); %there are different queues for the cluster depending on the number and length of the jobs you want to submit 
clusterconfig('slot', 1); %slot in the queu

%%

fif_list = dir('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/*.fif'); %creates a list with the .fif files

for ii = 340:341%:length(fif_list) %over the .fif files
    S = []; %structure 'S'                   
    S.dataset = [fif_list(ii).folder '/' fif_list(ii).name];
    D = spm_eeg_convert(S);
%     D = job2cluster(@cluster_spmobject, S); %actual function for conversion
end

%% Removing bad segments using OSLVIEW

%checks data for potential bad segments (periods)
%marking is done by right-clicking in the proximity of the event and click on 'mark event'
%a first click (green dashed label) marks the beginning of a bad period
%a second click indicates the end of a bad period (red)
%this will mean that we are not using about half of the data, but with such bad artefacts this is the best we can do
%we can still obtain good results with what remains
%NB: Push the disk button to save to disk (no prefix will be added, same name is kept)

%OBS! remember to check for bad segments of the signal both at 'megplanar' and 'megmag' channels (you can change the channels in the OSLVIEW interface)
%OBS! remember to mark the trial within the bad segments as 'badtrials' and use the label for removing them from the Averaging (after Epoching) 

spm_list = dir('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/spmeeg*.mat'); %path to SPM objects

for ii = 1:3%:length(spm_list) %over experimental blocks %OBS!
    D = spm_eeg_load([spm_list(ii).folder '/' spm_list(ii).name]);
    D = oslview(D);
    D.save(); %save the selected bad segments and/or channels in OSLVIEW
    disp(ii)
end

%% AFRICA denoising (part I)

%setting up the cluster
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/Cluster_ParallelComputing') %add the path to the function that submits the jobs to the cluster
clusterconfig('scheduler', 'cluster');
clusterconfig('long_running', 1); %there are different queues for the cluster depending on the number and length of the jobs you want to submit 
clusterconfig('slot', 1); %slot in the queu

%%

%ICA calculation
spm_list = dir('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/spmeeg*.mat');

for ii = 339:length(spm_list) %OBS!
    S = [];
    D = spm_eeg_load([spm_list(ii).folder '/' spm_list(ii).name]);
    S.D = D;
    
    jobid = job2cluster(@cluster_africa,S);
%   D = osl_africa(D,'do_ica',true,'do_ident',false,'do_remove',false,'used_maxfilter',true); 
%   D.save();
end

%% AFRICA denoising (part II)

% v = [11 12 19 32];
%visual inspection and removal of artifacted components
%look for EOG and ECG channels (usually the most correlated ones, but check a few more just in case)
spm_list = dir('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/spmeeg*.mat');

for ii = 339:length(spm_list) %OBS!%38:41
    D = spm_eeg_load([spm_list(ii).folder '/' spm_list(ii).name]);
    D = osl_africa(D,'do_ident','manual','do_remove',false,'artefact_channels',{'EOG','ECG'});
    %hacking the function to manage to get around the OUT OF MEMORY problem..
    S = [];
    S.D = D;
    jobid = job2cluster(@cluster_rembadcomp,S);
%   D.save();
    disp(ii)
end

%% Epoching: one epoch per old/new excerpt (baseline = (-100ms)

prefix_tobeadded = 'e'; %adds this prefix to epoched files
spm_list = dir('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/spmeeg*.mat');

for ii = 339:length(spm_list) %over .mat files
    D = spm_eeg_load([spm_list(ii).folder '/' spm_list(ii).name]); %load spm_list .mat files
    dummy = D.fname; %OBS! D.fname does not work, so we need to use a 'dummy' variable instead
    %if strcmp(dummy(22:26), 'speed') %checks whether characters 22 to 26 are equal to 'speed'; the loop continues if this is true (1) and it stops if this is false (0)
    events = D.events; %look for triggers
    %takes the correct triggers sent during the recording
    clear trigcor
    count_evval = 0; %???
    for ieve = 1:length(events) %over triggers
        if strcmp(events(ieve).type,'STI101_up') %only triggers at the beginning of each stimuli
            if events(ieve).value ~= 103 && events(ieve).value ~= 104 && events(ieve).value ~= 128 && events(ieve).value ~= 8 && events(ieve).value ~= 132 && events(ieve).value ~= 48 && events(ieve).value ~= 32 && events(ieve).value ~= 64 %discard 104 and 128 for random triggers
                count_evval = count_evval + 1;
                trigcor(count_evval,1) = events(ieve).time; %+ 0.010; %this takes the correct triggers and add 10ms of delay of the sound travelling into the tubes
                %variable with all the triggers we need
            end
        end
    end
    trl_sam = zeros(length(trigcor),3); %prepare the samples matrix with 0's in all its cells
    trl_sec = zeros(length(trigcor),3); %prepare the seconds matrix with 0's in all its cells
    %deftrig = zeros(length(trigcor),1); %this is not useful
    for k = 1:length(trigcor) %over selected triggers
        %deftrig(k,1) = 0.012 + trigcor(k,1); %adding a 0.012 seconds delay to the triggers sent during the experiment (this delay was due to technical reasons related to the stimuli)
        trl_sec(k,1) = trigcor(k,1) - 0.1000; %beginning time-window epoch in s (please note that we computed this operation two times, obtaining two slightly different pre-stimulus times.
        %this was done because for some computations was convenient to have a slightly longer pre-stimulus time
        %remove 1000ms of baseline
        if strcmp(dummy(22:26), 'speed')
            trl_sec(k,2) = trigcor(k,1) + 5.5; %end time-window epoch in seconds
        else
            trl_sec(k,2) = trigcor(k,1) + 4.4; %end time-window epoch in seconds
        end
        trl_sec(k,3) = trl_sec(k,2) - trl_sec(k,1); %range time-windows in seconds
        trl_sam(k,1) = round(trl_sec(k,1) * 250) + 1; %beginning time-window epoch in samples %250Hz per second
        trl_sam(k,2) = round(trl_sec(k,2) * 250) + 1; %end time-window epoch in samples
        trl_sam(k,3) = -25; %sample before the onset of the stimulus (corresponds to 0.100ms)
    end
    dif = trl_sam(:,2) - trl_sam(:, 1); %difference between the end and the beginning of each sample (just to make sure that everything is fine)
    if ~all(dif == dif(1)) %checking if every element of the vector are the same (i.e. the length of the trials is the same; we may have 1 sample of difference sometimes because of different rounding operations..)
        trl_sam(:,2) = trl_sam(:,1) + dif(1);
    end
    %creates the epochinfo structure that is required for the source reconstruction later
    epochinfo.trl = trl_sam;
    epochinfo.time_continuous = D.time;
    %switch the montage to 0 because for some reason OSL people prefer to do the epoching with the not denoised data
    D = D.montage('switch',0);
    %build structure for spm_eeg_epochs
    S = [];
    S.D = D;
    S.trl = trl_sam;
    S.prefix = prefix_tobeadded;
    D = spm_eeg_epochs(S);
    
    %store the epochinfo structure inside the D object
    D.epochinfo = epochinfo;
    D.save();
    %take bad segments registered in OSLVIEW and check if they overlap with the trials. if so, it gives the number of overlapped trials that will be removed later   
    count = 0;
    Bad_trials = zeros(length(trigcor),1);
    for kkk = 1:length(events) %over events
        if strcmp(events(kkk).type,'artefact_OSL')
            for k = 1:length(trl_sec) %over trials
                if events(kkk).time - trl_sec(k,2) < 0 %if end of trial is > than beginning of artifact
                    if trl_sec(k,1) < (events(kkk).time + events(kkk).duration) %if beginning of trial is < than end of artifact
                        Bad_trials(k,1) = 1; %it is a bad trial (stored here)
                        count = count + 1;
                    end
                end                  
            end
        end
    end
    %if bad trials were detected, their indices are stored within D.badtrials field
    disp(spm_list(ii).name);
    if count == 0
        disp('there are no bad trials marked in oslview');
    else
        D = badtrials(D,find(Bad_trials),1); %get the indices of the badtrials marked as '1' (that means bad)
%         D = conditions(D,find(Bad_trials),1); %get the indices of the badtrials marked as '1' (that means bad)
        epochinfo = D.epochinfo;
        xcv = find(Bad_trials == 1);
        %this should be done only later.. in any case.. not a problem..
        for jhk = 1:length(xcv)
            D = D.conditions(xcv(jhk),'Bad');
            epochinfo.conditionlabels(xcv(jhk)) = {'Bad'};
            disp([num2str(ii) ' - ' num2str(jhk) ' / ' num2str(length(xcv))])
        end
        D.epochinfo = epochinfo;
        D.save(); %saving on disk
        disp('bad trials are ')
        length(D.badtrials)
    end
    D.save();
    disp(ii)
end


%% Defining the conditions - All blocks

%define conditions - 1 epoch for each old/new excerpt (baseline = (-)100ms)

xlsx_dir_behav = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/BehavioralTaskMEG/Version_2/Final_xlsx'; %dir to MEG behavioral results (.xlsx files)
epoch_list = dir('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/e*.mat'); %dir to epoched files

for ii = 339:length(epoch_list) %over epoched data
    D = spm_eeg_load([epoch_list(ii).folder '/' epoch_list(ii).name]);
    dummy = D.fname;
    %barbaric solution.. to build the name to be read for the excel files with the MEG behavioral tasks performance
    if strcmp(dummy(18:23),'recogm')
        dumbloc = 'Block_3.xlsx';
        bl = 3;
    elseif strcmp(dummy(18:23),'recogs')
        dumbloc = 'Block_4.xlsx';
        bl = 4;
    elseif strcmp(dummy(18:23),'sameme')
        dumbloc = 'Block_5.xlsx';
        bl = 5;
    elseif strcmp(dummy(18:23),'visual')
        dumbloc = 'Block_6.xlsx';
        bl = 6;
    elseif strcmp(dummy(18:19),'pd')
        dumbloc = 'Project_PD.xlsx';
        bl = 7;
    end
    if strcmp(dummy((end-14):(end-14)+2),'bis')
        dumls = ['Subj_' dummy(13:16) 'bis_' dumbloc]; %getting subject ID directly from the SPM object to reduce the probability to make mistakes
    else
        dumls = ['Subj_' dummy(13:16) '_' dumbloc];
    end
    [~,~,raw_recog] = xlsread([xlsx_dir_behav '/' dumls]); %excel files
    %picking the current block
    if bl == 3 %block 3 (THIS IS THE ACTUAL BLOCK FOR THIS PAPER)
        for k = 1:length(D.trialonset)
            if raw_recog{(k + 1),3} == 0 %if there was no response
                D = D.conditions(k,'No_response');
            elseif strcmp(raw_recog{(k + 1),2}(8:9),'ol') && raw_recog{(k + 1),3} == 1 %old correct
                D = D.conditions(k,'Old_Correct'); %assign old correct
            elseif strcmp(raw_recog{(k + 1),2}(8:9),'ol') && raw_recog{(k + 1),3} == 2 %old incorrect
                D = D.conditions(k,'Old_Incorrect'); %otherwise assign new correct
            elseif strcmp(raw_recog{(k + 1),2}(14:15),'t1') && raw_recog{(k + 1),3} == 2 %new t1 correct
                D = D.conditions(k,'New_T1_Correct');
            elseif strcmp(raw_recog{(k + 1),2}(14:15),'t1') && raw_recog{(k + 1),3} == 1 %new t1 incorrect
                D = D.conditions(k,'New_T1_Incorrect');
            elseif strcmp(raw_recog{(k + 1),2}(14:15),'t2') && raw_recog{(k + 1),3} == 2 %new t2 correct
                D = D.conditions(k,'New_T2_Correct');
            elseif strcmp(raw_recog{(k + 1),2}(14:15),'t2') && raw_recog{(k + 1),3} == 1 %new t2 incorrect
                D = D.conditions(k,'New_T2_Incorrect');
            elseif strcmp(raw_recog{(k + 1),2}(14:15),'t3') && raw_recog{(k + 1),3} == 2 %new t3 correct
                D = D.conditions(k,'New_T3_Correct');
            elseif strcmp(raw_recog{(k + 1),2}(14:15),'t3') && raw_recog{(k + 1),3} == 1 %new t3 incorrect
                D = D.conditions(k,'New_T3_Incorrect');
            elseif strcmp(raw_recog{(k + 1),2}(14:15),'t4') && raw_recog{(k + 1),3} == 2 %new t4 correct
                D = D.conditions(k,'New_T4_Correct');
            elseif strcmp(raw_recog{(k + 1),2}(14:15),'t4') && raw_recog{(k + 1),3} == 1 %new t4 incorrect
                D = D.conditions(k,'New_T4_Incorrect');
            end
        end
    elseif bl == 4 %block 4
        for k = 1:length(D.trialonset)
            if raw_recog{(k + 1),3} == 0 %if there was no response
                D = D.conditions(k,'No_response');
            elseif strcmp(raw_recog{(k + 1),2}(6:13),'fast_old') && raw_recog{(k + 1),3} == 1 %old fast correct
                D = D.conditions(k,'Old_Fast_Correct'); %assign old fast correct
            elseif strcmp(raw_recog{(k + 1),2}(6:13),'fast_old') && raw_recog{(k + 1),3} == 2 %old fast incorrect
                D = D.conditions(k,'Old_Fast_Incorrect'); %assign old fast incorrect
            elseif strcmp(raw_recog{(k + 1),2}(6:13),'slow_old') && raw_recog{(k + 1),3} == 1 %old slow correct
                D = D.conditions(k,'Old_Slow_Correct'); %assign old slow correct
            elseif strcmp(raw_recog{(k + 1),2}(6:13),'slow_old') && raw_recog{(k + 1),3} == 2 %old slow incorrect
                D = D.conditions(k,'Old_Slow_Incorrect'); %assign old slow incorrect
            elseif strcmp(raw_recog{(k + 1),2}(8:15),'fast_new') && raw_recog{(k + 1),3} == 2 %new fast correct
                D = D.conditions(k,'New_Fast_Correct'); %assign new fast correct
            elseif strcmp(raw_recog{(k + 1),2}(8:15),'fast_new') && raw_recog{(k + 1),3} == 1 %new fast incorrect
                D = D.conditions(k,'New_Fast_Incorrect'); %assign new fast incorrect
            elseif strcmp(raw_recog{(k + 1),2}(8:15),'slow_new') && raw_recog{(k + 1),3} == 2 %new slow correct
                D = D.conditions(k,'New_Slow_Correct'); %assign new slow correct
            else
                D = D.conditions(k,'New_Slow_Incorrect'); %assign new incorrect
            end
        end
    elseif bl == 5 %block 5
        for k = 1:20 %over the 20 encoding trials
            D = D.conditions(k,'Encoding');
        end
        for k = 1:(length(D.trialonset)-20) %over recognition trials (this would work even if the last trial was out of bonds and thus was not epoched..)
            if raw_recog{(k + 23),3} == 0 %if there was no response
                D = D.conditions(k + 20,'No_response');
            elseif strcmp(raw_recog{(k + 23),2}(9:11),'enc') && raw_recog{(k + 23),3} == 1 %encoding correct
                D = D.conditions(k + 20,'Old_Correct'); %assign encoding correct
            elseif strcmp(raw_recog{(k + 23),2}(9:11),'enc') && raw_recog{(k + 23),3} == 2 %encoding incorrect
                D = D.conditions(k + 20,'Old_Incorrect'); %otherwise assign encoding incorrect
            elseif strcmp(raw_recog{(k + 23),2}(9:11),'rec') && raw_recog{(k + 23),3} == 2 %recognition correct
                D = D.conditions(k + 20,'New_Correct'); %assign recognition correct
            else
                D = D.conditions(k + 20,'New_Incorrect'); %assign new incorrect
            end
        end
    elseif bl == 6 %block 6
        for k = 1:20 %over the 20 encoding trials
            D = D.conditions(k,'Encoding');
        end
        for k = 1:(length(D.trialonset)-20) %over recognition trials (this would work even if the last trial was out of bonds and thus was not epoched..)
            if raw_recog{(k + 23),3} == 0 %if there was no response
                D = D.conditions(k + 20,'No_response');
            elseif strcmp(raw_recog{(k + 23),2}(8:10),'old') && raw_recog{(k + 23),3} == 1 %encoding correct
                D = D.conditions(k + 20,'Old_Correct'); %assign encoding correct
            elseif strcmp(raw_recog{(k + 23),2}(8:10),'old') && raw_recog{(k + 23),3} == 2 %encoding incorrect
                D = D.conditions(k + 20,'Old_Incorrect'); %otherwise assign encoding incorrect
            elseif strcmp(raw_recog{(k + 23),2}(8:10),'new') && raw_recog{(k + 23),3} == 2 %recognition correct
                D = D.conditions(k + 20,'New_Correct'); %assign recognition correct
            else
                D = D.conditions(k + 20,'New_Incorrect'); %assign new incorrect
            end
        end
    elseif bl == 7
        for k = 1:length(D.trialonset) %over trials
            if strcmp(raw_recog{(k + 1),2}(3:7),'porco')
                D = D.conditions(k,'pd');
            else
                D = D.conditions(k,'pm');
            end
        end
    end
    %this is for every block
    if ~isempty(D.badtrials) %overwriting badtrials (if any) on condition labels
        BadTrials = D.badtrials;
        for badcount = 1:length(BadTrials) %over bad trials
            D = D.conditions(BadTrials(badcount),'Bad_trial');
        end
    end
    D = D.montage('switch',1);
    D.epochinfo.conditionlabels = D.conditions; %to add for later use in the source reconstruction
    D.save(); %saving data on disk
    disp(num2str(ii))
end

%% COMPUTING STATISTICS OF BEHAVIORAL TASKS IN MEG (Figure S3)

xlsx_dir_behav = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/BehavioralTaskMEG/Version_2/Final_xlsx/'; %dir to MEG behavioral results (.xlsx files)
bb = 3;
list_beh = dir([xlsx_dir_behav 'Subj*' num2str(bb) '.xlsx']); %normal blocks
bl = bb;
if bl == 3
    Block_3 = cell(length(list_beh)+1,13);
elseif bl == 4
    Block_4 = cell(length(list_beh)+1,11);
elseif bl == 5
    Block_5 = cell(length(list_beh)+1,7);
elseif bl == 6
    Block_6 = cell(length(list_beh)+1,7);
end
for ii = 1:length(list_beh) %over subjects for block bb
    %barbaric solution.. to build the name to be read for the excel files with the MEG behavioral tasks performance
    [~,~,raw_recog] = xlsread([list_beh(ii).folder '/' list_beh(ii).name]); %excel files
    %picking the current block
    %legend
    Block_3{1,1} = 'Subject'; Block_3{1,2} = 'OLD_Cor'; Block_3{1,3} = 'OLD_Cor %'; Block_3{1,4} = 'New_T1_Cor'; Block_3{1,5} = 'New_T1_Cor %'; %1st row
    Block_3{1,6} = 'New_T2_Cor'; Block_3{1,7} = 'New_T2_Cor %'; Block_3{1,8} = 'New_T3_Cor'; Block_3{1,9} = 'New_T3_Cor %'; Block_3{1,10} = 'New_T4_Cor'; Block_3{1,11} = 'New_T4_Cor %'; Block_3{1,12} = 'No response'; Block_3{1,13} = 'No response %'; %1st row
    Block_3{1,14} = 'OLD_RT'; Block_3{1,15} = 'NEWT1_RT'; Block_3{1,16} = 'NEWT2_RT'; Block_3{1,17} = 'NEWT3_RT'; Block_3{1,18} = 'NEWT4_RT';
    nr = 0; old = 0; n1 = 0; n2 = 0; n3 = 0; n4 = 0;
    ort = []; n1rt = []; n2rt = []; n3rt = []; n4rt = [];
    for k = 1:135
        if raw_recog{(k + 1),3} == 0 %if there was no response
            nr = nr + 1;
        elseif strcmp(raw_recog{(k + 1),2}(8:9),'ol') && raw_recog{(k + 1),3} == 1 %old correct
            old = old + 1;
            ort = cat(1,ort,raw_recog{(k + 1),4});
        elseif strcmp(raw_recog{(k + 1),2}(14:15),'t1') && raw_recog{(k + 1),3} == 2 %new t1 correct
            n1 = n1 + 1;
            n1rt = cat(1,n1rt,raw_recog{(k + 1),4});
        elseif strcmp(raw_recog{(k + 1),2}(14:15),'t2') && raw_recog{(k + 1),3} == 2 %new t2 correct
            n2 = n2 + 1;
            n2rt = cat(1,n2rt,raw_recog{(k + 1),4});
        elseif strcmp(raw_recog{(k + 1),2}(14:15),'t3') && raw_recog{(k + 1),3} == 2 %new t3 correct
            n3 = n3 + 1;
            n3rt = cat(1,n3rt,raw_recog{(k + 1),4});
        elseif strcmp(raw_recog{(k + 1),2}(14:15),'t4') && raw_recog{(k + 1),3} == 2 %new t4 correct
            n4 = n4 + 1;
            n4rt = cat(1,n4rt,raw_recog{(k + 1),4});
        end
    end
    disp(num2str(['Block ' num2str(bb) ' - Subject ' num2str(ii)]))
    Block_3{ii+1,1} = list_beh(ii).name(6:9); Block_3{ii+1,2} = old; Block_3{ii+1,3} = (old/27)*100; Block_3{ii+1,4} = n1; Block_3{ii+1,5} = (n1/27)*100;
    Block_3{ii+1,6} = n2; Block_3{ii+1,7} = (n2/27)*100; Block_3{ii+1,8} = n3; Block_3{ii+1,9} = (n3/27)*100; Block_3{ii+1,10} = n4; Block_3{ii+1,11} = (n4/27)*100; Block_3{ii+1,12} = nr; Block_3{ii+1,13} = (nr/135)*100;
    Block_3{ii+1,14} = mean(ort); Block_3{ii+1,15} = mean(n1rt); Block_3{ii+1,16} = mean(n2rt); Block_3{ii+1,17} = mean(n3rt); Block_3{ii+1,18} = mean(n4rt);
end
Block_3_t = cell2table(Block_3); %remove the possible empty cell

%statistics (ANOVAs)
data = cell2mat(Block_3(2:end,[2 4 6 8 10])); %correct responses
datart = cell2mat(Block_3(2:end,14:18)); %reaction times
%correct responses
% [p,t,stats] = anova1(data); %'off' for not showing the plot
[p,t,stats] = kruskalwallis(data); %'off' for not showing the plot
[c,m,h,nms] = multcompare(stats,'ctype','tukey-kramer'); %perform multiple comparison test based on anova1 output = post-hoc analysis (c needs to be saved for every time point and every channel)
%reaction times
% [p,t,stats] = anova1(datart); %'off' for not showing the plot
[p,t,stats] = kruskalwallis(datart); %'off' for not showing the plot
[c,m,h,nms] = multcompare(stats,'ctype','tukey-kramer'); %perform multiple comparison test based on anova1 output = post-hoc analysis (c needs to be saved for every time point and every channel)

%plotting correct responses
%some default color specifications for later plotting..
% cb = [0.5 0.8 0.9; 1 1 0.7; 0.7 0.8 0.9; 0.8 0.5 0.4; 0.5 0.7 0.8; 1 0.8 0.5; 0.7 1 0.4; 1 0.7 1; 0.6 0.6 0.6; 0.7 0.5 0.7; 0.8 0.9 0.8; 1 1 0.4];
cl = [0 0 0.8; 0 0 0.8; 0 0 0.8; 0 0 0.8; 0 0 0.8;];
datadum = cell(1,size(data,2));
for ii = 1:size(data,2)
    datadum{1,ii} = data(:,ii);
end
figure
rm_raincloud2(datadum',cl)
grid minor
set(gcf,'color','w')
set(gcf,'Position',[200,200,400,550])
xlabel('Correct responses'); %set(gca,'YDir','normal');
xlim([-5 38])
FigS3 = []; %providing source data
FigS3.response_accuracy = datadum;

%plotting reaction times
%some default color specifications for later plotting..
% cb = [0.5 0.8 0.9; 1 1 0.7; 0.7 0.8 0.9; 0.8 0.5 0.4; 0.5 0.7 0.8; 1 0.8 0.5; 0.7 1 0.4; 1 0.7 1; 0.6 0.6 0.6; 0.7 0.5 0.7; 0.8 0.9 0.8; 1 1 0.4];
cl = [0 0 0.8; 0 0 0.8; 0 0 0.8; 0 0 0.8; 0 0 0.8;];
datadum = cell(1,size(datart,2));
for ii = 1:size(datart,2)
    datadum{1,ii} = datart(:,ii);
end
figure
rm_raincloud2(datadum',cl)
grid minor
set(gcf,'color','w')
set(gcf,'Position',[200,200,400,550])
xlabel('Time (ms)'); %set(gca,'YDir','normal');
xlim([1500 3500])
FigS3.reaction_times = datadum;
save FigS3.mat FigS3

%% Averaging and Combining planar gradiometers

%settings for cluster (parallel computing)
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/Cluster_ParallelComputing') %add the path to the function that submits the jobs to the cluster
clusterconfig('scheduler', 'cluster'); %set automatically the long run queue
clusterconfig('long_running', 1); %set automatically the long run queue
clusterconfig('slot', 1); %set manually the job cluster slots
% between 1 and 12 (n x 8gb of ram)

%% averaging

output_prefix_to_be_set = 'm';

v = [1]; %only a selection of files

epoch_list = dir ('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/e*mat'); %dir to epoched files (encoding)
for ii = 1:length(v)%1:length(epoch_list) %over epoched files
    %distribute 
    input = [];
    input.D = [epoch_list(v(ii)).folder '/' epoch_list(v(ii)).name];
    input.prefix = output_prefix_to_be_set;
    jobid = job2cluster(@sensor_average, input); % this is the command for send the job to the cluster, in the brackets you can find the name on the function to run (afeter the @) and the variable for the input (in this case input)
    % look the script for more details about the function work
end

%% combining planar gradiometers

average_list = dir ('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/m*mat'); %dir to epoched files (encoding)
v = [339:342]; %only a selection of files

for ii = 1:length(v)%1:length(average_list) %over files
    input = [];
    input.D = [average_list(v(ii)).folder '/' average_list(v(ii)).name];
    D = spm_eeg_load(input.D);
    D = D.montage('switch',1);
    D.save();
    jobid = job2cluster(@combining_planar_cluster, input); % this is the command for send the job to the cluster, in the brackets you can find the name on the function to run (afeter the @) and the variable for the input (in this case input)
end

%% LBPD_startup_D

pathl = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD'; %path to stored functions
addpath(pathl);
LBPD_startup_D(pathl);
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/Cluster_ParallelComputing') %add the path where is the function for submit the jobs to the server

%% Extracting and plotting MEG sensor data (Figure S5)

%%% This is currently configured for plotting MEG channel MEG0211 (= channels_plot = 9) 

block = 3; % 3 = recogminor; 4 = recogspeed; 5 = samemel; 6 = visualpat; 7 = pd
channels_plot = []; %13;95;;9; %95, 13, 15, 11, 141, 101, 9 % empty for plotting single channels; otherwise number(s) of channels to be averaged and plotted (e.g. [13] or [13 18])
% channels_plot = []; % empty for plotting single channels; otherwise number(s) of channels to be averaged and plotted (e.g. [13] or [13 18])

%1321 1411
save_data = 1;
load_data = 0; %set 1 if you want to load the data instead of extracting it from SPM objects
% v = [1]; %subjects
%bad 8,9 and a bit 6 (recogminor)

S = [];
%computing data
if block == 3
    S.conditions = {'Old_Correct','New_T1_Correct','New_T2_Correct','New_T3_Correct','New_T4_Correct'};
%     S.conditions = {'Old_Correct','New_T3_Correct'};
elseif block == 4
    % S.conditions = {'Old_Fast_Correct','Old_Slow_Correct','New_Fast_Correct','New_Slow_Correct'};
%     S.conditions = {'Old_Fast_Correct','New_Fast_Correct'};
    S.conditions = {'Old_Slow_Correct','New_Slow_Correct'};
elseif block == 5
    S.conditions = {'Old_Correct','New_Correct'};
%     S.conditions = {'Old_Correct','Encoding'};
elseif block == 6
    S.conditions = {'Old_Correct','New_Correct'};
%     S.conditions = {'Old_Correct','Encoding'};
elseif block == 7
    S.conditions = {'pd','pm'};
end
if block == 3
    list = dir ('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/P*_recogminor*_tsssdsm.mat'); %dir to epoched files (encoding)
elseif block == 4
    list = dir ('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/P*recogspeed*_tsssdsm.mat'); %dir to epoched files (encoding)
elseif block == 5
    list = dir ('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/P*samemel*_tsssdsm.mat'); %dir to epoched files (encoding)
elseif block == 6
    list = dir ('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/P*visualpat*_tsssdsm.mat'); %dir to epoched files (encoding)
elseif block == 7
    list = dir ('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/P*pd*_tsssdsm.mat'); %dir to epoched files (encoding)
end
v = 1:length(list); %subjects
% v = [];
if ~exist('chanlabels','var')
    load('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter/MEG_sensors/recogminor_all_conditions.mat', 'chanlabels')
end
outdir = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/MEG_sensors'; %path to write output in
S.outdir = outdir;
S.data = [];
if load_data == 1 %if you already computed and saved on disk the t-tests you can load them here
%     load([outdir '/Block_' num2str(block) '.mat']);
    S.data = data_mat(:,:,v,:);
    S.chanlabels = chanlabels;
    S.time_real = time_sel;
else %otherwise you can extract the data from SPM MEEG objects (one for each subject)
%     S.spm_list = cell(1,length(list));
% v = 7;
    S.spm_list = cell(1,length(v));
    for ii = 1:length(v)
        S.spm_list(ii) = {[list(v(ii)).folder '/' list(v(ii)).name]};
    end
end

S.timeextract = []; %time-points to be extracted
S.centerdata0 = 0; %1 to make data starting at 0
S.save_data = save_data; %only meaningfull if you read data from SPM objects saved on disk
S.save_name_data = ['Block_' num2str(block)];

%individual waveform plotting
if isempty(channels_plot)
    S.waveform_singlechannels_label = 1; %1 to plot single channel waveforms
else
    S.waveform_singlechannels_label = 0; %1 to plot single channel waveforms
end
S.wave_plot_conditions_together = 0; %1 for plotting the average of all
S.mag_lab = 1; %1 for magnetometers; 2 for gradiometers
S.x_lim_temp_wave = [-0.1 3.4]; %limits for time (in secs) (E.g. [-0.1 3.4])
S.y_lim_ampl_wave = []; %limit for amplitude (E.g. [0 120] magnetometes, [0 6] gradiometers)

%averaged waveform plotting
if isempty(channels_plot)
    S.waveform_average_label = 0; %average of some channels
    S.left_mag = 95; %13 %37 (visual) %43 (visual) %199 (visual) %203 (visual) %channels for averaging
else
    S.waveform_average_label = 1; %average of some channels
    S.left_mag = channels_plot; %13 %37 (visual) %43 (visual) %199 (visual) %203 (visual) %channels for averaging
end
% S.left_mag = [2:2:204];
S.legc = 1; %set 1 for legend
% S.left_mag = 99;
S.signtp = {[]};
% S.sr = 150; %sampling rate (Hz)
S.avewave_contrast = 0; %1 to plot the contrast between conditions (averaged waveform)
S.save_label_waveaverage = 0;
S.label_plot = 'c';
%t-tests
S.t_test_for_permutations = 0;
S.cond_ttests_tobeplotted_topoplot = [1 2]; %this is for both topoplot and t-tests!! (here [1 2] means cond1 vs cond2!!!!!!!)

%topoplotting
S.topoplot_label = 0;
S.fieldtrip_mask = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External';
S.topocontr = 0;
S.topocondsing = [1]; %condition for topoplot
% S.xlim = [0.75 0.85]; %time topolot
% S.xlim = [1.1 1.2]; %time topolot
S.xlim = [0.25 0.25]; 
S.zlimmag = []; %magnetometers amplitude topoplot limits
S.zlimgrad = []; %gradiometers amplitude topoplot limits
S.colormap_spec = 0;
% x = []; x.bottom = [0 0 1]; x.botmiddle = [0 0.5 1]; x.middle = [1 1 1]; x.topmiddle = [1 1 0.5]; x.top = [1 0.95 0]; %yellow - blue
x = []; x.bottom = [0 0 0.5]; x.botmiddle = [0 0.5 1]; x.middle = [1 1 1]; x.topmiddle = [1 0 0]; x.top = [0.6 0 0]; %red - blue
S.colormap_spec_x = x;
S.topoplot_save_label = 0;
color_line = colormap(lines(5)); %extracting some colours from a colormap
S.color_line = color_line;
S.color_line(1,:) = color_line(2,:);
S.color_line(2,:) = color_line(1,:);
S.color_line(5,:) = [0.4 0.4 0.4];

[out] = MEG_sensors_plotting_ttest_LBPD_D2(S);

%%

%% *** DECODING (MULTIVARIATE PATTERN ANALYSIS) ***

% Figures 2a, 2b, 2c and S4

%%

%% decoding - support vector machine (SVM) - multivariate pattern analysis - FIGURE 3A

%The decoding consisted of support vector machines (SVM) implemented in
%external functions that can be found at the following link:
%http://www.csie.ntu.edu.tw/~cjlin/libsvm/
%Those functions can be also used in Matlab with a few adjustments done to achieve a proper imlplementation
%as described in the paper and including crucial steps such as, for example, cross validation.
%In our case, we implemented them by using AU server, therefore the codes immediately below show how we provided the data to the AU server.
%After that, we provide the output of the SVM algorithm so that you can reproduce
%the statistics and the plotting solutions that we presented in the paper.
%If you wish to know more about this SVM implementation, you are
%very welcome to contact us.


%% settings for AU server and SVM algorithm

cd /projects/MINDLAB2017_MEG-LearningBach/scripts/MITDecodingToServer/Decoding
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/MITDecodingToServer/Decoding')
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/MITDecodingToServer/Decoding/scilearnlab')
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/MITDecodingToServer/Decoding/scilearnlab/external');
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/MITDecodingToServer/Decoding/scilearnlab/external/libsvm-3.21')
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/MITDecodingToServer/Decoding/scilearnlab/external/libsvm-3.21/matlab')
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/MITDecodingToServer/Decoding/SpatialLocation')
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/MITDecodingToServer/Decoding/SpatialLocation/plotchannel')
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/MITDecodingToServer/Decoding/permutationlab')
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/MITDecodingToServer/Decoding/permutationlab/private1')
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/MITDecodingToServer/Decoding/permutationlab/developer')
make

addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/Cluster_ParallelComputing')
% config of the cluster server
clusterconfig('slot', 1); % set manually the job cluster slots
% between 1 and 12 (n x 8gb of ram)
clusterconfig('scheduler','cluster'); % set automatically the long run queue
clusterconfig('long_running', 1); % set automatically the long run queue

%% ORIGINAL SUBMISSION

%% preparing data for decoding functions (for AU server)

time = [1 1126]; %defining time of interest (this refers to 0.100sec before the onset of the stimuli)
numperm = 100; %number of permutations
kfold = 5; %number of categories for k-fold classification
pairl = 1; %1 = pairwise decoding; 2 = multiclass (confusion matrix)

list = dir ('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/e*_recogminor*_tsssdsm.mat'); %dir to epoched files (encoding)
savedirdum = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Decoding';
D1 = cell(length(list),2);
subjs_real = cell(1,length(list));
for ii = 1:length(list) %length(subjnum)
        %loading data
        D = spm_eeg_load([list(ii).folder '/' list(ii).name]);
        old = find(strcmp('Old_Correct',D.conditions)); %getting condition old indices
        for nn = 1:4 %over NewTx conditions
            new = find(strcmp(['New_T' num2str(nn) '_Correct'],D.conditions)); %getting condition newTx indices
            savedir = [savedirdum '/Old_vs_NewT' num2str(nn)];
            mkdir(savedir);
            clear data condid ind
            condid(1:length(old)) = {'O'}; %creating condition labels (first old)
            condid(length(old):(length(old)+length(new))) = {'N'}; %(then new)
            ind = [old new]; %rearranging indices (first old and then new)
            idxMEGc1 = find(strcmp(D.chanlabels,'MEG0111')); %getting extreme channels indexes
            idxMEGc2 = find(strcmp(D.chanlabels,'MEG2643')); %same
            data = D(idxMEGc1:idxMEGc2,time(1):time(2),ind); %extracting data
            %structure with inputs for AU server
            S = [];
            S.pairl = pairl;
            S.data = data;
            S.condid = condid;
            S.numperm = numperm;
            S.kfold = kfold;
            S.savedir = savedir;
            S.ii = ii;
            disp(ii)
            %job to cluster
            jobid = job2cluster(@decoding, S);
        end
    if ii == 1 %saving D.time (time in seconds), only once
        time_sel = D.time(time(1):time(2));
        save([savedir 'time.mat'],'time_sel');
    end
end


%% permutation test on decoding accuracy time-series (pairwise-decoding and temporal generalization)

T_cond = 0; %1 = Old vs NewT1; 2 = Old vs NewT2; 3 = Old vs NewT3; 4 = Old vs NewT4; 0 = all together (only waveform plot);
tempgen_l = 1; %1 computing (and plotting) statistics for temporal generalization; 0 = not doing it (it takes time..)
topoplot_l = 0; %1 for topoplot; 0 for no topoplot (remember that you also have to specify the time for the topo-plot)


%directory where decoding results are stored in
datadirdum = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Decoding';
%loading and reshpaping data outputted by AU server (SVM algorithm) for statistics
load([datadirdum '/time.mat']);
if T_cond ~= 0
    datadir = [datadirdum '/Old_vs_NewT' num2str(T_cond)];
    listPD = dir([datadir '/PD*']);
    listGT = dir([datadir '/TG*']);
    %loading data
    dd1 = zeros(length(time_sel),length(listPD));
    ddpatt = zeros(306,length(time_sel),length(listPD));
    ddTG = zeros(length(time_sel),length(time_sel),length(listPD));
    cnt = 0;
    for ii = 1:length(listPD)
        cnt = cnt + 1;
        load([listPD(ii).folder '/' listPD(ii).name]);
        dd1(:,cnt) = d.d;
        ddpatt(:,:,cnt) = d.pattern;
        load([listGT(ii).folder '/' listGT(ii).name])
        ddTG(:,:,cnt) = d.d;
        disp(ii)
    end
    %average over subjects of pairwise decoding accuracy time-series
    dd12 = mean(dd1,2);
    dd12std = std(dd1,0,2)./sqrt(size(dd1,2)); %standard errors
    %reshaping pairwise decoding acuracy time-series for submitting them to statistical testing
    dd2 = dd1 - 50; %subtracting 50 since the function tests significance against 0 (here represented by 50% chance level)
    dat = cell(1,size(dd2,2));
    for ii = 1:size(dd2,2)
        dat(1,ii) = {dd2(:,ii)'};
    end
    stat = pl_permtest(dat,'alpha',0.05); %actual function
    %plotting p-values after FDR correction for multiple comparisons
    xx = stat.FDR.criticalmap;
    CCt = bwconncomp(xx); %getting significant 'islands'
    %plotting average over subjects (pairwise decoding accuracy) and the corresponding significance
    figure;hold on;grid on;grid minor;
    plot(time_sel,dd12,'k','linewidth',2,'DisplayName',['M_vs_NewT' num2str(T_cond)]);
    hold on
    plot(time_sel,dd12 + dd12std,':','color','k','linewidth',0.5);
    hold on
    plot(time_sel,dd12 - dd12std,':','color','k','linewidth',0.5);
    xlabel('Time (s)');ylabel('Decoding accuracy (%)');
    hold on
    % Show significant time window
    patch_color = [.85 .85 .85]; % Color of grey box marking time range
    ylims = get(gca,'YLim');
    sgf = CCt.PixelIdxList;
    for ii = 1:length(sgf)
        dumbum = sgf{ii};
        sgf2 = time_sel(dumbum);
%         patch([sgf2(1) sgf2(end) sgf2(end) sgf2(1)],[ylims(1) ylims(1) ylims(2) ylims(2)],patch_color,'EdgeColor','none','FaceAlpha',.7)
        patch([sgf2(1) sgf2(end) sgf2(end) sgf2(1)],[40 40 70 70],patch_color,'EdgeColor','none','FaceAlpha',.7)
        hold on
    end
    plot(time_sel,dd12,'k','linewidth',2,'DisplayName',['M_vs_NewT' num2str(T_cond)]);
    hold on
    plot(time_sel,dd12 + dd12std,':','color','k','linewidth',0.5);
    hold on
    plot(time_sel,dd12 - dd12std,':','color','k','linewidth',0.5);
    xlabel('Time (s)');ylabel('Decoding accuracy (%)'); %(slightly stupid..) trick to get the time-series in front..
    set(gcf,'Color','w')
    legend('show')
    xlim([time_sel(1) time_sel(end)])
    ylim([40 70]);
    if topoplot_l == 1
        %plotting patterns (derived from weights) for pairwise decoding
        warning('have you selected the proper time-window for topoplot..?')
        signt = [0.5 2.0]; %significant time-point(s) you want to plot
        %extracting magnetometers for plotting purposes
        dp = mean(ddpatt,3); %average over subjects
        avg = dp(1:3:end,:); %extracting magnetometers only
        %plotting topoplot (fieldtrip function)
        %creating the mask for the data
        fieldtrip_example = load('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/fieldmask.mat');
        label2 = fieldtrip_example.M_timb_lock_gradplanarComb.label;
        label = label2(103:end);
        %setting labels according to the ones in the layout
        for ii = 1:102
            label{ii,1} = [label{ii,1}(1:3) label{ii,1}(5:end)];
        end
        cfgdummy = fieldtrip_example.M_timb_lock_gradplanarComb.cfg;
        cfgdummy.previous = [];
        data = [];
        data.cfg = cfgdummy;
        data.time(1,:) = time_sel(:);
        data.label = label;
        data.dimord = fieldtrip_example.M_timb_lock_gradplanarComb.dimord;
        data.grad = fieldtrip_example.M_timb_lock_gradplanarComb.grad;
        data.avg = avg;
        %creating the cfg for the actual plotting (magnetometers)
        cfg = [];
        cfg.layout = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/neuromag306mag.lay';
        cfg.colorbar = 'yes';
        cfg.xlim = signt; %set temporal limits (in seconds)
        cfg.zlim = [];
        cfg.colormap = 'jet';
        figure
        ft_topoplotER(cfg,data);
        set(gcf,'Color','w')
        %colormap with white for 0 values
        x = []; x.bottom = [0 0 0.5]; x.botmiddle = [0 0.5 1]; x.middle = [1 1 1]; x.topmiddle = [1 0 0]; x.top = [0.6 0 0]; %red - blue
        colormap(bluewhitered_PD(0,x))
    end
    %plotting temoral generalization
    ddTGm = mean(ddTG,3);
    figure; imagesc(time_sel,time_sel,ddTGm); xlabel('Train time (s)'); ylabel('Test time (s)'); set(gca,'YDir','normal');
    colorbar
    set(gcf,'Color','w')
    %reshaping temporal generalization decoding acuracy for submitting it to statistical testing
    ddTG2 = ddTG(1:776,1:776,:) - 50; %subtracting 50 since the function tests significance against 0 (here represented by 50% chance level)
    dat = cell(1,size(ddTG2,3));
    for ii = 1:size(ddTG2,3)
        dat(1,ii) = {ddTG2(:,:,ii)};
    end
    if tempgen_l == 1
        %computing statistics and plotting results
        stat = pl_permtest(dat,'alpha',0.05); %actual function
        STAT{T_cond} = stat; %storing statistics to compare images later and decide the limits for the colorbars
        time_selj = time_sel(1:776);
        figure; imagesc(time_selj,time_selj,stat.statmap); xlabel('Train time (s)'); ylabel('Test time (s)'); set(gca,'YDir','normal');
        colorbar
        set(gcf,'Color','w')
        caxis([-5 10])

        %testing contrast results (correcting for multiple comparison) by performing Monte Carlo simulations
        P2 = zeros(776,776);
        statt = stat.statmap;
        P2(statt>2.6) = 1; %threshold t-val = 2.6 corresponding to p-val < 0.01 (obtained by dividing 0.05 by the 4 comparisons employed here)
        thresh = 0;
        permut = 1000;
        threshMC = 0.001;
        perm_max = 1;
        t1 = time_sel; t2 = t1;
        
        [ OUT ] = twoD_MCS_LBPD_D( P2, thresh, permut, threshMC, perm_max, t1 , t2 )
        
        outdir = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/Beam_abs_1_sens_1_freq_broadband/Sources_Decoding_MEG0211/MCS/FinalImagesToWorkBench/Statistics';
        PDn = cell2table(OUT); %table
        writetable(PDn,[outdir '/TemporalGeneralization_OldvsNewT' num2str(T_cond) 'Testing_Training.xlsx'],'Sheet',1); %printing excel file
    end
else %plotting together time series of decoding accuracy for all 4 decoding analyses (i.e. Old vs NewT1, Old vs NewT2, Old vs NewT3, Old vs NewT4)
    figure
    %     COL{1} = 'b'; COL{2} = 'r'; COL{3} = 'k'; COL{4} = 'g'; COL{5} = 'c'; COL{6} = 'm'; %colors (TO BE MADE THE SAME AS FOR WAVEFORMS OF DIFFERENT CONDITIONS!!)
    color_line = colormap(lines(5)); %extracting some colours from a colormap
    color_line2 = color_line;
    color_line2(1,:) = color_line(2,:);
    color_line2(2,:) = color_line(1,:);
    color_line2(5,:) = [0.4 0.4 0.4];
    FigS4data = zeros(83,1126,4);
    FigS4 = [];
    for nn = 1:4 %over NewTX conditions
        datadir = [datadirdum '/Old_vs_NewT' num2str(nn)];
        listPD = dir([datadir '/PD*']);
        listGT = dir([datadir '/TG*']);
        %loading data
        dd1 = zeros(length(time_sel),length(listPD));
        ddpatt = zeros(306,length(time_sel),length(listPD));
        ddTG = zeros(length(time_sel),length(time_sel),length(listPD));
        cnt = 0;
        for ii = 1:length(listPD)
            cnt = cnt + 1;
            load([listPD(ii).folder '/' listPD(ii).name]);
            dd1(:,cnt) = d.d;
            ddpatt(:,:,cnt) = d.pattern;
            load([listGT(ii).folder '/' listGT(ii).name])
            ddTG(:,:,cnt) = d.d;
            disp([num2str(nn) ' - ' num2str(ii)])
        end
        %average over subjects of pairwise decoding accuracy time-series
        dd12 = mean(dd1,2);
        dd12std = std(dd1,0,2)./sqrt(size(dd1,2)); %standard errors
        %reshaping pairwise decoding acuracy time-series for submitting them to statistical testing
        dd2 = dd1 - 50; %subtracting 50 since the function tests significance against 0 (here represented by 50% chance level)
        dat = cell(1,size(dd2,2));
        for ii = 1:size(dd2,2)
            dat(1,ii) = {dd2(:,ii)'};
        end
        if nn == 4
            FigS4data(1:82,:,nn) = dd2';
            FigS4data(83,:,nn) = NaN;
            for nnn = 1:4 %over NewTX conditions
                datadir = [datadirdum '/Old_vs_NewT' num2str(nnn)];
                listPD = dir([datadir '/PD*']);
                FigS4(nnn).list = listPD;
            end
        else
            FigS4data(:,:,nn) = dd2';
        end
        %plotting average over subjects (pairwise decoding accuracy)
        grid on;grid minor;
        plot(time_sel,dd12,'color',color_line2(nn+1,:),'linewidth',2,'DisplayName',['NewT' num2str(nn)]);
        hold on
        plot(time_sel,dd12 + dd12std,':','color',color_line2(nn+1,:),'linewidth',0.5,'DisplayName',['NewT' num2str(nn)]);
        hold on
        plot(time_sel,dd12 - dd12std,':','color',color_line2(nn+1,:),'linewidth',0.5,'DisplayName',['NewT' num2str(nn)]);
        xlabel('Time (s)');ylabel('Decoding accuracy (%)');
        hold on
    end
    legend('show')
    set(gcf,'Color','w')
    xlim([time_sel(1) time_sel(end)])
    FigS4b.data = FigS4data;
    FigS4b.subjects = FigS4;
    FigS4 = FigS4b;
    FigS4.time = time_sel;
    save FigS4.mat FigS4
end
% %trick to export results in excel sheet.. I initialized dumbo in the command window
% dumbo(:,T_cond) = stat.FDR.criticalmap;
% dumbo = dumbo';
% dumbo2 = zeros(5,1126);
% dumbo2(1,1:1126) = time_sel;
% dumbo2(2:5,1:1126) = dumbo;
% outdir = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/Beam_abs_1_sens_1_freq_broadband/Sources_Decoding_MEG0211/MCS/FinalImagesToWorkBench/Statistics';
% xlswrite([outdir '/PairwiseDecoding.xlsx'],dumbo2)

%% REVISION I

%% A FEW ADJUSTMENTS OF THE DECODING CODES (e.g. to produce temporal generalisation matrix showing only the significant values after correction for multiple comparisons

%% permutation test on decoding accuracy time-series (pairwise-decoding and temporal generalization)

clear
close all
T_cond = 4; %1 = Old vs NewT1; 2 = Old vs NewT2; 3 = Old vs NewT3; 4 = Old vs NewT4; 0 = all together (only waveform plot);
tempgen_l = 1; %1 computing (and plotting) statistics for temporal generalization; 0 = not doing it (it takes time..)
topoplot_l = 0; %1 for topoplot; 0 for no topoplot (remember that you also have to specify the time for the topo-plot)

%directory where decoding results are stored in
datadirdum = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Decoding';
%loading and reshpaping data outputted by AU server (SVM algorithm) for statistics
load([datadirdum '/time.mat']);
if T_cond ~= 0
    datadir = [datadirdum '/Old_vs_NewT' num2str(T_cond)];
    listPD = dir([datadir '/PD*']);
    listGT = dir([datadir '/TG*']);
    %loading data
    dd1 = zeros(length(time_sel),length(listPD));
    ddpatt = zeros(306,length(time_sel),length(listPD));
    ddTG = zeros(length(time_sel),length(time_sel),length(listPD));
    cnt = 0;
    for ii = 1:length(listPD)
        cnt = cnt + 1;
        load([listPD(ii).folder '/' listPD(ii).name]);
        dd1(:,cnt) = d.d;
        ddpatt(:,:,cnt) = d.pattern;
        load([listGT(ii).folder '/' listGT(ii).name])
        ddTG(:,:,cnt) = d.d;
        disp(ii)
    end
    %average over subjects of pairwise decoding accuracy time-series
    dd12 = mean(dd1,2);
    dd12std = std(dd1,0,2)./sqrt(size(dd1,2)); %standard errors
    %reshaping pairwise decoding acuracy time-series for submitting them to statistical testing
    dd2 = dd1 - 50; %subtracting 50 since the function tests significance against 0 (here represented by 50% chance level)
    dat = cell(1,size(dd2,2));
    for ii = 1:size(dd2,2)
        dat(1,ii) = {dd2(:,ii)'};
    end
    stat = pl_permtest(dat,'alpha',0.05); %actual function
    %plotting p-values after FDR correction for multiple comparisons
    xx = stat.FDR.criticalmap;
    CCt = bwconncomp(xx); %getting significant 'islands'
    %plotting average over subjects (pairwise decoding accuracy) and the corresponding significance
    figure;hold on;grid on;grid minor;
    plot(time_sel,dd12,'k','linewidth',2,'DisplayName',['M_vs_NewT' num2str(T_cond)]);
    hold on
    plot(time_sel,dd12 + dd12std,':','color','k','linewidth',0.5);
    hold on
    plot(time_sel,dd12 - dd12std,':','color','k','linewidth',0.5);
    xlabel('Time (s)');ylabel('Decoding accuracy (%)');
    hold on
    % Show significant time window
    patch_color = [.85 .85 .85]; % Color of grey box marking time range
    ylims = get(gca,'YLim');
    sgf = CCt.PixelIdxList;
    for ii = 1:length(sgf)
        dumbum = sgf{ii};
        sgf2 = time_sel(dumbum);
%         patch([sgf2(1) sgf2(end) sgf2(end) sgf2(1)],[ylims(1) ylims(1) ylims(2) ylims(2)],patch_color,'EdgeColor','none','FaceAlpha',.7)
        patch([sgf2(1) sgf2(end) sgf2(end) sgf2(1)],[40 40 70 70],patch_color,'EdgeColor','none','FaceAlpha',.7)
        hold on
    end
    plot(time_sel,dd12,'k','linewidth',2,'DisplayName',['M_vs_NewT' num2str(T_cond)]);
    hold on
    plot(time_sel,dd12 + dd12std,':','color','k','linewidth',0.5);
    hold on
    plot(time_sel,dd12 - dd12std,':','color','k','linewidth',0.5);
    xlabel('Time (s)');ylabel('Decoding accuracy (%)'); %(slightly stupid..) trick to get the time-series in front..
    set(gcf,'Color','w')
    legend('show')
    xlim([time_sel(1) time_sel(end)])
    ylim([40 70]);
    if topoplot_l == 1
        %plotting patterns (derived from weights) for pairwise decoding
        warning('have you selected the proper time-window for topoplot..?')
        signt = [1.6 1.8]; %significant time-point(s) you want to plot
        %extracting magnetometers for plotting purposes
        dp = mean(ddpatt,3); %average over subjects
        avg = dp(1:3:end,:); %extracting magnetometers only
        %plotting topoplot (fieldtrip function)
        %creating the mask for the data
        fieldtrip_example = load('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/fieldmask.mat');
        label2 = fieldtrip_example.M_timb_lock_gradplanarComb.label;
        label = label2(103:end);
        %setting labels according to the ones in the layout
        for ii = 1:102
            label{ii,1} = [label{ii,1}(1:3) label{ii,1}(5:end)];
        end
        cfgdummy = fieldtrip_example.M_timb_lock_gradplanarComb.cfg;
        cfgdummy.previous = [];
        data = [];
        data.cfg = cfgdummy;
        data.time(1,:) = time_sel(:);
        data.label = label;
        data.dimord = fieldtrip_example.M_timb_lock_gradplanarComb.dimord;
        data.grad = fieldtrip_example.M_timb_lock_gradplanarComb.grad;
        data.avg = avg;
        %creating the cfg for the actual plotting (magnetometers)
        cfg = [];
        cfg.layout = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/neuromag306mag.lay';
        cfg.colorbar = 'yes';
        cfg.xlim = signt; %set temporal limits (in seconds)
        cfg.zlim = [];
        cfg.colormap = 'jet';
        figure
        ft_topoplotER(cfg,data);
        set(gcf,'Color','w')
        %colormap with white for 0 values
        x = []; x.bottom = [0 0 0.5]; x.botmiddle = [0 0.5 1]; x.middle = [1 1 1]; x.topmiddle = [1 0 0]; x.top = [0.6 0 0]; %red - blue
        colormap(bluewhitered_PD(0,x))
    end
    %plotting temoral generalization
    ddTGm = mean(ddTG,3);
    figure; imagesc(time_sel,time_sel,ddTGm); xlabel('Train time (s)'); ylabel('Test time (s)'); set(gca,'YDir','normal');
    colorbar
    set(gcf,'Color','w')
    %reshaping temporal generalization decoding acuracy for submitting it to statistical testing
    ddTG2 = ddTG(1:776,1:776,:) - 50; %subtracting 50 since the function tests significance against 0 (here represented by 50% chance level)
    dat = cell(1,size(ddTG2,3));
    for ii = 1:size(ddTG2,3)
        dat(1,ii) = {ddTG2(:,:,ii)};
    end
    if tempgen_l == 1
        %computing statistics and plotting results
        stat = pl_permtest(dat,'alpha',0.05); %actual function
        STAT{T_cond} = stat; %storing statistics to compare images later and decide the limits for the colorbars
        time_selj = time_sel(1:776);
        figure; imagesc(time_selj,time_selj,stat.statmap); xlabel('Train time (s)'); ylabel('Test time (s)'); set(gca,'YDir','normal');
        colorbar
        set(gcf,'Color','w')
        caxis([-5 10])
        %testing contrast results (correcting for multiple comparison) by performing Monte Carlo simulations
        P2 = zeros(776,776);
%         statt = stat.statmap;
        P2(stat.statmappv<0.05) = 1; %threshold t-val = 2.6 corresponding to p-val < 0.01 (obtained by dividing 0.05 by the 4 comparisons employed here)
        thresh = 0;
        permut = 1000;
        threshMC = 0.001;
        perm_max = 1;
        t1 = time_sel; t2 = t1;
        [ OUT ] = twoD_MCS_LBPD_D( P2, thresh, permut, threshMC, perm_max, t1 , t2 )
        %plotting only significant clusters
        bumbum = zeros(776,776);
        bumbum(OUT{1,8}==1) = stat.statmap(OUT{1,8}==1);
        figure; imagesc(time_selj,time_selj,bumbum); xlabel('Train time (s)'); ylabel('Test time (s)'); set(gca,'YDir','normal');
        colorbar
        set(gcf,'Color','w')
        caxis([-5 10])
        %exporting figures
        outdir = '/aux/MINDLAB2021_MEG-TempSeqAges/12_10_2023';
        export_fig([outdir '/TemporalGeneralization_OldvsNewT' num2str(T_cond) 'Testing_Training.png'])
        export_fig([outdir '/TemporalGeneralization_OldvsNewT' num2str(T_cond) 'Testing_Training.eps'])
        %exporting table
        PDn = cell2table(OUT(:,1:7)); %table
        writetable(PDn,[outdir '/TemporalGeneralization_OldvsNewT' num2str(T_cond) 'Testing_Training.xlsx'],'Sheet',1); %printing excel file
        %exporting data for NC source data
        outdir = '/aux/MINDLAB2021_MEG-TempSeqAges/09_11_2023';
        PDn = table(bumbum); %table
        writetable(PDn,[outdir '/Temporalgen_Cond' num2str(T_cond) '.xlsx'],'Sheet',1); %printing excel file
    end
else %plotting together time series of decoding accuracy for all 4 decoding analyses (i.e. Old vs NewT1, Old vs NewT2, Old vs NewT3, Old vs NewT4)
    figure
    %     COL{1} = 'b'; COL{2} = 'r'; COL{3} = 'k'; COL{4} = 'g'; COL{5} = 'c'; COL{6} = 'm'; %colors (TO BE MADE THE SAME AS FOR WAVEFORMS OF DIFFERENT CONDITIONS!!)
    color_line = colormap(lines(5)); %extracting some colours from a colormap
    color_line2 = color_line;
    color_line2(1,:) = color_line(2,:);
    color_line2(2,:) = color_line(1,:);
    color_line2(5,:) = [0.4 0.4 0.4];
    LOC = zeros(1,length(time_sel),length(listPD),4);
    for nn = 1:4 %over NewTX conditions
        datadir = [datadirdum '/Old_vs_NewT' num2str(nn)];
        listPD = dir([datadir '/PD*']);
        listGT = dir([datadir '/TG*']);
        %loading data
        dd1 = zeros(length(time_sel),length(listPD));
        ddpatt = zeros(306,length(time_sel),length(listPD));
        ddTG = zeros(length(time_sel),length(time_sel),length(listPD));
        cnt = 0;
        for ii = 1:length(listPD)
            cnt = cnt + 1;
            load([listPD(ii).folder '/' listPD(ii).name]);
            dd1(:,cnt) = d.d;
            LOC(1,:,cnt,nn) = d.d; %to be exported for local (on laptop) purposes (e.g. plotting)
            ddpatt(:,:,cnt) = d.pattern;
            load([listGT(ii).folder '/' listGT(ii).name])
            ddTG(:,:,cnt) = d.d;
            disp([num2str(nn) ' - ' num2str(ii)])
        end
        save('/aux/MINDLAB2021_MEG-TempSeqAges/12_10_2023/decpack.mat','LOC')
        %average over subjects of pairwise decoding accuracy time-series
        dd12 = mean(dd1,2);
        dd12std = std(dd1,0,2)./sqrt(size(dd1,2)); %standard errors
        %reshaping pairwise decoding acuracy time-series for submitting them to statistical testing
        dd2 = dd1 - 50; %subtracting 50 since the function tests significance against 0 (here represented by 50% chance level)
        dat = cell(1,size(dd2,2));
        for ii = 1:size(dd2,2)
            dat(1,ii) = {dd2(:,ii)'};
        end
        %plotting average over subjects (pairwise decoding accuracy)
        grid on;grid minor;
        plot(time_sel,dd12,'color',color_line2(nn+1,:),'linewidth',2,'DisplayName',['NewT' num2str(nn)]);
        hold on
        plot(time_sel,dd12 + dd12std,':','color',color_line2(nn+1,:),'linewidth',0.5,'DisplayName',['NewT' num2str(nn)]);
        hold on
        plot(time_sel,dd12 - dd12std,':','color',color_line2(nn+1,:),'linewidth',0.5,'DisplayName',['NewT' num2str(nn)]);
        xlabel('Time (s)');ylabel('Decoding accuracy (%)');
        hold on
    end
    legend('show')
    set(gcf,'Color','w')
    xlim([time_sel(1) time_sel(end)])
end
% %trick to export results in excel sheet.. I initialized dumbo in the command window
% dumbo(:,T_cond) = stat.FDR.criticalmap;
% dumbo = dumbo';
% dumbo2 = zeros(5,1126);
% dumbo2(1,1:1126) = time_sel;
% dumbo2(2:5,1:1126) = dumbo;
% outdir = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/Beam_abs_1_sens_1_freq_broadband/Sources_Decoding_MEG0211/MCS/FinalImagesToWorkBench/Statistics';
% xlswrite([outdir '/PairwiseDecoding.xlsx'],dumbo2)

%% 

%% *** MEG SENSORS ANALYSIS ***

%% Step 1: Getting average weights on the MEG channels for each of the significant time-window of the decoding and averaging them (Figure 2b)

%directory where decoding results are stored in
datadirdum = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Decoding';
%loading and reshpaping data outputted by AU server (SVM algorithm) for statistics
load([datadirdum '/time.mat']);

%computing and extracting significant time-points for decoding (all conditions)
% dumbum = [];
% avg = zeros(102,4); %magnetometers x comparisons in decoding (i.e. Old vs each category of New)
avg = zeros(306,4); %magnetometers x comparisons in decoding (i.e. Old vs each category of New)
for nn = 1:4 %over NewTX conditions
    dumbum = [];
    datadir = [datadirdum '/Old_vs_NewT' num2str(nn)];
    listPD = dir([datadir '/PD*']);
    listGT = dir([datadir '/TG*']);
    %loading data
    dd1 = zeros(length(time_sel),length(listPD));
    ddpatt = zeros(306,length(time_sel),length(listPD));
    ddTG = zeros(length(time_sel),length(time_sel),length(listPD));
    cnt = 0;
    for ii = 1:length(listPD)
        cnt = cnt + 1;
        load([listPD(ii).folder '/' listPD(ii).name]);
        dd1(:,cnt) = d.d;
        ddpatt(:,:,cnt) = d.pattern;
        load([listGT(ii).folder '/' listGT(ii).name])
        ddTG(:,:,cnt) = d.d;
        disp([num2str(nn) ' - ' num2str(ii)])
    end
    %average over subjects of pairwise decoding accuracy time-series
    dd12 = mean(dd1,2);
    dd12std = std(dd1,0,2)./sqrt(size(dd1,2)); %standard errors
    %reshaping pairwise decoding acuracy time-series for submitting them to statistical testing
    dd2 = dd1 - 50; %subtracting 50 since the function tests significance against 0 (here represented by 50% chance level)
    dat = cell(1,size(dd2,2));
    for ii = 1:size(dd2,2)
        dat(1,ii) = {dd2(:,ii)'};
    end
    stat = pl_permtest(dat,'alpha',0.05); %actual function
    %plotting p-values after FDR correction for multiple comparisons
    xx = stat.FDR.criticalmap;
    CCt = bwconncomp(xx); %getting significant 'islands'
    sgf = CCt.PixelIdxList;
    for ii = 1:length(sgf)
        dumbum = cat(1,dumbum,sgf{ii});
    end
    dp = mean(ddpatt,3); %average over subjects
%     avg(:,nn) = mean(dp(1:3:end,dumbum),2); %extracting magnetometers and significant time-points only and then averaging over time
    avg(:,nn) = mean(dp(:,dumbum),2); %extracting significant time-points only and then averaging over time
end
avg2 = mean(avg,2); %average over comparisons in decoding
%plotting topoplot (fieldtrip function)
%creating the mask for the data
fieldtrip_example = load('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/fieldmask.mat');
load('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/plot3D_LF.mat','label')
cfgdummy = fieldtrip_example.M_timb_lock_gradplanarComb.cfg;
cfgdummy.previous = [];
data = [];
data.cfg = cfgdummy;
data.time = 1;
data.label = label;
data.dimord = fieldtrip_example.M_timb_lock_gradplanarComb.dimord;
data.grad = fieldtrip_example.M_timb_lock_gradplanarComb.grad;
data.avg = avg2;
%creating the cfg for the actual plotting (magnetometers)
cfg = [];
cfg.layout = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/neuromag306all.lay';
cfg.colorbar = 'yes';
%     cfg.xlim = [1 1]; %set temporal limits (in seconds)
cfg.zlim = [-1100 1100];
cfg.colormap = 'jet';
figure
ft_topoplotER(cfg,data);
set(gcf,'Color','w')
%colormap with white for 0 values
x = []; x.bottom = [0 0 0.5]; x.botmiddle = [0 0.5 1]; x.middle = [1 1 1]; x.topmiddle = [1 0 0]; x.top = [0.6 0 0]; %red - blue
colormap(bluewhitered_PD(0,x))

%only one-channel topoplot
fieldtrip_example = load('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/fieldmask.mat');
load('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/plot3D_LF.mat','label')
cfgdummy = fieldtrip_example.M_timb_lock_gradplanarComb.cfg;
cfgdummy.previous = [];
data = [];
data.cfg = cfgdummy;
data.time = 1;
data.label = label;
data.dimord = fieldtrip_example.M_timb_lock_gradplanarComb.dimord;
data.grad = fieldtrip_example.M_timb_lock_gradplanarComb.grad;
data.avg = zeros(306,1);
%creating the cfg for the actual plotting (magnetometers)
cfg = [];
cfg.layout = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/neuromag306all.lay';
cfg.colorbar = 'yes';
%     cfg.xlim = [1 1]; %set temporal limits (in seconds)
cfg.zlim = [];
cfg.colormap = 'jet';
figure
ft_topoplotER(cfg,data);
set(gcf,'Color','w')
%colormap with white for 0 values
x = []; x.bottom = [0 0 0.5]; x.botmiddle = [0 0.5 1]; x.middle = [1 1 1]; x.topmiddle = [1 0 0]; x.top = [0.6 0 0]; %red - blue
colormap(bluewhitered_PD(0,x))

%magnetometers only
fieldtrip_example = load('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/fieldmask.mat');
load('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/plot3D_LF.mat','label')
cfgdummy = fieldtrip_example.M_timb_lock_gradplanarComb.cfg;
cfgdummy.previous = [];
data = [];
data.cfg = cfgdummy;
data.time = 1;
data.label = label(1:3:end);
data.dimord = fieldtrip_example.M_timb_lock_gradplanarComb.dimord;
data.grad = fieldtrip_example.M_timb_lock_gradplanarComb.grad;
data.avg = avg2(1:3:end);
cfg = [];% label2 = fieldtrip_example.M_timb_lock_gradplanarComb.label;
% label = label2(103:end);
% %setting labels according to the ones in the layout
% for ii = 1:102
%     label{ii,1} = [label{ii,1}(1:3) label{ii,1}(5:end)];
% end
cfg.layout = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/neuromag306mag.lay';
cfg.colorbar = 'yes';
cfg.colormap = 'jet';
cfg.zlim = [-1100 1100];
figure
ft_topoplotER(cfg,data);
set(gcf,'Color','w')
%colormap with white for 0 values
x = []; x.bottom = [0 0 0.5]; x.botmiddle = [0 0.5 1]; x.middle = [1 1 1]; x.topmiddle = [1 0 0]; x.top = [0.6 0 0]; %red - blue
colormap(bluewhitered_PD(0,x))
%exporting data for NC source data
bum = zeros(204,2); bum(1:102,1) = avg2(1:3:end); %sotring magnetometers

%gradiometers only
avggrad = avg2; avggrad(1:3:end) = [];
labelg = label; labelg(1:3:end) = [];
data.label = labelg;
cfg = [];
cfg.layout = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/neuromag306planar.lay';
cfg.colorbar = 'yes';
cfg.colormap = 'jet';
cfg.zlim = [-22 22];
data.avg = avggrad;
figure
ft_topoplotER(cfg,data);
set(gcf,'Color','w')
%colormap with white for 0 values
x = []; x.bottom = [0 0 0.5]; x.botmiddle = [0 0.5 1]; x.middle = [1 1 1]; x.topmiddle = [1 0 0]; x.top = [0.6 0 0]; %red - blue
colormap(bluewhitered_PD(0,x))
%exporting data for NC source data
bum(:,2) = data.avg; %sotring gradiometers
outdir = '/aux/MINDLAB2021_MEG-TempSeqAges/09_11_2023';
PDn = table(bum); %table
writetable(PDn,[outdir '/DecodingWeights.xlsx'],'Sheet',1); %printing excel file

%% Step 2: Final averaged weights for the decoding and extracting the of highest values in absolute terms (higher than the average plus one standard deviation), independently for magnetometers and gradiometers (Figure 2d)

%mean + 1 standard deviation solution (formally better)
%magnetometers
avg2mag = avg2(1:3:end);
ave = mean(abs(avg2mag));
stand = std((abs(avg2mag)),0,1); 
thr = ave + stand;
[sorted , b] = sort(abs(avg2(1:3:end)),'descend');
bmag = b(sorted>thr);
labelmag = label(1:3:end); labelmag = labelmag(b); %sorting labels according to magnetometers value
labelmag = labelmag(sorted>thr); %showing main channels contributing to the decoding output
%gradiometers
avg2grad = avg2; avg2grad(1:3:end) = []; %erasing magnetometers
labelgrad = label; labelgrad(1:3:end) = []; %erasing magnetometers label
ave = mean(abs(avg2grad));
stand = std((abs(avg2grad)),0,1); 
thr = ave + stand;
[sorted , b] = sort(abs(avg2grad),'descend');
bgrad = b(sorted>thr); labelgrad = labelgrad(b);
labelgrad = labelgrad(sorted>thr); %showing main channels contributing to the decoding output

save('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/Papers/Systematic_Variations_DCM/NatureCommunications_Revision/Revision_I/MEGSensors_FunctionalParcellation/Decoding_Weights_MEGchannels.mat','labelmag','labelgrad','bmag','bgrad');

%% Step 3a - magnetometers: They shows two main clusters (one per hemisphere) which are analysed independently, as follows (only magnetometers here, gradiometers are below, after the source reconstruction) (Figure 2d)

%loading previously computed labels of main contributing MEG channels
load('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/Papers/Systematic_Variations_DCM/NatureCommunications_Revision/Revision_I/MEGSensors_FunctionalParcellation/Decoding_Weights_MEGchannels.mat');
%loading already pre-processed MEG sensors (magnetometers and combined planar gradiometers) 
load('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/MEG_sensors/Block_3.mat');

%selecting MEG channels for left (1) and right (2) hemisphere
%magnetometers
bmag2{1} = [3 5 2 12 8 55 1]; bmag2{2} = [48 54 51 44 49 47 99 50]; %hard coding for selecting the proper channel indices
% 1)magnetometers
data_mat2 = data_mat(1:2:end,:,:,:);
dum2 = zeros(2,size(data_mat2,2),size(data_mat2,3),size(data_mat2,4));
for ii = 1:2
    dum2(ii,:,:,:) = mean(data_mat2(bmag2{ii},:,:,:),1);
end
%t-tests
P = zeros(size(dum2,1),size(dum2,2),(size(dum2,4)-1)); %clusters x time-points x contrasts (every NewTX versus Old)
T = zeros(size(dum2,1),size(dum2,2),(size(dum2,4)-1)); %clusters x time-points x contrasts (every NewTX versus Old)
for ii = 1:size(dum2,1) %over clusters (2)
    for jj = 1:size(dum2,2) %ove time-points
        for cc = 1:(size(dum2,4)-1) %over the 4 NewTX
            %codes t-test
            [~,p,~,stats] = ttest(squeeze(dum2(ii,jj,:,1)),squeeze(dum2(ii,jj,:,cc+1))); %contrasting cond1 (Old) with cond X (depending on vectcond it can be either NewT1 or NewT2 or NewT3 or NewT4
            P(ii,jj,cc) = p;
            T(ii,jj,cc) = stats.tstat;
        end
        disp([num2str(ii) ' - ' num2str(jj)])
    end
end
%MCS
p_thresh = 0.05; %threshold for binarising p-values vector
load('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/time_normal.mat'); %loading time
outdir = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/Papers/Systematic_Variations_DCM/NatureCommunications_Revision/Revision_I/MEGSensors_FunctionalParcellation/MEG_sensors_stats';
mkdir(outdir);
SIGN = cell(size(dum2,1),(size(dum2,4)-1));
for ii = 1:size(dum2,1) %over clusters
    for cc = 1:(size(dum2,4)-1) %over the 4 NewTX
        Pbin = zeros(1,size(dum2,2));
        Pbin(P(ii,:,cc)<p_thresh) = 1;
        tvals = T(ii,:,cc);
        [ sign_clust ] = oneD_MCS_LBPD_D( Pbin, 0, 1, 1000, 0.001, time(1:size(dum2,2)), tvals ) %Monte Carlo simulation function to correct for multiple comparisons
        PDn = cell2table(sign_clust); %table
        writetable(PDn,[outdir '/Mag_Hem_' num2str(ii) '_OldvsNewT' num2str(cc) '.xlsx'],'Sheet',1); %printing excel file
        SIGN{ii,cc} = sign_clust;
    end
end
save('/aux/MINDLAB2021_MEG-TempSeqAges/14_09_2023/Sign_2clustMEGsensors_MAG.mat','SIGN');
%plotting
%defining colors
color_line = colormap(lines(5)); %extracting some colours from a colormap
color_line2 = color_line;
color_line2(1,:) = color_line(2,:);
color_line2(2,:) = color_line(1,:);
color_line2(5,:) = [0.4 0.4 0.4];
conds{1} = ' Mem '; conds{2} = ' NewT1 '; conds{3} = ' NewT2 '; conds{4} = ' NewT3 '; conds{5} = ' NewT4 ';
load('/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/Source_LBPD/time_345.mat')
% COL{1} = 'b'; COL{2} = 'r'; COL{3} = 'k'; COL{4} = 'g'; COL{5} = 'c'; COL{6} = 'm'; %colors
ROIN{1} = 'Mag LH'; ROIN{2} = 'Mag RH';
% load('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/time_normal.mat')
for ii = 1:length(ROIN) %over clusters
    figure
    for cc = 1:length(conds) %over conditions
        plot(squeeze(time(1:size(dum2,2))),squeeze(nanmean(dum2(ii,:,:,cc),3)),'Color',color_line2(cc,:),'LineWidth',2,'DisplayName',[conds{cc}])
        hold on
        plot(squeeze(time(1:size(dum2,2))),squeeze(nanmean(dum2(ii,:,:,cc),3)) + (nanstd(dum2(ii,:,:,cc),0,3)./sqrt(size(dum2,3))),':','Color',color_line2(cc,:),'LineWidth',0.5,'HandleVisibility','off')
        hold on
        plot(squeeze(time(1:size(dum2,2))),squeeze(nanmean(dum2(ii,:,:,cc),3)) - (nanstd(dum2(ii,:,:,cc),0,3)./sqrt(size(dum2,3))),':','Color',color_line2(cc,:),'LineWidth',0.5,'HandleVisibility','off')
        hold on
    end
    grid minor
    title(ROIN{ii})
    xlim([-0.1 3.4])
%     ylim(limmy)
    set(gcf,'color','w')
    legend('show')
end

%identifying the peaks of positive and negative activity for all conditions within the time-window of response to each musical tone
%left hemisphere
%average over subjects
duml = squeeze(nanmean(dum2(1,:,:,:),3));
%getting maximum values for each condition, and then overall maximum value
[maxv, imax] = max(duml(139:238,:)); [maxv2, imax2] = max(maxv); disp(['max, tone 2, left H, time sec : ' num2str(time_sel(imax(imax2)+139)) ' - time sam : ' num2str(imax(imax2)+139)]) %getting maximum value
%minimum value (discard the names of the variables..)
[maxv, imax] = min(duml(139:238,:)); [maxv2, imax2] = min(maxv); disp(['min, tone 2, left H, time sec : ' num2str(time_sel(imax(imax2)+139)) ' - time sam : ' num2str(imax(imax2)+139)]) %getting maximum value
%tone 3
%getting maximum values for each condition, and then overall maximum value
[maxv, imax] = max(duml(239:338,:)); [maxv2, imax2] = max(maxv); disp(['max, tone 3, left H, time sec : ' num2str(time_sel(imax(imax2)+239)) ' - time sam : ' num2str(imax(imax2)+239)]) %getting maximum value
%minimum value (discard the names of the variables..)
[maxv, imax] = min(duml(239:338,:)); [maxv2, imax2] = min(maxv); disp(['min, tone 3, left H, time sec : ' num2str(time_sel(imax(imax2)+239)) ' - time sam : ' num2str(imax(imax2)+239)]) %getting maximum value
%tone 4
%getting maximum values for each condition, and then overall maximum value
[maxv, imax] = max(duml(339:438,:)); [maxv2, imax2] = max(maxv); disp(['max, tone 4, left H, time sec : ' num2str(time_sel(imax(imax2)+339)) ' - time sam : ' num2str(imax(imax2)+339)]) %getting maximum value
%minimum value (discard the names of the variables..)
[maxv, imax] = min(duml(339:438,:)); [maxv2, imax2] = min(maxv); disp(['min, tone 4, left H, time sec : ' num2str(time_sel(imax(imax2)+339)) ' - time sam : ' num2str(imax(imax2)+339)]) %getting maximum value
%tone 5
%getting maximum values for each condition, and then overall maximum value
[maxv, imax] = max(duml(439:538,:)); [maxv2, imax2] = max(maxv); disp(['max, tone 5, left H, time sec : ' num2str(time_sel(imax(imax2)+439)) ' - time sam : ' num2str(imax(imax2)+439)]) %getting maximum value
%minimum value (discard the names of the variables..)
[maxv, imax] = min(duml(439:538,:)); [maxv2, imax2] = min(maxv); disp(['min, tone 5, left H, time sec : ' num2str(time_sel(imax(imax2)+439)) ' - time sam : ' num2str(imax(imax2)+439)]) %getting maximum value

%right hemisphere
%average over subjects
duml = squeeze(nanmean(dum2(2,:,:,:),3));
%tone 2
%getting maximum values for each condition, and then overall maximum value
[maxv, imax] = max(duml(139:238,:)); [maxv2, imax2] = max(maxv); disp(['max, tone 2, right H, time sec : ' num2str(time_sel(imax(imax2)+139)) ' - time sam : ' num2str(imax(imax2)+139)]) %getting maximum value
%minimum value (discard the names of the variables..)
[maxv, imax] = min(duml(139:238,:)); [maxv2, imax2] = min(maxv); disp(['min, tone 2, right H, time sec : ' num2str(time_sel(imax(imax2)+139)) ' - time sam : ' num2str(imax(imax2)+139)]) %getting maximum value
%tone 3
%getting maximum values for each condition, and then overall maximum value
[maxv, imax] = max(duml(239:338,:)); [maxv2, imax2] = max(maxv); disp(['max, tone 3, right H, time sec : ' num2str(time_sel(imax(imax2)+239)) ' - time sam : ' num2str(imax(imax2)+239)]) %getting maximum value
%minimum value (discard the names of the variables..)
[maxv, imax] = min(duml(239:338,:)); [maxv2, imax2] = min(maxv); disp(['min, tone 3, right H, time sec : ' num2str(time_sel(imax(imax2)+239)) ' - time sam : ' num2str(imax(imax2)+239)]) %getting maximum value
%tone 4
%getting maximum values for each condition, and then overall maximum value
[maxv, imax] = max(duml(339:438,:)); [maxv2, imax2] = max(maxv); disp(['max, tone 4, right H, time sec : ' num2str(time_sel(imax(imax2)+339)) ' - time sam : ' num2str(imax(imax2)+339)]) %getting maximum value
%minimum value (discard the names of the variables..)
[maxv, imax] = min(duml(339:438,:)); [maxv2, imax2] = min(maxv); disp(['min, tone 4, right H, time sec : ' num2str(time_sel(imax(imax2)+339)) ' - time sam : ' num2str(imax(imax2)+339)]) %getting maximum value
%tone 5
%getting maximum values for each condition, and then overall maximum value
[maxv, imax] = max(duml(439:538,:)); [maxv2, imax2] = max(maxv); disp(['max, tone 5, right H, time sec : ' num2str(time_sel(imax(imax2)+439)) ' - time sam : ' num2str(imax(imax2)+439)]) %getting maximum value
%minimum value (discard the names of the variables..)
[maxv, imax] = min(duml(439:538,:)); [maxv2, imax2] = min(maxv); disp(['min, tone 5, right H, time sec : ' num2str(time_sel(imax(imax2)+439)) ' - time sam : ' num2str(imax(imax2)+439)]) %getting maximum value

%% Step 3b - digression - musical expertise and neural responses

% OBS!! to run this cell you first need to run the cell above (so you can compute dum2, which is the matrix with the neural data) 

[~,~,raw] = xlsread('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/Papers/Systematic_Variations_DCM/NatureCommunications_Revision/Revision_I/MusicalTrainingPracticeWM.xlsx');
beha = cell2mat(raw(2:end,3)); %extracting inforation about musical training
vect = [203 179 299 274 389 361 477 460; 181 234 273 298 364 390 459 506]; %time-samples of peaks (cluster 1 and cluster 2)
vectcond = [2,2,2,3,2,4,2,5;2,2,3,2,4,2,5,2]; %condition X to be contrasted against condition 1 (Old)
rho = zeros(size(vect,1),size(vect,2));
p = zeros(size(vect,1),size(vect,2));
for cc = 1:2 %over clusters
   for ii = 1:size(vect,2) %oer peaks     
       datadum = squeeze(mean(dum2(cc,vect(cc,ii)-5:vect(cc,ii)+5,~isnan(beha),vectcond(cc,ii)),2)); %extracting data: cluster 1 or 2, peak plus/minus 20ms, subjects with behavioral information about musical training, condition to be associated to the peaks
       [rho(cc,ii),p(cc,ii)] = corr(datadum,beha(~isnan(beha))); %correlation
   end
end
%FDR correction
[pthr,pcor,padj] = fdr(p(1,:))
[pthr,pcor,padj] = fdr(p(2,:))

%% Step 3c - gradiometers (same analysis as for the magnetometers but without source reconstruction) (Figure 2e)

%NOT USING COMBINED PLANAR GRADIOMETERS

epoch_list = dir ('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/mes*recogminor*mat'); %dir to epoched files (encoding)
% bgrad2{1} = [11 13 17 23 4 10 104];
% bgrad2{2} = [108 98 94 198 100 88 184 96 86 2 110 15]; %hard coding for selecting the proper channel indices
bgrad2{1} = [10 13 17 23 5 11 104 85 80 26 102 90 20 198]; %hard coding for selecting the proper channel indices (polarity cluster 1)
bgrad2{2} = [107 96 93 197 100 183 99 2 110 15 7 8 117 154 87 4 97 101]; %hard coding for selecting the proper channel indices (polarity cluster 2)


conditions = {'Old_Correct','New_T1_Correct','New_T2_Correct','New_T3_Correct','New_T4_Correct'};
dataa = zeros(2,1126,length(conditions),length(epoch_list));
% dataa = zeros(5,1126,length(conditions),length(epoch_list));
% azz = [1:306];
for iii = 1:length(epoch_list)%1:length(epoch_list) %over epoched files 
    D = spm_eeg_load([epoch_list(iii).folder '/' epoch_list(iii).name]); %loading data
    %getting conditions index with the right order (according to variable 'conditions') 
    conds = zeros(1,length(conditions));
    for ii = 1:length(conditions)
        conds(ii) = sum(find(strcmp((D.conditions),conditions{ii}))); %getting the condition numbers ordered according to the specification of the user in S.conditions (sum to avoid the problem of zero-vector occurring if it does not find the match between the specified label and the data label.. of course if you specify two times the same label in the structure S it crashes as well.. but you would need to be quite peculiar to do that so I assume you will not do that..)
    end
    %extracting only
    chans = D.chanlabels;
    idxMEGc1 = find(strcmp(chans,'MEG0111')); %getting extreme channels indexes (mag) (the first MEG channel may have a different index..)
    dum = D(idxMEGc1:305+idxMEGc1,:,conds); %extracting MEG channels only
    dum(1:3:end,:,:) = []; %removing magnetometers
    dataa(1,:,:,iii) = nanmean(dum(bgrad2{1},:,:),1); %cluster 1 (left hemisphere)
    dataa(2,:,:,iii) = nanmean(dum(bgrad2{2},:,:),1); %cluster 2 (right hemisphere)
% for pp= 1:306%length(azz)
%     dataa(pp,:,:,iii) = dum(azz(pp),:,:); %cluster 1 (left hemisphere)
% end
    disp(iii)
end
%t-tests
P = zeros(size(dataa,1),size(dataa,2),(size(dataa,3)-1)); %clusters x time-points x contrasts (every NewTX versus Old)
T = zeros(size(dataa,1),size(dataa,2),(size(dataa,3)-1)); %clusters x time-points x contrasts (every NewTX versus Old)
for ii = 1:size(dataa,1) %over clusters (2)
    for jj = 1:size(dataa,2) %ove time-points
        for cc = 1:(size(dataa,3)-1) %over the 4 NewTX
            %codes t-test
            [~,p,~,stats] = ttest(squeeze(dataa(ii,jj,1,:)),squeeze(dataa(ii,jj,cc+1,:))); %contrasting cond1 (Old) with cond X (depending on vectcond it can be either NewT1 or NewT2 or NewT3 or NewT4
            P(ii,jj,cc) = p;
            T(ii,jj,cc) = stats.tstat;
        end
        disp([num2str(ii) ' - ' num2str(jj)])
    end
end
%MCS
p_thresh = 0.05; %threshold for binarising p-values vector
load('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/time_normal.mat'); %loading time
outdir = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/Papers/Systematic_Variations_DCM/NatureCommunications_Revision/Revision_I/MEGSensors_FunctionalParcellation/MEG_sensors_stats';
mkdir(outdir);
SIGN = cell(size(dataa,1),(size(dataa,3)-1));
for ii = 1:size(dataa,1) %over clusters
    for cc = 1:(size(dataa,3)-1) %over the 4 NewTX
        Pbin = zeros(1,size(dataa,2));
        Pbin(P(ii,:,cc)<p_thresh) = 1;
        tvals = T(ii,:,cc);
        [ sign_clust ] = oneD_MCS_LBPD_D( Pbin, 0, 1, 1000, 0.001, time(1:size(dataa,2)), tvals ) %Monte Carlo simulation function to correct for multiple comparisons
        PDn = cell2table(sign_clust); %table
        writetable(PDn,[outdir '/Grad_Clust_' num2str(ii) '_OldvsNewT' num2str(cc) '.xlsx'],'Sheet',1); %printing excel file
        SIGN{ii,cc} = sign_clust;
    end
end
save('/aux/MINDLAB2021_MEG-TempSeqAges/14_09_2023/Sign_2clustMEGsensors_GRAD.mat','SIGN');
%plotting
%defining colors
color_line = colormap(lines(5)); %extracting some colours from a colormap
color_line2 = color_line;
color_line2(1,:) = color_line(2,:);
color_line2(2,:) = color_line(1,:);
color_line2(5,:) = [0.4 0.4 0.4];
clear conds
conds{1} = ' Mem '; conds{2} = ' NewT1 '; conds{3} = ' NewT2 '; conds{4} = ' NewT3 '; conds{5} = ' NewT4 ';
load('/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/Source_LBPD/time_345.mat')
ROIN{1} = 'Grad Cluster 1'; ROIN{2} = 'Grad Cluster 2';
% ROIN = 1:length(azz)
for ii = 1:length(ROIN) %over clusters
    figure
    for cc = 1:length(conds) %over conditions
        plot(squeeze(time(1:size(dataa,2))),squeeze(nanmean(dataa(ii,:,cc,:),4)),'Color',color_line2(cc,:),'LineWidth',2,'DisplayName',[conds{cc}])
        hold on
        plot(squeeze(time(1:size(dataa,2))),squeeze(nanmean(dataa(ii,:,cc,:),4)) + (nanstd(dataa(ii,:,cc,:),0,4)./sqrt(size(dataa,4))),':','Color',color_line2(cc,:),'LineWidth',0.5,'HandleVisibility','off')
        hold on
        plot(squeeze(time(1:size(dataa,2))),squeeze(nanmean(dataa(ii,:,cc,:),4)) - (nanstd(dataa(ii,:,cc,:),0,4)./sqrt(size(dataa,4))),':','Color',color_line2(cc,:),'LineWidth',0.5,'HandleVisibility','off')
        hold on
    end
    grid minor
    title(ROIN{ii})
% title(num2str(azz(ii)))
    xlim([-0.1 3.4])
%     ylim(limmy)
    set(gcf,'color','w')
    legend('show')
end

%% Step 3d - Digression - Topoplots for the peaks of differential activity between Old and New conditions (for supplementary material) (Figure S6)

%% barbaric (manual) solution for exporting the topoplots.. automatic solutions could not work.. anyway, this is correct

export_fig(['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/Papers/Systematic_Variations_DCM/NatureCommunications_Revision/Revision_I/Images_Data/' names(pp,:) '_cond' num2str(cc) '_mag.eps'])
export_fig(['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/Papers/Systematic_Variations_DCM/NatureCommunications_Revision/Revision_I/Images_Data/' names(pp,:) '_cond' num2str(cc) '_mag.png'])
disp('barabba')

%%

export_fig(['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/Papers/Systematic_Variations_DCM/NatureCommunications_Revision/Revision_I/Images_Data/' names(pp,:) '_cond' num2str(cc) '_grad.eps'])
export_fig(['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/Papers/Systematic_Variations_DCM/NatureCommunications_Revision/Revision_I/Images_Data/' names(pp,:) '_cond' num2str(cc) '_grad.png'])
disp('barabba')

%% actual computations

peaks = [203 179 299 274 389 361 477 460 181 234 273 298 364 390 459 506];
names = ['L_t2_max'; 'L_t2_min'; 'L_t3_max'; 'L_t3_min'; 'L_t4_max'; 'L_t4_min'; 'L_t5_max'; 'L_t5_min'; 'R_t2_max'; 'R_t2_min'; 'R_t3_max'; 'R_t3_min'; 'R_t4_max'; 'R_t4_min'; 'R_t5_max'; 'R_t5_min';'N100____'];
load_data = 1;
save_data = 0;
datadirdum = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Decoding';
%loading and reshpaping data outputted by AU server (SVM algorithm) for statistics
load([datadirdum '/time.mat'])
vectcond = [2,2,2,3,2,4,2,5,2,2,3,2,4,2,5,2,1]; %condition X to be contrasted against condition 1 (Old)
gradl = [1 3.5];
magl = [-110 110];
export_l = 1;
jkl = 4;
ccc = 1;
if export_l == 0
    gradl = [];
    magl = [];
end

for pp = jkl%length(peaks) + 1
    for cc = ccc
        if export_l == 1
            close all
        end
        if cc == 1
            condd = 1;
        else
            condd = vectcond(pp);
        end
        block = 3; % 3 = recogminor; 4 = recogspeed; 5 = samemel; 6 = visualpat; 7 = pd
        channels_plot = []; % empty for plotting single channels; otherwise number(s) of channels to be averaged and plotted (e.g. [13] or [13 18])
        if pp == length(peaks)+1
            barb = [0.1 0.102]; %N100
        else
            barb = [time_sel(peaks(pp)) - 0.020 time_sel(peaks(pp)) + 0.020]; %all other peaks
        end
        S = [];
        %computing data
        if block == 3
            S.conditions = {'Old_Correct','New_T1_Correct','New_T2_Correct','New_T3_Correct','New_T4_Correct'};
            %     S.conditions = {'Old_Correct','New_T1_Correct'};
        elseif block == 4
            % S.conditions = {'Old_Fast_Correct','Old_Slow_Correct','New_Fast_Correct','New_Slow_Correct'};
            %     S.conditions = {'Old_Fast_Correct','New_Fast_Correct'};
            S.conditions = {'Old_Slow_Correct','New_Slow_Correct'};
        elseif block == 5
            S.conditions = {'Old_Correct','New_Correct'};
            %     S.conditions = {'Old_Correct','Encoding'};
        elseif block == 6
            S.conditions = {'Old_Correct','New_Correct'};
            %     S.conditions = {'Old_Correct','Encoding'};
        elseif block == 7
            S.conditions = {'pd','pm'};
        end
        if block == 3
            list = dir ('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/P*_recogminor*_tsssdsm.mat'); %dir to epoched files (encoding)
        elseif block == 4
            list = dir ('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/P*recogspeed*_tsssdsm.mat'); %dir to epoched files (encoding)
        elseif block == 5
            list = dir ('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/P*samemel*_tsssdsm.mat'); %dir to epoched files (encoding)
        elseif block == 6
            list = dir ('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/P*visualpat*_tsssdsm.mat'); %dir to epoched files (encoding)
        elseif block == 7
            list = dir ('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/P*pd*_tsssdsm.mat'); %dir to epoched files (encoding)
        end
        v = 1:length(list); %subjects
        % v = [2];
        if ~exist('chanlabels','var')
            load('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter/MEG_sensors/recogminor_all_conditions.mat', 'chanlabels')
        end
        outdir = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/MEG_sensors'; %path to write output in
        S.outdir = outdir;
        S.data = [];
        if load_data == 1 %if you already computed and saved on disk the t-tests you can load them here
            if ~exist('data_mat','var')
                load([outdir '/Block_' num2str(block) '.mat']);
            end
            S.data = data_mat(:,:,v,:);
            S.chanlabels = chanlabels;
            S.time_real = time_sel;
        else %otherwise you can extract the data from SPM MEEG objects (one for each subject)
            %     S.spm_list = cell(1,length(list));
            % v = 7;
            S.spm_list = cell(1,length(v));
            for ii = 1:length(v)
                S.spm_list(ii) = {[list(v(ii)).folder '/' list(v(ii)).name]};
            end
        end
        S.timeextract = []; %time-points to be extracted
        S.centerdata0 = 0; %1 to make data starting at 0
        S.save_data = save_data; %only meaningfull if you read data from SPM objects saved on disk
        S.save_name_data = ['Block_' num2str(block)];
        %individual waveform plotting
        if isempty(channels_plot)
            S.waveform_singlechannels_label = 1; %1 to plot single channel waveforms
        else
            S.waveform_singlechannels_label = 0; %1 to plot single channel waveforms
        end
        S.wave_plot_conditions_together = 0; %1 for plotting the average of all
        S.mag_lab = 2; %1 for magnetometers; 2 for gradiometers
        S.x_lim_temp_wave = [-0.1 3.4]; %limits for time (in secs) (E.g. [-0.1 3.4])
        S.y_lim_ampl_wave = []; %limit for amplitude (E.g. [0 120] magnetometes, [0 6] gradiometers)
        %averaged waveform plotting
        if isempty(channels_plot)
            S.waveform_average_label = 0; %average of some channels
            S.left_mag = 95; %13 %37 (visual) %43 (visual) %199 (visual) %203 (visual) %channels for averaging
        else
            S.waveform_average_label = 1; %average of some channels
            S.left_mag = channels_plot; %13 %37 (visual) %43 (visual) %199 (visual) %203 (visual) %channels for averaging
        end
        % S.left_mag = [2:2:204];
        S.legc = 1; %set 1 for legend
        % S.left_mag = 99;
        S.signtp = {[]};
        % S.sr = 150; %sampling rate (Hz)
        S.avewave_contrast = 0; %1 to plot the contrast between conditions (averaged waveform)
        S.save_label_waveaverage = 0;
        S.label_plot = 'c';
        %t-tests
        S.t_test_for_permutations = 0;
        S.cond_ttests_tobeplotted_topoplot = [1 2]; %this is for both topoplot and t-tests!! (here [1 2] means cond1 vs cond2!!!!!!!)
        %topoplotting
        S.topoplot_label = 1;
        S.fieldtrip_mask = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External';
        S.topocontr = 0;
        S.topocondsing = [condd]; %condition for topoplot
        % S.xlim = [0.75 0.85]; %time topolot
        % S.xlim = [1.1 1.2]; %time topolot
        S.xlim = barb;
        S.zlimmag = magl; %magnetometers amplitude topoplot limits
        S.zlimgrad = gradl; %gradiometers amplitude topoplot limits
        S.colormap_spec = 0;
        % x = []; x.bottom = [0 0 1]; x.botmiddle = [0 0.5 1]; x.middle = [1 1 1]; x.topmiddle = [1 1 0.5]; x.top = [1 0.95 0]; %yellow - blue
        x = []; x.bottom = [0 0 0.5]; x.botmiddle = [0 0.5 1]; x.middle = [1 1 1]; x.topmiddle = [1 0 0]; x.top = [0.6 0 0]; %red - blue
        S.colormap_spec_x = x;
        S.topoplot_save_label = 0;
        color_line = colormap(lines(5)); %extracting some colours from a colormap
        S.color_line = color_line;
        S.color_line(1,:) = color_line(2,:);
        S.color_line(2,:) = color_line(1,:);
        S.color_line(5,:) = [0.4 0.4 0.4];
        
        [out] = MEG_sensors_plotting_ttest_LBPD_D2(S);
%         if export_l == 1
%             figure(2)
%             export_fig(['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/Papers/Systematic_Variations_DCM/NatureCommunications_Revision/Revision_I/Images_Data/' names(pp,:) '_cond' num2str(cc) '_mag.eps'])
%             export_fig(['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/Papers/Systematic_Variations_DCM/NatureCommunications_Revision/Revision_I/Images_Data/' names(pp,:) '_cond' num2str(cc) '_mag.png'])
%             figure(3)
%             export_fig(['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/Papers/Systematic_Variations_DCM/NatureCommunications_Revision/Revision_I/Images_Data/' names(pp,:) '_cond' num2str(cc) '_grad.eps'])
%             export_fig(['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/Papers/Systematic_Variations_DCM/NatureCommunications_Revision/Revision_I/Images_Data/' names(pp,:) '_cond' num2str(cc) '_grad.png'])
%         end
%         disp('barabba')
    end
end

%%

%% *** SOURCE RECONSTRUCTION (LBPD) ***

%%

%% CREATING 8mm PARCELLATION FOR EASIER INSPECTION IN FSLEYES
%OBS!! This section is done only for better handling of some visualization purposes, but it does not affect any of the beamforming algorithm;
% it is just important not to mix up the MNI coordinates, thus I would recommend to use the following lines

%1) USE load_nii TO LOAD A PREVIOUS NIFTI IMAGE
imag_8mm = load_nii('/scratch7/MINDLAB2017_MEG-LearningBach/DTI_Portis/Templates/MNI152_T1_8mm_brain.nii.gz');
Minfo = size(imag_8mm.img); %get info about the size of the original image
M8 = zeros(Minfo(1), Minfo(2), Minfo(3)); %Initialize an empty matrix with the same dimensions as the original .nii image
cc = 0; %set a counter
M1 = imag_8mm.img;
for ii = 1:Minfo(1) %loop across each voxel of every dimension
    for jj = 1:Minfo(2)
        for zz = 1:Minfo(3)
            if M1(ii,jj,zz) ~= 0 %if we have an actual brain voxel
                cc = cc+1;
                M8(ii,jj,zz) = cc;
            end
        end
    end
end
%2) PUT YOUR MATRIX IN THE FIELD ".img"
imag_8mm.img = M8; %assign values to new matrix 
%3) SAVE NIFTI IMAGE USING save_nii
save_nii(imag_8mm,'/scratch7/MINDLAB2017_MEG-LearningBach/DTI_Portis/Templates/MNI152_8mm_brain_diy.nii.gz');
%4) USE FSLEYES TO LOOK AT THE FIGURE
%Create parcellation on the 8mm template
for ii = 1:3559 %for each 8mm voxel
    cmd = ['fslmaths /scratch7/MINDLAB2017_MEG-LearningBach/DTI_Portis/Templates/MNI152_8mm_brain_diy.nii.nii.gz -thr ' num2str(ii) ' -uthr ' num2str(ii) ' -bin /scratch7/MINDLAB2017_MEG-LearningBach/DTI_Portis/Templates/AAL_80mm_3559ROIs/' num2str(ii) '.nii.gz'];
    system(cmd)
    disp(ii)
end
%5) GET MNI COORDINATES OF THE NEW FIGURE AND SAVE THEM ON DISK
MNI8 = zeros(3559,3);
for mm = 1:3559 %over brain voxel
    path_8mm = ['/scratch7/MINDLAB2017_MEG-LearningBach/DTI_Portis/Templates/parcel_80mm_3559ROIs/' num2str(mm) '.nii.gz']; %path for each of the 3559 parcels
    [mni_coord,pkfo] = osl_mnimask2mnicoords(path_8mm);  %getting MNI coordinates
    MNI8(mm,:) = mni_coord; %storing MNI coordinates
end
%saving on disk
save('/scratch7/MINDLAB2017_MEG-LearningBach/DTI_Portis/Templates/MNI152_8mm_coord_dyi.mat', 'MNI8');

%% CONVERSION T1 - DICOM TO NIFTI

addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/osl/dicm2nii'); %adds path to the dcm2nii folder in osl
MRIsubj = dir('/projects/MINDLAB2020_MEG-AuditoryPatternRecognition/raw/0*');
MRIoutput = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/MRI_nifti';
MRIout_block{1} = 'Block_3'; MRIout_block{2} = 'Block_4'; MRIout_block{3} = 'Block_5'; MRIout_block{4} = 'Block_6'; MRIout_block{5} = 'Block_7';

for bb = 1:5 %over experimental blocks
    for ii = 1:length(MRIsubj) %over subjects
        asd = ['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/MRI_nifti/Block_' num2str(bb+2) '/' MRIsubj(ii).name];
        if ~exist(asd,'dir') %checking whether the directory exists
            mkdir(asd); %if not, creating it
        end
        if isempty(dir([asd '/*.nii'])) %if there are no nifti images.. I need to convert them
            flagg = 0;
            MRIMEGdate = dir([MRIsubj(ii).folder '/' MRIsubj(ii).name '/20*']);
            niiFolder = [MRIoutput '/' MRIout_block{bb} '/' MRIsubj(ii).name];
            for jj = 1:length(MRIMEGdate) %over dates of recording
                if ~isempty(dir([MRIMEGdate(jj).folder '/' MRIMEGdate(jj).name '/MR*'])) %if we get an MRI recording
                    MRI2 = dir([MRIMEGdate(jj).folder '/' MRIMEGdate(jj).name '/MR/*fatsat']); %looking for T1
                    if ~isempty(MRI2) %if we have it
                        flagg = 1; %determining that I could convert MRI T1
                        dcmSource = [MRI2(1).folder '/' MRI2(1).name '/files/'];
                        if ii ~= 68 || jj ~= 3 %this is because subject 0068 got two MRIs stored.. but the second one (indexed by jj = 3) is of another subject (0086); in this moment, subject 0086 is perfectly fine, but in subject 0068 there are still the two MRIs (for 0068 (jj = 2) and for 0086 (jj = 3))
                            dicm2nii(dcmSource, niiFolder, '.nii');
                        end
                    end
                end
            end
            if flagg == 0
                warning(['subject ' MRIsubj(ii).name ' has no MRI T1']);
            end
        end
        disp(ii)
    end
end

%% SETTING FOR CLUSTER (PARALLEL COMPUTING)

% clusterconfig('scheduler', 'none'); %If you do not want to submit to the cluster, but simply want to test the script on the hyades computer, you can instead of 'cluster', write 'none'
clusterconfig('scheduler', 'cluster'); %If you do not want to submit to the cluster, but simply want to test the script on the hyades computer, you can instead of 'cluster', write 'none'
clusterconfig('long_running', 1); % This is the cue we want to use for the clsuter. There are 3 different cues. Cue 0 is the short one, which should be enough for us
clusterconfig('slot', 1); %slot is memory, and 1 memory slot is about 8 GB. Hence, set to 2 = 16 GB
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/Cluster_ParallelComputing')

%% RHINO coregistration

%block to be run RHINO coregistrartion on
block = 6; % 3 = recogminor; 4 = recogspeed; 5 = samemel; 6 = visualpat; 7 = pd

if block == 3
    list = dir ('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/e*recogminor*_tsssdsm.mat'); %dir to epoched files (encoding)
elseif block == 4
    list = dir ('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/e*recogspeed*_tsssdsm.mat'); %dir to epoched files (encoding)
elseif block == 5
    list = dir ('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/e*samemel*_tsssdsm.mat'); %dir to epoched files (encoding)
elseif block == 6
    list = dir ('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/e*visualpat*_tsssdsm.mat'); %dir to epoched files (encoding)
elseif block == 7
    list = dir ('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/e*pd*_tsssdsm.mat'); %dir to epoched files (encoding)
end

%running rhino
%OBS! check that all MEG data are in the same order and number as MRI nifti files!
a = ['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/MRI_nifti/Block_' num2str(block)]; %set path to MRI subjects' folders
for ii = 35%1:length(list) %OBS! change this depending on atonal vs. major
    S = [];
    S.ii = ii;
    S.D = [list(ii).folder '/' list(ii).name]; %path to major files
    D = spm_eeg_load(S.D);
    if ~isfield(D,'inv') %checking if the coregistration was already run
        dummyname = D.fname;
        if 7 == exist([a '/' dummyname(13:16)],'dir') %if you have the MRI folder
            dummymri = dir([a '/' dummyname(13:16) '/*.nii']); %path to nifti files (ending with .nii)
            if ~isempty(dummymri)
                S.mri = [dummymri(1).folder '/' dummymri(1).name];
                %standard parameters
                S.useheadshape = 1;
                S.use_rhino = 1; %set 1 for rhino, 0 for no rhino
                %         S.forward_meg = 'MEG Local Spheres';
                S.forward_meg = 'Single Shell'; %CHECK WHY IT SEEMS TO WORK ONLY WITH SINGLE SHELL!!
                S.fid.label.nasion = 'Nasion';
                S.fid.label.lpa = 'LPA';
                S.fid.label.rpa = 'RPA';
                jobid = job2cluster(@coregfunc,S); %running with parallel computing
            else
                warning(['subject ' dummyname(13:16) ' does not have the MRI'])
            end
        end
    else
        if isempty(D.inv{1}) %checking whether the coregistration was run but now it is empty..
            warning(['subject ' D.fname ' has an empty rhino..']);
        end
    end
    disp(ii)
end

%% checking (or copying) RHINO

copy_label = 0; % 1 = pasting inv RHINO from epoched data (where it was computed) to continuous data; 0 = simply showing RHINO coregistration
%block to be run RHINO coregistration on
block = 3; % 3 = recogminor; 4 = recogspeed; 5 = samemel; 6 = visualpat; 7 = pd

if block == 3
    list = dir ('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/e*recogminor*_tsssdsm.mat'); %dir to epoched files (encoding)
elseif block == 4
    list = dir ('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/e*recogspeed*_tsssdsm.mat'); %dir to epoched files (encoding)
elseif block == 5
    list = dir ('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/e*samemel*_tsssdsm.mat'); %dir to epoched files (encoding)
elseif block == 6
    list = dir ('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/e*visualpat*_tsssdsm.mat'); %dir to epoched files (encoding)
elseif block == 7
    list = dir ('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/e*pd*_tsssdsm.mat'); %dir to epoched files (encoding)
end

for ii = 1:length(list)
    D = spm_eeg_load([list(ii).folder '/' list(ii).name]);
    if isfield(D,'inv')
        if copy_label == 0 %simply displaying RHINO coregistration
            if isfield(D,'inv') %checking if the coregistration was already run
                rhino_display(D)
            end
        else %pasting inv RHINO from epoched data (where it was computed) to continuous data
            inv_rhino = D.inv;
            D2 = spm_eeg_load([list(ii).folder '/' list(ii).name(2:end)]);
            D2.inv = inv_rhino;
            D2.save();
        end
    end
    disp(['Block ' num2str(block) ' - Subject ' num2str(ii)])
end

%%

%% BEAMFORMING

%% SETTING FOR CLUSTER (PARALLEL COMPUTING)

% clusterconfig('scheduler', 'none'); %If you do not want to submit to the cluster, but simply want to test the script on the hyades computer, you can instead of 'cluster', write 'none'
clusterconfig('scheduler', 'cluster'); %If you do not want to submit to the cluster, but simply want to test the script on the hyades computer, you can instead of 'cluster', write 'none'
clusterconfig('long_running', 1); % This is the cue we want to use for the clsuter. There are 3 different cues. Cue 0 is the short one, which should be enough for us
clusterconfig('slot', 1); %slot is memory, and 1 memory slot is about 8 GB. Hence, set to 2 = 16 GB
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/Cluster_ParallelComputing')

%% FUNCTION FOR SOURCE RECONSTRUCTION

%user settings
clust_l = 1; %1 = using cluster of computers (CFIN-MIB, Aarhus University); 0 = running locally
timek = 1:1026; %time-points
freqq = []; %frequency range (empty [] for broad band)
% freqq = [0.1 1]; %frequency range (empty [] for broad band)
% freqq = [2 8]; %frequency range (empty [] for broad band)
sensl = 1; %1 = magnetometers only; 2 = gradiometers only; 3 = both magnetometers and gradiometers (SUGGESTED 1!)
workingdir2 = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD'; %high-order working directory (a subfolder for each analysis with information about frequency, time and absolute value will be created)
block = 3; % 3 = recogminor; 4 = recogspeed; 5 = samemel; 6 = visualpat; 7 = pd
invers = 1; %1-4 = different ways (e.g. mean, t-values, etc.) to aggregate trials and then source reconstruct only one trial; 5 for single trial independent source reconstruction

if isempty(freqq)
    absl = 0; % 1 = absolute value of sources; 0 = not
else
    absl = 0;
end
%actual computation
%list of subjects with coregistration (RHINO - OSL/FSL) - epoched
if block == 3
    list = dir ('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/e*recogminor*_tsssdsm.mat'); %dir to epoched files (encoding)
    condss = {'Old_Correct','New_T1_Correct','New_T2_Correct','New_T3_Correct','New_T4_Correct'};
    load('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/time_normal.mat');
elseif block == 4
    list = dir ('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/e*recogspeed*_tsssdsm.mat'); %dir to epoched files (encoding)
    condss = {'Old_Fast_Correct','Old_Slow_Correct','New_Fast_Correct','New_Slow_Correct'};
elseif block == 5
    list = dir ('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/e*samemel*_tsssdsm.mat'); %dir to epoched files (encoding)
    condss = {'Encoding','Old_Correct','New_Correct'};
elseif block == 6
    list = dir ('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/e*visualpat*_tsssdsm.mat'); %dir to epoched files (encoding)
    condss = {'Encoding','Old_Correct','New_Correct'};
elseif block == 7
    list = dir ('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/e*pd*_tsssdsm.mat'); %dir to epoched files (encoding)
    condss = {'pd','pm'};
end
if isempty(freqq)
    workingdir = [workingdir2 '/Block_' num2str(block) '/Beam_abs_' num2str(absl) '_sens_' num2str(sensl) '_freq_broadband_invers_' num2str(invers)];
else
    workingdir = [workingdir2 '/Block_' num2str(block) '/Beam_abs_' num2str(absl) '_sens_' num2str(sensl) '_freq_' num2str(freqq(1)) '_' num2str(freqq(2)) '_invers_' num2str(invers)];
end
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/Cluster_ParallelComputing');
if ~exist(workingdir,'dir') %creating working folder if it does not exist
    mkdir(workingdir)
end
for ii = 1:length(list) %over subjects
    S = [];
    if ~isempty(freqq) %if you want to apply the bandpass filter, you need to provide continuous data
        %             disp(['copying continuous data for subj ' num2str(ii)])
        %thus pasting it here
        %             copyfile([list_c(ii).folder '/' list_c(ii).name],[workingdir '/' list_c(ii).name]); %.mat file
        %             copyfile([list_c(ii).folder '/' list_c(ii).name(1:end-3) 'dat'],[workingdir '/' list_c(ii).name(1:end-3) 'dat']); %.dat file
        %and assigning the path to the structure S
        S.norm_megsensors.MEGdata_c = [list(ii).folder '/' list(ii).name(2:end)];
    end
    %copy-pasting epoched files
    %         disp(['copying epoched data for subj ' num2str(ii)])
    %         copyfile([list(ii).folder '/' list(ii).name],[workingdir '/' list(ii).name]); %.mat file
    %         copyfile([list(ii).folder '/' list(ii).name(1:end-3) 'dat'],[workingdir '/' list(ii).name(1:end-3) 'dat']); %.dat file
    
    S.Aarhus_cluster = clust_l; %1 for parallel computing; 0 for local computation
    
    S.norm_megsensors.zscorel_cov = 1; % 1 for zscore normalization; 0 otherwise
    S.norm_megsensors.workdir = workingdir;
    S.norm_megsensors.MEGdata_e = [list(ii).folder '/' list(ii).name];
    S.norm_megsensors.freq = freqq; %frequency range
    S.norm_megsensors.forward = 'Single Shell'; %forward solution (for now better to stick to 'Single Shell')
    
    S.beamfilters.sensl = sensl; %1 = magnetometers; 2 = gradiometers; 3 = both MEG sensors (mag and grad) (SUGGESTED 3!)
    S.beamfilters.maskfname = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/MNI152_T1_8mm_brain.nii.gz'; % path to brain mask: (e.g. 8mm MNI152-T1: '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/MNI152_T1_8mm_brain.nii.gz')
    
    S.inversion.znorml = 0; % 1 for inverting MEG data using the zscored normalized one; (SUGGESTED 0 IN BOTH CASES!)
    %                                 0 to normalize the original data with respect to maximum and minimum of the experimental conditions if you have both magnetometers and gradiometers.
    %                                 0 to use original data in the inversion if you have only mag or grad (while e.g. you may have used zscored-data for covariance matrix)
    %
    S.inversion.timef = timek; %data-points to be extracted (e.g. 1:300); leave it empty [] for working on the full length of the epoch
    S.inversion.conditions = condss; %cell with characters for the labels of the experimental conditions (e.g. {'Old_Correct','New_Correct'})
    S.inversion.bc = [1 26]; %extreme time-samples for baseline correction (leave empty [] if you do not want to apply it)
    S.inversion.abs = absl; %1 for absolute values of sources time-series (recommendnded 1!)
    S.inversion.effects = invers;
    
    S.smoothing.spatsmootl = 0; %1 for spatial smoothing; 0 otherwise
    S.smoothing.spat_fwhm = 100; %spatial smoothing fwhm (suggested = 100)
    S.smoothing.tempsmootl = 0; %1 for temporal smoothing; 0 otherwise
    S.smoothing.temp_param = 0.01; %temporal smoothing parameter (suggested = 0.01)
    S.smoothing.tempplot = [1 2030 3269]; %vector with sources indices to be plotted (original vs temporally smoothed timeseries; e.g. [1 2030 3269]). Leave empty [] for not having any plot.
    
    S.nifti = 1; %1 for plotting nifti images of the reconstructed sources of the experimental conditions
    S.out_name = ['SUBJ_' list(ii).name(13:16)]; %name (character) for output nifti images (conditions name is automatically detected and added)
    
    if clust_l ~= 1 %useful  mainly for begugging purposes
        MEG_SR_Beam_LBPD(S);
    else
        jobid = job2cluster(@MEG_SR_Beam_LBPD,S); %running with parallel computing
    end
end

%% SETTING FOR CLUSTER (PARALLEL COMPUTING)

% clusterconfig('scheduler', 'none'); %If you do not want to submit to the cluster, but simply want to test the script on the hyades computer, you can instead of 'cluster', write 'none'
clusterconfig('scheduler', 'cluster'); %If you do not want to submit to the cluster, but simply want to test the script on the hyades computer, you can instead of 'cluster', write 'none'
clusterconfig('long_running', 0); % This is the cue we want to use for the clsuter. There are 3 different cues. Cue 0 is the short one, which should be enough for us
clusterconfig('slot', 3); %slot is memory, and 1 memory slot is about 8 GB. Hence, set to 2 = 16 GB
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/Cluster_ParallelComputing')

%% STATISTICS OVER PARTICIPANTS (USING PARALLEL COMPUTING, PARTICULARLY USEFUL IF YOU HAVE SEVERAL CONTRASTS)

block = 3; % 3 = recogminor; 4 = recogspeed; 5 = samemel; 6 = visualpat; 7 = pd
clust = 1; % 1 = using Aarhus cluster (parallel computing); 0 = run locally
analys_n = 2; %analysis number (in the list indexed below)

%building structure
asd = dir(['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_' num2str(block) '/Beam*']);
S = [];
% dumm = dir([asd(analys_n).folder '/' asd(analys_n).name '/SUBJ*mat']);
% S.list = dumm(1:10);
S.workingdir = [asd(analys_n).folder '/' asd(analys_n).name '/test']; %path where the data from MEG_SR_Beam_LBPD.m is stored
S.sensl = 1; % 1 = magnetometers only; 2 = gradiometers only; 3 = both magnetometers and gradiometers.
S.plot_nifti = 1; %1 to plot nifti images; 0 otherwise
S.plot_nifti_name = []; %character with name for nifti files (it may be useful if you run separate analysis); Leave empty [] to not  specify any name
% S.contrast = [1 0 0 0 0 0 -1; 0 1 0 0 0 0 -1; 0 0 1 0 0 0 -1; 0 0 0 1 0 0 -1; 0 0 0 0 1 0 -1; 0 0 0 0 0 1 -1; 1 1 1 1 1 1 -1]; %one row per contrast (e.g. having 3 conditions, [1 -1 0; 1 -1 -1; 0 1 -1]; two or more -1 or 1 are interpreted as the mean over them first and then the contrast. Leave empty [] for no contrasts. 
if block == 3
    S.contrast = [1 -1 0 0 0; 1 0 -1 0 0; 1 0 0 -1 0; 1 0 0 0 -1];
elseif block == 4
        S.contrast = [1 0 -1 0; 0 1 0 -1];
elseif block == 5 || block == 6
    S.contrast = [0 1 -1; -1 1 0];
else
    S.contrast = [1 -1];
end
S.effects = 1; %mean over subjects for now
if clust == 1
    S.Aarhus_clust = 1; %1 to use paralle computing (Aarhus University, contact me, Leonardo Bonetti, for more information; leonardo.bonetti@clin.au.dk)
    %actual function
    jobid = job2cluster(@MEG_SR_Stats1_Fast_LBPD,S); %running with parallel computing
else
    S.Aarhus_clust = 0; %1 to use paralle computing (Aarhus University, contact me, Leonardo Bonetti, for more information; leonardo.bonetti@clin.au.dk)
    MEG_SR_Stats1_Fast_LBPD(S)
end

%% Step 4: Source reconstruction for the peaks of differential activity between Old and New conditions (Figure 3)

% Using the above defined peaks, arranged as hemisphere (l-r), tone (2-5), max-min (thus 16 peaks in total)
% The peaks are reported here in time-samples
peaks = [203 179 299 274 389 361 477 460 181 234 273 298 364 390 459 506];
names = ['L_t2_max'; 'L_t2_min'; 'L_t3_max'; 'L_t3_min'; 'L_t4_max'; 'L_t4_min'; 'L_t5_max'; 'L_t5_min'; 'R_t2_max'; 'R_t2_min'; 'R_t3_max'; 'R_t3_min'; 'R_t4_max'; 'R_t4_min'; 'R_t5_max'; 'R_t5_min'];
vectcond = [2,2,2,3,2,4,2,5,2,2,3,2,4,2,5,2]; %condition X to be contrasted against condition 1 (Old)
list = dir(['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/Beam_abs_1_sens_1_freq_broadband/SUBJ*.mat']);
data = zeros(3559,length(peaks),5,length(list)); %sources x time-windows of intereste x conditions x subjects
for ii = 1:length(list) %over subjects
    load([list(ii).folder '/' list(ii).name])
    for jj = 1:length(peaks) %over time-windows
        data(:,jj,:,ii) = mean(OUT.sources_ERFs(:,peaks(jj)-5:peaks(jj)+5,:),2); %averge over the selected jj time-window
    end
    disp(ii)
end
%computing t-tests
P = zeros(3559,length(vectcond));
T = zeros(3559,length(vectcond));
for ii = 1:size(data,1) %over sources
    for jj = 1:length(vectcond) %over time-windows
        [~,p,~,stats] = ttest(squeeze(data(ii,jj,1,:)),squeeze(data(ii,jj,vectcond(jj),:))); %contrasting cond1 (Old) with cond X (depending on vectcond it can be either NewT1 or NewT2 or NewT3 or NewT4
        P(ii,jj) = p;
        T(ii,jj) = stats.tstat;
    end
    disp(ii)
end
%creating nifti images
outdir = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/Papers/Systematic_Variations_DCM/NatureCommunications_Revision/Revision_I/MEGSensors_FunctionalParcellation/MEG_sources_stats_peaks';
maskk = load_nii('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/MNI152_8mm_brain_diy.nii.gz'); %getting the mask for creating the figure
for iii = 1:size(T,2)
    fnameniit = [outdir '/' names(iii,:) '_T.nii.gz']; %path and name of the image to be saved
    fnameniip = [outdir '/' names(iii,:) '_P.nii.gz']; %path and name of the image to be saved
    SO = T(:,iii);
    SOp = 1 - P(:,iii);
    %building nifti image
    SS = size(maskk.img);
    dumimg = zeros(SS(1),SS(2),SS(3),1); %no time here so size of 4D = 1
    dumimgp = zeros(SS(1),SS(2),SS(3),1); %no time here so size of 4D = 1
    for ii = 1:size(SO,1) %over brain sources
        dum = find(maskk.img == ii); %finding index of sources ii in mask image (MNI152_8mm_brain_diy.nii.gz)
        [i1,i2,i3] = ind2sub([SS(1),SS(2),SS(3)],dum); %getting subscript in 3D from index
        dumimg(i1,i2,i3,:) = SO(ii,:); %storing values for all time-points in the image matrix
        dumimgp(i1,i2,i3,:) = SOp(ii,:); %storing values for all time-points in the image matrix
    end
    %T
    nii = make_nii(dumimg,[8 8 8]); %making nifti image from 3D data matrix (and specifying the 8 mm of resolution)
    nii.img = dumimg; %storing matrix within image structure
    nii.hdr.hist = maskk.hdr.hist; %copying some information from maskk
    disp(['saving nifti image - ' num2str(iii)])
    save_nii(nii,fnameniit); %printing image
    %P
    nii = make_nii(dumimg,[8 8 8]); %making nifti image from 3D data matrix (and specifying the 8 mm of resolution)
    nii.img = dumimgp; %storing matrix within image structure
    nii.hdr.hist = maskk.hdr.hist; %copying some information from maskk
    disp(['saving nifti image - ' num2str(iii)])
    save_nii(nii,fnameniip); %printing image
end
%FDR correction and creation of T-values images keeping only the FDR-significant voxels
pthr = zeros(size(P,2),1);
for pp = 1:size(P,2) %over peaks
    [pthr(pp,1)] = fdr(P(:,pp));
end
%creating FDR-corrected nifti images
outdir = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/Papers/Systematic_Variations_DCM/NatureCommunications_Revision/Revision_I/MEGSensors_FunctionalParcellation/MEG_sources_stats_peaks/FDR_Corrected_Tvalue';
maskk = load_nii('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/MNI152_8mm_brain_diy.nii.gz'); %getting the mask for creating the figure
DUMM = zeros(3559,size(T,2));
for iii = 1:size(T,2)
    fnameniit = [outdir '/' names(iii,:) '_FDR_T.nii.gz']; %path and name of the image to be saved
    SO = T(:,iii);
    Pdum = P(:,iii);
    SO(Pdum>pthr(iii)) = 0; %thresholding the T-values with FDR-corrected p-values
    %building nifti image
    SS = size(maskk.img);
    dumimg = zeros(SS(1),SS(2),SS(3),1); %no time here so size of 4D = 1
    dumimgp = zeros(SS(1),SS(2),SS(3),1); %no time here so size of 4D = 1
    for ii = 1:size(SO,1) %over brain sources
        dum = find(maskk.img == ii); %finding index of sources ii in mask image (MNI152_8mm_brain_diy.nii.gz)
        [i1,i2,i3] = ind2sub([SS(1),SS(2),SS(3)],dum); %getting subscript in 3D from index
        dumimg(i1,i2,i3,:) = SO(ii,:); %storing values for all time-points in the image matrix
%         dumimgp(i1,i2,i3,:) = SOp(ii,:); %storing values for all time-points in the image matrix
    end
    %T
    nii = make_nii(dumimg,[8 8 8]); %making nifti image from 3D data matrix (and specifying the 8 mm of resolution)
    nii.img = dumimg; %storing matrix within image structure
    nii.hdr.hist = maskk.hdr.hist; %copying some information from maskk
    disp(['saving nifti image - ' num2str(iii)])
    save_nii(nii,fnameniit); %printing image
    %exporting data for NC source data
    DUMM(:,iii) = SO;
end
PDn = table(DUMM);
writetable(PDn,['/aux/MINDLAB2021_MEG-TempSeqAges/09_11_2023/FDRCorrected_Sources.xlsx'],'Sheet',1); %printing excel file

%after this, we created the final images using Workbench

% Extracting information about the clusters at source level and reporting it in xlsx files
%here we obtain information about the brain regions forming the clusters
%the tables can be found in SUPPLEMENTARY MATERIALS

path = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/Papers/Systematic_Variations_DCM/NatureCommunications_Revision/Revision_I/MEGSensors_FunctionalParcellation/MEG_sources_stats_peaks/FDR_Corrected_Tvalue';
list = dir([path '/*gz']);
MNI = [];
for ii = 1:length(list)
    fname = [path '/' list(ii).name]; %tone x cluster 1
    S = [];
    S.fname = fname;
    %actual function
    PDn = Extract_BrainCluster_Information_3D_LBPD_D(S);
    tvall = cell2mat(table2array(PDn(3:end,3))); %extracting t-values
    idxx = find(tvall>mean(tvall)+(1*std(tvall))); %finding voxels with values bigger than the average and one standard deviation
    MNI = cat(1,MNI,cell2mat(table2array(PDn(idxx + 2,4))));
    writetable(PDn,[path '/' list(ii).name(1:end-7) '.xlsx'],'Sheet',1); %printing excel file
end
save('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/Papers/Systematic_Variations_DCM/NatureCommunications_Revision/Revision_I/MEGSensors_FunctionalParcellation/MNI_tval.mat','MNI');

%%

%% *** STATISTICS ON ALL MEG SENSORS (each sensor independently) ***

% Figure S5

load('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/Papers/Systematic_Variations_DCM/NatureCommunications_Revision/Revision_I/Images_Data/data_MEGsensors_all.mat');
%t-tests
P = zeros(size(dataa,1),size(dataa,2),(size(dataa,3)-1)); %clusters x time-points x contrasts (every NewTX versus Old)
T = zeros(size(dataa,1),size(dataa,2),(size(dataa,3)-1)); %clusters x time-points x contrasts (every NewTX versus Old)
for ii = 1:size(dataa,1) %over MEG sensors
    for jj = 1:size(dataa,2) %ove time-points
        for cc = 1:(size(dataa,3)-1) %over the 4 NewTX
            %codes t-test
            [~,p,~,stats] = ttest(squeeze(dataa(ii,jj,1,:)),squeeze(dataa(ii,jj,cc+1,:))); %contrasting cond1 (Old) with cond X (depending on vectcond it can be either NewT1 or NewT2 or NewT3 or NewT4
            P(ii,jj,cc) = p;
            T(ii,jj,cc) = stats.tstat;
        end
        disp([num2str(ii) ' - ' num2str(jj)])
    end
end
%MCS
p_thresh = 0.05; %threshold for binarising p-values vector
load('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/time_normal.mat'); %loading time
outdir = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/Papers/Systematic_Variations_DCM/NatureCommunications_Revision/Revision_I/MEGSensors_FunctionalParcellation/All_MEG_sensors_stats';
mkdir(outdir);
SIGN = cell(size(dataa,1),(size(dataa,3)-1));
for ii = 1:size(dataa,1) %over MEG sensors
    for cc = 1:(size(dataa,3)-1) %over the 4 NewTX
        Pbin = zeros(1,size(dataa,2));
        Pbin(P(ii,:,cc)<p_thresh) = 1;
        tvals = T(ii,:,cc);
        [ sign_clust ] = oneD_MCS_LBPD_D( Pbin, 0, 1, 1000, 0.001, time(1:size(dataa,2)), tvals ) %Monte Carlo simulation function to correct for multiple comparisons
        PDn = cell2table(sign_clust); %table
        writetable(PDn,[outdir '/MEG_sensor_'  chanlabs{ii} '_OldvsNewT' num2str(cc) '.xlsx'],'Sheet',1); %printing excel file
        SIGN{ii,cc} = sign_clust;
    end
end
save('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/Papers/Systematic_Variations_DCM/NatureCommunications_Revision/Revision_I/MEGSensors_FunctionalParcellation/All_MEG_sensors_stats/Sign.mat','SIGN')

%%

%% *** AAL PARCELLATION ***

% Figures 4 and S7

%loading data in 3559-voxel space (8mm)
%loading AAL parcels label
AAL_lab = importdata('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/aal_ROIs_cerebellum.txt');
%list of source reconstructed data (for each subject)
list = dir('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/Beam_abs_0_sens_1_freq_broadband/SUBJ*.mat');
%getting sign of the voxels based on the aggregated source reconstruction (over participants) - AVERAGED TRIALS
load('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/Beam_abs_0_sens_1_freq_broadband/sources_main_effects.mat');

%actual computation
%adjusting polarity
timex = 45:52;
vect = zeros(3559,1);
for jj = 1:3559 %over brain voxels
    if squeeze(mean(t_val_s(jj,timex,1),2)) > 0 %if the data in voxel jj is positive during N100 time
        vect(jj,1) = -1; %storing a vector with 1 and -1 to be used for later statistics
    else
        vect(jj,1) = 1;
    end
end
load([list(1).folder '/' list(1).name])
dum2 = zeros(116,size(OUT.sources_ERFs,2),size(OUT.sources_ERFs,3),length(list));
for ii = 1:length(list) %over subjects
    disp(['loading source reconstructed data for subject ' num2str(ii) ' / ' num2str(length(list))])
    load([list(ii).folder '/' list(ii).name])
    dum = zeros(size(OUT.sources_ERFs,1),size(OUT.sources_ERFs,2),size(OUT.sources_ERFs,3));
    for cc = 1:size(dum,3) %over conditions
        for jj = 1:size(OUT.sources_ERFs,1) %over brain voxels
            dum(jj,:,cc) = OUT.sources_ERFs(jj,:,cc) .* vect(jj,1); %reversing (or not)..
            disp(['subject - ' num2str(ii) ' - condition ' num2str(cc) ' - source ' num2str(jj)])
        end
    end
    for pp = 1:length(AAL_lab) %over ROIs (with possibility of averaging together voxels of different ROIs (e.g. left and right auditory cortex (AC))
        %exracting voxels for pp parcel of AAL
        S = [];
        S.input = 2; %1 = MNI coordinates; 2 = AAL ROIs; 3 = general image with non-zero values
        S.coordd = []; %coordinates in MNI space (x,y,z)
        S.AAL_ROIs = [pp]; %AAL ROIs numbers you want to use
        % S.image = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/Beam_abs_0_sens_1_freq_broadband/Contr_1_abs_0.nii.gz';
        S.image = [];
        %actual function
        idx_LBPD = From3DNifti_OrMNICoords_2_CoordMatrix_8mm_LBPD_D(S);
        for cc = 1:size(dum,3) %over conditions
            dum2(pp,:,cc,ii) = mean(dum(idx_LBPD,:,cc),1);
        end
    end
end
dataa = dum2;
save(['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/Papers/Systematic_Variations_DCM/NatureCommunications_Revision/Revision_I/MEGSensors_FunctionalParcellation/AAL/AAL_data.mat'],'dataa')

%% Plotting

v = [91:116]; %AAL parcels to plot

% load('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/Papers/Systematic_Variations_DCM/NatureCommunications_Revision/Revision_I/MEGSensors_FunctionalParcellation/AAL/data.mat');
%loading AAL parcels label
AAL_lab = importdata('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/aal_ROIs_cerebellum.txt');
%defining colors
color_line = colormap(lines(5)); %extracting some colours from a colormap
color_line2 = color_line;
color_line2(1,:) = color_line(2,:);
color_line2(2,:) = color_line(1,:);
color_line2(5,:) = [0.4 0.4 0.4];
clear conds
conds{1} = ' Mem '; conds{2} = ' NewT1 '; conds{3} = ' NewT2 '; conds{4} = ' NewT3 '; conds{5} = ' NewT4 ';
load('/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/Source_LBPD/time_345.mat')
for ii = 1:length(v) %over clusters
    figure
    for cc = 1:length(conds) %over conditions
        plot(squeeze(time(1:size(dataa,2))),squeeze(nanmean(dataa(v(ii),:,cc,:),4)),'Color',color_line2(cc,:),'LineWidth',2,'DisplayName',[conds{cc}])
        hold on
        plot(squeeze(time(1:size(dataa,2))),squeeze(nanmean(dataa(v(ii),:,cc,:),4)) + (nanstd(dataa(v(ii),:,cc,:),0,4)./sqrt(size(dataa,4))),':','Color',color_line2(cc,:),'LineWidth',0.5,'HandleVisibility','off')
        hold on
        plot(squeeze(time(1:size(dataa,2))),squeeze(nanmean(dataa(v(ii),:,cc,:),4)) - (nanstd(dataa(v(ii),:,cc,:),0,4)./sqrt(size(dataa,4))),':','Color',color_line2(cc,:),'LineWidth',0.5,'HandleVisibility','off')
        hold on
    end
    grid minor
    title(AAL_lab{v(ii)})
    xlim([-0.1 3.4])
%     ylim(limmy)
    set(gcf,'color','w')
    legend('show')
end

%% Statistics

load('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/Papers/Systematic_Variations_DCM/NatureCommunications_Revision/Revision_I/MEGSensors_FunctionalParcellation/AAL/AAL_data.mat');
%loading AAL parcels label
AAL_lab = importdata('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/aal_ROIs_cerebellum.txt');
%t-tests
P = zeros(size(dataa,1),size(dataa,2),(size(dataa,3)-1)); %clusters x time-points x contrasts (every NewTX versus Old)
T = zeros(size(dataa,1),size(dataa,2),(size(dataa,3)-1)); %clusters x time-points x contrasts (every NewTX versus Old)
for ii = 1:size(dataa,1) %over AAL parcels
    for jj = 1:size(dataa,2) %ove time-points
        for cc = 1:(size(dataa,3)-1) %over the 4 NewTX
            %codes t-test
            [~,p,~,stats] = ttest(squeeze(dataa(ii,jj,1,:)),squeeze(dataa(ii,jj,cc+1,:))); %contrasting cond1 (Old) with cond X (depending on vectcond it can be either NewT1 or NewT2 or NewT3 or NewT4
            P(ii,jj,cc) = p;
            T(ii,jj,cc) = stats.tstat;
        end
        disp([num2str(ii) ' - ' num2str(jj)])
    end
end
%MCS
p_thresh = 0.05; %threshold for binarising p-values vector
load('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/time_normal.mat'); %loading time
outdir = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/Papers/Systematic_Variations_DCM/NatureCommunications_Revision/Revision_I/MEGSensors_FunctionalParcellation/AAL/Stats';
mkdir(outdir);
SIGN = cell(size(dataa,1),(size(dataa,3)-1));
for ii = 1:size(dataa,1) %over AAL parcels
    for cc = 1:(size(dataa,3)-1) %over the 4 NewTX
        Pbin = zeros(1,size(dataa,2));
        Pbin(P(ii,:,cc)<p_thresh) = 1;
        tvals = T(ii,:,cc);
        [ sign_clust ] = oneD_MCS_LBPD_D( Pbin, 0, 1, 1000, 0.001, time(1:size(dataa,2)), tvals ) %Monte Carlo simulation function to correct for multiple comparisons
        PDn = cell2table(sign_clust); %table
        writetable(PDn,[outdir '/AAL_'  AAL_lab{ii} '_OldvsNewT' num2str(cc) '.xlsx'],'Sheet',1); %printing excel file
        SIGN{ii,cc} = sign_clust;
    end
end
save([outdir '/Sign.mat'],'SIGN')

%% SOURCE LEAKAGE CORRECTION - AAL

cond = 1; %condition that you want; 1 = Old; 2-5 = NewTX (X = 1-4)
ROIs = [1]; %selected ROIs: 1 = ACL, 2 = ACR, 3 = HIPPL, 4 = HIPPR, 5 = ACC, 6 = MC. They correspond to AAL ID 79, 80 (Heschl's L, R), 37, 38 (Hipp L, R), 31-32 (anterior cingulate) 33-34 (medial cingulate)
limmy = []; %limit for y-axis (amplitude of the signal)

if ~exist('dum2','var')
    load('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/Papers/Systematic_Variations_DCM/NatureCommunications_Revision/Revision_I/MEGSensors_FunctionalParcellation/AAL/AAL_data.mat');
    load('/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/Source_LBPD/time_345.mat')
    %loading AAL parcels label
    AAL_lab = importdata('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/aal_ROIs_cerebellum.txt');
    dum2 = dataa;
end
conds{1} = ' Mem '; conds{2} = ' NewT1 '; conds{3} = ' NewT2 '; conds{4} = ' NewT3 '; conds{5} = ' NewT4 ';
tso2 = zeros(6,size(dum2,2),5,size(dum2,4));
tso = zeros(6,size(dum2,2),size(dum2,4));
ROIN{1} = 'ACL'; ROIN{2} = 'ACR'; ROIN{3} = 'HIPPL'; ROIN{4} = 'HIPPR'; ROIN{5} = 'ACC'; ROIN{6} = 'MC';
vectt = [79 80 37 38 31 32 33 34]; %8 selected ROIs
% tso2 = zeros(size(tso,1),size(tso,2),size(tso,3),size(tso,4));
tso2(1,:,:,:) = dum2(vectt(1),:,:,:); %ACL
tso2(2,:,:,:) = dum2(vectt(2),:,:,:); %ACR
tso2(3,:,:,:) = dum2(vectt(3),:,:,:); %HIPPL
tso2(4,:,:,:) = dum2(vectt(4),:,:,:); %HIPPR
tso2(5,:,:,:) = mean(dum2(vectt(5:6),:,:,:),1); %ACC
tso2(6,:,:,:) = mean(dum2(vectt(7:8),:,:,:),1); %MC
% tso = tso2;
for ii = 1:size(dum2,4) %over subjects
    tso(:,:,ii) = remove_source_leakage_b(squeeze(tso2(:,:,cond,ii)),'symmetric'); %source leakage correction
    disp(ii)
end
%implementation for all conditions to export to local computer and plotting there
tso = zeros(6,size(dum2,2),5,size(dum2,4));
for ii = 1:size(dum2,4) %over subjects
    for cc = 1:size(dum2,3) %over conditions
        tso(:,:,cc,ii) = remove_source_leakage_b(squeeze(tso2(:,:,cc,ii)),'symmetric'); %source leakage correction
        disp(ii)
    end
end
save sources_leakagecorrected_AAL.mat tso ROIN conds
% close all
%defining colors
color_line = colormap(lines(5)); %extracting some colours from a colormap
color_line2 = color_line;
color_line2(1,:) = color_line(2,:);
color_line2(2,:) = color_line(1,:);
color_line2(5,:) = [0.4 0.4 0.4];
LIN{1} = 3.2; LIN{2} = 2.8; LIN{3} = 2.3; LIN{4} = 1.8; LIN{5} = 1;
% load('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/time_normal.mat')
figure
clear CC
for ii = 1:length(ROIs) %over ROIs
    plot(time(1:size(tso,2)),nanmean(tso(ROIs(ii),:,:),3),'Color',color_line2(ii,:),'LineWidth',2,'DisplayName',[ROIN{ROIs(ii)}])
    hold on
    plot(time(1:size(tso,2)),nanmean(tso(ROIs(ii),:,:),3) + (nanstd(tso(ROIs(ii),:,:),0,3)./sqrt(size(tso,3))),':','Color',color_line2(ii,:),'LineWidth',0.5,'HandleVisibility','off')
    hold on
    plot(time(1:size(tso,2)),nanmean(tso(ROIs(ii),:,:),3) - (nanstd(tso(ROIs(ii),:,:),0,3)./sqrt(size(tso,3))),':','Color',color_line2(ii,:),'LineWidth',0.5,'HandleVisibility','off')
    hold on
end
grid minor
% xlim([-0.1 time(size(dum2,2))])
xlim([-0.1 3.4])
if ~isempty(limmy)
    ylim(limmy)
end
set(gcf,'color','w')
legend('show')

%% STATISTICS AFTER SOURCE LEAKAGE CORRECTION (Figure S9)

load('/aux/MINDLAB2021_MEG-TempSeqAges/19_10_2023/sources_leakagecorrected_AAL.mat');
ROIN{1} = 'ACL'; ROIN{2} = 'ACR'; ROIN{3} = 'HIPPL'; ROIN{4} = 'HIPPR'; ROIN{5} = 'ACC'; ROIN{6} = 'MC';
dataa = tso;
%t-tests
P = zeros(size(dataa,1),size(dataa,2),(size(dataa,3)-1)); %clusters x time-points x contrasts (every NewTX versus Old)
T = zeros(size(dataa,1),size(dataa,2),(size(dataa,3)-1)); %clusters x time-points x contrasts (every NewTX versus Old)
for ii = 1:size(dataa,1) %over AAL parcels
    for jj = 1:size(dataa,2) %ove time-points
        for cc = 1:(size(dataa,3)-1) %over the 4 NewTX
            %codes t-test
            [~,p,~,stats] = ttest(squeeze(dataa(ii,jj,1,:)),squeeze(dataa(ii,jj,cc+1,:))); %contrasting cond1 (Old) with cond X (depending on vectcond it can be either NewT1 or NewT2 or NewT3 or NewT4
            P(ii,jj,cc) = p;
            T(ii,jj,cc) = stats.tstat;
        end
        disp([num2str(ii) ' - ' num2str(jj)])
    end
end
%MCS
p_thresh = 0.05; %threshold for binarising p-values vector
load('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/time_normal.mat'); %loading time
outdir = '/aux/MINDLAB2021_MEG-TempSeqAges/19_10_2023/Stats';
mkdir(outdir);
SIGN = cell(size(dataa,1),(size(dataa,3)-1));
for ii = 1:size(dataa,1) %over AAL parcels
    for cc = 1:(size(dataa,3)-1) %over the 4 NewTX
        Pbin = zeros(1,size(dataa,2));
        Pbin(P(ii,:,cc)<p_thresh) = 1;
        tvals = T(ii,:,cc);
        [ sign_clust ] = oneD_MCS_LBPD_D( Pbin, 0, 1, 1000, 0.001, time(1:size(dataa,2)), tvals ) %Monte Carlo simulation function to correct for multiple comparisons
        PDn = cell2table(sign_clust); %table
        writetable(PDn,[outdir '/'  ROIN{ii} '_OldvsNewT' num2str(cc) '.xlsx'],'Sheet',1); %printing excel file
        SIGN{ii,cc} = sign_clust;
    end
end
save([outdir '/Sign.mat'],'SIGN')

%%

%% *** FOCUS ON PREDICTION ERROR ***

% Figure 5

%% STRONGER RESPONSES TO FIRST SOUND WHERE VARIATION IS INSERTED IN HIPPOCAMPUS AND ACC, BUT NOT IN HESCHL'S GYRUS

ROIN{1} = 'LHG'; ROIN{2} = 'RHG'; ROIN{3} = 'LHP'; ROIN{4} = 'RHP'; ROIN{5} = 'ACC'; ROIN{6} = 'MC'; % ROIs in dum2
HG = {[159:169] [246:256] [344:354] [434:444]}; %time-points (already with plus/minus 20ms) of the peaks for ACL
HPCC = {[176:186] [268:278] [363:373] [454:464]}; %time-points (already with plus/minus 20ms) of the peaks for VMPFC

if ~exist('dum2','var')
    load('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/Papers/Systematic_Variations_DCM/NatureCommunications_Revision/Revision_I/MEGSensors_FunctionalParcellation/AAL/AAL_data.mat');
    load('/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/Source_LBPD/time_345.mat')
    %loading AAL parcels label
    AAL_lab = importdata('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/aal_ROIs_cerebellum.txt');
    dum2 = dataa;
    clear dataa
end
conds{1} = ' Mem '; conds{2} = ' NewT1 '; conds{3} = ' NewT2 '; conds{4} = ' NewT3 '; conds{5} = ' NewT4 ';
tso2 = zeros(6,size(dum2,2),5,size(dum2,4));
tso = zeros(6,size(dum2,2),size(dum2,4));
vectt = [79 80 37 38 31 32 33 34]; %8 selected ROIs
% tso2 = zeros(size(tso,1),size(tso,2),size(tso,3),size(tso,4));
tso2(1,:,:,:) = dum2(vectt(1),:,:,:); %ACL
tso2(2,:,:,:) = dum2(vectt(2),:,:,:); %ACR
tso2(3,:,:,:) = dum2(vectt(3),:,:,:); %HIPPL
tso2(4,:,:,:) = dum2(vectt(4),:,:,:); %HIPPR
tso2(5,:,:,:) = mean(dum2(vectt(5:6),:,:,:),1); %ACC
tso2(6,:,:,:) = mean(dum2(vectt(7:8),:,:,:),1); %MC
SIGN = cell(length(ROIN),3,2); %ROIs x conditions - 1 x 2 (ANOVA statistics and post-hoc tests)
cntp = 0;
PP = zeros(18,1);
for ss = 1:size(tso2,1) %over ROIs
    for cc = 1:4 %over NTX conditions
        cnt = 0;
        dumbo = zeros(83,5-cc); %subjects x peaks
        for pp = cc:4 %over the peaks (that are 4 for cc = 2 (NT1), 3 for cc = 3 (NT2), etc.)
            cnt = cnt + 1;
            if ss < 3 %if Heschl's gyrus (either left or right)
                dumbo(:,cnt) = squeeze(mean(tso2(ss,HG{pp},cc+1,:),2)); %extracting ROI ss, time-window of Heschl's gyrus, condition NTcc, all participants
            else
                dumbo(:,cnt) = squeeze(mean(tso2(ss,HPCC{pp},cc+1,:),2)); %extracting ROI ss, time-window of higher-order ROIs, condition NTcc, all participants
            end
        end
        if cc < 4
            cntp = cntp + 1;
            %ANOVA
            [p,t,stats] = anova1(dumbo); %'off' for not showing the plot
            [c,m,h,nms] = multcompare(stats,'ctype','tukey-kramer'); %perform multiple comparison test based on anova1 output = post-hoc analysis (c needs to be saved for every time point and every channel)
            SIGN{ss,cc,1} = t; %storing ANOVA statistics
            SIGN{ss,cc,2} = c; %storing post-hoc tests
            PP(cntp) = p;
            ccell = num2cell(c);
            PDn = cell2table(ccell); %table
            writetable(PDn,['/aux/MINDLAB2021_MEG-TempSeqAges/03_11_2023/NT' num2str(cc) '_ROI_' ROIN{ss} '.xlsx'],'Sheet',1); %printing excel file
            PDn = table(dumbo);
            writetable(PDn,['/aux/MINDLAB2021_MEG-TempSeqAges/09_11_2023/Violin_NT' num2str(cc) '_ROI_' ROIN{ss} '.xlsx'],'Sheet',1); %printing excel file
        end
        %plotting
        cl = [0 0 0.8; 0 0 0.8; 0 0 0.8; 0 0 0.8; 0 0 0.8;];
        %         dumbo = dumbo';
        datadum = cell(1,size(dumbo,2));
        for ii = 1:size(dumbo,2)
            datadum{1,ii} = dumbo(:,ii);
        end
        figure
        rm_raincloud2(datadum',cl)
        grid minor
        set(gcf,'color','w')
        set(gcf,'Position',[200,200,400,550])
        xlabel('Amplitude'); %set(gca,'YDir','normal');
        title(['ROI ' ROIN{ss} ' - NT' num2str(cc)])
%         export_fig(['/aux/MINDLAB2021_MEG-TempSeqAges/03_11_2023/NT' num2str(cc) '_ROI_' ROIN{ss} '.eps'])
%         export_fig(['/aux/MINDLAB2021_MEG-TempSeqAges/03_11_2023/NT' num2str(cc) '_ROI_' ROIN{ss} '.png'])
    end
end
%FDR correction for multiple comparisons
[pthr,pcor,padj] = fdr(PP);

%%

%% *** DEFINING FUNCTIONAL PARCELS ***

%%% Originally, I tried the functional k-means clustering, then I decided
%%% to proceed with a simpler solution, as described below

%%% To be noted, this is only a supplementary and confirmatory analysis for
%%% this study and the results and figures related to it are entirely 
%%% reported in the supplementary information

%% ORIGINAL SOLUTION

%% *** SIMPLER SOLUTION TO GET THE FUNCTIONAL PARCELS (ROIs) ***

%%% 1) Getting key time-points of different activity between Old and New (and a short time-window around them).
%%%    This were based on the results obtained from the decdding and from the inspection of the prototypical MEG0211 channel.
%%% 2) Thresholding them (e.g. with a high t-value of 3 or even higher)
%%% 3) Running cluster-based MCS to remove spurious brain voxels who were randomly involved in steps 1) and 2)
%%% 4) In some cases, when the activity was left or right lateralized, getting the mirror ROI in the opposite hemisphere 

%%

%% ACTUAL COMPUTATION

% 1) Getting peaks
%%% The peaks of the different activity between Old and New corresponded to:
% - slow peak after the second (and third and fourth) tone of the sequence
% - slow peak after the last tone
% - later prediction error peak (NewT2) - hippocampus and VMPFC
% - first prediction error peak (NewT1) - auditory cortex

%0)output directory
outdir = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/Beam_abs_1_sens_1_freq_broadband/DCM/ROIs_MainActivity';

%1)computing ROIs
%loading data
load('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/Beam_abs_1_sens_1_freq_broadband/sources_contrast_1.mat');
Vstat2 = t_val_s;
names{1} = 'MC'; names{2} = 'HITR'; names{3} = 'HITL'; names{4} = 'VMPFC'; names{5} = 'ACL'; names{6} = 'ACR';
%1-MC
%defining peaks (using time-windows of + and - 5 time-points (+ and - 20ms)) and getting voxels above threshold in those peaks, independently for each time-point (taking all voxels with at least one significant time-point)
MC = zeros(size(Vstat2,1),11); %preparing dummy variable
MC(Vstat2(:,200:210) > 3) = 1; %getting ones when values are above the threshold (t-value = 3 here)
MC = sum(MC,2); %summing voxels to get how many times they were significant
MC(MC~=0) = 1; %if a voxel was significant at least one time, we take it
ROIs{1} = MC;
%2-HITR
%defining peaks (using time-windows of + and - 5 time-points (+ and - 20ms)) and getting voxels above threshold in those peaks, independently for each time-point (taking all voxels with at least one significant time-point)
HITR = zeros(size(Vstat2,1),11); %preparing dummy variable
HITR(Vstat2(:,483:493) > 3.2) = 1; %getting ones when values are above the threshold (t-value = 3 here)
HITR = sum(HITR,2); %summing voxels to get how many times they were significant
HITR(HITR~=0) = 1; %if a voxel was significant at least one time, we take it
ROIs{2} = HITR;
%5-ACL
%defining peaks (using time-windows of + and - 5 time-points (+ and - 20ms)) and getting voxels above threshold in those peaks, independently for each time-point (taking all voxels with at least one significant time-point)
ACL = zeros(size(Vstat2,1),11); %preparing dummy variable
ACL(Vstat2(:,159:169) < -3) = 1; %getting ones when values are above the threshold (t-value = 3 here)
ACL = sum(ACL,2); %summing voxels to get how many times they were significant
ACL(ACL~=0) = 1; %if a voxel was significant at least one time, we take it
ROIs{5} = ACL;
%4-VMPFC
%loading data
load('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/Beam_abs_1_sens_1_freq_broadband/sources_contrast_2.mat');
Vstat2 = t_val_s;
%defining peaks (using time-windows of + and - 5 time-points (+ and - 20ms)) and getting voxels above threshold in those peaks, independently for each time-point (taking all voxels with at least one significant time-point)
VMPFC = zeros(size(Vstat2,1),11); %preparing dummy variable
VMPFC(Vstat2(:,269:279) < -3.7) = 1; %getting ones when values are above the threshold (t-value = 3 here)
VMPFC = sum(VMPFC,2); %summing voxels to get how many times they were significant
VMPFC(VMPFC~=0) = 1; %if a voxel was significant at least one time, we take it
ROIs{4} = VMPFC;

%2)creating nifti images of computed ROIs
maskk = load_nii('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/MNI152_8mm_brain_diy.nii.gz'); %getting the mask for creating the figure
for iii = 1:length(ROIs)
    fnamenii = [outdir '/' names{iii} '.nii.gz']; %path and name of the image to be saved
    SO = ROIs{iii};
    if ~isempty(SO)
        %building nifti image
        SS = size(maskk.img);
        dumimg = zeros(SS(1),SS(2),SS(3),1); %no time here so size of 4D = 1
        for ii = 1:size(SO,1) %over brain sources
            dum = find(maskk.img == ii); %finding index of sources ii in mask image (MNI152_8mm_brain_diy.nii.gz)
            [i1,i2,i3] = ind2sub([SS(1),SS(2),SS(3)],dum); %getting subscript in 3D from index
            dumimg(i1,i2,i3,:) = SO(ii,:); %storing values for all time-points in the image matrix
        end
        nii = make_nii(dumimg,[8 8 8]); %making nifti image from 3D data matrix (and specifying the 8 mm of resolution)
        nii.img = dumimg; %storing matrix within image structure
        nii.hdr.hist = maskk.hdr.hist; %copying some information from maskk
        disp(['saving nifti image - ROI ' num2str(iii)])
        save_nii(nii,fnamenii); %printing image
    end
end

%3)cluster-based MCS to clean things up (either by removing random voxels automatically or by selecting only the biggest cluster)
%Cluster permutations test
%loading main image 
list = dir('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/Beam_abs_1_sens_1_freq_broadband/DCM/ROIs_MainActivity/*gz');
%OBS! you may get a warning since the skeletonized image is not exactly in MNI space, but close enough
[ mni_coords, xform ] = osl_mnimask2mnicoords('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/MNI152_T1_8mm_brain.nii.gz');
MASK = load_nii('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/MNI152_T1_8mm_brain.nii.gz');
%preparing function for MCS
for ii = 1:length(list) %over time-windows
    %getting MNI coordinates
    T = load_nii([list(ii).folder '/' list(ii).name]);
    %extracting matrix with statistics
    T2 = T.img(:,:,:); %extracting time-window ii
    %mask for non-0 voxels in brain imges (basically building a layout for actual brain voxels)
    mask = zeros(size(T2,1),size(T2,2),size(T2,3));
    data = T2;
    mask(MASK.img~=0) = 1; %assigning 1 when you have real brain voxels
    %preparation of information and data for the actual function
    S = [];
    S.T = T2;
    S.outdir = ['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/Beam_abs_1_sens_1_freq_broadband/DCM/ROIs_MainActivity/MCS']; %output path
    S.parcelfile = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/aal_8mm_try5.nii.gz';
    S.labels = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/AAL_labels.mat';
    S.MNIcoords = mni_coords; %MNI coordinates of 8mm MNI152T1 brain
    S.data = data;
    S.mask = mask; %mask of the brain layout you have your results in
    %THE NUMBERS FOR THE PERMUTATIONS HAVE BEEN CHANGED IN DIFFERENT BRAIN AREAS.. DOESN'T REALLY MATTER SINCE HERE I JUST WANT TO GET THE PROPER ROIs AND NOT DOING ANY STATISTICAL TEST
    S.permut = 1000; %number of permutations for Monte Carlo simulation
    S.clustmax = 1; %set 1 for only max cluster size of each permutation MCS (more strict); set 0 for every size of each cluster detected for each permutation MCS (less strict).
    S.permthresh = 0.001; %threshold for MCS
    %UNTIL HERE
    S.anal_name = [list(ii).name(1:end-7) '_MCS']; %name for the analysis (used to identify and save image and results)
    
    %actual function
    PP = BrainSources_MonteCarlosim_3D_LBPD_D(S);
    
    disp(ii)
end

%%

%% *** ADDITIONAL ANALYSES ON/FOR THE FUNCTIONAL PARCELLATION (in supplementary materials) ***

%% USING K-MEANS CLUSTERING ON SPATIAL MNI COORDINATES TO ESTABLISH HOW MANY MAIN BROAD BRAIN REGIONS WERE IDENTIFIED BY THE PREVIOUS ANALYSIS

load('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/Papers/Systematic_Variations_DCM/NatureCommunications_Revision/Revision_I/MEGSensors_FunctionalParcellation/MNI_tval.mat');

SF = zeros(1,20);
for ff = 1:length(SF) %over spatial cluster solutions
    [idk2f,~,smdf] = kmeans(MNI,ff,'Replicates',100); %kmeans
    SF(1,ff) = sum(smdf); %saving sum of square distances
    disp(ff)
end

%silhouette to establish best clustering solution
eva = evalclusters(MNI,'kmeans','silhouette','KList',1:length(SF));
eva.OptimalK
figure
plot(eva)
grid minor
set(gcf,'Color','w')
%best clustering solution is 4 clusters
[idk,Cf,smdf] = kmeans(MNI,4,'Replicates',100); %kmeans
%please, note that every time you run the k-means clustering you shoudl get identical or almost identical results.
%However, the ordinal number of the clusters changes on everyrun since the start of the k-means clustering is random by definition.

%plotting
% 1) plotting all selected brain voxels with the same color
openfig('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/BrainTemplate_GT.fig')
hold on
for ii = 1:size(MNI,1)
    plot3(MNI(ii, 1), MNI(ii, 2), MNI(ii, 3), '.','Color','k', 'MarkerSize', 10); %plotting each voxel
    hold on
end
rotate3d on;
axis off
axis vis3d
axis equal
set(gcf,'Color','w')

view([0 90])
export_fig('black_top.png')
export_fig('black_top.eps')
view([-90 0])
export_fig('black_left.png')
export_fig('black_left.eps')
view([90 0])
export_fig('black_right.png')
export_fig('black_right.eps')
view([0 -90])
export_fig('black_top2.png')
export_fig('black_top2.eps')
view([-180 0])
export_fig('black_back.png')
export_fig('black_back.eps')
view([0 0])
export_fig('black_front.png')
export_fig('black_front.eps')

% 2) plotting clusters in 3D
openfig('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/BrainTemplate_GT.fig')
hold on
Col = ['r','b','c','y'];
% for ii = 1:size(Cf,1) %over clusters
%     plot3(Cf(ii, 1), Cf(ii, 2), Cf(ii, 3), '.','Color',Col(ii), 'MarkerSize', 40); %plotting centroids
%     hold on
% end
for ii = 1:size(MNI,1)
    plot3(MNI(ii, 1), MNI(ii, 2), MNI(ii, 3), '.','Color',Col(idk(ii)), 'MarkerSize', 10); %plotting each voxel
    hold on
end
rotate3d on;
axis off
axis vis3d
axis equal
set(gcf,'Color','w')

view([0 90])
export_fig('top.png')
export_fig('top.eps')
view([-90 0])
export_fig('left.png')
export_fig('left.eps')
view([90 0])
export_fig('right.png')
export_fig('right.eps')
view([0 -90])
export_fig('top2.png')
export_fig('top2.eps')
view([-180 0])
export_fig('back.png')
export_fig('back.eps')
view([0 0])
export_fig('front.png')
export_fig('front.eps')

%not reported in the paper
%producing four nifti images of the clusterised voxels
%getting proper coordinates
S = [];
S.input = 1; %1 = MNI coordinates; 2 = AAL ROIs; 3 = general image with non-zero values
S.coordd = MNI; %coordinates in MNI space (x,y,z)
S.AAL_ROIs = []; %AAL ROIs numbers you want to use
% S.image = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/Beam_abs_0_sens_1_freq_broadband/Contr_1_abs_0.nii.gz';
S.image = '/scratch7/MINDLAB2017_MEG-LearningBach/DTI_Portis/Templates/parcel_80mm_3559ROIs/3000.nii.gz';
%actual function
idx_LBPD = From3DNifti_OrMNICoords_2_CoordMatrix_8mm_LBPD_D(S);
MNII = zeros(3559,4);
for ii = 1:length(idx_LBPD)
    if idk(ii) == 1
        MNII(idx_LBPD(ii),1) = 1;
    elseif idk(ii) == 2
        MNII(idx_LBPD(ii),2) = 1;
    elseif idk(ii) == 3
        MNII(idx_LBPD(ii),3) = 1;
    elseif idk(ii) == 4
        MNII(idx_LBPD(ii),4) = 1;
    end
end
%producing nifti image
S = [];
S.data = MNII; %data (voxels x ROIs (time-points))
S.fname = '/aux/MINDLAB2021_MEG-TempSeqAges/14_09_2023/Kmeans'; %path and name for the image to be saved
S.names{1} = 'MC_k'; S.names{2} = 'RHIT_k'; S.names{3} = 'VMPFC_k'; S.names{4} = 'LAC_k'; %names for the different images
names = S.names;
FromCoordMatrix_2_3DNifti_8mm_LBPD_D(S);

%% 3D plotting of 4 broad ROIs

% Figure S14

list = dir('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/MainROIs/*GM.nii.gz');
ROIs_4 = [5 4 6 1 2 3]; %ROIs from 4-broad-ROIs

% %loading original parcellation
% load('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/Beam_abs_1_sens_1_freq_broadband/DCM/ROIs_MainActivity/FinalROIs/ROIsDCM.mat')
% %loading MNI coordinates for LBPD indices
load('/scratch7/MINDLAB2017_MEG-LearningBach/DTI_Portis/Templates/MNI152_8mm_coord_dyi.mat');
% ROIs_4 = [1 2 4 5]; %ROIs from 4-broad-ROIs
Col = ['r','b','c','y'];
for iii = 1:4
    
    %FROM 3D NIFTI IMAGE (OR ALL ROIs OR MNI COORDINATES) TO LBPD COORDINATES
    S = [];
    S.input = 3; %1 = MNI coordinates; 2 = AAL ROIs; 3 = general image with non-zero values
    S.coordd = []; %coordinates in MNI space (x,y,z)
    S.AAL_ROIs = []; %AAL ROIs numbers you want to use
    % S.image = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/Beam_abs_0_sens_1_freq_broadband/Contr_1_abs_0.nii.gz';
    S.image = [list(ROIs_4(iii)).folder '/' list(ROIs_4(iii)).name];
    %actual function
    idx_LBPD = From3DNifti_OrMNICoords_2_CoordMatrix_8mm_LBPD_D(S);
    
    MNI = MNI8(idx_LBPD,:);
    openfig('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/BrainTemplate_GT.fig')
    hold on
    for ii = 1:size(MNI,1)
        plot3(MNI(ii, 1), MNI(ii, 2), MNI(ii, 3), '.','Color',Col(iii), 'MarkerSize', 10); %plotting each voxel
        hold on
    end
    rotate3d on;
    axis off
    axis vis3d
    axis equal
    set(gcf,'Color','w')
    view([0 90])
    % export_fig([names{ROIs_4(ii)} '_4ROIs_top.png'])
    % export_fig([names{ROIs_4(ii)} '_4ROIs_top.eps'])
    % view([-90 0])
    % export_fig([names{ROIs_4(ii)} '_4ROIs_left.png'])
    % export_fig([names{ROIs_4(ii)} '_4ROIs_left.eps'])
    % view([90 0])
    % export_fig([names{ROIs_4(ii)} '_4ROIs_right.png'])
    % export_fig([names{ROIs_4(ii)} '_4ROIs_right.eps'])
    % view([0 -90])
    % export_fig([names{ROIs_4(ii)} '_4ROIs_top2.png'])
    % export_fig([names{ROIs_4(ii)} '_4ROIs_top2.eps'])
    % view([-180 0])
    % export_fig([names{ROIs_4(ii)} '_4ROIs_back.png'])
    % export_fig([names{ROIs_4(ii)} '_4ROIs_back.eps'])
    % view([0 0])
    % export_fig([names{ROIs_4(ii)} '_4ROIs_front.png'])
    % export_fig([names{ROIs_4(ii)} '_4ROIs_front.eps'])
end

%%

%% GETTING FINAL ROIs (MAIN CLUSTER) AND ROIs MIRRORED IN THE OTHER HEMISPHERE WHEN THE ORIGINAL ROIs WERE LEFT- OR RIGHT-LATERALIZED

list = dir('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/Beam_abs_1_sens_1_freq_broadband/DCM/ROIs_MainActivity/MCS/*gz');
load('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/MNI152_8mm_coord_dyi.mat'); %loading new cooridnate system..
maskk = load_nii('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/MNI152_8mm_brain_diy.nii.gz'); %getting the mask for creating the figure

ROIs = zeros(3559,length(list)); %initializing ROIs
for mm = 1:length(list) %over ROIs
    roi = load_nii([list(mm).folder '/' list(mm).name]); %loading parcel (ROI) mm
    for ii = 1:size(maskk.img,1) %over dimension 1 of the 3D image
        for jj = 1:size(maskk.img,2) %over dimension 2 of the 3D image
            for zz = 1:size(maskk.img,3) %over dimension 3 of the 3D image
                if roi.img(ii,jj,zz) ~= 0 %if there is a voxel forming the ROI
                    ROIs(maskk.img(ii,jj,zz),mm) = 1; %assigning one to that voxel in ROIs (by taking number of voxel from the maskk.img)
                end
            end
        end
    end
end
%list of nifti images for the previously computed ROIs: 1)ACL; 2)HITR; 3)MC; 4)VMPFC
%3-HITL
ab = find(ROIs(:,2)==1);
HITL = zeros(3559,1);
for ii = 1:length(ab) %over voxels of HITR
    %[MIN,IDX] = min(abs(sum(MNI8-[MNI8(ab(ii),1)*(-1),MNI8(ab(ii),2),MNI8(ab(ii),3)],2)))
    dummy = zeros(3559,3);
    dummy((MNI8(ab(ii),1)*(-1)+4)==MNI8(:,1),1) = 1; %getting values in the other hemisphere (multiplying by -1 o e the sign and adding 4 since there is a very small (negligible) imperfection of the coordinate system (symmetric across hemispheres voxels are slightly misalligned.. nothing that really matters here in MEG)
    dummy((MNI8(ab(ii),2)==MNI8(:,2)),2) = 1;
    dummy((MNI8(ab(ii),3)==MNI8(:,3)),3) = 1;
    HITL((sum(dummy,2)==3),1) = 1;
end
%6-ACR
ab = find(ROIs(:,1)==1);
ACR = zeros(3559,1);
for ii = 1:length(ab) %over voxels of HITR
    %[MIN,IDX] = min(abs(sum(MNI8-[MNI8(ab(ii),1)*(-1),MNI8(ab(ii),2),MNI8(ab(ii),3)],2)))
    dummy = zeros(3559,3);
    dummy((MNI8(ab(ii),1)*(-1)+4)==MNI8(:,1),1) = 1; %getting values in the other hemisphere (multiplying by -1 o reverse the sign and adding 4 since there is a very small (negligible) imperfection of the coordinate system (symmetric across hemispheres voxels are slightly misalligned.. nothing that really matters here in MEG)
    dummy((MNI8(ab(ii),2)==MNI8(:,2)),2) = 1;
    dummy((MNI8(ab(ii),3)==MNI8(:,3)),3) = 1;
    ACR((sum(dummy,2)==3),1) = 1;
end

%reshaping everything in order
%order in ROIs: 1)ACL; 2)HITR; 3)MC; 4)VMPFC
%final order
names{1} = 'MC'; names{2} = 'HITR'; names{3} = 'HITL'; names{4} = 'VMPFC'; names{5} = 'ACL'; names{6} = 'ACR';
ROIs_DCM = zeros(3559,6);
ROIs_DCM(:,1) = ROIs(:,3); %MC
ROIs_DCM(:,2) = ROIs(:,2); %HITR
ROIs_DCM(:,3) = HITL; %HITL
ROIs_DCM(:,4) = ROIs(:,4); %VMPFC
ROIs_DCM(:,5) = ROIs(:,1); %ACL
ROIs_DCM(:,6) = ACR; %ACR
%saving on disk
save('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/Beam_abs_1_sens_1_freq_broadband/DCM/ROIs_MainActivity/FinalROIs/ROIsDCM','ROIs_DCM','names')

%% creating nifti images for the final ROIs (Figure S10)

maskk = load_nii('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/MNI152_8mm_brain_diy.nii.gz'); %getting the mask for creating the figure
load('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/Beam_abs_1_sens_1_freq_broadband/DCM/ROIs_MainActivity/FinalROIs/ROIsDCM.mat');
for iii = 1:size(ROIs_DCM,2)
    fnamenii = ['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/MainROIs/' names{iii} '.nii.gz']; %path and name of the image to be saved
    SO = ROIs_DCM(:,iii);
    if ~isempty(SO)
        %building nifti image
        SS = size(maskk.img);
        dumimg = zeros(SS(1),SS(2),SS(3),1); %no time here so size of 4D = 1
        for ii = 1:size(SO,1) %over brain sources
            dum = find(maskk.img == ii); %finding index of sources ii in mask image (MNI152_8mm_brain_diy.nii.gz)
            [i1,i2,i3] = ind2sub([SS(1),SS(2),SS(3)],dum); %getting subscript in 3D from index
            dumimg(i1,i2,i3,:) = SO(ii,:); %storing values for all time-points in the image matrix
        end
        nii = make_nii(dumimg,[8 8 8]); %making nifti image from 3D data matrix (and specifying the 8 mm of resolution)
        nii.img = dumimg; %storing matrix within image structure
        nii.hdr.hist = maskk.hdr.hist; %copying some information from maskk
        disp(['saving nifti image - condition ' num2str(iii)])
        save_nii(nii,fnamenii); %printing image
    end
end

%% FROM LBPD COORDINATES TO 3D NIFTI IMAGE
%same solution provided with a function

load('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/Beam_abs_1_sens_1_freq_broadband/DCM/ROIs_MainActivity/FinalROIs/ROIsDCM.mat');
S.data = ROIs_DCM; %data (voxels x ROIs (time-points))
S.fname = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/MainROIs/Test'; %path and name for the image to be saved
names{1} = 'MC'; names{2} = 'HITR'; names{3} = 'HITL'; names{4} = 'VMPFC'; names{5} = 'ACL'; names{6} = 'ACR';
S.names = names; %names for the different images
%actual function
FromCoordMatrix_2_3DNifti_8mm_LBPD_D(S);

%% FROM 3D NIFTI IMAGE (OR ALL ROIs OR MNI COORDINATES) TO LBPD COORDINATES
%not needed here but useful for future references

S = [];
S.input = 3; %1 = MNI coordinates; 2 = AAL ROIs; 3 = general image with non-zero values
S.coordd = [38 18 -16; -22 58 -16]; %coordinates in MNI space (x,y,z)
S.AAL_ROIs = [80]; %AAL ROIs numbers you want to use
% S.image = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/Beam_abs_0_sens_1_freq_broadband/Contr_1_abs_0.nii.gz';
S.image = '/scratch7/MINDLAB2017_MEG-LearningBach/DTI_Portis/Templates/parcel_80mm_3559ROIs/3000.nii.gz';
%actual function
idx_LBPD = From3DNifti_OrMNICoords_2_CoordMatrix_8mm_LBPD_D(S);

%% 

%% *** QUANTIFYING HIERARCHY IN THE BRAIN - DYNAMIC CAUSAL MODELLING ***

%%

%% PREPARING DATA I - 6 PARCELS FROM AAL PARCELLATION

tonel = 6; %time-window related to the tone that you want (from 1 to 5; 6 means the whole time-window)

AALROIs = []; %ROIs from AAL
ROIII = {5,6,2,3,4,1}; %selected ROIs; 1 = MC, 2 = HITR, 3 = HITL, 4 = VMPFC, 5 = ACL, 6 = ACR
%list for epoched SPM objects (to be cloned afterwards)
list_spm = dir('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/e*recogminor*.mat');
%list of source reconstructed data (for each subject): SINGLE TRIALS!!
list = dir('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/Beam_abs_0_sens_1_freq_broadband_invers_5/SUBJ*.mat');
% data = zeros(3559,1026,length(list));
load('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/Beam_abs_1_sens_1_freq_broadband/DCM/ROIs_MainActivity/FinalROIs/ROIsDCM.mat')
%getting sign of the voxels based on the aggregated source reconstruction (over participants) - AVERAGED TRIALS
load('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/Beam_abs_0_sens_1_freq_broadband/sources_main_effects.mat')
vect = zeros(3559,1);
for jj = 1:3559 %over brain voxels
    if squeeze(mean(t_val_s(jj,45:52,1),2)) > 0 %if the data in voxel jj is positive during N100 time
        vect(jj,1) = -1; %storing a vector with 1 and -1 to be used for later statistics
    else
        vect(jj,1) = 1;
    end
end
%bulding the SPM objects
for ii = 1:length(list) %over subjects
    disp(['loading source reconstructed data for subject ' num2str(ii) ' / ' num2str(length(list))])
    clear OUT
    load([list(ii).folder '/' list(ii).name])
    duMM = [];
    condds = zeros(5,1);
    for cc = 1:length(OUT.sources_ERFs) %over conditions
        dum = OUT.sources_ERFs{cc};
        for jj = 1:size(dum,1) %over brain voxels
            for zz = 1:size(dum,3) %over trials
                dum(jj,:,zz) = dum(jj,:,zz) .* vect(jj,1); %reversing (or not)..
            end
            disp(['subject - ' num2str(ii) ' - condition ' num2str(cc) ' - source ' num2str(jj)])
        end
        if tonel == 2 || tonel == 3 || tonel == 4 || tonel == 5 %if we want a subset of the data (e.g. only time-window for tone 2)
            if tonel == 2 %time-window for tone x
                extr = [89 226];
            elseif tonel == 3
                extr = [176 314];
            elseif tonel == 4
                extr = [264 401];
            elseif tonel == 5
                extr = [351 489];
            end
            dumdum = dum(:,extr(1):extr(2),:); %extracting time-points
%             dumudum2 = dumdum - mean(dum(:,1:XXXXXX,:),2); %baseline correction %for now no additional baseline correction
            dum = dumdum; %dum is the new data in the requested time-window
        end
        dum2 = zeros(length(ROIII),size(dum,2),size(dum,3));
        for pp = 1:length(ROIII) %over ROIs
            S = [];
            S.input = 2; %1 = MNI coordinates; 2 = AAL ROIs; 3 = general image with non-zero values
            S.coordd = []; %coordinates in MNI space (x,y,z)
            S.AAL_ROIs = [AALROIs(pp)]; %AAL ROIs numbers you want to use
            % S.image = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/Beam_abs_0_sens_1_freq_broadband/Contr_1_abs_0.nii.gz';
            S.image = [];
            %actual function
            idx_LBPD = From3DNifti_OrMNICoords_2_CoordMatrix_8mm_LBPD_D(S);
            dum2(pp,:,:) = mean(dum(idx_LBPD,:,:),1);
%             dum2(pp,:,:) = mean(dum(sum(ROIs_DCM(:,ROIII{pp}),2)~=0,:,:),1); %getting the indices of voxels of both ROIs (if you want two ROIs)            
        end
        duMM = cat(3,duMM,dum2);
        condds(cc,1) = size(dum,3); %storing number of correct trials for each conditions;
    end
    disp(['building new SPM object for subject ' num2str(ii) ' / ' num2str(length(list))])
    D = spm_eeg_load([list_spm(ii).folder '/' list_spm(ii).name]);
    %hard coding here..
    if tonel == 6
        Dnew = clone(D, ['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/Papers/Systematic_Variations_DCM/NatureCommunications_Revision/Revision_I/DCM_AAL/DataforDCM_' fname(D)], [length(ROIII), [], 1026, size(duMM,3)]);
    else
        Dnew = clone(D, ['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/Papers/Systematic_Variations_DCM/NatureCommunications_Revision/Revision_I/DCM_AAL/DataforDCM_tone' num2str(tonel) '_' fname(D)], [length(ROIII), [], size(dum,2), size(duMM,3)]);
    end
    Dnew = chanlabels(Dnew, 1, 'ACL'); Dnew = chanlabels(Dnew, 2, 'ACR'); Dnew = chanlabels(Dnew, 3, 'HITL'); Dnew = chanlabels(Dnew, 4, 'HITR'); Dnew = chanlabels(Dnew, 5, 'VMPFC'); Dnew = chanlabels(Dnew, 6, 'MC');
    Dnew = chantype(Dnew, 1:5, 'LFP');
    %barbaric solution to assign proper conditions.. since it seemed that other apparently more reasonable solutions were not working
    for qq = 1:condds(1)
        Dnew = Dnew.conditions(qq,'Old_Correct');
    end
    dumo = cumsum(condds(1:2),1);
    for qq = dumo(1)+1:dumo(end)
        Dnew = Dnew.conditions(qq,'NewT1_Correct');
    end
    dumo = cumsum(condds(1:3),1);
    for qq = dumo(2)+1:dumo(end)
        Dnew = Dnew.conditions(qq,'NewT2_Correct');
    end
    dumo = cumsum(condds(1:4),1);
    for qq = dumo(3)+1:dumo(end)
        Dnew = Dnew.conditions(qq,'NewT3_Correct');
    end
    dumo = cumsum(condds(1:5),1);
    for qq = dumo(4)+1:dumo(end)
        Dnew = Dnew.conditions(qq,'NewT4_Correct');
    end
    Dnew(:,:,:) = duMM; %actual data
    Dnew.save();
    clear Dnew D duMM dum dum2
end

%% PREPARING DATA II (FUNCTIONAL PARCELLATION)

tonel = 2; %time-window related to the tone that you want (from 1 to 5; 6 means the whole time-window)

ROIII = {5,6,3,2,4,1}; %selected ROIs; 1 = MC, 2 = HITR, 3 = HITL, 4 = VMPFC, 5 = ACL, 6 = ACR
%list for epoched SPM objects (to be cloned afterwards)
list_spm = dir('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/e*recogminor*.mat');
%list of source reconstructed data (for each subject): SINGLE TRIALS!!
list = dir('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/Beam_abs_0_sens_1_freq_broadband_invers_5/SUBJ*.mat');
% data = zeros(3559,1026,length(list));
load('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/Beam_abs_1_sens_1_freq_broadband/DCM/ROIs_MainActivity/FinalROIs/ROIsDCM.mat')
%getting sign of the voxels based on the aggregated source reconstruction (over participants) - AVERAGED TRIALS
load('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/Beam_abs_0_sens_1_freq_broadband/sources_main_effects.mat')
vect = zeros(3559,1);
for jj = 1:3559 %over brain voxels
    if squeeze(mean(t_val_s(jj,45:52,1),2)) > 0 %if the data in voxel jj is positive during N100 time
        vect(jj,1) = -1; %storing a vector with 1 and -1 to be used for later statistics
    else
        vect(jj,1) = 1;
    end
end
%bulding the SPM objects
for ii = 1:length(list) %over subjects
    disp(['loading source reconstructed data for subject ' num2str(ii) ' / ' num2str(length(list))])
    clear OUT
    load([list(ii).folder '/' list(ii).name])
    duMM = [];
    condds = zeros(5,1);
    for cc = 1:length(OUT.sources_ERFs) %over conditions
        dum = OUT.sources_ERFs{cc};
        for jj = 1:size(dum,1) %over brain voxels
            for zz = 1:size(dum,3) %over trials
                dum(jj,:,zz) = dum(jj,:,zz) .* vect(jj,1); %reversing (or not)..
            end
            disp(['subject - ' num2str(ii) ' - condition ' num2str(cc) ' - source ' num2str(jj)])
        end
        if tonel == 2 || tonel == 3 || tonel == 4 || tonel == 5 %if we want a subset of the data (e.g. only time-window for tone 2)
            if tonel == 2 %time-window for tone x
                extr = [89 226];
            elseif tonel == 3
                extr = [176 314];
            elseif tonel == 4
                extr = [264 401];
            elseif tonel == 5
                extr = [351 489];
            end
            dumdum = dum(:,extr(1):extr(2),:); %extracting time-points
%             dumudum2 = dumdum - mean(dum(:,1:XXXXXX,:),2); %baseline correction %for now no additional baseline correction
            dum = dumdum; %dum is the new data in the requested time-window
        end
        dum2 = zeros(length(ROIII),size(dum,2),size(dum,3));
        for pp = 1:length(ROIII) %over ROIs
            dum2(pp,:,:) = mean(dum(sum(ROIs_DCM(:,ROIII{pp}),2)~=0,:,:),1); %getting the indices of voxels of both ROIs (if you want two ROIs)
        end
        duMM = cat(3,duMM,dum2);
        condds(cc,1) = size(dum,3); %storing number of correct trials for each conditions;
    end
    disp(['building new SPM object for subject ' num2str(ii) ' / ' num2str(length(list))])
    D = spm_eeg_load([list_spm(ii).folder '/' list_spm(ii).name]);
    %hard coding here..
    if tonel == 6
        Dnew = clone(D, ['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/Papers/Systematic_Variations_DCM/NatureCommunications_Revision/Revision_I/DCM_swap/DataforDCM_' fname(D)], [length(ROIII), [], 1026, size(duMM,3)]);
    else
        Dnew = clone(D, ['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/Papers/Systematic_Variations_DCM/NatureCommunications_Revision/Revision_I/DCM_swap/DataforDCM_tone' num2str(tonel) '_' fname(D)], [length(ROIII), [], size(dum,2), size(duMM,3)]);
    end
    Dnew = chanlabels(Dnew, 1, 'ACL'); Dnew = chanlabels(Dnew, 2, 'ACR'); Dnew = chanlabels(Dnew, 3, 'HITL'); Dnew = chanlabels(Dnew, 4, 'HITR'); Dnew = chanlabels(Dnew, 5, 'VMPFC'); Dnew = chanlabels(Dnew, 6, 'MC');
    Dnew = chantype(Dnew, 1:5, 'LFP');
    %barbaric solution to assign proper conditions.. since it seemed that other apparently more reasonable solutions were not working
    for qq = 1:condds(1)
        Dnew = Dnew.conditions(qq,'Old_Correct');
    end
    dumo = cumsum(condds(1:2),1);
    for qq = dumo(1)+1:dumo(end)
        Dnew = Dnew.conditions(qq,'NewT1_Correct');
    end
    dumo = cumsum(condds(1:3),1);
    for qq = dumo(2)+1:dumo(end)
        Dnew = Dnew.conditions(qq,'NewT2_Correct');
    end
    dumo = cumsum(condds(1:4),1);
    for qq = dumo(3)+1:dumo(end)
        Dnew = Dnew.conditions(qq,'NewT3_Correct');
    end
    dumo = cumsum(condds(1:5),1);
    for qq = dumo(4)+1:dumo(end)
        Dnew = Dnew.conditions(qq,'NewT4_Correct');
    end
    Dnew(:,:,:) = duMM; %actual data
    Dnew.save();
    clear Dnew D duMM dum dum2
end

%% DOING THE SAME AS ABOVE BUT ON THE CLUSTER OF COMPUTERS IN AARHUS

% config of the cluster server
clusterconfig('slot', 1); % set manually the job cluster slots
% between 1 and 12 (n x 8gb of ram)
clusterconfig('scheduler','none'); % set automatically the long run queue
clusterconfig('long_running', 1); % set automatically the long run queue

%% CASE I - AAL

addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/Cluster_ParallelComputing');
S = [];
S.tonel = 5; %time-window related to the tone that you want (from 1 to 5; 6 means the whole time-window)
AALROIsdum = [79 80 37 38 31 32 33 34]; %ROIs from AAL (L Hesch, R Hesch, L Hipp, R Hipp, L Cing Ant, R Cing Ant (to be combined), L Cing Mid, R Cing Mid (to be combined))

ROIII = {5,6,2,3,4,1}; %selected ROIs; 1 = MC, 2 = HITR, 3 = HITL, 4 = VMPFC, 5 = ACL, 6 = ACR
%list for epoched SPM objects (to be cloned afterwards)
list_spm = dir('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/e*recogminor*.mat');
%list of source reconstructed data (for each subject): SINGLE TRIALS!!
list = dir('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/Beam_abs_0_sens_1_freq_broadband_invers_5/SUBJ*.mat');
% data = zeros(3559,1026,length(list));
load('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/Beam_abs_1_sens_1_freq_broadband/DCM/ROIs_MainActivity/FinalROIs/ROIsDCM.mat')
%getting sign of the voxels based on the aggregated source reconstruction (over participants) - AVERAGED TRIALS
load('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/Beam_abs_0_sens_1_freq_broadband/sources_main_effects.mat')
vect = zeros(3559,1);
for jj = 1:3559 %over brain voxels
    if squeeze(mean(t_val_s(jj,45:52,1),2)) > 0 %if the data in voxel jj is positive during N100 time
        vect(jj,1) = -1; %storing a vector with 1 and -1 to be used for later statistics
    else
        vect(jj,1) = 1;
    end
end
%building structure with inputs
S.vect = vect;
S.ROIII = ROIII;
S.list_spm = list_spm;
S.list = list;

AALROIs = zeros(length(vect),6);
cnt = 4;
for pp = 1:6 %over ROIs
    S2 = [];
    S2.input = 2; %1 = MNI coordinates; 2 = AAL ROIs; 3 = general image with non-zero values
    S2.coordd = []; %coordinates in MNI space (x,y,z)
    S2.image = [];
    if pp > 4 %if we have medial cingulate and anterior cingulate/VMPFC
        dumdum = [];
        for ii = 1:2 %over L and R hemisphere
            cnt = cnt + 1;
            S2.AAL_ROIs = [AALROIsdum(cnt)]; %AAL ROIs numbers you want to use
            %actual function
            dumdum = cat(1,dumdum,From3DNifti_OrMNICoords_2_CoordMatrix_8mm_LBPD_D(S2));
        end
        dumbo = dumdum;
    else
        S2.AAL_ROIs = [AALROIsdum(pp)]; %AAL ROIs numbers you want to use
        dumbo = From3DNifti_OrMNICoords_2_CoordMatrix_8mm_LBPD_D(S2); %otherwise only single AAL ROI voxels
    end
    AALROIs(dumbo,pp) = 1; %assigning 1 to the proper voxels
end
S.AALROIs = AALROIs;

%bulding the SPM objects
for ii = 1:length(list) %over subjects
    S.ii = ii;
    jobid = job2cluster(@BuildingDataforDCM_NotGeneral_AAL, S);
end

%% CASE II - FUNCTIONAL PARCELLATION

addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/Cluster_ParallelComputing');

S = [];
S.tonel = 5; %time-window related to the tone that you want (from 1 to 5; 6 means the whole time-window)

ROIII = {5,6,3,2,4,1}; %selected ROIs; 1 = MC, 2 = HITR, 3 = HITL, 4 = VMPFC, 5 = ACL, 6 = ACR
%list for epoched SPM objects (to be cloned afterwards)
list_spm = dir('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/e*recogminor*.mat');
%list of source reconstructed data (for each subject): SINGLE TRIALS!!
list = dir('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/Beam_abs_0_sens_1_freq_broadband_invers_5/SUBJ*.mat');
% data = zeros(3559,1026,length(list));
load('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/Beam_abs_1_sens_1_freq_broadband/DCM/ROIs_MainActivity/FinalROIs/ROIsDCM.mat')
%getting sign of the voxels based on the aggregated source reconstruction (over participants) - AVERAGED TRIALS
load('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/Beam_abs_0_sens_1_freq_broadband/sources_main_effects.mat')
vect = zeros(3559,1);
for jj = 1:3559 %over brain voxels
    if squeeze(mean(t_val_s(jj,45:52,1),2)) > 0 %if the data in voxel jj is positive during N100 time
        vect(jj,1) = -1; %storing a vector with 1 and -1 to be used for later statistics
    else
        vect(jj,1) = 1;
    end
end
%building structure with inputs
S.ROIs_DCM = ROIs_DCM;
S.vect = vect;
S.ROIII = ROIII;
S.list_spm = list_spm;
S.list = list;

%bulding the SPM objects
for ii = 4%1:length(list) %over subjects
    S.ii = ii;
    jobid = job2cluster(@BuildingDataforDCM_NotGeneral_swap, S);
end

%% ACTUAL INVERSION OF THE DCM MODELS

%% setting up the cluster

addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/Cluster_ParallelComputing') %add the path to the function that submits the jobs to the cluster
clusterconfig('scheduler', 'cluster'); %'none' or 'cluster'
clusterconfig('long_running', 1); %there are different queues for the cluster depending on the number and length of the jobs you want to submit 
clusterconfig('slot', 1); %slot in the queu
addpath('/projects/MINDLAB2020_MEG-AuditoryPatternRecognition/scripts/leonardo/DCM')

%% RESHAPING CODES FROM MARTIN DIETZ FOR DCM FOR EVOKED RESPONSES

%% DCM specifications - CASE I - AAL

model_n = 6; %total number of models
modvec = [1:6]; % %models you want to run
subjvect = [1:83]; %selected subjects (from the overall list of subjects prepared for DCM analysis)
tonel = 3; %time-window related to the tone that you want (from 1 to 5; 6 means the whole time-window)
condd = 3; %1 = old > newt1; 2 = newt1 > old; 3 = newt2 > old; 4 = newt3 > old; 5 = newt4 > old
contrl = 1; %1 for old> newt1 or newtx > old; 0 for only old (without contrast/baseline)

% 1 = ACL; 2 = ACR; 3 = HITL; 4 = HITR; 5 = VMPFC; 6 = MC

LC = cell(1,model_n); NLC = cell(1,model_n); INP = cell(1,model_n); BB = cell(1,model_n);

%model specification (model 1) - AC -> HIT and VMPFC and MC (FEEDFORWARD AND BACKWARD) (input to AC)
%1st row = from the four ROIs to the 1st ROIads; 2nd row = from the four ROIs to the 2nd ROI
LC{1} = [0 0 1 1 1 1; 0 0 1 1 1 1; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0]; %feedforward connections
NLC{1} = [0 0 0 0 0 0; 0 0 0 0 0 0; 1 1 0 0 0 0; 1 1 0 0 0 0; 1 1 0 0 0 0; 1 1 0 0 0 0]; %backward connections
INP{1} = [1 1 0 0 0 0]'; %inputs (to ROIs)

%model specification (model 2) - AC -> HIT and VMPFC and MC (FEEDFORWARD ONLY)
%1st row = from the four ROIs to the 1st ROIads; 2nd row = from the four ROIs to the 2nd ROI
LC{2} = [0 0 1 1 1 1; 0 0 1 1 1 1; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0]; %feedforward connections
NLC{2} = [0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0]; %backward connections
INP{2} = [1 1 0 0 0 0]'; %inputs (to ROIs)

%model specification (model 3) - AC -> HIT -> VMPFC and MC (FEEDFORWARD AND BACKWARD)
%1st row = from the four ROIs to the 1st ROIads; 2nd row = from the four ROIs to the 2nd ROI
LC{3} = [0 0 1 1 0 0; 0 0 1 1 0 0; 0 0 0 0 1 1; 0 0 0 0 1 1; 0 0 0 0 0 0; 0 0 0 0 0 0]; %feedforward connections
NLC{3} = [0 0 0 0 0 0; 0 0 0 0 0 0; 1 1 0 0 0 0; 1 1 0 0 0 0; 0 0 1 1 0 0; 0 0 1 1 0 0]; %backward connections
INP{3} = [1 1 0 0 0 0]'; %inputs (to ROIs)

%model specification (model 4 - AC -> VMPFC and MC -> HIT (FEEDFORWARD AND BACKWARD)
%1st row = from the four ROIs to the 1st ROIads; 2nd row = from the four ROIs to the 2nd ROI
LC{4} = [0 0 0 0 1 1; 0 0 0 0 1 1; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 1 1 0 0; 0 0 1 1 0 0]; %feedforward connections
NLC{4} = [0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 1 1; 0 0 0 0 1 1; 1 1 0 0 0 0; 1 1 0 0 0 0]; %backward connections
INP{4} = [1 1 0 0 0 0]'; %inputs (to ROIs)

%model specification (model 5 - AC -> MC -> HIT and VMPFC (FEEDFORWARD AND BACKWARD)
%1st row = from the four ROIs to the 1st ROIads; 2nd row = from the four ROIs to the 2nd ROI
LC{5} = [0 0 0 0 0 1; 0 0 0 0 0 1; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 1 1 1 0]; %feedforward connections
NLC{5} = [0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 1; 0 0 0 0 0 1; 0 0 0 0 0 1; 1 1 0 0 0 0]; %backward connections
INP{5} = [1 1 0 0 0 0]'; %inputs (to ROIs)

%model specification (model 6) - AC -> HIT and VMPFC -> MC (FEEDFORWARD AND BACKWARD)
%1st row = from the four ROIs to the 1st ROIads; 2nd row = from the four ROIs to the 2nd ROI
LC{6} = [0 0 1 1 1 0; 0 0 1 1 1 0; 0 0 0 0 0 1; 0 0 0 0 0 1; 0 0 0 0 0 1; 0 0 0 0 0 0]; %feedforward connections
NLC{6} = [0 0 0 0 0 0; 0 0 0 0 0 0; 1 1 0 0 0 0; 1 1 0 0 0 0; 1 1 0 0 0 0; 0 0 1 1 1 0]; %backward connections
INP{6} = [1 1 0 0 0 0]'; %inputs (to ROIs)

%specifying B
clear B
for ii = 1:model_n %over different models
    B = LC{ii} + NLC{ii}; %getting both feedforward and feedback connections and assuming that you expect them to be different in different conditions
    B(B>0) = ones;
    BB{ii} = B;
end
%output directory
outdir = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/Papers/Systematic_Variations_DCM/NatureCommunications_Revision/Revision_I/DCM_AAL/Out';
mkdir(outdir);
%loading only one subject at the moment
if tonel == 1 || tonel == 6 %the data for only first tone or for all five tones in a row is the same
    list = dir(['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/Papers/Systematic_Variations_DCM/NatureCommunications_Revision/Revision_I/DCM_AAL/DataforDCM_espmeeg*.mat']);
elseif tonel == 2 %otherwise you get the data that was specifically prepared for tone 2
    list = dir(['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/Papers/Systematic_Variations_DCM/NatureCommunications_Revision/Revision_I/DCM_AAL/DataforDCM_tone2*.mat']);
elseif tonel == 3 %or tone 3
    list = dir(['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/Papers/Systematic_Variations_DCM/NatureCommunications_Revision/Revision_I/DCM_AAL/DataforDCM_tone3*.mat']);
elseif tonel == 4 %tone 4
    list = dir(['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/Papers/Systematic_Variations_DCM/NatureCommunications_Revision/Revision_I/DCM_AAL/DataforDCM_tone4*.mat']);
elseif tonel == 5 %tone 5
    list = dir(['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/Papers/Systematic_Variations_DCM/NatureCommunications_Revision/Revision_I/DCM_AAL/DataforDCM_tone5*.mat']);
end

%time window (ms)
if tonel < 6 %because for tones 2,3,4,5 we start 100 ms before the tone of interest (e.g. 2,3,4,5)
    toi = [-100 350];
elseif tonel == 6
    toi = [-100 1750]; 
end

%sources
sour = {'HeschL' 'HeschR' 'HippL' 'HippR' 'VMPFC','MC'};
for subj = 1:length(subjvect)%length(list) %over subjects
    %loading subject
    D = spm_eeg_load([list(subjvect(subj)).folder '/' list(subjvect(subj)).name]);
    dummy = D.fname;
    if tonel == 1 || tonel == 6
        subjJ = dummy(24:27); %subject's name
    else
        subjJ = dummy(30:33); %subject's name
    end
    %call for SPM
    spm('defaults','eeg');
    spm_get_defaults('cmdline',true);
    for ii = 1:length(modvec) %over requested models
        %DCM settings
        clear DCM
        DCM.xY.Dfile    = fullfile(D.path,D.fname); %path to file (SPM object)
        DCM.xY.modality = 'LFP'; %modality
        DCM.xY.name     = sour; %name of the sources
        DCM.xY.Ic       = indchannel(D,sour)'; %channel indices of the sources
        DCM.Lpos        = []; %meaningful only if you run the source reconstruction together with the DCM analysis
        DCM.m = [];
        %DCM options
        DCM.options.analysis = 'ERP'; %evoked responses
        DCM.options.model = 'CMC';
        DCM.options.spatial = 'LFP'; %LFP since I alredy source reconstructed the data
        if contrl == 0
            DCM.options.trials   = [1]; %condition (Old)
        else
            if condd == 1
                DCM.options.trials   = [1 2]; %conditions (here Old and Newt1)
            elseif condd == 2
                DCM.options.trials   = [2 1]; %conditions (here NewT1 and Old)
            elseif condd == 3
                DCM.options.trials   = [3 1]; %conditions (here NewT2 and Old)
            elseif condd == 4
                DCM.options.trials   = [4 1]; %conditions (here NewT3 and Old)
            elseif condd == 5
                DCM.options.trials   = [5 1]; %conditions (here NewT4 and Old)
            end
        end
        DCM.options.Tdcm = toi; %time-window
%         DCM.options.Fdcm = foi; %frequency spectrum
        
%         DCM.options.Rft = 5; %Morlet wavelet coefficient
        DCM.options.onset = 50; %input prior mean (ms) (when the first event of interest in the brain should happen)
        DCM.options.dur = 16; %input prior std (ms) (flexibility of that event to happen..)
        DCM.options.Nmodes = length(DCM.xY.Ic); %SVD frequency modes
        DCM.options.h = 1; %DCT terms
        DCM.options.han = 1; %hanning
        DCM.options.D = 1; %down-sampling
        
        %(I guess) some options for source reconstruction (not needed in this case so set to 0)
%         DCM.options.lock = 0;
%         DCM.options.multiC = 0;
%         DCM.options.location = 0;
%         DCM.options.symmetry = 0;
%         DCM.options.Fmodes = 8;

        % design matrix
        if contrl == 1
            DCM.xU.X = [1 0]'; % condition 1 (Old) > condition 2 (NewT1) or condition 3 (NewT3) > condition 1 (Old)
            %DCM name for output file
            DCM.name = [outdir '/DCM_tone' num2str(tonel) '_EVK_SUBJ' subjJ '_model_' num2str(modvec(ii)) '_cond' num2str(DCM.options.trials(1)) '_vs_cond' num2str(DCM.options.trials(2))]; % subject subj, model j
        else
            DCM.xU.X = []'; % condition 1 (Old) > condition 2 (NewT1) or condition 3 (NewT3) > condition 1 (Old)
            %DCM name for output file
            DCM.name = [outdir '/DCM_tone' num2str(tonel) '_EVK_SUBJ' subjJ '_model_' num2str(modvec(ii)) '_cond' num2str(DCM.options.trials(1)) 'only']; % subject subj, model j
        end
%         DCM.xU.X = [];
        %region names
        DCM.Sname = sour;
        DCM.xU.name = 'name';
        DCM.M.nograph = 1;
        %neuronal model (as specified above)
        DCM.A{1} = LC{modvec(ii)}; %feedforward connections
        DCM.A{2} = NLC{modvec(ii)}; %feedback connections
        DCM.A{3} = zeros(size(DCM.A{1},1));
        DCM.C = INP{modvec(ii)}; %inputs (to ROIs)
        DCM.B{1} = BB{modvec(ii)}; %change in coupling due to condition 1 (in {} parenthesis)
%         DCM.B{1} = []; %change in coupling due to condition 1 (in {} parenthesis)
        
        %actual function for running time-frequency decomposition, model inversion, etc. on the Aarhus cluster
        jobid = job2cluster(@DCM_cluster_evk,DCM); %running with parallel computing
    end
end

%bigger F wins

%% DCM SPECIFICATIONS - CASE II - FUNCTIONAL PARCELLATION

model_n = 6; %total number of models
modvec = [1:6]; % %models you want to run
subjvect = [1:83]; %selected subjects (from the overall list of subjects prepared for DCM analysis)
tonel = 5; %time-window related to the tone that you want (from 1 to 5; 6 means the whole time-window)
condd = 1; %1 = old > newt1; 2 = newt1 > old; 3 = newt2 > old; 4 = newt3 > old; 5 = newt4 > old
contrl = 1; %1 for old> newt1 or newtx > old; 0 for only old (without contrast/baseline)

% 1 = ACL; 2 = ACR; 3 = HITL; 4 = HITR; 5 = VMPFC; 6 = MC

LC = cell(1,model_n); NLC = cell(1,model_n); INP = cell(1,model_n); BB = cell(1,model_n);

%model specification (model 1) - AC -> HIT and VMPFC and MC (FEEDFORWARD AND BACKWARD) (input to AC)
%1st row = from the four ROIs to the 1st ROIads; 2nd row = from the four ROIs to the 2nd ROI
LC{1} = [0 0 1 1 1 1; 0 0 1 1 1 1; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0]; %feedforward connections
NLC{1} = [0 0 0 0 0 0; 0 0 0 0 0 0; 1 1 0 0 0 0; 1 1 0 0 0 0; 1 1 0 0 0 0; 1 1 0 0 0 0]; %backward connections
INP{1} = [1 1 0 0 0 0]'; %inputs (to ROIs)

%model specification (model 2) - AC -> HIT and VMPFC and MC (FEEDFORWARD ONLY)
%1st row = from the four ROIs to the 1st ROIads; 2nd row = from the four ROIs to the 2nd ROI
LC{2} = [0 0 1 1 1 1; 0 0 1 1 1 1; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0]; %feedforward connections
NLC{2} = [0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0]; %backward connections
INP{2} = [1 1 0 0 0 0]'; %inputs (to ROIs)

%model specification (model 3) - AC -> HIT -> VMPFC and MC (FEEDFORWARD AND BACKWARD)
%1st row = from the four ROIs to the 1st ROIads; 2nd row = from the four ROIs to the 2nd ROI
LC{3} = [0 0 1 1 0 0; 0 0 1 1 0 0; 0 0 0 0 1 1; 0 0 0 0 1 1; 0 0 0 0 0 0; 0 0 0 0 0 0]; %feedforward connections
NLC{3} = [0 0 0 0 0 0; 0 0 0 0 0 0; 1 1 0 0 0 0; 1 1 0 0 0 0; 0 0 1 1 0 0; 0 0 1 1 0 0]; %backward connections
INP{3} = [1 1 0 0 0 0]'; %inputs (to ROIs)

%model specification (model 4 - AC -> VMPFC and MC -> HIT (FEEDFORWARD AND BACKWARD)
%1st row = from the four ROIs to the 1st ROIads; 2nd row = from the four ROIs to the 2nd ROI
LC{4} = [0 0 0 0 1 1; 0 0 0 0 1 1; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 1 1 0 0; 0 0 1 1 0 0]; %feedforward connections
NLC{4} = [0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 1 1; 0 0 0 0 1 1; 1 1 0 0 0 0; 1 1 0 0 0 0]; %backward connections
INP{4} = [1 1 0 0 0 0]'; %inputs (to ROIs)

%model specification (model 5 - AC -> MC -> HIT and VMPFC (FEEDFORWARD AND BACKWARD)
%1st row = from the four ROIs to the 1st ROIads; 2nd row = from the four ROIs to the 2nd ROI
LC{5} = [0 0 0 0 0 1; 0 0 0 0 0 1; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 1 1 1 0]; %feedforward connections
NLC{5} = [0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 1; 0 0 0 0 0 1; 0 0 0 0 0 1; 1 1 0 0 0 0]; %backward connections
INP{5} = [1 1 0 0 0 0]'; %inputs (to ROIs)

%model specification (model 6) - AC -> HIT and VMPFC -> MC (FEEDFORWARD AND BACKWARD)
%1st row = from the four ROIs to the 1st ROIads; 2nd row = from the four ROIs to the 2nd ROI
LC{6} = [0 0 1 1 1 0; 0 0 1 1 1 0; 0 0 0 0 0 1; 0 0 0 0 0 1; 0 0 0 0 0 1; 0 0 0 0 0 0]; %feedforward connections
NLC{6} = [0 0 0 0 0 0; 0 0 0 0 0 0; 1 1 0 0 0 0; 1 1 0 0 0 0; 1 1 0 0 0 0; 0 0 1 1 1 0]; %backward connections
INP{6} = [1 1 0 0 0 0]'; %inputs (to ROIs)


%specifying B (THINK ABOUT THAT!!)
for ii = 1:model_n %over different models
    B = LC{ii} + NLC{ii}; %getting both feedforward and feedback connections and assuming that you expect them to be different in different conditions
    B(B>0) = ones;
    BB{ii} = B;
end
%output directory
outdir = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/Papers/Systematic_Variations_DCM/NatureCommunications_Revision/Revision_I/DCM_swap/Out';
mkdir(outdir);
%loading only one subject at the moment
if tonel == 1 || tonel == 6 %the data for only first tone or for all five tones in a row is the same
    list = dir(['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/Papers/Systematic_Variations_DCM/NatureCommunications_Revision/Revision_I/DCM_swap/DataforDCM_espmeeg*.mat']);
elseif tonel == 2 %otherwise you get the data that was specifically prepared for tone 2
    list = dir(['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/Papers/Systematic_Variations_DCM/NatureCommunications_Revision/Revision_I/DCM_swap/DataforDCM_tone2*.mat']);
elseif tonel == 3 %or tone 3
    list = dir(['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/Papers/Systematic_Variations_DCM/NatureCommunications_Revision/Revision_I/DCM_swap/DataforDCM_tone3*.mat']);
elseif tonel == 4 %tone 4
    list = dir(['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/Papers/Systematic_Variations_DCM/NatureCommunications_Revision/Revision_I/DCM_swap/DataforDCM_tone4*.mat']);
elseif tonel == 5 %tone 5
    list = dir(['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/Papers/Systematic_Variations_DCM/NatureCommunications_Revision/Revision_I/DCM_swap/DataforDCM_tone5*.mat']);
end

%time window (ms)
if tonel < 6 %because for tones 2,3,4,5 we start 100 ms before the tone of interest (e.g. 2,3,4,5)
    toi = [-100 350];
elseif tonel == 6
    toi = [-100 1750]; 
end

%sources
sour = {'ACL' 'ACR' 'HITL' 'HITR' 'VMPFC','MC'};
for subj = 1:length(subjvect)%length(list) %over subjects
    %loading subject
    D = spm_eeg_load([list(subjvect(subj)).folder '/' list(subjvect(subj)).name]);
    dummy = D.fname;
    if tonel == 1 || tonel == 6
        subjJ = dummy(24:27); %subject's name
    else
        subjJ = dummy(30:33); %subject's name
    end
    %call for SPM
    spm('defaults','eeg');
    spm_get_defaults('cmdline',true);
    for ii = 1:length(modvec) %over requested models
        %DCM settings
        clear DCM
        DCM.xY.Dfile    = fullfile(D.path,D.fname); %path to file (SPM object)
        DCM.xY.modality = 'LFP'; %modality
        DCM.xY.name     = sour; %name of the sources
        DCM.xY.Ic       = indchannel(D,sour)'; %channel indices of the sources
        DCM.Lpos        = []; %meaningful only if you run the source reconstruction together with the DCM analysis
        DCM.m = [];
        %DCM options
        DCM.options.analysis = 'ERP'; %evoked responses
        DCM.options.model = 'CMC';
        DCM.options.spatial = 'LFP'; %LFP since I alredy source reconstructed the data
        if contrl == 0
            DCM.options.trials   = [1]; %condition (Old)
        else
            if condd == 1
                DCM.options.trials   = [1 2]; %conditions (here Old and Newt1)
            elseif condd == 2
                DCM.options.trials   = [2 1]; %conditions (here NewT1 and Old)
            elseif condd == 3
                DCM.options.trials   = [3 1]; %conditions (here NewT2 and Old)
            elseif condd == 4
                DCM.options.trials   = [4 1]; %conditions (here NewT3 and Old)
            elseif condd == 5
                DCM.options.trials   = [5 1]; %conditions (here NewT4 and Old)
            end
        end
        DCM.options.Tdcm = toi; %time-window
%         DCM.options.Fdcm = foi; %frequency spectrum
        
%         DCM.options.Rft = 5; %Morlet wavelet coefficient
        DCM.options.onset = 50; %input prior mean (ms) (when the first event of interest in the brain should happen)
        DCM.options.dur = 16; %input prior std (ms) (flexibility of that event to happen..)
        DCM.options.Nmodes = length(DCM.xY.Ic); %SVD frequency modes
        DCM.options.h = 1; %DCT terms
        DCM.options.han = 1; %hanning
        DCM.options.D = 1; %down-sampling
        
        %(I guess) some options for source reconstruction (not needed in this case so set to 0)
%         DCM.options.lock = 0;
%         DCM.options.multiC = 0;
%         DCM.options.location = 0;
%         DCM.options.symmetry = 0;
%         DCM.options.Fmodes = 8;

        % design matrix
        if contrl == 1
            DCM.xU.X = [1 0]'; % condition 1 (Old) > condition 2 (NewT1) or condition 3 (NewT3) > condition 1 (Old)
            %DCM name for output file
            DCM.name = [outdir '/DCM_tone' num2str(tonel) '_EVK_SUBJ' subjJ '_model_' num2str(modvec(ii)) '_cond' num2str(DCM.options.trials(1)) '_vs_cond' num2str(DCM.options.trials(2))]; % subject subj, model j
        else
            DCM.xU.X = []'; % condition 1 (Old) > condition 2 (NewT1) or condition 3 (NewT3) > condition 1 (Old)
            %DCM name for output file
            DCM.name = [outdir '/DCM_tone' num2str(tonel) '_EVK_SUBJ' subjJ '_model_' num2str(modvec(ii)) '_cond' num2str(DCM.options.trials(1)) 'only']; % subject subj, model j
        end
%         DCM.xU.X = [];
        %region names
        DCM.Sname = sour;
        DCM.xU.name = 'name';
        DCM.M.nograph = 1;
        %neuronal model (as specified above)
        DCM.A{1} = LC{modvec(ii)}; %feedforward connections
        DCM.A{2} = NLC{modvec(ii)}; %feedback connections
        DCM.A{3} = zeros(size(DCM.A{1},1));
        DCM.C = INP{modvec(ii)}; %inputs (to ROIs)
        DCM.B{1} = BB{modvec(ii)}; %change in coupling due to condition 1 (in {} parenthesis)
%         DCM.B{1} = []; %change in coupling due to condition 1 (in {} parenthesis)
        
        %actual function for running time-frequency decomposition, model inversion, etc. on the Aarhus cluster
        jobid = job2cluster(@DCM_cluster_evk,DCM); %running with parallel computing
    end
end

%bigger F wins

%% MODEL COMPARISON - ACROSS SUBJECTS - CASE I - AAL

% Figure 6

tonel = 5; %time-window related to the tone that you want (from 1 to 5)
model_n = [1;6]; %number of models
ppp = [1:83]; %number of participants (in ascendent order)
condd = 5; %1 = old > newt1; 2 = newt1 > old; 3 = newt2 > old; 4 = newt3 > old; 5 = newt4 > old

%for NTX against M only some specific tones were computed (e.g. tone 2 for NT1 vs M, tone 3 for NT2 vs M, etc.) 
if condd == 2
    tonel = 2;
elseif condd == 3
    tonel = 3;
elseif condd == 4
    tonel = 4;
elseif condd == 5
    tonel = 5;
end
list = dir('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/Papers/Systematic_Variations_DCM/NatureCommunications_Revision/Revision_I/DCM_AAL/Out/DCM_tone2_EVK_S*_4_cond1_vs_cond2.mat');
if condd == 1
    conddd = '_cond1_vs_cond2';
elseif condd == 2
    conddd = '_cond2_vs_cond1';
    tonel = 2;
elseif condd == 3
    conddd = '_cond3_vs_cond1';
    tonel = 3;
elseif condd == 4
    conddd = '_cond4_vs_cond1';
    tonel = 4;
elseif condd == 5
    conddd = '_cond5_vs_cond1';
    tonel = 5;
end
F = zeros(length(ppp),length(model_n)); %to be updated here as well
sbarb = cell(length(ppp),length(model_n));
for ii = 1:length(ppp) %over subjects
    for jj = 1:length(model_n) %over models
        barbb = [list(1).folder '/' list(ppp(ii)).name(1:8) num2str(tonel) list(ppp(ii)).name(10:29) num2str(model_n(jj)) conddd '.mat'];
        load(barbb); %loading data
        sbarb{ii,jj} = barbb; %checking that you loaded the right data
        F(ii,jj) = DCM.F;
        disp(ii)
    end
end
addpath('/projects/MINDLAB2020_MEG-AuditoryPatternRecognition/scripts/leonardo/DCM')
BMS = bms_rfx(F); %computing statistics over participants
mean(F,1);
% [p,tbl,stats] = anova1(F)
BMS
bum = zeros(6,2);
bum(:,1) = BMS.pp'; bum(:,2) = BMS.pep';
outdir = '/aux/MINDLAB2021_MEG-TempSeqAges/09_11_2023';
PDn = table(bum); %table
writetable(PDn,[outdir '/DCM_tone_' num2str(tonel) '_cond_' num2str(condd) '.xlsx'],'Sheet',1); %printing excel file
%plotting for the manuscript
%bar plotting
%posterior (expected) probability
[m,i] = max(BMS.pp);
W = zeros(1,length(BMS.pp));
W(i) = m;
figure
bar(1:6,BMS.pp,0.7,'b')
hold on
bar(W,'k')
ylim([0 1])
xlabel('Models')
ylabel('Probability')
set(gcf,'Color','w')
title('Posterior probability')
export_fig(['/aux/MINDLAB2021_MEG-TempSeqAges/11_10_2023/PP_Tone' num2str(tonel) '_contrast' num2str(condd) '.eps'])
export_fig(['/aux/MINDLAB2021_MEG-TempSeqAges/11_10_2023/PP_Tone' num2str(tonel) '_contrast' num2str(condd) '.png'])
%protected exceedance probability
[m,i] = max(BMS.pep);
W = zeros(1,length(BMS.pep));
W(i) = m;
figure
bar(1:6,BMS.pep,0.7,'b')
hold on
bar(W,'k')
ylim([0 1])
xlabel('Models')
ylabel('Probability')
set(gcf,'Color','w')
title('Protected exceedance probability')
export_fig(['/aux/MINDLAB2021_MEG-TempSeqAges/11_10_2023/PEP_Tone' num2str(tonel) '_contrast' num2str(condd) '.eps'])
export_fig(['/aux/MINDLAB2021_MEG-TempSeqAges/11_10_2023/PEP_Tone' num2str(tonel) '_contrast' num2str(condd) '.png'])

%% MODEL COMPARISON - ACROSS SUBJECTS - CASE II - ORIGINAL PARCELLATION

% Figure S16

tonel = 5; %time-window related to the tone that you want (from 1 to 5)
% model_n = [1:10]; %number of models
% model_n = [1:6]; %number of models
model_n = [7 8 9 4 5 10]; %number of models
ppp = [1:83]; %number of participants (in ascendent order)
condd = 5; %1 = old > newt1; 2 = newt1 > old; 3 = newt2 > old; 4 = newt3 > old; 5 = newt4 > old

close all
%for NTX against M only some specific tones were computed (e.g. tone 2 for NT1 vs M, tone 3 for NT2 vs M, etc.) 
if condd == 2
    tonel = 2;
elseif condd == 3
    tonel = 3;
elseif condd == 4
    tonel = 4;
elseif condd == 5
    tonel = 5;
end
list = dir('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/Papers/Systematic_Variations_DCM/NatureCommunications_Revision/Revision_I/DCM_swap/Out/DCM_tone2_EVK_S*_4_cond1_vs_cond2.mat');
if condd == 1
    conddd = '_cond1_vs_cond2';
elseif condd == 2
    conddd = '_cond2_vs_cond1';
    tonel = 2;
elseif condd == 3
    conddd = '_cond3_vs_cond1';
    tonel = 3;
elseif condd == 4
    conddd = '_cond4_vs_cond1';
    tonel = 4;
elseif condd == 5
    conddd = '_cond5_vs_cond1';
    tonel = 5;
end
F = zeros(length(ppp),length(model_n)); %to be updated here as well
sbarb = cell(length(ppp),length(model_n));
for ii = 1:length(ppp) %over subjects
    for jj = 1:length(model_n) %over models
        barbb = [list(1).folder '/' list(ppp(ii)).name(1:8) num2str(tonel) list(ppp(ii)).name(10:29) num2str(model_n(jj)) conddd '.mat'];
        load(barbb); %loading data
        sbarb{ii,jj} = barbb; %checking that you loaded the right data
        F(ii,jj) = DCM.F;
        disp(ii)
    end
end
addpath('/projects/MINDLAB2020_MEG-AuditoryPatternRecognition/scripts/leonardo/DCM')
BMS = bms_rfx(F); %computing statistics over participants
mean(F,1);
BMS
%plotting for the manuscript
% bar plotting
%posterior (expected) probability
[m,i] = max(BMS.pp);
W = zeros(1,length(BMS.pp));
W(i) = m;
figure
bar(1:6,BMS.pp,0.7,'b')
hold on
bar(W,'k')
ylim([0 1])
xlabel('Models')
ylabel('Probability')
set(gcf,'Color','w')
title('Posterior probability')
% export_fig(['/aux/MINDLAB2021_MEG-TempSeqAges/19_10_2023/PP_Tone' num2str(tonel) '_contrast' num2str(condd) '.eps'])
% export_fig(['/aux/MINDLAB2021_MEG-TempSeqAges/19_10_2023/PP_Tone' num2str(tonel) '_contrast' num2str(condd) '.png'])
%protected exceedance probability
[m,i] = max(BMS.pep);
W = zeros(1,length(BMS.pep));
W(i) = m;
figure
bar(1:6,BMS.pep,0.7,'b')
hold on
bar(W,'k')
ylim([0 1])
xlabel('Models')
ylabel('Probability')
set(gcf,'Color','w')
title('Protected exceedance probability')
% export_fig(['/aux/MINDLAB2021_MEG-TempSeqAges/19_10_2023/PEP_Tone' num2str(tonel) '_contrast' num2str(condd) '.eps'])
% export_fig(['/aux/MINDLAB2021_MEG-TempSeqAges/19_10_2023/PEP_Tone' num2str(tonel) '_contrast' num2str(condd) '.png'])

%% codes for additional brain figure (for plotting the models within a brain template) - Figures 6 and S16

for ppp = 1:7
    modeln = ppp; %model number (with reference to the models reported above for DCM
    
    name_model = ['Model_' num2str(modeln -1)];
    %loading MNI coordinates of AAL 2-mm centroids
    load('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/MNI_RC_AAL_2mm.mat')
    %coordinates for centroids of ROIs
    vcoord = [-42.3467 -20.0356 8.6756; 42.3467 -20.0356 8.6756; 0 -16.1010 40.2143; -25.2554 -21.9635 -11.3841; 25.2554 -21.9635 -11.3841; 0 52.5174 -8.8567];
    if modeln == 1
        error('not supported for model 1')
    elseif modeln == 2 || modeln == 3
        vcol = [0 0.4470 0.7410; 0 0.4470 0.7410; 0.8500 0.0250 0.1980; 0.8500 0.0250 0.1980; 0.8500 0.0250 0.1980; 0.8500 0.0250 0.1980]; %color codes for ROIs
    elseif modeln == 4
        vcol = [0 0.4470 0.7410; 0 0.4470 0.7410; 0 0 0; 0.8500 0.0250 0.1980; 0.8500 0.0250 0.1980; 0 0 0;]; %color codes for ROIs
    elseif modeln == 5
        vcol = [0 0.4470 0.7410; 0 0.4470 0.7410; 0.8500 0.0250 0.1980; 0 0 0; 0 0 0; 0.8500 0.0250 0.1980]; %color codes for ROIs
    elseif modeln == 6
        vcol = [0 0.4470 0.7410; 0 0.4470 0.7410; 0.8500 0.0250 0.1980; 0 0 0; 0 0 0; 0 0 0]; %color codes for ROIs
    elseif modeln == 7
        vcol = [0 0.4470 0.7410; 0 0.4470 0.7410; 0 0 0; 0.8500 0.0250 0.1980; 0.8500 0.0250 0.1980; 0.8500 0.0250 0.1980]; %color codes for ROIs
    end
    %loading a brain template
    openfig('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/BrainTemplate_GT.fig')
    hold on
    for ii = 1:length(vcoord) %over ROIs
        plot3(vcoord(ii, 1), vcoord(ii, 2), vcoord(ii, 3), ['.'], 'Color', vcol(ii,:), 'MarkerSize', 25); %centroid of ROI ii
        hold on
    end
    rotate3d on; axis off; axis vis3d; axis equal
    view([0 90])
    export_fig([name_model '_top.eps'])
    export_fig([name_model '_top.png'])
    view([-90 0])
    export_fig([name_model '_left.eps'])
    export_fig([name_model '_left.png'])
    view([90 0])
    export_fig([name_model '_right.eps'])
    export_fig([name_model '_right.png'])
    view([0 -90])
    export_fig([name_model '_top2.eps'])
    export_fig([name_model '_top2.png'])
    view([0 0])
    export_fig([name_model '_back.eps'])
    export_fig([name_model '_back.png'])
    view([-180 0])
    export_fig([name_model '_front.eps'])
    export_fig([name_model '_front.png'])
end

%%

%% *** TIME-FREQUENCY ANALYSIS ***

%% setting up the cluster

addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/Cluster_ParallelComputing') %add the path to the function that submits the jobs to the cluster
clusterconfig('scheduler', 'cluster');
clusterconfig('long_running', 1); %there are different queues for the cluster depending on the number and length of the jobs you want to submit 
clusterconfig('slot', 6); %slot in the queu

%% PROPER INDUCED RESPONSES - MEG CHANNELS (8 SELECTED)

basel = 1; % 1 = baseline correction (subtraction); 0 = no baseline correction

addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/MITDecodingToServer/BrainstormDemo_CodesTimeFrequencyDecomposition/functions') %path to Dimitrios function for Morlet transform
f = 1:1:60; %frequencies used
CHAN{1} = 'MEG211'; CHAN{2} = 'MEG1311'; CHAN{3} = 'MEG241'; CHAN{4} = 'MEG1331'; CHAN{5} = 'MEG1631'; CHAN{6} = 'MEG2441'; CHAN{7} = 'MEG192'; CHAN{8} = 'MEG2341';
chans_idx = [5 47 8 49 61 94 72 90];
list = dir(['//scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/espmeeg_*recogminor*.mat']);
P2 = zeros(length(CHAN),length(f),1026,5,length(list));
for ii = 1:length(list)
    D = spm_eeg_load([list(ii).folder '/' list(ii).name]);
    chans = D.chanlabels;
    idxMEGc1 = find(strcmp(chans,'MEG0111')); %getting extreme channels indexes (mag)
    chanlab = chans(idxMEGc1:3:idxMEGc1+305);
    IND{1} = find(contains(D.conditions,'Old_Correct')); %getting indices of old condition trials
    IND{2} = find(contains(D.conditions,'New_T1_Correct')); %same for NewT1
    IND{3} = find(contains(D.conditions,'New_T2_Correct')); %same
    IND{4} = find(contains(D.conditions,'New_T3_Correct'));
    IND{5} = find(contains(D.conditions,'New_T4_Correct'));
    for jj = 1:5 %over conditions
        data = D(idxMEGc1:3:idxMEGc1+305,1:1026,IND{jj}); %extracting magnetometers, proper time and proper trials
        data2 = data(chans_idx,:,:); %extracting proper channels
        P2dum = zeros(length(CHAN),length(f),1026,size(data2,3)); %allocating space for computation oversingle trials
        for tt = 1:size(data2,3) %over trials
            if basel ~= 1
                %no baseline correction (revision I)
                P2dum(:,:,:,tt) = morlet_transform(data2(:,:,tt),D.time(1:1026),f); %frequency decomposition over single trials
            else
                %baseline correction (revision II)
                dumbo = morlet_transform(data2(:,:,tt),D.time(1:1026),f); %frequency decomposition over single trials
                P2dum(:,:,:,tt) = dumbo - mean(dumbo(:,:,1:25),3);
            end
        end
        P2(:,:,:,jj,ii) = mean(P2dum,4); %averaging TF results over trials
        disp(['subj ' num2str(ii) ' - cond ' num2str(jj)])
    end
    disp(ii)
end
time = D.time;
if basel == 1
    %baseline correction (revision II)
    save('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/Papers/Systematic_Variations_DCM/NatureCommunications_Revision/Revision_I/Time_Frequency_Analysis/TF_ERFs/BaselineCorr_Time_Frequency_AllSubjects_AveragedTrials_MEGsensors.mat','P2','time','f');
else
    %no baseline correction (revision I)
    save('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/Papers/Systematic_Variations_DCM/NatureCommunications_Revision/Revision_I/Time_Frequency_Analysis/TF_ERFs/Time_Frequency_AllSubjects_AveragedTrials_MEGsensors.mat','P2','time','f');
end

%% INDUCED RESPONSES - MEG SOURCES - AAL

%(COMPUTING TIME-FREQUENCY ANALYSIS ON SIGNLE VOXELS AND THEN AVERAGING THE RESULTS)

%setting up the cluster
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/Cluster_ParallelComputing') %add the path to the function that submits the jobs to the cluster
clusterconfig('scheduler', 'cluster');
clusterconfig('long_running', 1); %there are different queues for the cluster depending on the number and length of the jobs you want to submit 
clusterconfig('slot', 1); %slot in the queu

%% AAL - here you provide either coordinates in MNI space or the AAL ROIs (following the AAL order)

ROIn = 8;
basel = 1; % 1 = baseline correction (subtraction); 0 = no baseline correction

vecttname = {'Occ_Sup_L';'Occ_Sup_R';'HeschlL';'HeschlR';'HippL';'HippR';'ACC';'MC'};
namee = vecttname{ROIn}; %if you use the AAL ROI(s), here you can specify the name to be used when saving the data
vectt = {[49];[50];[79];[80];[37];[38];[31 32];[33 34]}; %8 selected ROIs

%setting up the cluster (this one should be local.. otherwise it does not work for some unknown reasons.. 
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/Cluster_ParallelComputing') %add the path to the function that submits the jobs to the cluster
clusterconfig('scheduler', 'none'); % 'none' or 'cluster'
clusterconfig('long_running', 1); %there are different queues for the cluster depending on the number and length of the jobs you want to submit 
clusterconfig('slot', 1); %slot in the queu

if basel == 1
    outdir = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/Papers/Systematic_Variations_DCM/NatureCommunications_Revision/Revision_I/Time_Frequency_Analysis/LBPD_Parcels/FinalImages/BaselCorr';
else
    outdir = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/Papers/Systematic_Variations_DCM/NatureCommunications_Revision/Revision_I/Time_Frequency_Analysis/LBPD_Parcels/FinalImages';
end
mkdir(outdir)
S = [];
S.Aarhus_clust = 2; %0 = working locally; integer number (i.e. 1) = sending one job for each subject to the Aarhus cluster (the number you insert here corresponds to the slots of memory that you allocate for each job.
S.average_trials = 0; %1 = frequency decomposition after averaging the trials (evoked responses); 0 = frequency decomposition for single trials and then average of the time-frequency results (induced responses)
S.f = 1:1:60; %frequencies used
load('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/time_normal.mat');
S.time = time(1:1026);
if basel == 1
    S.baselcorr = [1 25];
end
S.conds = [1 2 3 4 5]; %experimental conditions
% S.coordd = [-38 18 -16; -22 58 -16; -30 34 -16]; %coordinates (empty for AAL ROIs)
S.coordd = []; %coordinates (empty for AAL ROIs)
% S.AAL_ROIs =[5,6,9,10,15,16,25,26]; %AAL ROIs indices (meaningful only if S.coordd is empty)
S.AAL_ROIs = [vectt{ROIn}]; %AAL ROIs indices (meaningful only if S.coordd is empty)
if S.average_trials == 1
    S.subjlist = dir(['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/Beam_abs_0_sens_1_freq_broadband/SUBJ*.mat']);
else
    S.subjlist = dir(['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/Beam_abs_0_sens_1_freq_broadband_invers_5/SUBJ*.mat']); 
end
if isempty(S.coordd)
    S.outdir = [outdir '/Time_Freq_AllSubjs_AAL_' namee]; %remember to specify the proper AAL ROI name that you want!!
else
    S.outdir = [outdir '/Time_Freq_AllSubjs_singlevoxels'];
end

jobid = job2cluster(@InducedResponses_Morlet_Coords_AALROIs_LBPD_D,S);

%% PLOTTING TF RESULTS (AND T-TESTS)

%% INDUCED RESPONSES - MEG SENSORS

% Figure S13

condition = 1; % 1 = Old; 2 = NewT1; 3 = NewT2; 4 = NewT3; 5 = NewT4
fff = 1:60; %frequencies to be plotted
exportl = 1; % 1 = export figures
ttestl = 1; %1 = t-tests Old vs NewT1; 0 = single condition
basel = 1; % 1 = baseline correction (subtraction); 0 = no baseline correction

CHAN{1} = 'MEG211'; CHAN{2} = 'MEG1311'; CHAN{3} = 'MEG241'; CHAN{4} = 'MEG1331'; CHAN{5} = 'MEG1631'; CHAN{6} = 'MEG2441'; CHAN{7} = 'MEG1921'; CHAN{8} = 'MEG2341'; %channels name
if ttestl ~= 1
%     loading data
%                 if ~exist('P2','var')
        if basel == 1
            load('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/Papers/Systematic_Variations_DCM/NatureCommunications_Revision/Revision_I/Time_Frequency_Analysis/TF_ERFs/BaselineCorr_Time_Frequency_AllSubjects_AveragedTrials_MEGsensors.mat');
        else
            load('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/Papers/Systematic_Variations_DCM/NatureCommunications_Revision/Revision_I/Time_Frequency_Analysis/TF_ERFs/Time_Frequency_AllSubjects_AveragedTrials_MEGsensors.mat');
        end
%                 end
    for ii = 1:8
        chani = ii; % 1 = MEG211; 2 = MEG1311; 3 = MEG241; 4 = MEG1331; 5 = MEG1631; 6 = MEG2441; 7 = MEG192; 8 = MEG2341
        %mean over subjects (for plotting purposes) and removing pre-stimulus time
        Pold = mean(P2,5); %mean over subjects
        %plotting
        figure
        imagesc(time,f(fff),squeeze(Pold(chani,fff,:,condition)))
        set(gca,'YDir','normal') %plotting frequencies in descending order
        xlabel('time (s)'); ylabel('f (Hz)');
        colorbar
        caxis([-3500 3500])
        set(gcf,'color','w')
        %colormap with white for 0 values
        x = []; x.bottom = [0 0 0.5]; x.botmiddle = [0 0.5 1]; x.middle = [1 1 1]; x.topmiddle = [1 0 0]; x.top = [0.6 0 0]; %red - blue
        colormap(bluewhitered_PD(0,x))
        title(['Chan ' CHAN{chani} ' - cond ' num2str(condition)])
        if exportl == 1
            export_fig(['/aux/MINDLAB2021_MEG-TempSeqAges/09_11_2023/SourceData/Figure_S13_Gemma/' CHAN{chani} '_Cond' num2str(condition) '.png'])
            export_fig(['/aux/MINDLAB2021_MEG-TempSeqAges/09_11_2023/SourceData/Figure_S13_Gemma/' CHAN{chani} '_Cond' num2str(condition) '.eps'])
%             export_fig(['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/Papers/Systematic_Variations_DCM/NatureCommunications_Revision/Revision_I/Time_Frequency_Analysis/LBPD_Parcels/FinalImages/MEGsensors_Induced/' CHAN{chani} '_Cond' num2str(condition) '.png'])
%             export_fig(['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/Papers/Systematic_Variations_DCM/NatureCommunications_Revision/Revision_I/Time_Frequency_Analysis/LBPD_Parcels/FinalImages/MEGsensors_Induced/' CHAN{chani} '_Cond' num2str(condition) '.eps'])
        end
    end
else
    %loading data
    if basel == 1
        load('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/Papers/Systematic_Variations_DCM/NatureCommunications_Revision/Revision_I/Time_Frequency_Analysis/TF_ERFs/BaselineCorr_Time_Frequency_AllSubjects_AveragedTrials_MEGsensors.mat');
    else
        load('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/Papers/Systematic_Variations_DCM/NatureCommunications_Revision/Revision_I/Time_Frequency_Analysis/TF_ERFs/Time_Frequency_AllSubjects_AveragedTrials_MEGsensors.mat');
    end
    Tallthresh = zeros(size(P2,2),size(P2,3),8);
    %computing t-tests
    for cc = 1:8 %over MEG channels
        P = zeros(size(P2,2),size(P2,3));
        T = zeros(size(P2,2),size(P2,3));
        for ii = 1:size(P2,2) %over frequencies
            for jj = 1:size(P2,3) %over time-points
                [~,p,~,stats] = ttest(squeeze(P2(cc,ii,jj,1,:)),squeeze(P2(cc,ii,jj,2,:))); %contrasting Old vs New
                P(ii,jj) = p;
                T(ii,jj) = stats.tstat;
            end
            disp(ii)
        end
        %testing contrast results (correcting for multiple comparison) by performing Monte Carlo simulations
        P3 = zeros(size(T,1),size(T,2));
        P3(abs(T)>2) = 1; %threshold
        thresh = 0;
        permut = 1000;
        threshMC = 0.001;
        perm_max = 1;
        t1 = f(fff); t2 = time;
        [ OUT ] = twoD_MCS_LBPD_D( P3, thresh, permut, threshMC, perm_max, t1 , t2 )
        PDn = cell2table(OUT(:,1:7)); %table (not getting the last column (8) because it contains a large matrix useful for plotting purposes)
        writetable(PDn,['/aux/MINDLAB2021_MEG-TempSeqAges/24_02_2024/MEGSensors/TF' CHAN{cc} '_cond1_vs_cond2.xlsx'],'Sheet',1); %printing excel file
        %plotting
        figure
        T2 = T;
        T2(OUT{1,8}==0) = 0;
        T2(:,1:25) = 0;
        imagesc(time,f(fff),squeeze(T2(fff,:)))
        set(gca,'YDir','normal') %plotting frequencies in descending order
        xlabel('time (s)'); ylabel('f (Hz)');
        colorbar
        caxis([-6 6])
        set(gcf,'color','w')
        %colormap with white for 0 values
        x = []; x.bottom = [0 0 0.5]; x.botmiddle = [0 0.5 1]; x.middle = [1 1 1]; x.topmiddle = [1 0 0]; x.top = [0.6 0 0]; %red - blue
        colormap(bluewhitered_PD(0,x))
        title([CHAN{cc} ' - Cond 1 vs Cond 2'])
        export_fig(['/aux/MINDLAB2021_MEG-TempSeqAges/24_02_2024/MEGSensors/' CHAN{cc} '_Cond1_vs_Cond2.png'])
        export_fig(['/aux/MINDLAB2021_MEG-TempSeqAges/24_02_2024/MEGSensors/' CHAN{cc} '_Cond1_vs_Cond2.eps'])
        Tallthresh(:,:,cc) = T2;
    end
end

%% plotting SELECTED AAL ROIs PLUS OCCIPITAL SUPERIOR L or R (ALSO FROM AAL)

% Figures 7 and S12 

condition = 0; % 1 = Old; 2 = NewT1; 3 = NewT2; 4 = NewT3; 5 = NewT4; 0 = Old vs NewTcc (all categories of New)
fff = 1:60; %frequencies to be plotted
listn = 8; %index of the analysis to be run (with reference to the variable "list")
loadl = 1; %1 for loading; 0 for not loading, in case you already loaded the data and want to plot different conditions
basel = 1; % 1 = baseline correction (subtraction); 0 = no baseline correction
exportl = 1; %1 for exporting figures
% jes = [10 1400];
% jes = [-1400 1400];
jes = [-600 600];

if loadl == 1
    %loading data
    if basel == 1
        list = dir(['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/Papers/Systematic_Variations_DCM/NatureCommunications_Revision/Revision_I/Time_Frequency_Analysis/LBPD_Parcels/FinalImages/BaselCorr/Time_*']);
    else
        list = dir(['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/Papers/Systematic_Variations_DCM/NatureCommunications_Revision/Revision_I/Time_Frequency_Analysis/LBPD_Parcels/FinalImages/Time_*']);
    end
    list = list([list(:).isdir]); %keeping only the folders (otherwise you would get .mat files as well..)
    list2 = dir([list(listn).folder '/' list(listn).name '/SUBJ*.mat']);
    load([list2(1).folder '/' list2(1).name]); %to get some information
    %         load([list(listn).folder '/' list(listn).name]);
    Pold = zeros(size(Psubj,2),size(Psubj,3),size(Psubj,4),length(list2));
    for ii = 1:length(list2)
        load([list2(ii).folder '/' list2(ii).name]);
        Pold(:,:,:,ii) = squeeze(mean(Psubj,1));
        disp(['loading subject ' num2str(ii)])
    end
end
%plotting
if condition ~= 0
    %mean over subjects (for plotting purposes) and removing pre-stimulus time
    %     Pold = squeeze(mean(mean(P2,5),1)); %mean over subjects and over voxels
    figure
    imagesc(time,f(fff),squeeze(Pold(fff,:,condition)))
    set(gca,'YDir','normal') %plotting frequencies in descending order
    xlabel('time (s)'); ylabel('f (Hz)');
    colorbar
    if ~isempty(jes)
        caxis(jes)
    end
    set(gcf,'color','w')
    %colormap with white for 0 values
    x = []; x.bottom = [0 0 0.5]; x.botmiddle = [0 0.5 1]; x.middle = [1 1 1]; x.topmiddle = [1 0 0]; x.top = [0.6 0 0]; %red - blue
    colormap(bluewhitered_PD(0,x))
    title([list(listn).name(end-10:end) ' - condition ' num2str(condition)])
    if exportl == 1
        export_fig(['/aux/MINDLAB2021_MEG-TempSeqAges/09_11_2023/SourceData/Figure_7_Gemma/TF_' list(listn).name(24:end) '_Cond' num2str(condition) '.png'])
        export_fig(['/aux/MINDLAB2021_MEG-TempSeqAges/09_11_2023/SourceData/Figure_7_Gemma/TF_' list(listn).name(24:end) '_Cond' num2str(condition) '.eps'])
    end
    %export data for NC data source
    if condition < 3 %only M and NT1
       bumbum = squeeze(Pold(fff,:,condition));
       PDn = table(bumbum); %table (not getting the last column (8) because it contains a large matrix useful for plotting purposes)
       writetable(PDn,['/aux/MINDLAB2021_MEG-TempSeqAges/09_11_2023/SourceData/Figure_7_Gemma/TF_' list(listn).name(24:end) '_Cond' num2str(condition) '.xlsx'],'Sheet',1); %printing excel file
    end
else
    Tall = zeros(size(Pold,1),size(Pold,2),4);
    Tallthresh = zeros(size(Pold,1),size(Pold,2),4);
    for cc = 1:4 %over contrasts (i.e. 'old' versus 'newtcc')
%         Pold = squeeze(mean(P2,1)); %mean over voxels
        %t-tests
        %computing t-tests
        P = zeros(size(Pold,1),size(Pold,2));
        T = zeros(size(Pold,1),size(Pold,2));
        for ii = 1:size(P,1) %over frequencies
            for jj = 1:size(P,2) %over time-points
                [~,p,~,stats] = ttest(squeeze(Pold(ii,jj,1,:)),squeeze(Pold(ii,jj,cc+1,:))); %contrasting cond1 (Old) with cond X (depending on vectcond it can be either NewT1 or NewT2 or NewT3 or NewT4
                P(ii,jj) = p;
                T(ii,jj) = stats.tstat;
            end
            disp(ii)
        end
        %testing contrast results (correcting for multiple comparison) by performing Monte Carlo simulations
        P3 = zeros(size(T,1),size(T,2));
%         P3(abs(T)>2) = 1; %threshold
        P3(P<0.05) = 1; %threshold
        thresh = 0;
        permut = 1000;
        threshMC = 0.001;
        perm_max = 1;
        t1 = f(fff); t2 = time;
        [ OUT ] = twoD_MCS_LBPD_D( P3, thresh, permut, threshMC, perm_max, t1 , t2 )
        PDn = cell2table(OUT(:,1:7)); %table (not getting the last column (8) because it contains a large matrix useful for plotting purposes)
        writetable(PDn,['/aux/MINDLAB2021_MEG-TempSeqAges/24_02_2024/AAL/TF_' list(listn).name(24:end) '_cond1_vs_cond' num2str(cc+1) '.xlsx'],'Sheet',1); %printing excel file
        %plotting
        figure
        T2 = T;
        T2(OUT{1,8}==0) = 0;
        T2(:,1:25) = 0;
        imagesc(time,f(fff),squeeze(T2(fff,:)))
        set(gca,'YDir','normal') %plotting frequencies in descending order
        xlabel('time (s)'); ylabel('f (Hz)');
        colorbar
        caxis([-6 6])
        set(gcf,'color','w')
        %colormap with white for 0 values
        x = []; x.bottom = [0 0 0.5]; x.botmiddle = [0 0.5 1]; x.middle = [1 1 1]; x.topmiddle = [1 0 0]; x.top = [0.6 0 0]; %red - blue
        colormap(bluewhitered_PD(0,x))
        Tallthresh(:,:,cc) = T2;
        title([list(listn).name(end-10:end) ' - Cond 1 vs Cond ' num2str(cc+1)])
        if exportl == 1
            export_fig(['/aux/MINDLAB2021_MEG-TempSeqAges/24_02_2024/AAL/Thresh_TF_' list(listn).name(24:end) '_cond1_vs_cond' num2str(cc+1) '.png'])
            export_fig(['/aux/MINDLAB2021_MEG-TempSeqAges/24_02_2024/AAL/Thresh_TF_' list(listn).name(24:end) '_cond1_vs_cond' num2str(cc+1) '.eps'])
        end
        %export data for NC data source
        if cc == 1 %only M vs NT1
            bumbum = squeeze(T2(fff,:));
            PDn = table(bumbum); %table (not getting the last column (8) because it contains a large matrix useful for plotting purposes)
            writetable(PDn,['/aux/MINDLAB2021_MEG-TempSeqAges/24_02_2024/AAL/TF_SourceData_' list(listn).name(24:end) '_cond1_vs_cond' num2str(cc+1) '.xlsx'],'Sheet',1); %printing excel file
        end
        %plotting
        figure
        T2 = T;
        imagesc(time,f(fff),squeeze(T2(fff,:)))
        set(gca,'YDir','normal') %plotting frequencies in descending order
        xlabel('time (s)'); ylabel('f (Hz)');
        colorbar
        caxis([-6 6])
        set(gcf,'color','w')
        %colormap with white for 0 values
        x = []; x.bottom = [0 0 0.5]; x.botmiddle = [0 0.5 1]; x.middle = [1 1 1]; x.topmiddle = [1 0 0]; x.top = [0.6 0 0]; %red - blue
        colormap(bluewhitered_PD(0,x))
        title([list(listn).name(end-10:end) ' - Cond 1 vs Cond ' num2str(cc+1)])
        Tall(:,:,cc) = T2;
        if exportl == 1
            export_fig(['/aux/MINDLAB2021_MEG-TempSeqAges/24_02_2024/AAL/TF_' list(listn).name(24:end) '_cond1_vs_cond' num2str(cc+1) '.png'])
            export_fig(['/aux/MINDLAB2021_MEG-TempSeqAges/24_02_2024/AAL/TF_' list(listn).name(24:end) '_cond1_vs_cond' num2str(cc+1) '.eps'])
        end
        if listn > 6
            save(['/aux/MINDLAB2021_MEG-TempSeqAges/24_02_2024/AAL/' list(listn).name '.mat'],'Tallthresh','Tall','Pold')
        else
            save(['/aux/MINDLAB2021_MEG-TempSeqAges/24_02_2024/AAL/' list(listn).name '.mat'],'Tallthresh','Tall')
        end
    end
end

%% creating source data for Figure S12

list = dir(['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/Papers/Systematic_Variations_DCM/NatureCommunications_Revision/Revision_I/Time_Frequency_Analysis/LBPD_Parcels/FinalImages/BaselCorr/Time_*']);
TALL = zeros(60,1026,4,length(list));
TALLTHRESH = zeros(60,1026,4,length(list));
for ii = 1:length(list) %over ROIs
    load(['/aux/MINDLAB2021_MEG-TempSeqAges/24_02_2024/AAL/' list(ii).name])
    TALL(:,:,:,ii) = Tall;
    TALLTHRESH(:,:,:,ii) = Tallthresh;
    disp(ii)
end

%%

%% TIME-FREQUENCY ANALYSIS - FUNCTIONAL PARCELLATION

%% setting up the cluster
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/Cluster_ParallelComputing') %add the path to the function that submits the jobs to the cluster
clusterconfig('scheduler', 'cluster');
clusterconfig('long_running', 1); %there are different queues for the cluster depending on the number and length of the jobs you want to submit 
clusterconfig('slot', 1); %slot in the queu

%% here you provide a matrix in LBPD coordinates with 0s and 1s to index specific ROIs 

basel = 1; % 1 = baseline correction (subtraction); 0 = no baseline correction

ROII{1} = 'MC'; ROII{2} = 'HITR'; ROII{3} = 'HITL'; ROII{4} = 'VMPFC'; ROII{5} = 'ACL'; ROII{6} = 'ACR';
%loading matrix in LBPD coordinates with 0s and 1s to index specific ROIs 
load('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/Beam_abs_1_sens_1_freq_broadband/DCM/ROIs_MainActivity/FinalROIs/ROIsDCM.mat');
for ii = 1:6 %over ROIs
    S = [];
    S.ROI_n = ii; %1 = MC; 2 = HITLR; 3 = VMPFC; 4 = ACLR; 5 = ACL; 6 = ACR
    S.f = 1:1:60; %frequencies used
    S.mask = ROIs_DCM;
    list = dir(['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/Beam_abs_0_sens_1_freq_broadband_invers_5/SUBJ*.mat']);
%         S.subjlist = list(82:end);
    S.subjlist = list;
    load('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/time_normal.mat');
    time(1027:end) = [];
    S.time = time;
    S.Aarhus_clust = 1; %0 = working locally; integer number (i.e. 1) = sending one job for each subject to the Aarhus cluster (the number you insert here corresponds to the slots of memory that you allocate for each job.
    S.single_subj = 1; %1 for saving single subject data; 0 for not saving it
    if basel == 1
        S.baselcorr = [1 25];
        S.outdir = ['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/Papers/Systematic_Variations_DCM/NatureCommunications_Revision/Revision_I/Time_Frequency_Analysis/LBPD_Parcels/BaselCorr/' ROII{ii} '_SingleSubjects_SingleVoxels.mat'];
        mkdir(S.outdir)
    else
        S.outdir = ['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/Papers/Systematic_Variations_DCM/NatureCommunications_Revision/Revision_I/Time_Frequency_Analysis/LBPD_Parcels/' ROII{ii} '_SingleSubjects_SingleVoxels.mat'];
    end
    %actual function
    jobid = job2cluster(@InducedResponses_Morlet_ROIs_LBPD_D,S);
end

%% STATISTICS AND PLOTTING - TF - ORIGINAL PARCELLATION

%LOADING INDEPENDENT SUBJECTS (COMPUTED INDEPENDENT VOXELS AND THEN AVERAGED IN THE ROIs)

% Figures S17 and S18

condition = 0; % 1 = Old; 2 = NewT1; 3 = NewT2; 4 = NewT3; 5 = NewT5; 0 = Old - NewT1
contr = [1 2 3 4 5]; %conditions to be contrasted (NOW IT'S FROM 1 TO 5 BECAUSE I CREATED A LOOP FOR DOING ALL CONTRASTS), meaningul only if condition == 0
fff = 1:60; %frequencies to be plotted
ROIn = [1:6]; % 1 = MC; 2 = HITL; 3 = HITR; 4 = VMPFC; 5 = ACL; 6 = ACR with Morlet computed independently on each voxel and each trial and then the output was averaged (twice)
loadl = 1; %1 for loading; 0 for not loading, in case you already loaded the data and want to plot different conditions
basel = 1; % 1 = baseline correction (subtraction); 0 = no baseline correction
exportl = 1; %1 for exporting figures

if basel == 1
    outdir = '/aux/MINDLAB2021_MEG-TempSeqAges/24_02_2024/FunctionalROI';
    indir = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/Papers/Systematic_Variations_DCM/NatureCommunications_Revision/Revision_I/Time_Frequency_Analysis/LBPD_Parcels/BaselCorr';
else
    outdir = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/Papers/Systematic_Variations_DCM/NatureCommunications_Revision/Revision_I/Time_Frequency_Analysis/LBPD_Parcels/FinalImages';
    indir = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/Papers/Systematic_Variations_DCM/NatureCommunications_Revision/Revision_I/Time_Frequency_Analysis/LBPD_Parcels';
end
ROII{1} = 'MC'; ROII{2} = 'HITL'; ROII{3} = 'HITR'; ROII{4} = 'VMPFC'; ROII{5} = 'ACL'; ROII{6} = 'ACR';
BARB = zeros(60,1026,6);
for ss = 1:length(ROIn)
    ROI = ROIn(ss);
    if loadl == 1
        if ROI == 1
            list = dir([indir '/MC*.mat']);
        elseif ROI == 2
            list = dir([indir '/HITL*.mat']);
        elseif ROI == 3
            list = dir([indir '/HITR*.mat']);
        elseif ROI == 4
            list = dir([indir '/VMPFC*.mat']);
        elseif ROI == 5
            list = dir([indir '/ACL*.mat']);
        elseif ROI == 6
            list = dir([indir '/ACR*.mat']);
        end
        load([list(1).folder '/' list(1).name]);
        P2 = zeros(size(Psubj,1),size(Psubj,2),size(Psubj,3),length(list));
        for ii = 1:length(list) %over subjects
            load([list(ii).folder '/' list(ii).name]);
            P2(:,:,:,ii) = Psubj;
            disp(ii)
        end
    end
    if condition ~= 0
        %mean over subjects (for plotting purposes) and removing pre-stimulus time
        figure
        Pold = squeeze(mean(P2,4)); %mean over subjects
        imagesc(time,f(fff),squeeze(Pold(fff,:,condition)))
        set(gca,'YDir','normal') %plotting frequencies in descending order
        xlabel('time (s)'); ylabel('f (Hz)');
        colorbar
        set(gcf,'color','w')
        caxis([-600 600])
        %colormap with white for 0 values
        x = []; x.bottom = [0 0 0.5]; x.botmiddle = [0 0.5 1]; x.middle = [1 1 1]; x.topmiddle = [1 0 0]; x.top = [0.6 0 0]; %red - blue
        colormap(bluewhitered_PD(0,x))
        title([ROII{ROI} ' - cond ' num2str(condition)])
        if exportl == 1
            export_fig([outdir '/ROIs_Induced_' ROII{ROI} '_Cond' num2str(condition) '.png'])
            export_fig([outdir '/ROIs_Induced_' ROII{ROI} '_Cond' num2str(condition) '.eps'])
        end
    end
    if condition == 0
        Tall = zeros(60,1026,4,6); %this was temporary just to save source data files for Nature Communications
        Tallthresh = zeros(60,1026,4,6); %this was temporary just to save source data files for Nature Communications
        for cc = 1:4 %over contrasts (old versus NewTX)
            Pold = P2;
            %     Pold = Pold - mean(Pold(:,1:26,:,:),2); %baseline correction
            %t-tests
            %computing t-tests
            P = zeros(size(Pold,1),size(Pold,2));
            T = zeros(size(Pold,1),size(Pold,2));
            for ii = 1:size(P,1) %over frequencies
                for jj = 1:size(P,2) %over time-points
                    [~,p,~,stats] = ttest(squeeze(Pold(ii,jj,contr(1),:)),squeeze(Pold(ii,jj,contr(cc+1),:))); %contrasting Old vs New
                    P(ii,jj) = p;
                    T(ii,jj) = stats.tstat;
                end
                disp(ii)
            end
            Tall(:,:,cc,ss) = T; %this was temporary just to save source data for Nature Communications
%         end
            %testing contrast results (correcting for multiple comparison) by performing Monte Carlo simulations
            P3 = zeros(size(T,1),size(T,2));
            P3(abs(T)>2) = 1; %threshold t-val = 2.6 corresponding to p-val < 0.01 (obtained by dividing 0.05 by the 4 comparisons employed here)
            thresh = 0;
            permut = 1000;
            threshMC = 0.001;
            perm_max = 1;
            t1 = f(fff); t2 = time;
            [ OUT ] = twoD_MCS_LBPD_D( P3, thresh, permut, threshMC, perm_max, t1 , t2 )
            PDn = cell2table(OUT(:,1:7)); %table
            writetable(PDn,[outdir '/TF_cond' num2str(contr(1)) '_vs_cond' num2str(contr(cc+1)) '_' ROII{ROIn(ss)} '.xlsx'],'Sheet',1); %printing excel file
            %plotting
            figure
            T2 = T;
            T2(OUT{1,8}==0) = 0;
            T2(:,1:24) = 0;
            imagesc(time,f(fff),squeeze(T2(fff,:)))
            set(gca,'YDir','normal') %plotting frequencies in descending order
            xlabel('time (s)'); ylabel('f (Hz)');
            colorbar
            caxis([-6 6])
            set(gcf,'color','w')
            %colormap with white for 0 values
            x = []; x.bottom = [0 0 0.5]; x.botmiddle = [0 0.5 1]; x.middle = [1 1 1]; x.topmiddle = [1 0 0]; x.top = [0.6 0 0]; %red - blue
            colormap(bluewhitered_PD(0,x))
            title([ROII{ROIn(ss)} ' - Cond' num2str(contr(1)) ' vs cond' num2str(contr(cc+1))])
            Tallthresh(:,:,cc,ss) = T;
            if exportl == 1
                export_fig([outdir '/TF_cond' num2str(contr(1)) '_vs_cond' num2str(contr(cc+1)) '_' ROII{ROIn(ss)} '.eps'])
                export_fig([outdir '/TF_cond' num2str(contr(1)) '_vs_cond' num2str(contr(cc+1)) '_' ROII{ROIn(ss)} '.png'])
            end
            %thresholded
            figure
            T2 = T;
            T2(OUT{1,8}==0) = 0;
            T2(:,1:24) = 0;
            imagesc(time,f(fff),squeeze(T2(fff,:)))
            set(gca,'YDir','normal') %plotting frequencies in descending order
            xlabel('time (s)'); ylabel('f (Hz)');
            colorbar
            caxis([-6 6])
            set(gcf,'color','w')
            %colormap with white for 0 values
            x = []; x.bottom = [0 0 0.5]; x.botmiddle = [0 0.5 1]; x.middle = [1 1 1]; x.topmiddle = [1 0 0]; x.top = [0.6 0 0]; %red - blue
            colormap(bluewhitered_PD(0,x))
            title([ROII{ROIn(ss)} ' - Cond' num2str(contr(1)) ' vs cond' num2str(contr(cc+1))])
            Tallthresh(:,:,cc,ss) = T;
            if exportl == 1
                export_fig([outdir '/TF_SignOnly_cond' num2str(contr(1)) '_vs_cond' num2str(contr(cc+1)) '_' ROII{ROIn(ss)} '.eps'])
                export_fig([outdir '/TF_SignOnly_cond' num2str(contr(1)) '_vs_cond' num2str(contr(cc+1)) '_' ROII{ROIn(ss)} '.png'])
            end
        end            
    end
end

%%

%% *** PLOTTING (AND DOING STATISTICS) ON THE TIME SERIES OF THE FUNCTIONALLY DERIVED ROIs ***

%% PREPARING DATA FROM EACH SUBJECT USING THE ABOVE SELECTED ROIs (THIS SECTION WOULD ALSO PLOT DATA FROM OTHER DATASETS COLLECTED AT THE SAME TIME AS THE ONE USED IN THIS PAPER)

%task = 6 is the task reported in this paper
task = 6; %1 = elderly; 2 = fast/slow; 3 = learningbach; 4 = numbers (auditory); 5 = encoding block 5; 6 = Block 3 APR2020; 7 = MMN

%loading ROIs coordinates
load('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/Beam_abs_1_sens_1_freq_broadband/DCM/ROIs_MainActivity/FinalROIs/ROIsDCM.mat')
if task == 1
    %list of source reconstructed data (for each subject)
    list = dir('/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/Source_LBPD/Block_3/Beam_abs_0_sens_1_freq_broadband_invers_1/SUBJ*.mat');
    %getting sign of the voxels based on the aggregated source reconstruction (over participants) - AVERAGED TRIALS
    load('/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/Source_LBPD/Block_3/Beam_abs_0_sens_1_freq_broadband_invers_1/sources_main_effects.mat');
elseif task == 2
    %list of source reconstructed data (for each subject)
    list = dir('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_4/Beam_abs_0_sens_1_freq_broadband_invers_1/SUBJ*.mat');
    %getting sign of the voxels based on the aggregated source reconstruction (over participants) - AVERAGED TRIALS
    load('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_4/Beam_abs_0_sens_1_freq_broadband_invers_1/sources_main_effects.mat');
elseif task == 3
    %list of source reconstructed data (for each subject)
    list = dir('/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/source_LBPD/LearningBach/Beam_abs_0_sens_1_freq_broadband_time_1_466/SUBJ*.mat');
    %getting sign of the voxels based on the aggregated source reconstruction (over participants) - AVERAGED TRIALS
    load('/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/source_LBPD/LearningBach/Beam_abs_0_sens_1_freq_broadband_time_1_466/sources_main_effects.mat');
elseif task == 4
    %list of source reconstructed data (for each subject)
    list = dir('/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/Source_LBPD/Block_4/Beam_abs_0_sens_1_freq_broadband_invers_1/SUBJ*.mat');
    %getting sign of the voxels based on the aggregated source reconstruction (over participants) - AVERAGED TRIALS
    load('/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/Source_LBPD/Block_4/Beam_abs_0_sens_1_freq_broadband_invers_1/sources_main_effects.mat');
elseif task == 5
    %list of source reconstructed data (for each subject)
    list = dir('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_5/Beam_abs_0_sens_1_freq_broadband_invers_1/SUBJ*.mat');
    %getting sign of the voxels based on the aggregated source reconstruction (over participants) - AVERAGED TRIALS
    load('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_5/Beam_abs_0_sens_1_freq_broadband_invers_1/sources_main_effects.mat');
elseif task == 6
    %list of source reconstructed data (for each subject)
    list = dir('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/Beam_abs_0_sens_1_freq_broadband/SUBJ*.mat');
    %getting sign of the voxels based on the aggregated source reconstruction (over participants) - AVERAGED TRIALS
    load('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/Beam_abs_0_sens_1_freq_broadband/sources_main_effects.mat');
elseif task == 7
    list = dir('/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/Source_LBPD/Block_2/Beam_abs_0_sens_1_freq_broadband_invers_1/SUBJ*.mat');
    load('/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/Source_LBPD/Block_2/Beam_abs_0_sens_1_freq_broadband_invers_1/MMNSubtracted_Average.mat');
end
if task == 3
    timex = 32:37;
elseif task == 5 || task == 4
    timex = 49:55;
elseif task == 7
    timex = 150:160;
else
    timex = 45:52;
end
vect = zeros(3559,1);
for jj = 1:3559 %over brain voxels
    if task == 7
        if squeeze(mean(mean(data(jj,timex,1,:),4),2)) > 0 %if the data in voxel jj is positive during MMN time
            vect(jj,1) = -1; %storing a vector with 1 and -1 to be used for later statistics
        else
            vect(jj,1) = 1;
        end
    else
        if squeeze(mean(t_val_s(jj,timex,1),2)) > 0 %if the data in voxel jj is positive during N100 time
            vect(jj,1) = -1; %storing a vector with 1 and -1 to be used for later statistics
        else
            vect(jj,1) = 1;
        end
    end
end
ROIII = {1,2,3,4,5,6}; %selected ROIs; 1 = MC, 2 = HITR, 3 = HITL, 4 = VMPFC, 5 = ACL, 6 = ACR
load([list(1).folder '/' list(1).name])
if task == 7
    dum2 = zeros(length(ROIII),size(OUT.sources_ERFs,2),2,length(list));
else
    dum2 = zeros(length(ROIII),size(OUT.sources_ERFs,2),size(OUT.sources_ERFs,3),length(list));
end
for ii = 1:length(list) %over subjects
    disp(['loading source reconstructed data for subject ' num2str(ii) ' / ' num2str(length(list))])
    load([list(ii).folder '/' list(ii).name])
    if task == 7
        dum = zeros(size(OUT.sources_ERFs,1),size(OUT.sources_ERFs,2),2);
    else
        dum = zeros(size(OUT.sources_ERFs,1),size(OUT.sources_ERFs,2),size(OUT.sources_ERFs,3));
    end
    if task == 7
        dum3 = (OUT.sources_ERFs(:,:,3) - OUT.sources_ERFs(:,:,1));
        dum4 = (OUT.sources_ERFs(:,:,4) - OUT.sources_ERFs(:,:,2));
    end
    for cc = 1:size(dum,3) %over conditions
        for jj = 1:size(OUT.sources_ERFs,1) %over brain voxels
            if task == 7
                if cc == 1 %here I subtract standard from deviant to obtain MMN (either global)
                    dum(jj,:,cc) = dum3(jj,:) .* vect(jj,1); %reversing (or not)..
                else %or local
                    dum(jj,:,cc) = dum4(jj,:) .* vect(jj,1); %reversing (or not)..
                end
            else
                dum(jj,:,cc) = OUT.sources_ERFs(jj,:,cc) .* vect(jj,1); %reversing (or not)..
            end
            disp(['block ' num2str(task) ' - subject - ' num2str(ii) ' - condition ' num2str(cc) ' - source ' num2str(jj)])
        end
    end
    for cc = 1:size(dum,3) %over conditions
        for pp = 1:length(ROIII) %over ROIs (with possibility of averaging together voxels of different ROIs (e.g. left and right auditory cortex (AC))
            dum2(pp,:,cc,ii) = mean(dum(sum(ROIs_DCM(:,ROIII{pp}),2)~=0,:,cc),1); %getting the indices of voxels of both ROIs (if you want two ROIs)
        end
    end
end
if task == 7
    condds = cell(1,2); condds{1} = 'GlobalMMN'; condds{2} = 'LocalMMN';
else
    condds = OUT.S.inversion.conditions;
end
save([list(1).folder '/ROIs_6.mat'],'dum2','condds')

%% PLOTTING MAIN ROIs AND CONDITIONS TOGETHER (ALL SUBJECTS AT THE SAME TIME)

conditions = [1 2]; %vector with conditions; tasks 1 and 4: 1 = Old; 2 = NewT1; 3 = NewT2; 4 = NewT3; 5 = NewT4
%                      task 5: 1 = encoding; 2 = Old; 3 = NewT1); task 3: 1 = Old; 2 = New
%                      task 2: 1 = OldF; 2 = OldS; 3 = NewF; 4 = NewS
ROII = {6}; %selected ROIs; 1 = MC, 2 = HITR, 3 = HITL, 4 = VMPFC, 5 = ACL, 6 = ACR; ([2 3] averages together ROII 2 and 3; [2:4] averages together ROIs 2,3,4)
task = 5; %1 = block 3 TSA2021; 2 = fast/slow; 3 = learningbach; 4 = numbers (auditory); 5 = encoding block 5; 6 = block 3 APR2020
limmy = []; %limit for y-axis (amplitude of the signal)
% limmy = [];
leg_l = 1; %1 for legend; 0 otherwise
figexp = 0; %1 = export figures; 0 = not export figures

% close all
%defining colors
color_line = colormap(lines(5)); %extracting some colours from a colormap
color_line2 = color_line;
color_line2(1,:) = color_line(2,:);
color_line2(2,:) = color_line(1,:);
color_line2(5,:) = [0.4 0.4 0.4];
if task == 1
    load('/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/Source_LBPD/Block_3/Beam_abs_0_sens_1_freq_broadband_invers_1/ROIs_6.mat');
    conds{1} = ' Mem '; conds{2} = ' NewT1 '; conds{3} = ' NewT3 ';
elseif task == 2
    load('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_4/Beam_abs_0_sens_1_freq_broadband_invers_1/ROIs_6.mat');
    conds{1} = ' MemF '; conds{2} = ' MemS '; conds{3} = ' NewF '; conds{4} = ' NewS ';
elseif task == 3
    load('/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/source_LBPD/LearningBach/Beam_abs_0_sens_1_freq_broadband_time_1_466/ROIs_6.mat');
    conds{1} = ' Mem '; conds{2} = ' New ';
elseif task == 4
    load('/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/Source_LBPD/Block_4/Beam_abs_0_sens_1_freq_broadband_invers_1/ROIs_6.mat');
    conds{1} = ' Mem '; conds{2} = ' NewT1 '; conds{3} = ' NewT2 '; conds{4} = ' NewT3 '; conds{5} = ' NewT4 ';
elseif task == 5
    load('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_5/Beam_abs_0_sens_1_freq_broadband_invers_1/ROIs_6.mat');
    conds{1} = ' Enc '; conds{2} = ' Mem '; conds{3} = ' NewT1 ';
elseif task == 6
    load('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/Beam_abs_0_sens_1_freq_broadband/ROIs_6.mat');
    conds{1} = ' Mem '; conds{2} = ' NewT1 '; conds{3} = ' NewT2 '; conds{4} = ' NewT3 '; conds{5} = ' NewT4 ';
end
if task == 3
    load('/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/source_LBPD/LearningBach/time.mat')
else
    load('/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/Source_LBPD/time_345.mat')
end
% COL{1} = 'b'; COL{2} = 'r'; COL{3} = 'k'; COL{4} = 'g'; COL{5} = 'c'; COL{6} = 'm'; %colors
if length(conditions) == 2 %line width
    LIN{1} = 2; LIN{2} = 1;
elseif length(conditions) == 3 %line width
    LIN{1} = 2.3; LIN{2} = 1.8; LIN{3} = 1;
elseif length(conditions) == 4 %line width
    LIN{1} = 2.8; LIN{2} = 2.3; LIN{3} = 1.8; LIN{4} = 1;
else
    LIN{1} = 3.2; LIN{2} = 2.8; LIN{3} = 2.3; LIN{4} = 1.8; LIN{5} = 1;
end
ROIN{1} = 'MC'; ROIN{2} = 'HITR'; ROIN{3} = 'HITL'; ROIN{4} = 'VMPFC'; ROIN{5} = 'ACL'; ROIN{6} = 'ACR';
% load('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/time_normal.mat')
figure
clear CC
for cc = 1:length(conditions) %over conditions
    for ii = 1:length(ROII) %over ROIs
        if length(ROII) == 1
            plot(time(1:size(dum2,2)),nanmean(dum2(ROII{ii},:,conditions(cc),:),4),'Color',color_line2(cc,:),'LineWidth',2,'DisplayName',[conds{conditions(cc)}])
            hold on
            plot(time(1:size(dum2,2)),nanmean(dum2(ROII{ii},:,conditions(cc),:),4) + (nanstd(dum2(ROII{ii},:,conditions(cc),:),0,4)./sqrt(size(dum2,4))),':','Color',color_line2(cc,:),'LineWidth',0.5,'DisplayName',[conds{conditions(cc)}])
            hold on
            plot(time(1:size(dum2,2)),nanmean(dum2(ROII{ii},:,conditions(cc),:),4) - (nanstd(dum2(ROII{ii},:,conditions(cc),:),0,4)./sqrt(size(dum2,4))),':','Color',color_line2(cc,:),'LineWidth',0.5,'DisplayName',[conds{conditions(cc)}])
        else
            if length(ROII{ii}) == 1
                plot(time(1:size(dum2,2)),nanmean(dum2(ROII{ii},:,conditions(cc),:),4),'Color',color_line(ii,:),'LineWidth',LIN{cc})
            else
                plot(time(1:size(dum2,2)),nanmean(nanmean(dum2(ROII{ii},:,conditions(cc),:),4),1),'Color',color_line(ii,:),'LineWidth',LIN{cc})
            end
        end
        hold on
    end
end
grid minor
if length(ROII) == 1
    if leg_l == 1
        legend('show')
    end
    title(ROIN{ROII{1}})
else %elaborated way (and probably barbaric) to create legends..
    dumR = cell(1,length(ROII)); %trick to get the legend..
    for ii = 1:length(ROII)
        if length(ROII{ii}) == 1
            dumR(ii) = ROIN(ROII{ii});
        else
            dumm = [];
            for pp = 1:length(ROII{ii})
                dumm = strcat(dumm,ROIN{ROII{ii}(pp)});
            end
            dumR(ii) = {dumm};
        end
    end
    if leg_l == 1
        legend(dumR)
    end
    s = [];
    for ii = 1:length(conditions) %trick to get the title..
        s = strcat(s,conds{conditions(ii)});
    end
    title(s)
end
% xlim([-0.1 time(size(dum2,2))])
xlim([-0.1 3.4])
if ~isempty(limmy)
    ylim(limmy)
end
set(gcf,'color','w')
if figexp == 1
    dumm = [];
    for pp = 1:length(ROII)
        dumm = strcat(dumm,ROIN{ROII{pp}});
    end
    if ~isempty(limmy)
        export_fig(['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/Beam_abs_1_sens_1_freq_broadband/Sources_Decoding_MEG0211/MCS/FinalImagesToWorkBench/Figures_3_4_5/Cond_task' num2str(task) '_' conds{conditions} '_' dumm '_Scaled.png'])
        export_fig(['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/Beam_abs_1_sens_1_freq_broadband/Sources_Decoding_MEG0211/MCS/FinalImagesToWorkBench/Figures_3_4_5/Cond_task' num2str(task) '_' conds{conditions} '_' dumm '_Scaled.eps'])
    else
        export_fig(['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/Beam_abs_1_sens_1_freq_broadband/Sources_Decoding_MEG0211/MCS/FinalImagesToWorkBench/Figures_3_4_5/Cond_task' num2str(task) '_' conds{conditions} '_' dumm '.png'])
        export_fig(['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/Beam_abs_1_sens_1_freq_broadband/Sources_Decoding_MEG0211/MCS/FinalImagesToWorkBench/Figures_3_4_5/Cond_task' num2str(task) '_' conds{conditions} '_' dumm '.eps'])
    end
end

%% COMPUTING STATISTICS CONTRASTING OLD VERSUS THE FOUR CATEGORIES OF NEW, ONE AT A TIME (THIS IS DONE FOR EACH OF THE SIX ROIs)

% Figure S11

%this is the file that you have to load and use to run the statistics
load('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/Beam_abs_0_sens_1_freq_broadband/ROIs_6.mat');
clear ROIN
ROIN{1} = 'MC'; ROIN{2} = 'HITR'; ROIN{3} = 'HITL'; ROIN{4} = 'VMPFC'; ROIN{5} = 'ACL'; ROIN{6} = 'ACR'; % ROIs in dum2

%t-tests
P = zeros(size(dum2,1),size(dum2,2),(size(dum2,3)-1)); %ROIs x time-points x contrasts (every NewTX versus Old)
T = zeros(size(dum2,1),size(dum2,2),(size(dum2,3)-1)); %ROIs x time-points x contrasts (every NewTX versus Old)
for ii = 1:size(dum2,1) %over ROIs
    for jj = 1:size(dum2,2) %ove time-points
        for cc = 1:(size(dum2,3)-1) %over the 4 NewTX
            %codes t-test
            [~,p,~,stats] = ttest(squeeze(dum2(ii,jj,1,:)),squeeze(dum2(ii,jj,cc+1,:))); %contrasting cond1 (Old) with cond X (depending on vectcond it can be either NewT1 or NewT2 or NewT3 or NewT4
            P(ii,jj,cc) = p;
            T(ii,jj,cc) = stats.tstat;
        end
        disp([num2str(ii) ' - ' num2str(jj)])
    end
end
%MCS
p_thresh = 0.05; %threshold for binarising p-values vector
load('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/time_normal.mat'); %loading time
outdir = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/Beam_abs_1_sens_1_freq_broadband/Sources_Decoding_MEG0211/MCS/FinalImagesToWorkBench/Statistics';
for ii = 1:size(dum2,1) %over ROIs
    for cc = 1:(size(dum2,3)-1) %over the 4 NewTX
        Pbin = zeros(1,size(dum2,2));
        Pbin(P(ii,:,cc)<p_thresh) = 1;
        tvals = T(ii,:,cc);
        [ sign_clust ] = oneD_MCS_LBPD_D( Pbin, 0, 1, 1000, 0.001, time(1:size(dum2,2)), tvals ) %Monte Carlo simulation function to correct for multiple comparisons
        PDn = cell2table(sign_clust); %table
        writetable(PDn,[outdir '/' ROIN{ii} '_OldvsNewT' num2str(cc) '.xlsx'],'Sheet',1); %printing excel file
    end
end

%%

%% SOURCE LEAKAGE CORRECTION - ORIGINAL PARCELLATION

cond = 5; %condition that you want; 1 = Old; 2-5 = NewTX (X = 1-4)
ROIs = [2]; %selected ROIs; 1 = MC, 2 = HITR, 3 = HITL, 4 = VMPFC, 5 = ACL, 6 = ACR
limmy = []; %limit for y-axis (amplitude of the signal)

if ~exist('dum2','var')
    load('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/Beam_abs_0_sens_1_freq_broadband/ROIs_6.mat');
    load('/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/Source_LBPD/time_345.mat')
end
conds{1} = ' Mem '; conds{2} = ' NewT1 '; conds{3} = ' NewT2 '; conds{4} = ' NewT3 '; conds{5} = ' NewT4 ';
tso = zeros(size(dum2,1),size(dum2,2),size(dum2,4));
for ii = 1:size(dum2,4) %over subjects
    tso(:,:,ii) = remove_source_leakage_b(squeeze(dum2(:,:,cond,ii)),'symmetric'); %source leakage correction
    disp(ii)
end

% %implementation for all conditions to export to local computer and plotting there
% tso = zeros(size(dum2,1),size(dum2,2),size(dum2,3),size(dum2,4));
% for ii = 1:size(dum2,4) %over subjects
%     for cc = 1:size(dum2,3) %over conditions
%         tso(:,:,cc,ii) = remove_source_leakage_b(squeeze(dum2(:,:,cc,ii)),'symmetric'); %source leakage correction
%         disp(ii)
%     end
% end
% save sources_leakagecorrected.mat tso ROIN conds
close all
%defining colors
color_line = colormap(lines(5)); %extracting some colours from a colormap
color_line2 = color_line;
color_line2(1,:) = color_line(2,:);
color_line2(2,:) = color_line(1,:);
color_line2(5,:) = [0.4 0.4 0.4];
LIN{1} = 3.2; LIN{2} = 2.8; LIN{3} = 2.3; LIN{4} = 1.8; LIN{5} = 1;
ROIN{1} = 'MC'; ROIN{2} = 'HITR'; ROIN{3} = 'HITL'; ROIN{4} = 'VMPFC'; ROIN{5} = 'ACL'; ROIN{6} = 'ACR';
% load('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/time_normal.mat')
figure
clear CC
for ii = 1:length(ROIs) %over ROIs
    plot(time(1:size(tso,2)),nanmean(tso(ROIs(ii),:,:),3),'Color',color_line2(ii,:),'LineWidth',2,'DisplayName',[ROIN{ROIs(ii)}])
    hold on
    plot(time(1:size(tso,2)),nanmean(tso(ROIs(ii),:,:),3) + (nanstd(tso(ROIs(ii),:,:),0,3)./sqrt(size(tso,3))),':','Color',color_line2(ii,:),'LineWidth',0.5,'HandleVisibility','off')
    hold on
    plot(time(1:size(tso,2)),nanmean(tso(ROIs(ii),:,:),3) - (nanstd(tso(ROIs(ii),:,:),0,3)./sqrt(size(tso,3))),':','Color',color_line2(ii,:),'LineWidth',0.5,'HandleVisibility','off')
    hold on
end
grid minor

% xlim([-0.1 time(size(dum2,2))])
xlim([-0.1 3.4])
if ~isempty(limmy)
    ylim(limmy)
end
set(gcf,'color','w')
legend('show')

%% STATISTICS AFTER SOURCE LEAKAGE CORRECTION

% Figure S15

load('/aux/MINDLAB2021_MEG-TempSeqAges/26_09_2023/sources_leakagecorrected.mat');
ROIN{1} = 'MC'; ROIN{2} = 'HITR'; ROIN{3} = 'HITL'; ROIN{4} = 'VMPFC'; ROIN{5} = 'ACL'; ROIN{6} = 'ACR';
dataa = tso;
%t-tests
P = zeros(size(dataa,1),size(dataa,2),(size(dataa,3)-1)); %clusters x time-points x contrasts (every NewTX versus Old)
T = zeros(size(dataa,1),size(dataa,2),(size(dataa,3)-1)); %clusters x time-points x contrasts (every NewTX versus Old)
for ii = 1:size(dataa,1) %over AAL parcels
    for jj = 1:size(dataa,2) %ove time-points
        for cc = 1:(size(dataa,3)-1) %over the 4 NewTX
            %codes t-test
            [~,p,~,stats] = ttest(squeeze(dataa(ii,jj,1,:)),squeeze(dataa(ii,jj,cc+1,:))); %contrasting cond1 (Old) with cond X (depending on vectcond it can be either NewT1 or NewT2 or NewT3 or NewT4
            P(ii,jj,cc) = p;
            T(ii,jj,cc) = stats.tstat;
        end
        disp([num2str(ii) ' - ' num2str(jj)])
    end
end
%MCS
p_thresh = 0.05; %threshold for binarising p-values vector
load('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/time_normal.mat'); %loading time
outdir = '/aux/MINDLAB2021_MEG-TempSeqAges/26_09_2023/Stats';
mkdir(outdir);
SIGN = cell(size(dataa,1),(size(dataa,3)-1));
for ii = 1:size(dataa,1) %over AAL parcels
    for cc = 1:(size(dataa,3)-1) %over the 4 NewTX
        Pbin = zeros(1,size(dataa,2));
        Pbin(P(ii,:,cc)<p_thresh) = 1;
        tvals = T(ii,:,cc);
        [ sign_clust ] = oneD_MCS_LBPD_D( Pbin, 0, 1, 1000, 0.001, time(1:size(dataa,2)), tvals ) %Monte Carlo simulation function to correct for multiple comparisons
        PDn = cell2table(sign_clust); %table
        writetable(PDn,[outdir '/'  ROIN{ii} '_OldvsNewT' num2str(cc) '.xlsx'],'Sheet',1); %printing excel file
        SIGN{ii,cc} = sign_clust;
    end
end
save([outdir '/Sign.mat'],'SIGN')

%%

%% *** ADDITIONAL IMAGES FOR METHODS FIGURE ***

%% METHOD FIGURE - DEPICTION OF 8-MM PARCELLATION

% Figure 1e

load('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/MNI152_8mm_coord_dyi.mat')
name_model = 'Parc_8mm';

%preparing color specifications
vcol = rand(size(MNI8,1),size(MNI8,2));
%loading a brain template
openfig('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/BrainTemplate_GT.fig')
hold on
cnt = 0;
for aa = 1:90 %over AAL ROIs
    %getting LBPD coordinates for AAL ROIs aa
    S = [];
    S.input = 2; %1 = MNI coordinates; 2 = AAL ROIs; 3 = general image with non-zero values
    S.coordd = []; %coordinates in MNI space (x,y,z)
    S.AAL_ROIs = [aa]; %AAL ROIs numbers you want to use
    % S.image = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/Beam_abs_0_sens_1_freq_broadband/Contr_1_abs_0.nii.gz';
    S.image = '/scratch7/MINDLAB2017_MEG-LearningBach/DTI_Portis/Templates/parcel_80mm_3559ROIs/3000.nii.gz';
    %actual function
    idx_LBPD = From3DNifti_OrMNICoords_2_CoordMatrix_8mm_LBPD_D(S);
    for ii = 1:size(idx_LBPD,1) %over brain voxels
        cnt = cnt + 1;
        plot3(MNI8(idx_LBPD(ii), 1), MNI8(idx_LBPD(ii), 2), MNI8(idx_LBPD(ii), 3), ['.'], 'Color', vcol(cnt,:), 'MarkerSize', 10); %centroid of ROI ii
        hold on
    end
end
rotate3d on; axis off; axis vis3d; axis equal

% view([0 90])
% export_fig([name_model '_top.eps'])
% export_fig([name_model '_top.png'])
% view([-90 0])
% export_fig([name_model '_left.eps'])
% export_fig([name_model '_left.png'])
view([90 0])
export_fig([name_model '_right.eps'])
export_fig([name_model '_right.png'])
% view([0 -90])
% export_fig([name_model '_top2.eps'])
% export_fig([name_model '_top2.png'])
% view([0 0])
% export_fig([name_model '_back.eps'])
% export_fig([name_model '_back.png'])
% view([-180 0])
% export_fig([name_model '_front.eps'])
% export_fig([name_model '_front.png'])

%% METHOD FIGURE - DEPICTION OF AAL PARCELLATION

% Figure 1e

load('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/MNI152_8mm_coord_dyi.mat')
name_model = 'Parc_AAL';

%preparing color specificationst
vcol = rand(90,3);
%loading a brain template
openfig('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/BrainTemplate_GT.fig')
hold on
for aa = 1:90 %over AAL ROIs
    %getting LBPD coordinates for AAL ROIs aa
    S = [];
    S.input = 2; %1 = MNI coordinates; 2 = AAL ROIs; 3 = general image with non-zero values
    S.coordd = []; %coordinates in MNI space (x,y,z)
    S.AAL_ROIs = [aa]; %AAL ROIs numbers you want to use
    % S.image = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/Beam_abs_0_sens_1_freq_broadband/Contr_1_abs_0.nii.gz';
    S.image = '/scratch7/MINDLAB2017_MEG-LearningBach/DTI_Portis/Templates/parcel_80mm_3559ROIs/3000.nii.gz';
    %actual function
    idx_LBPD = From3DNifti_OrMNICoords_2_CoordMatrix_8mm_LBPD_D(S);
    for ii = 1:size(idx_LBPD,1) %over brain voxels
        plot3(MNI8(idx_LBPD(ii), 1), MNI8(idx_LBPD(ii), 2), MNI8(idx_LBPD(ii), 3), ['.'], 'Color', vcol(aa,:), 'MarkerSize', 10); %centroid of ROI ii
        hold on
    end
end
rotate3d on; axis off; axis vis3d; axis equal

% view([0 90])
% export_fig([name_model '_top.eps'])
% export_fig([name_model '_top.png'])
% view([-90 0])
% export_fig([name_model '_left.eps'])
% export_fig([name_model '_left.png'])
view([90 0])
export_fig([name_model '3_right.eps'])
export_fig([name_model '3_right.png'])
% view([0 -90])
% export_fig([name_model '_top2.eps'])
% export_fig([name_model '_top2.png'])
% view([0 0])
% export_fig([name_model '_back.eps'])
% export_fig([name_model '_back.png'])
% view([-180 0])
% export_fig([name_model '_front.eps'])
% export_fig([name_model '_front.png'])

%%

%%

%% *** PLOTTING ON LOCAL COMPUTER FOR BETTER QUALITY IMAGES ***

%%

%% LOCAL COMPUTER - REVISION I - NATURE COMMUNICATIONS

% This code was used locally to prepare plots with higher resolution compared to the solutions offered by the Aarhus cluster of computers.
% It is mostly the same as reported above, with minor modifications to adapt to the local computer.
% I decided to report it integrally for enhanced readability

%%

%% AAL (single ROIs) and relationship with original parcellation

addpath('/Users/au550322/Documents/AarhusUniversitet/Postdoc/DataCollection_AuditoryPatternRecognition2020/Papers/SystematicVariationSequence/Codes_Images_Local');

%loading data, labels, significant time-windows, etc.
load('/Users/au550322/Documents/AarhusUniversitet/Postdoc/DataCollection_AuditoryPatternRecognition2020/Papers/SystematicVariationSequence/NatureCommunications/FirstRevision/14_09_2023/AAL_Stats/Sign.mat');
load('/Users/au550322/Documents/AarhusUniversitet/Postdoc/DataCollection_AuditoryPatternRecognition2020/Papers/SystematicVariationSequence/NatureCommunications/FirstRevision/14_09_2023/AAL_parcels_nifti/ROIs_to_AAL.mat');
load('/Users/au550322/Documents/AarhusUniversitet/Postdoc/DataCollection_AuditoryPatternRecognition2020/Papers/SystematicVariationSequence/NatureCommunications/FirstRevision/14_09_2023/AAL_data.mat');
load('/Users/au550322/Documents/AarhusUniversitet/Postdoc/DataCollection_AuditoryPatternRecognition2020/Papers/SystematicVariationSequence/NatureCommunications/FirstRevision/14_09_2023/AAL_label_cerebellum.mat');
load('/Users/au550322/Documents/AarhusUniversitet/Postdoc/DataCollection_AuditoryPatternRecognition2020/Papers/SystematicVariationSequence/Codes_Images_Local/time.mat'); %loading time        

%% actual code for plotting

col_l = 1; %1 for significant time-windows with different colors; 1 = only grey color
ylimm = [-35 15]; %amplitude limits; leave empty [] for automatic adjustment

ROI = 1; %ROIs (they are 6)
condition = [1:5]; %1 = Old; 2 = NewT1; 3 = NewT2; 4 = NewT3; 5 = NewT4;
export_l = 1; %1 = export images; 0 = not
cnt = 0;
close all
for ii = 7%:size(ROIs_to_AAL,1) %over original parcels (+1 which is the voxels which did not belong to any parcel)
    lineplotdum = [20 1]; %number instructing where to place the lines showing significant time-windows; leave empty [] if you want to have shaded colors instead
    lineplot = repmat(lineplotdum,length(ROIs_to_AAL{ii,1}),1);
    for pp = 1:length(ROIs_to_AAL{ii,1}) %over AAL parcels within parcel ii from the original parcellation
        cnt = cnt + 1;
%         figure(ii)
        clear ROIN
        S = [];
        S.ii = pp;
        S.conds = {'Old ','NewT1','NewT2','NewT3','NewT4'};
        ROIN{1} = lab(ROIs_to_AAL{ii,1}(pp),:);
        % load('/Users/au550322/Documents/AarhusUniversitet/Postdoc/DataCollection_AuditoryPatternRecognition2020/Papers/SystematicVariationSequence/Codes_Images_Local/ROIs_6.mat'); %loading data
        %structure for the function
        data2 = dataa(ROIs_to_AAL{ii,1}(pp),:,:,:);
        S.data = permute(data2,[1 2 4 3]);
%         S.data = data2(:,1:1026,:,:);
        S.STE = 2; %1 = dot lines for standard error; 2 = shadows
        S.transp = 0.3; %transparency for standard errors shadow
        S.time_real = time_sel(1:1026);
        S.colorline = [1 0 0; 0.3686 0.6314 0.7412; 0.1882 0.4902 0.8118; 0.0784 0.1569 0.5255; 0 0 0];
        if export_l == 1
            S.legendl = 0;
        else
            S.legendl = 0;
        end
        S.x_lim = [-0.1 3.4]; % Set x limits
        S.y_lim = ylimm; %Set y limits
        S.ROI_n = ROI;
        S.condition_n = condition;
        S.ROIs_labels = ROIN(ROI);
        S.lineplot = lineplot(cnt,:);
        %         S.subplot = [];
        %rearranging the significant time-windows
        signtp_col = [];
        bum = [];
        for ss = 1:4 %over contrasts between conditions
            %             signtp_col = cat(1,signtp_col,ones(size(SIGN{ROIs_to_AAL{ii,1}(pp),ss},1),1)*(ss+1)); %getting the color code (number of the significant time-windows for each contrast)
            %             for ll = 1:sbam %over the number of the significant time-windows for contrast ss
            %                 bum = cat(2,bum,SIGN{ROIs_to_AAL{ii,1}(pp),ss}(ll,3)); %concatenating all significant time-windows, for one contrast at a time
            %             end
            sbam = size(SIGN{ROIs_to_AAL{ii,1}(pp),ss},1);
            for ll = 1:sbam %over the number of the significant time-windows for contrast ss
                if SIGN{ROIs_to_AAL{ii,1}(pp),ss}{ll,3}(1) > 0.350 && SIGN{ROIs_to_AAL{ii,1}(pp),ss}{ll,3}(1) < 2.5  %we do not want to plot very small spurious random differences between conditions; they are available though in the supplementary tables
                    bum = cat(2,bum,SIGN{ROIs_to_AAL{ii,1}(pp),ss}(ll,3)); %concatenating all significant time-windows, for one contrast at a time
                    signtp_col = cat(1,signtp_col,(ss+1)); %getting the color code (number of the significant time-windows for each contrast)
                end
            end
        end
        S.signtp = bum;
        if col_l == 1
            S.signtp_col = signtp_col;
        else
            S.signtp_col = [];
        end
        
        waveform_plotting_local_v2(S) %actual function
        title(lab(ROIs_to_AAL{ii,1}(pp),:))
        
        if export_l == 1
            exportgraphics(gcf,['/Users/au550322/Documents/AarhusUniversitet/Postdoc/DataCollection_AuditoryPatternRecognition2020/Papers/SystematicVariationSequence/NatureCommunications/FirstRevision/Images_ToGemma/Timeseries/AAL_ROIs/ROI_' ROIN{ROI} '.pdf'],'Resolution',300)
        end
    end
end

%% CORRELATIONS BETWEEN VOXELS WITHIN EACH PARCEL COMPARED TO VOXELS FROM OTHER PARCELS

data = mean(dataa,4); %mean over subjects
CORR = zeros(size(ROIs_to_AAL,1),size(ROIs_to_AAL,1),size(data,3));
for cc = 1:size(data,3) %over conditions
    for ii = 1:size(ROIs_to_AAL,1) %over original parcels (+1 which is the voxels which did not belong to any parcel)
        for qq = 1:size(ROIs_to_AAL,1) %over original parcels (+1 which is the voxels which did not belong to any parcel)
            dum = zeros(length(ROIs_to_AAL{ii,1}),length(ROIs_to_AAL{qq,1}));
            for pp = 1:length(ROIs_to_AAL{ii,1}) %over AAL parcels within parcel ii from the original parcellation
                for tt = 1:length(ROIs_to_AAL{qq,1}) %over AAL parcels within parcel ii from the original parcellation
                    if pp ~= tt %not computing correlations of voxels with themselves
                        dum(pp,tt) = corr(data(ROIs_to_AAL{ii,1}(pp),:,cc)',data(ROIs_to_AAL{qq,1}(tt),:,cc)'); %correlation
                    end
                end
            end
            dum2 = dum(dum~=0); %getting rid of the diagonal
            CORR(ii,qq,cc) = mean(mean(dum2)); %mean over the correlations for all voxels in ROI ii
        end
    end
%     figure
%     imagesc(CORR(:,:,cc))
%     caxis([-0.01 0.7])
%     colorbar
end
corm = mean(CORR,3);
figure
imagesc(triu(corm))
% caxis([-0.01 0.7])
colorbar
set(gcf,'color','w')

%% ALL MEG SENSORS

addpath('/Users/au550322/Documents/AarhusUniversitet/Postdoc/DataCollection_AuditoryPatternRecognition2020/Papers/SystematicVariationSequence/Codes_Images_Local');

%loading data, labels, significant time-windows, etc.
load('/Users/au550322/Documents/AarhusUniversitet/Postdoc/DataCollection_AuditoryPatternRecognition2020/Papers/SystematicVariationSequence/NatureCommunications/FirstRevision/14_09_2023/All_MEG_sensors_stats/Sign.mat');
load('/Users/au550322/Documents/AarhusUniversitet/Postdoc/DataCollection_AuditoryPatternRecognition2020/Papers/SystematicVariationSequence/NatureCommunications/FirstRevision/14_09_2023/data_MEGsensors_all.mat');
load('/Users/au550322/Documents/AarhusUniversitet/Postdoc/DataCollection_AuditoryPatternRecognition2020/Papers/SystematicVariationSequence/Codes_Images_Local/time.mat'); %loading time        

%% actual code for plotting

col_l = 1; %1 for significant time-windows with different colors; 1 = only grey color
lineplotdum = [20 1]; %number instructing where to place the lines showing significant time-windows; leave empty [] if you want to have shaded colors instead
MEGn = [2:3:306 3:3:306]; %MEG channels you want to plot
ylimm = [-4.5 7.2]; %amplitude limits; leave empty [] for automatic adjustmentc
export_l = 1; %1 = export images; 0 = not

ROI = 1; %ROIs (they are 6)
condition = [1:5]; %1 = Old; 2 = NewT1; 3 = NewT2; 4 = NewT3; 5 = NewT4;
for ii = 1:length(MEGn)
    close all
    lineplot = repmat(lineplotdum,length(MEGn),1);
    clear ROIN
    S = [];
    S.ii = ii;
    S.conds = {'Old ','NewT1','NewT2','NewT3','NewT4'};
    data2 = dataa(MEGn(ii),1:1026,:,:);
    S.data = permute(data2,[1 2 4 3]);
    %         S.data = data2(:,1:1026,:,:);
    S.STE = 2; %1 = dot lines for standard error; 2 = shadows
    S.transp = 0.3; %transparency for standard errors shadow
    S.time_real = time_sel(1:1026);
    S.colorline = [1 0 0; 0.3686 0.6314 0.7412; 0.1882 0.4902 0.8118; 0.0784 0.1569 0.5255; 0 0 0];
    if export_l == 1
        S.legendl = 0;
    else
        S.legendl = 0;
    end
    S.x_lim = [-0.1 3.4]; % Set x limits
    S.y_lim = ylimm; %Set y limits
    S.ROI_n = ROI;
    S.condition_n = condition;
    S.ROIs_labels = chanlabs(MEGn(ii));
    %rearranging the significant time-windows
    signtp_col = [];
    bum = [];    
    for ss = 1:4 %over contrasts between conditions
        %     signtp_col = cat(1,signtp_col,ones(size(SIGN{ROI,ss},1),1)*(ss+1)); %getting the color code (number of the significant time-windows for each contrast)
        sbam = size(SIGN{MEGn(ii),ss},1);
        for ll = 1:sbam %over the number of the significant time-windows for contrast ss
            if SIGN{MEGn(ii),ss}{ll,3}(1) > 0.350 && SIGN{MEGn(ii),ss}{ll,3}(1) < 2.5  %we do not want to plot very small spurious random differences between conditions; they are available though in the supplementary tables
                bum = cat(2,bum,SIGN{MEGn(ii),ss}(ll,3)); %concatenating all significant time-windows, for one contrast at a time
                signtp_col = cat(1,signtp_col,(ss+1)); %getting the color code (number of the significant time-windows for each contrast)
            end
        end
    end
    S.signtp = bum;
    if col_l == 1
        S.signtp_col = signtp_col;
    else
        S.signtp_col = [];
    end
    S.lineplot = lineplot(ii,:);
    waveform_plotting_local_v2(S) %actual function
    title(chanlabs{MEGn(ii)})

    if export_l == 1
        exportgraphics(gcf,['/Users/au550322/Documents/AarhusUniversitet/Postdoc/DataCollection_AuditoryPatternRecognition2020/Papers/SystematicVariationSequence/NatureCommunications/FirstRevision/Images_ToGemma/Timeseries/MEG_Channels/' chanlabs{MEGn(ii)} '.pdf'],'Resolution',300)
    end
end

%% SELECTED CLUSTERS OF MEG CHANNELS

col_l = 1; %1 for significant time-windows with different colors; 1 = only grey color
lineplot = [20 1; 20 1; 20 1; 20 1]; %number instructing where to place the lines showing significant time-windows; leave empty [] if you want to have shaded colors instead
ylimm = [-125 83;-125 83;-1.2 2.2;-1.2 2.2]; %amplitude limits; leave empty [] for automatic adjustment
ROI = 1; %ROIs (they are 6)
condition = [1:5]; %1 = Old; 2 = NewT1; 3 = NewT2; 4 = NewT3; 5 = NewT4;
export_l = 0; %1 = export images; 0 = not

addpath('/Users/au550322/Documents/AarhusUniversitet/Postdoc/DataCollection_AuditoryPatternRecognition2020/Papers/SystematicVariationSequence/Codes_Images_Local');
%loading time
load('/Users/au550322/Documents/AarhusUniversitet/Postdoc/DataCollection_AuditoryPatternRecognition2020/Papers/SystematicVariationSequence/Codes_Images_Local/time.mat'); %loading time        
close all
cnt = 0;
DUMM = zeros(5,1026);
cnt2 = 0;
for ii = 1:2 %over mag and grad
    if ii == 1
        load('/Users/au550322/Documents/AarhusUniversitet/Postdoc/DataCollection_AuditoryPatternRecognition2020/Papers/SystematicVariationSequence/NatureCommunications/FirstRevision/14_09_2023/data_MEGsensors_2clust_mag.mat');
        namee = 'mag';
        load('/Users/au550322/Documents/AarhusUniversitet/Postdoc/DataCollection_AuditoryPatternRecognition2020/Papers/SystematicVariationSequence/NatureCommunications/FirstRevision/14_09_2023/Sign_2clustMEGsensors_MAG.mat');
    else
        load('/Users/au550322/Documents/AarhusUniversitet/Postdoc/DataCollection_AuditoryPatternRecognition2020/Papers/SystematicVariationSequence/NatureCommunications/FirstRevision/14_09_2023/data_MEGsensors_2clust_grad.mat');
        namee = 'grad';
        load('/Users/au550322/Documents/AarhusUniversitet/Postdoc/DataCollection_AuditoryPatternRecognition2020/Papers/SystematicVariationSequence/NatureCommunications/FirstRevision/14_09_2023/Sign_2clustMEGsensors_GRAD.mat');
    end
    for pp = 1:2 %over cluster
        cnt = cnt + 1;
        clear ROIN
        S = [];
        S.ii = cnt;
        S.conds = {'Old ','NewT1','NewT2','NewT3','NewT4'};
        S.lineplot = lineplot(cnt,:);
        
        data2 = dataa(pp,1:1026,:,:);
        if ii == 1
            S.data = data2;
        else
            S.data = permute(data2,[1 2 4 3]);
        end
        %         S.data = data2(:,1:1026,:,:);
        S.STE = 2; %1 = dot lines for standard error; 2 = shadows
        S.transp = 0.3; %transparency for standard errors shadow
        S.time_real = time_sel(1:1026);
        S.colorline = [1 0 0; 0.3686 0.6314 0.7412; 0.1882 0.4902 0.8118; 0.0784 0.1569 0.5255; 0 0 0];
        if export_l == 1
            S.legendl = 0;
        else
            S.legendl = 0;
        end
        S.x_lim = [-0.1 3.4]; % Set x limits
        S.y_lim = ylimm(cnt,:); %Set y limits
        S.ROI_n = ROI;
        S.condition_n = condition;
        S.ROIs_labels = {[namee ' clust ' num2str(pp)]};
        %rearranging the significant time-windows
        signtp_col = [];
        bum = [];
        for ss = 1:4 %over contrasts between conditions
            signtp_col = cat(1,signtp_col,ones(size(SIGN{pp,ss},1),1)*(ss+1)); %getting the color code (number of the significant time-windows for each contrast)
            sbam = size(SIGN{pp,ss},1);
            for ll = 1:sbam %over the number of the significant time-windows for contrast ss
                bum = cat(2,bum,SIGN{pp,ss}(ll,3)); %concatenating all significant time-windows, for one contrast at a time
            end
        end
        S.signtp = bum;
        if col_l == 1
            S.signtp_col = signtp_col;
        else
            S.signtp_col = [];
        end
        waveform_plotting_local_v2(S) %actual function
        
        if export_l == 1
            exportgraphics(gcf,['/Users/au550322/Documents/AarhusUniversitet/Postdoc/DataCollection_AuditoryPatternRecognition2020/Papers/SystematicVariationSequence/NatureCommunications/FirstRevision/Images_ToGemma/Timeseries/MEG_Channels_2clust/' namee '_Clust_' num2str(pp) '.pdf'],'Resolution',300)
        end
        title([namee '_Clust_' num2str(pp)])
        
        %exporting data for NC source data
        DUMM(1,:) = S.time_real;
        for iii = 1:5 %over conditions
            cnt2 = cnt2 + 1;
            bum = S.data(:,:,:,iii);
            DUMM(cnt2 + 1,:) = squeeze(mean(bum,3));
        end
    end
end
outdir = '/Users/au550322/Documents/AarhusUniversitet/Postdoc/DataCollection_AuditoryPatternRecognition2020/Papers/SystematicVariationSequence/NatureCommunications/FirstRevision/SourceData';
PDn = table(DUMM); %table
writetable(PDn,[outdir '/MEGClustersSelected.xlsx'],'Sheet',1); %printing excel file

%% PREVIOUS MEG SOURCES - NEW FUNCTION

col_l = 1; %1 for significant time-windows with different colors; 1 = only grey color
lineplot = [20 1]; %number instructing where to place the lines showing significant time-windows; leave empty [] if you want to have shaded colors instead
ylimm = []; %amplitude limits; leave empty [] for automatic adjustment
source_leak_corr = 1; %1 = source leakage corrected time series; 0 = original time series

ROI = 6; %ROIs (they are 6)
condition = [1:5]; %1 = Old; 2 = NewT1; 3 = NewT2; 4 = NewT3; 5 = NewT4;
export_l = 0; %1 = export images; 0 = not

addpath('/Users/au550322/Documents/AarhusUniversitet/Postdoc/DataCollection_AuditoryPatternRecognition2020/Papers/SystematicVariationSequence/Codes_Images_Local');
close all
clear ROIN
if ROI > 4
    if source_leak_corr == 1
        ylimm = [-0.085 0.04]; %amplitude limits; leave empty [] for automatic adjustment
    else
        ylimm = [-60 20]; %amplitude limits; leave empty [] for automatic adjustment
    end
else
    if source_leak_corr == 1
        ylimm = [-0.055 0.04]; %amplitude limits; leave empty [] for automatic adjustment
    else
        ylimm = [-38 23]; %amplitude limits; leave empty [] for automatic adjustment
    end
end
S = [];
S.conds = {'Old ','NewT1','NewT2','NewT3','NewT4'};
ROIN{1} = 'MC'; ROIN{2} = 'HITR'; ROIN{3} = 'HITL'; ROIN{4} = 'VMPFC'; ROIN{5} = 'ACL'; ROIN{6} = 'ACR';

if source_leak_corr == 1
    load('/Users/au550322/Documents/AarhusUniversitet/Postdoc/DataCollection_AuditoryPatternRecognition2020/Papers/SystematicVariationSequence/NatureCommunications/FirstRevision/14_09_2023/sources_leakagecorrected.mat');
    data2 = tso;
    load('/Users/au550322/Documents/AarhusUniversitet/Postdoc/DataCollection_AuditoryPatternRecognition2020/Papers/SystematicVariationSequence/NatureCommunications/FirstRevision/14_09_2023/Stats_sourceleakage/Sign.mat');
else
    load('/Users/au550322/Documents/AarhusUniversitet/Postdoc/DataCollection_AuditoryPatternRecognition2020/Papers/SystematicVariationSequence/NatureCommunications/FirstRevision/14_09_2023/Sign_OriginalROIs.mat');
    load('/Users/au550322/Documents/AarhusUniversitet/Postdoc/DataCollection_AuditoryPatternRecognition2020/Papers/SystematicVariationSequence/Codes_Images_Local/ROIs_6.mat'); %loading data
    data2 = dum2;
end
load('/Users/au550322/Documents/AarhusUniversitet/Postdoc/DataCollection_AuditoryPatternRecognition2020/Papers/SystematicVariationSequence/Codes_Images_Local/time.mat'); %loading time
%structure for the function
data2 = permute(data2,[1 2 4 3]);
S.data = data2(:,1:1026,:,:);
S.STE = 2; %1 = dot lines for standard error; 2 = shadows
S.transp = 0.3; %transparency for standard errors shadow
S.time_real = time_sel(1:1026);
S.colorline = [1 0 0; 0.3686 0.6314 0.7412; 0.1882 0.4902 0.8118; 0.0784 0.1569 0.5255; 0 0 0];
if export_l == 1
    S.legendl = 0;
else
    S.legendl = 1;
end
S.x_lim = [-0.1 3.4]; % Set x limits
S.y_lim = ylimm; %Set y limits
S.lineplot = lineplot;
S.ROI_n = ROI;
S.condition_n = condition;
S.ROIs_labels = ROIN(ROI);
%rearranging the significant time-windows
signtp_col = [];
bum = [];

for ss = 1:4 %over contrasts between conditions
    %     signtp_col = cat(1,signtp_col,ones(size(SIGN{ROI,ss},1),1)*(ss+1)); %getting the color code (number of the significant time-windows for each contrast)
    sbam = size(SIGN{ROI,ss},1);
    for ll = 1:sbam %over the number of the significant time-windows for contrast ss
        if SIGN{ROI,ss}{ll,3}(1) > 0.350 && SIGN{ROI,ss}{ll,3}(1) < 2.5  %we do not want to plot very small spurious random differences between conditions; they are available though in the supplementary tables
            bum = cat(2,bum,SIGN{ROI,ss}(ll,3)); %concatenating all significant time-windows, for one contrast at a time
            signtp_col = cat(1,signtp_col,(ss+1)); %getting the color code (number of the significant time-windows for each contrast)
        end
    end
end
S.signtp = bum;
if col_l == 1
    S.signtp_col = signtp_col;
else
    S.signtp_col = [];
end
waveform_plotting_local_v2(S) %actual function

if export_l == 1
    if source_leak_corr == 1
        exportgraphics(gcf,['/Users/au550322/Documents/AarhusUniversitet/Postdoc/DataCollection_AuditoryPatternRecognition2020/Papers/SystematicVariationSequence/NatureCommunications/FirstRevision/Images_ToGemma/FinalFigures/Figure_SX_SourceLeakage/ROI_' ROIN{ROI} '_SL.pdf'],'Resolution',300)
    else
        exportgraphics(gcf,['/Users/au550322/Documents/AarhusUniversitet/Postdoc/DataCollection_AuditoryPatternRecognition2020/Papers/SystematicVariationSequence/NatureCommunications/FirstRevision/Images_ToGemma/Timeseries/Previous_ROIs/ROI_' ROIN{ROI} '.pdf'],'Resolution',300)
    end
end

%% AAL - WITH AND WITHOUT SOURCE LEAKAGE CORRECTION

col_l = 1; %1 for significant time-windows with different colors; 1 = only grey color
lineplot = [20 1]; %number instructing where to place the lines showing significant time-windows; leave empty [] if you want to have shaded colors instead

SLL = 0; %1 = source leakage correciton; 0 = no source leakage correction
ROI = 1; %ROIs (they are 6)
condition = [1:5]; %1 = Old; 2 = NewT1; 3 = NewT2; 4 = NewT3; 5 = NewT4;
export_l = 0; %1 = export images; 0 = not

addpath('/Users/au550322/Documents/AarhusUniversitet/Postdoc/DataCollection_AuditoryPatternRecognition2020/Papers/SystematicVariationSequence/Codes_Images_Local');
close all
clear ROIN
if ROI < 3
    if SLL == 1
        ylimm = [-0.085 0.04]; %amplitude limits; leave empty [] for automatic adjustment
    else
        ylimm = [-130 50];
    end
else
    if SLL == 1
        ylimm = [-0.055 0.038]; %amplitude limits; leave empty [] for automatic adjustment
    else
        ylimm = [-60 40];
    end
end
% ylimm = [];
S = [];
S.conds = {'Old ','NewT1','NewT2','NewT3','NewT4'};
ROIN{1} = 'HESCHLL'; ROIN{2} = 'HESCHLR'; ROIN{3} = 'HIPPL'; ROIN{4} = 'HIPPR'; ROIN{5} = 'ACC'; ROIN{6} = 'MC';

if SLL == 1
    load('/Users/au550322/Documents/AarhusUniversitet/Postdoc/DataCollection_AuditoryPatternRecognition2020/Papers/SystematicVariationSequence/NatureCommunications/FirstRevision/14_09_2023/SL_AAL/sources_leakagecorrected_AAL.mat');
    load('/Users/au550322/Documents/AarhusUniversitet/Postdoc/DataCollection_AuditoryPatternRecognition2020/Papers/SystematicVariationSequence/NatureCommunications/FirstRevision/14_09_2023/SL_AAL/Sign.mat');
    data2 = tso;
else
    load('/Users/au550322/Documents/AarhusUniversitet/Postdoc/DataCollection_AuditoryPatternRecognition2020/Papers/SystematicVariationSequence/NatureCommunications/FirstRevision/14_09_2023/AAL_noSL/dataAALfinal.mat');
    load('/Users/au550322/Documents/AarhusUniversitet/Postdoc/DataCollection_AuditoryPatternRecognition2020/Papers/SystematicVariationSequence/NatureCommunications/FirstRevision/14_09_2023/AAL_noSL/Sign.mat');
    data2 = tso2;
end
load('/Users/au550322/Documents/AarhusUniversitet/Postdoc/DataCollection_AuditoryPatternRecognition2020/Papers/SystematicVariationSequence/Codes_Images_Local/time.mat'); %loading time
%structure for the function
data2 = permute(data2,[1 2 4 3]);
S.data = data2(:,1:1026,:,:);
S.STE = 2; %1 = dot lines for standard error; 2 = shadows
S.transp = 0.3; %transparency for standard errors shadow
S.time_real = time_sel(1:1026);
S.colorline = [1 0 0; 0.3686 0.6314 0.7412; 0.1882 0.4902 0.8118; 0.0784 0.1569 0.5255; 0 0 0];
if export_l == 1
    S.legendl = 0;
else
    S.legendl = 1;
end
S.x_lim = [-0.1 3.4]; % Set x limits
S.y_lim = ylimm; %Set y limits
S.lineplot = lineplot;
S.ROI_n = ROI;
S.condition_n = condition;
S.ROIs_labels = ROIN(ROI);
%rearranging the significant time-windows
signtp_col = [];
bum = [];
for ss = 1:4 %over contrasts between conditions
    %     signtp_col = cat(1,signtp_col,ones(size(SIGN{ROI,ss},1),1)*(ss+1)); %getting the color code (number of the significant time-windows for each contrast)
    sbam = size(SIGN{ROI,ss},1);
    for ll = 1:sbam %over the number of the significant time-windows for contrast ss
        if SIGN{ROI,ss}{ll,3}(1) > 0.350 && SIGN{ROI,ss}{ll,3}(1) < 2.5  %we do not want to plot very small spurious random differences between conditions; they are available though in the supplementary tables
            bum = cat(2,bum,SIGN{ROI,ss}(ll,3)); %concatenating all significant time-windows, for one contrast at a time
            signtp_col = cat(1,signtp_col,(ss+1)); %getting the color code (number of the significant time-windows for each contrast)
        end
    end
end
S.signtp = bum;
if col_l == 1
    S.signtp_col = signtp_col;
else
    S.signtp_col = [];
end
waveform_plotting_local_v2(S) %actual function

if export_l == 1
    if SLL == 1
        exportgraphics(gcf,['/Users/au550322/Documents/AarhusUniversitet/Postdoc/DataCollection_AuditoryPatternRecognition2020/Papers/SystematicVariationSequence/NatureCommunications/FirstRevision/Images_ToGemma/FinalFigures/Figure_SX_SourceLeakage_AAL/ROI_' ROIN{ROI} '_SL.pdf'],'Resolution',300)
    else
        exportgraphics(gcf,['/Users/au550322/Documents/AarhusUniversitet/Postdoc/DataCollection_AuditoryPatternRecognition2020/Papers/SystematicVariationSequence/NatureCommunications/FirstRevision/Images_ToGemma/FinalFigures/Figure_SX_AAL_6ROIs/ROI_' ROIN{ROI} '_SL.pdf'],'Resolution',300)
    end
end

%exporting data for NC source data
DUMM = zeros(31,1026);
cnt2 = 0;
DUMM(1,:) = S.time_real;
for ss = 1:6 %over ROIs
    for iii = 1:5 %over conditions
        cnt2 = cnt2 + 1;
        bum = S.data(ss,:,:,iii);
        DUMM(cnt2 + 1,:) = squeeze(mean(bum,3));
    end
end
outdir = '/Users/au550322/Documents/AarhusUniversitet/Postdoc/DataCollection_AuditoryPatternRecognition2020/Papers/SystematicVariationSequence/NatureCommunications/FirstRevision/SourceData';
PDn = table(DUMM); %table
writetable(PDn,[outdir '/MEGsurces_SelectedAALROIs.xlsx'],'Sheet',1); %printing excel file

%% PREVIOUS DECODING ANALYSIS - NEW PLOTTING FUNCTION

col_l = 1; %1 for significant time-windows with different colors; 1 = only grey color
lineplot = [30 1]; %number instructing where to place the lines showing significant time-windows; leave empty [] if you want to have shaded colors instead
ylimm = []; %amplitude limits; leave empty [] for automatic adjustment
export_l = 1; %1 = export images; 0 = not

addpath('/Users/au550322/Documents/AarhusUniversitet/Postdoc/DataCollection_AuditoryPatternRecognition2020/Papers/SystematicVariationSequence/Codes_Images_Local');
close all
clear ROIN
S = [];
S.conds = {'NewT1','NewT2','NewT3','NewT4'};

load('/Users/au550322/Documents/AarhusUniversitet/Postdoc/DataCollection_AuditoryPatternRecognition2020/Papers/SystematicVariationSequence/NatureCommunications/FirstRevision/14_09_2023/Decoding/SIGN.mat');
load('/Users/au550322/Documents/AarhusUniversitet/Postdoc/DataCollection_AuditoryPatternRecognition2020/Papers/SystematicVariationSequence/NatureCommunications/FirstRevision/14_09_2023/Decoding/decpack.mat');
load('/Users/au550322/Documents/AarhusUniversitet/Postdoc/DataCollection_AuditoryPatternRecognition2020/Papers/SystematicVariationSequence/Codes_Images_Local/time.mat'); %loading time
%structure for the function
S.data = LOC(:,1:1026,:,:);
S.STE = 2; %1 = dot lines for standard error; 2 = shadows
S.transp = 0.3; %transparency for standard errors shadow
S.time_real = time_sel(1:1026);
S.colorline = [0.3686 0.6314 0.7412; 0.1882 0.4902 0.8118; 0.0784 0.1569 0.5255; 0 0 0];
if export_l == 1
    S.legendl = 0;
else
    S.legendl = 1;
end
S.x_lim = [-0.1 3.4]; % Set x limits
S.y_lim = ylimm; %Set y limits
S.lineplot = lineplot;
S.ROI_n = 1;
S.condition_n = [1:4];
S.ROIs_labels = {''};
%rearranging the significant time-windows
signtp_col = [];
bum = [];
bigSIGN{1,1} = signcontr1; bigSIGN{1,2} = signcontr2; bigSIGN{1,3} = signcontr3; bigSIGN{1,4} = signcontr4;
for ss = 1:4 %over contrasts between conditions
    %     signtp_col = cat(1,signtp_col,ones(size(SIGN{ROI,ss},1),1)*(ss+1)); %getting the color code (number of the significant time-windows for each contrast)
    for ll = 1:length(bigSIGN{ss}) %over the number of the significant time-windows for contrast ss
        if time_sel(bigSIGN{1,ss}{ll}(1)) > 0.350 && time_sel(bigSIGN{1,ss}{ll}(end)) < 2.5 %we do not want to plot very small spurious random differences between conditions; they are available though in the supplementary tables
            dumbum = {[time_sel(bigSIGN{1,ss}{ll}(1)) time_sel(bigSIGN{1,ss}{ll}(end))]};
            bum = cat(2,bum,dumbum); %concatenating all significant time-windows, for one contrast at a time
            signtp_col = cat(1,signtp_col,(ss)); %getting the color code (number of the significant time-windows for each contrast)
        end
    end
end
S.signtp = bum;
if col_l == 1
    S.signtp_col = signtp_col;
else
    S.signtp_col = [];
end
waveform_plotting_local_v2(S) %actual function

if export_l == 1
    exportgraphics(gcf,['/Users/au550322/Documents/AarhusUniversitet/Postdoc/DataCollection_AuditoryPatternRecognition2020/Papers/SystematicVariationSequence/NatureCommunications/FirstRevision/Previous/Images_ToGemma/FinalFigures/Decoding/Decodingnew.pdf'],'Resolution',300)
end

%exporting data for NC source data
DUMM = zeros(5,1026);
DUMM(1,:) = S.time_real;
for ii = 1:4
    if ii == 4
    bum = S.data(:,:,1:82,ii);
    else
        bum = S.data(:,:,:,ii);
    end
    DUMM(ii+1,:) = squeeze(mean(bum,3));
end
outdir = '/Users/au550322/Documents/AarhusUniversitet/Postdoc/DataCollection_AuditoryPatternRecognition2020/Papers/SystematicVariationSequence/NatureCommunications/FirstRevision/SourceData';
PDn = table(DUMM); %table
writetable(PDn,[outdir '/Decoding_Timeseries.xlsx'],'Sheet',1); %printing excel file

%% renaming AAL parcels images

load('/Users/au550322/Documents/AarhusUniversitet/Postdoc/DataCollection_AuditoryPatternRecognition2020/Papers/SystematicVariationSequence/NatureCommunications/FirstRevision/14_09_2023/AAL_labels');
list = dir('/Users/au550322/Documents/AarhusUniversitet/Postdoc/DataCollection_AuditoryPatternRecognition2020/Papers/SystematicVariationSequence/NatureCommunications/FirstRevision/Images_ToGemma/FinalFigures/Figure_SX_AAL/AAL_Parcels_Workbench/*gii');
for ii = 1:size(lab,1)
    % Define the old file name and the new file name
    if mod(ii, 2) == 1
        oldFileName = [list(ii).folder '/' num2str(ii) '_L' list(ii).name(end-8:end)];
        newFileName = [list(ii).folder '/' lab(ii,:) list(ii).name(end-8:end)];
    else
        oldFileName = [list(ii).folder '/' num2str(ii) '_R' list(ii).name(end-8:end)];
        newFileName = [list(ii).folder '/' lab(ii,:) list(ii).name(end-8:end)];
    end

    % Use the movefile function to rename the file
    success = movefile(oldFileName, newFileName);
end

%% AAL - FOCUS ON PREDICTION ERROR IMAGE

col_l = 1; %1 for significant time-windows with different colors; 1 = only grey color
lineplot = [20 1]; %number instructing where to place the lines showing significant time-windows; leave empty [] if you want to have shaded colors instead
allconds = 0; %1 = all NT conditions together; 0 = 1 NT condition at a time
ROI = 6; %ROIs (they are 6)
condition = [1:5]; %1 = Old; 2 = NewT1; 3 = NewT2; 4 = NewT3; 5 = NewT4;
export_l = 1; %1 = export images; 0 = not
COL = [0.3686 0.6314 0.7412; 0.1882 0.4902 0.8118; 0.0784 0.1569 0.5255; 0 0 0];

addpath('/Users/au550322/Documents/AarhusUniversitet/Postdoc/DataCollection_AuditoryPatternRecognition2020/Papers/SystematicVariationSequence/Codes_Images_Local');
close all
clear ROIN
if ROI < 3
    ylimm = [-130 50];
else
    ylimm = [-60 40];
end
ROIN{1} = 'HESCHLL'; ROIN{2} = 'HESCHLR'; ROIN{3} = 'HIPPL'; ROIN{4} = 'HIPPR'; ROIN{5} = 'ACC'; ROIN{6} = 'MC';
load('/Users/au550322/Documents/AarhusUniversitet/Postdoc/DataCollection_AuditoryPatternRecognition2020/Papers/SystematicVariationSequence/NatureCommunications/FirstRevision/14_09_2023/AAL_noSL/dataAALfinal.mat');
load('/Users/au550322/Documents/AarhusUniversitet/Postdoc/DataCollection_AuditoryPatternRecognition2020/Papers/SystematicVariationSequence/NatureCommunications/FirstRevision/14_09_2023/AAL_noSL/Sign.mat');
data2 = tso2;
load('/Users/au550322/Documents/AarhusUniversitet/Postdoc/DataCollection_AuditoryPatternRecognition2020/Papers/SystematicVariationSequence/Codes_Images_Local/time.mat'); %loading time
%structure for the function
data2 = permute(data2,[1 2 4 3]);
data3 = data2(:,:,:,2:end);
% ylimm = [];

if allconds == 1 %all NT conditions together
    cc = 1;
    S = [];
    S.conds = {'NewT1','NewT2','NewT3','NewT4'};
    
    S.data = data3(:,1:1026,:,:);
    S.STE = 2; %1 = dot lines for standard error; 2 = shadows
    S.transp = 0.3; %transparency for standard errors shadow
    S.time_real = time_sel(1:1026);
    S.colorline = COL(cc:4,:);
    if export_l == 1
        S.legendl = 0;
    else
        S.legendl = 1;
    end
    S.x_lim = [-0.1 3.4]; % Set x limits
    S.y_lim = ylimm; %Set y limits
    S.lineplot = lineplot;
    S.ROI_n = ROI;
    S.condition_n = [cc:4];
    S.ROIs_labels = ROIN(ROI);
    %rearranging the significant time-windows
    signtp_col = [];
    bum = [];
    S.signtp = [];
    waveform_plotting_local_v2(S) %actual function
    
    if export_l == 1
        exportgraphics(gcf,['/Users/au550322/Documents/AarhusUniversitet/Postdoc/DataCollection_AuditoryPatternRecognition2020/Papers/SystematicVariationSequence/NatureCommunications/FirstRevision/Images_ToGemma/FinalFigures/Figure_XX_FocusPredictionError/ROI_' ROIN{ROI} '_All_NTs.pdf'],'Resolution',300)
    end
else %otherwise, one condition at a time
    for cc = 1:4
        S = [];
        S.conds = {'NewT1','NewT2','NewT3','NewT4'};
        
        S.data = data3(:,1:1026,:,:);
        S.STE = 2; %1 = dot lines for standard error; 2 = shadows
        S.transp = 0.3; %transparency for standard errors shadow
        S.time_real = time_sel(1:1026);
        S.colorline = COL(cc,:);
        if export_l == 1
            S.legendl = 0;
        else
            S.legendl = 1;
        end
        S.x_lim = [-0.1 3.4]; % Set x limits
        S.y_lim = ylimm; %Set y limits
        S.lineplot = lineplot;
        S.ROI_n = ROI;
        S.condition_n = [cc];
        S.ROIs_labels = ROIN(ROI);
        %rearranging the significant time-windows
        signtp_col = [];
        bum = [];
        S.signtp = [];
       
        waveform_plotting_local_v2(S) %actual function
        
        if export_l == 1
            exportgraphics(gcf,['/Users/au550322/Documents/AarhusUniversitet/Postdoc/DataCollection_AuditoryPatternRecognition2020/Papers/SystematicVariationSequence/NatureCommunications/FirstRevision/Images_ToGemma/FinalFigures/Figure_XX_FocusPredictionError/ROI_' ROIN{ROI} '_NT' num2str(cc) '.pdf'],'Resolution',300)
        end
        close all
    end
end

%%

%% *** FOCUS ON PREDICTION ERROR ***

%% STRONGER RESPONSES TO FIRST SOUND WHERE VARIATION IS INSERTED IN HIPPOCAMPUS AND ACC, BUT NOT IN HESCHL'S GYRUS

ROIN{1} = 'LHG'; ROIN{2} = 'RHG'; ROIN{3} = 'LHP'; ROIN{4} = 'RHP'; ROIN{5} = 'ACC'; ROIN{6} = 'MC'; % ROIs in dum2
HG = {[159:169] [246:256] [344:354] [434:444]}; %time-points (already with plus/minus 20ms) of the peaks for ACL
HPCC = {[176:186] [268:278] [363:373] [454:464]}; %time-points (already with plus/minus 20ms) of the peaks for VMPFC

addpath('/Users/au550322/Documents/AarhusUniversitet/Postdoc/DataCollection_AuditoryPatternRecognition2020/Papers/SystematicVariationSequence/NatureCommunications/SecondRevision_AcceptedInPrinciple/Gemma_Figures_TEMPORARY_TOBEDISCARDED')
if ~exist('dum2','var')
    load('/Users/au550322/Documents/AarhusUniversitet/Postdoc/DataCollection_AuditoryPatternRecognition2020/Papers/SystematicVariationSequence/NatureCommunications/FirstRevision/14_09_2023/AAL_data.mat');
    load('/Users/au550322/Documents/AarhusUniversitet/Postdoc/DataCollection_AuditoryPatternRecognition2020/Papers/SystematicVariationSequence/Codes_Images_Local/time.mat'); %loading time
    %loading AAL parcels label
    %     AAL_lab = importdata('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/aal_ROIs_cerebellum.txt');
    dum2 = dataa;
    clear dataa
end
conds{1} = ' Mem '; conds{2} = ' NewT1 '; conds{3} = ' NewT2 '; conds{4} = ' NewT3 '; conds{5} = ' NewT4 ';
tso2 = zeros(6,size(dum2,2),5,size(dum2,4));
tso = zeros(6,size(dum2,2),size(dum2,4));
vectt = [79 80 37 38 31 32 33 34]; %8 selected ROIs
% tso2 = zeros(size(tso,1),size(tso,2),size(tso,3),size(tso,4));
tso2(1,:,:,:) = dum2(vectt(1),:,:,:); %ACL
tso2(2,:,:,:) = dum2(vectt(2),:,:,:); %ACR
tso2(3,:,:,:) = dum2(vectt(3),:,:,:); %HIPPL
tso2(4,:,:,:) = dum2(vectt(4),:,:,:); %HIPPR
tso2(5,:,:,:) = mean(dum2(vectt(5:6),:,:,:),1); %ACC
tso2(6,:,:,:) = mean(dum2(vectt(7:8),:,:,:),1); %MC
SIGN = cell(length(ROIN),3,2); %ROIs x conditions - 1 x 2 (ANOVA statistics and post-hoc tests)
cntp = 0;
PP = zeros(18,1);
for ss = 1:size(tso2,1) %over ROIs
    for cc = 1:4 %over NTX conditions
        cnt = 0;
        dumbo = zeros(83,5-cc); %subjects x peaks
        for pp = cc:4 %over the peaks (that are 4 for cc = 2 (NT1), 3 for cc = 3 (NT2), etc.)
            cnt = cnt + 1;
            if ss < 3 %if Heschl's gyrus (either left or right)
                dumbo(:,cnt) = squeeze(mean(tso2(ss,HG{pp},cc+1,:),2)); %extracting ROI ss, time-window of Heschl's gyrus, condition NTcc, all participants
            else
                dumbo(:,cnt) = squeeze(mean(tso2(ss,HPCC{pp},cc+1,:),2)); %extracting ROI ss, time-window of higher-order ROIs, condition NTcc, all participants
            end
        end
        %plotting
        cl = [0 0 0.8; 0 0 0.8; 0 0 0.8; 0 0 0.8; 0 0 0.8;];
        %         dumbo = dumbo';
        datadum = cell(1,size(dumbo,2));
        for ii = 1:size(dumbo,2)
            datadum{1,ii} = dumbo(:,ii);
        end
        figure
        rm_raincloud2(datadum',cl)
        grid minor
        set(gcf,'color','w')
        set(gcf,'Position',[200,200,400,550])
        box on
        xlabel('Amplitude'); %set(gca,'YDir','normal');
        title(['ROI ' ROIN{ss} ' - NT' num2str(cc)])
       
        exportgraphics(gcf,['/Users/au550322/Documents/AarhusUniversitet/Postdoc/DataCollection_AuditoryPatternRecognition2020/Papers/SystematicVariationSequence/NatureCommunications/SecondRevision_AcceptedInPrinciple/Gemma_Figures_TEMPORARY_TOBEDISCARDED/NT' num2str(cc) '_ROI_' ROIN{ss} '.pdf'],'Resolution',300)
    end
end

%%
