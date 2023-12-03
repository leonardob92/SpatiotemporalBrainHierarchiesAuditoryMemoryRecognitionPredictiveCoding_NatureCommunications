%%

%% SPATIOTEMPORAL BRAIN HIERARCHIES OF AUDITORY MEMORY RECOGNITION AND PREDICTIVE CODING (DATA COLLECTION CODE: AUDITORY PATTERN RECOGNITION 2020) - LEONARDO BONETTI

%%

%% NATURE COMMUNICATIONS - REVISION I

%%

%% *** START UP FUNCTIONS.. (LBPD_startup_D) ***

%starting up some functions for LBPD toolbox.

%starting up some of the functions that I wrote for the following analysis
pathl = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD'; %path to stored functions (THIS MUST BECOME A FOLDER ON YOUR OWN COMPUTER)
addpath(pathl);
LBPD_startup_D(pathl);

%%

%%

%% *** MEG SENSORS ANALYSIS ***

%% Step 1: Getting average weights on the MEG channels for each of the significant time-window of the decoding and averaging them

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


%% Step 2: Final averaged weights for the decoding and extracting the of highest values in absolute terms (higher than the average plus one standard deviation), independently for magnetometers and gradiometers

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

%% Step 3a - magnetometers: They shows two main clusters (one per hemisphere) which are analysed independently, as follows (only magnetometers here, gradiometers are below, after the source reconstruction)

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

%% Step 3c - gradiometers (same analysis as for the magnetometers but without source reconstruction)

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

%% Step 3d - Digression - Topoplots for the peaks of differential activity between Old and New conditions (for supplementary material)

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
            S.waveform_singlechannels_label = 0; %1 to plot single channel waveforms
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

%% Step 4: Source reconstruction for the peaks of differential activity between Old and New conditions

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


%% Computing correlation between the AAL parcels time series (at group level)

load('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/Papers/Systematic_Variations_DCM/NatureCommunications_Revision/Revision_I/MEGSensors_FunctionalParcellation/AAL/AAL_data.mat');
dataa2 = nanmean(nanmean(dataa,4),3); %averaging oer subjects and conditions
corrm = corr(dataa2');
%correlation matrix
figure
imagesc(corrm)
axis square
colorbar
set(gcf,'Color','w')
set(gca,'CLim',[-1 1])
%correlation matrix L and R hemispheres
% order = [1:2:116 116:-2:1];
% dataa2b = dataa2(order,:);
% corrm = corr(dataa2b');
% figure
% imagesc(corrm)
% axis square
% colorbar
% set(gcf,'Color','w')

%% SOURCE LEAKAGE CORRECTION - AAL

cond = 5; %condition that you want; 1 = Old; 2-5 = NewTX (X = 1-4)
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

%% STATISTICS AFTER SOURCE LEAKAGE CORRECTION

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

%% QUANTIFYING HIERARCHY IN THE BRAIN (DYNAMIC CAUSAL MODELLING)

%% *** DYNAMIC CAUSAL MODELLING ***

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

%%

%% PREPARING DATA II (ORIGINAL PARCELLATION)

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

%%

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

%%

%% ACTUAL INVERSION OF THE DCM MODELS

%% setting up the cluster

addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/Cluster_ParallelComputing') %add the path to the function that submits the jobs to the cluster
clusterconfig('scheduler', 'cluster'); %'none' or 'cluster'
clusterconfig('long_running', 1); %there are different queues for the cluster depending on the number and length of the jobs you want to submit 
clusterconfig('slot', 1); %slot in the queu
addpath('/projects/MINDLAB2020_MEG-AuditoryPatternRecognition/scripts/leonardo/DCM')

%% RESHAPING CODES FROM MARTIN DIETZ'S FOR DCM FOR EVOKED RESPONSES

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
export_fig(['/aux/MINDLAB2021_MEG-TempSeqAges/19_10_2023/PP_Tone' num2str(tonel) '_contrast' num2str(condd) '.eps'])
export_fig(['/aux/MINDLAB2021_MEG-TempSeqAges/19_10_2023/PP_Tone' num2str(tonel) '_contrast' num2str(condd) '.png'])
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
export_fig(['/aux/MINDLAB2021_MEG-TempSeqAges/19_10_2023/PEP_Tone' num2str(tonel) '_contrast' num2str(condd) '.eps'])
export_fig(['/aux/MINDLAB2021_MEG-TempSeqAges/19_10_2023/PEP_Tone' num2str(tonel) '_contrast' num2str(condd) '.png'])

%% codes for additional brain figure (for plotting the models within a brain template)

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

%% METHOD FIGURE - DEPICTION OF 8-MM PARCELLATION

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

load('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/MNI152_8mm_coord_dyi.mat')
name_model = 'Parc_AAL';

%preparing color specifications
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

%% *** TIME-FREQUENCY ANALYSIS ***

%% setting up the cluster

addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/Cluster_ParallelComputing') %add the path to the function that submits the jobs to the cluster
clusterconfig('scheduler', 'cluster');
clusterconfig('long_running', 1); %there are different queues for the cluster depending on the number and length of the jobs you want to submit 
clusterconfig('slot', 6); %slot in the queu

%% PROPER INDUCED RESPONSES - MEG CHANNELS (8 SELECTED)

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
            P2dum(:,:,:,tt) = morlet_transform(data2(:,:,tt),D.time(1:1026),f); %frequency decomposition over single trials
        end
        P2(:,:,:,jj,ii) = mean(P2dum,4); %averaging TF results over trials
        disp(['subj ' num2str(ii) ' - cond ' num2str(jj)])
    end
    disp(ii)
end
time = D.time;
save('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/Papers/Systematic_Variations_DCM/NatureCommunications_Revision/Revision_I/Time_Frequency_Analysis/TF_ERFs/Time_Frequency_AllSubjects_AveragedTrials_MEGsensors.mat','P2','time','f');

%% INDUCED RESPONSES - MEG SOURCES - AAL
%(COMPUTING TIME-FREQUENCY ANALYSIS ON SIGNLE VOXELS AND THEN AVERAGING THE RESULTS)

%setting up the cluster
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/Cluster_ParallelComputing') %add the path to the function that submits the jobs to the cluster
clusterconfig('scheduler', 'cluster');
clusterconfig('long_running', 1); %there are different queues for the cluster depending on the number and length of the jobs you want to submit 
clusterconfig('slot', 1); %slot in the queu

%% AAL - here you provide either coordinates in MNI space or the AAL ROIs (following the AAL order)

ROIn = 7;

vecttname = {'Occ_Sup_L';'Occ_Sup_R';'HeschlL';'HeschlR';'HippL';'HippR';'ACC';'MC'};
namee = vecttname{ROIn}; %if you use the AAL ROI(s), here you can specify the name to be used when saving the data
vectt = {[49];[50];[79];[80];[37];[38];[31 32];[33 34]}; %8 selected ROIs

%setting up the cluster (this one should be local.. otherwise it does not work for some unknown reasons.. 
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/Cluster_ParallelComputing') %add the path to the function that submits the jobs to the cluster
clusterconfig('scheduler', 'none'); % 'none' or 'cluster'
clusterconfig('long_running', 1); %there are different queues for the cluster depending on the number and length of the jobs you want to submit 
clusterconfig('slot', 1); %slot in the queu

outdir = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/Papers/Systematic_Variations_DCM/NatureCommunications_Revision/Revision_I/Time_Frequency_Analysis/LBPD_Parcels/FinalImages';
S = [];
S.Aarhus_clust = 2; %0 = working locally; integer number (i.e. 1) = sending one job for each subject to the Aarhus cluster (the number you insert here corresponds to the slots of memory that you allocate for each job.
S.average_trials = 0; %1 = frequency decomposition after averaging the trials (evoked responses); 0 = frequency decomposition for single trials and then average of the time-frequency results (induced responses)
S.f = 1:1:60; %frequencies used
load('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/time_normal.mat');
S.time = time(1:1026);
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

condition = 1; % 1 = Old; 2 = NewT1; 3 = NewT2; 4 = NewT3; 5 = NewT4
fff = 1:60; %frequencies to be plotted
exportl = 1; % 1 = export figures
ttestl = 1; %1 = t-tests Old vs NewT1; 0 = single condition

CHAN{1} = 'MEG211'; CHAN{2} = 'MEG1311'; CHAN{3} = 'MEG241'; CHAN{4} = 'MEG1331'; CHAN{5} = 'MEG1631'; CHAN{6} = 'MEG2441'; CHAN{7} = 'MEG1921'; CHAN{8} = 'MEG2341'; %channels name
if ttestl ~= 1
    for ii = 1:8
        chani = ii; % 1 = MEG211; 2 = MEG1311; 3 = MEG241; 4 = MEG1331; 5 = MEG1631; 6 = MEG2441; 7 = MEG192; 8 = MEG2341
        %loading data
        if ~exist('P2','var')
            load('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/Papers/Systematic_Variations_DCM/NatureCommunications_Revision/Revision_I/Time_Frequency_Analysis/TF_ERFs/Time_Frequency_AllSubjects_AveragedTrials_MEGsensors.mat');
        end
        %mean over subjects (for plotting purposes) and removing pre-stimulus time
        Pold = mean(P2,5); %mean over subjects
        %plotting
        figure
        imagesc(time,f(fff),squeeze(Pold(chani,fff,:,condition)))
        set(gca,'YDir','normal') %plotting frequencies in descending order
        xlabel('time (s)'); ylabel('f (Hz)');
        colorbar
        caxis([0 6000])
        set(gcf,'color','w')
        %colormap with white for 0 values
        x = []; x.bottom = [0 0 0.5]; x.botmiddle = [0 0.5 1]; x.middle = [1 1 1]; x.topmiddle = [1 0 0]; x.top = [0.6 0 0]; %red - blue
        colormap(bluewhitered_PD(0,x))
        title(['Chan ' CHAN{chani} ' - cond ' num2str(condition)])
        if exportl == 1
            export_fig(['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/Papers/Systematic_Variations_DCM/NatureCommunications_Revision/Revision_I/Time_Frequency_Analysis/LBPD_Parcels/FinalImages/MEGsensors_Induced/' CHAN{chani} '_Cond' num2str(condition) '.png'])
            export_fig(['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/Papers/Systematic_Variations_DCM/NatureCommunications_Revision/Revision_I/Time_Frequency_Analysis/LBPD_Parcels/FinalImages/MEGsensors_Induced/' CHAN{chani} '_Cond' num2str(condition) '.eps'])
        end
    end
else
    %loading data
    if ~exist('P2','var')
        load('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/Papers/Systematic_Variations_DCM/NatureCommunications_Revision/Revision_I/Time_Frequency_Analysis/TF_ERFs/Time_Frequency_AllSubjects_AveragedTrials_MEGsensors.mat');
    end
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
        writetable(PDn,['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/Papers/Systematic_Variations_DCM/NatureCommunications_Revision/Revision_I/Time_Frequency_Analysis/LBPD_Parcels/FinalImages/MEGsensors_Induced/TF' CHAN{cc} '_cond1_vs_cond2.xlsx'],'Sheet',1); %printing excel file
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
        export_fig(['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/Papers/Systematic_Variations_DCM/NatureCommunications_Revision/Revision_I/Time_Frequency_Analysis/LBPD_Parcels/FinalImages/MEGsensors_Induced/' CHAN{cc} '_Cond1_vs_Cond2.png'])
        export_fig(['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/Papers/Systematic_Variations_DCM/NatureCommunications_Revision/Revision_I/Time_Frequency_Analysis/LBPD_Parcels/FinalImages/MEGsensors_Induced/' CHAN{cc} '_Cond1_vs_Cond2.eps'])
    end
end

%% plotting SELECTED AAL ROIs PLUS OCCIPITAL SUPERIOR L or R (ALSO FROM AAL)

condition = 2; % 1 = Old; 2 = NewT1; 3 = NewT2; 4 = NewT3; 5 = NewT4; 0 = Old vs NewTcc (all categories of New)
fff = 1:60; %frequencies to be plotted
listn = 8; %index of the analysis to be run (with reference to the variable "list")
loadl = 0; %1 for loading; 0 for not loading, in case you already loaded the data and want to plot different conditions
exportl = 0; %1 for exporting figures
jes = [10 1400];
% jes = [];

if loadl == 1
    %loading data
    list = dir(['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/Papers/Systematic_Variations_DCM/NatureCommunications_Revision/Revision_I/Time_Frequency_Analysis/LBPD_Parcels/FinalImages/Time_*']);
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
        export_fig(['/aux/MINDLAB2021_MEG-TempSeqAges/20_10_2023/TF_' list(listn).name(24:end) '_Cond' num2str(condition) '.png'])
        export_fig(['/aux/MINDLAB2021_MEG-TempSeqAges/20_10_2023/TF_' list(listn).name(24:end) '_Cond' num2str(condition) '.eps'])
    end
    %export data for NC data source
    if condition < 3 %only M and NT1
       bumbum = squeeze(Pold(fff,:,condition));
       PDn = table(bumbum); %table (not getting the last column (8) because it contains a large matrix useful for plotting purposes)
       writetable(PDn,['/aux/MINDLAB2021_MEG-TempSeqAges/09_11_2023/TF_' list(listn).name(24:end) '_Cond' num2str(condition) '.xlsx'],'Sheet',1); %printing excel file
    end
else
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
        writetable(PDn,['/aux/MINDLAB2021_MEG-TempSeqAges/20_10_2023/TF_' list(listn).name(24:end) '_cond1_vs_cond' num2str(cc+1) '.xlsx'],'Sheet',1); %printing excel file
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
        title([list(listn).name(end-10:end) ' - Cond 1 vs Cond ' num2str(cc+1)])
        if exportl == 1
            export_fig(['/aux/MINDLAB2021_MEG-TempSeqAges/20_10_2023/Thresh_TF_' list(listn).name(24:end) '_cond1_vs_cond' num2str(cc+1) '.png'])
            export_fig(['/aux/MINDLAB2021_MEG-TempSeqAges/20_10_2023/Thresh_TF_' list(listn).name(24:end) '_cond1_vs_cond' num2str(cc+1) '.eps'])
        end
        %export data for NC data source
        if cc == 1 %only M vs NT1
            bumbum = squeeze(T2(fff,:));
            PDn = table(bumbum); %table (not getting the last column (8) because it contains a large matrix useful for plotting purposes)
            writetable(PDn,['/aux/MINDLAB2021_MEG-TempSeqAges/09_11_2023/TF_' list(listn).name(24:end) '_cond1_vs_cond' num2str(cc+1) '.xlsx'],'Sheet',1); %printing excel file
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
        if exportl == 1
            export_fig(['/aux/MINDLAB2021_MEG-TempSeqAges/20_10_2023/TF_' list(listn).name(24:end) '_cond1_vs_cond' num2str(cc+1) '.png'])
            export_fig(['/aux/MINDLAB2021_MEG-TempSeqAges/20_10_2023/TF_' list(listn).name(24:end) '_cond1_vs_cond' num2str(cc+1) '.eps'])
        end
    end
end

%%

%%

%% *** *** *** 

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

list = dir('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/MainROIs/*GM.nii.gz');
ROIs_4 = [5 4 6 1 2 3]; %ROIs from 4-broad-ROIs

% %loading original parcellation
% load('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/Beam_abs_1_sens_1_freq_broadband/DCM/ROIs_MainActivity/FinalROIs/ROIsDCM.mat')
% %loading MNI coordinates for LBPD indices
% load('/scratch7/MINDLAB2017_MEG-LearningBach/DTI_Portis/Templates/MNI152_8mm_coord_dyi.mat');
% ROIs_4 = [1 2 4 5]; %ROIs from 4-broad-ROIs
Col = ['r','b','c','y'];
for iii = 2%:4
    
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

%% COMPARISON BETWEEN AAL ROIs AND THE FUNCTIONAL PARCELLATION

% 1)creating nifti images (in 8mm) for each AAL ROIs
% creating indepdent image files with AAL ROIs
% command to split one ROI (80) from AAL
% cmd = 'fslmaths aal_1mm.nii.gz -thr 80 -uthr 80 -bin 80.nii.gz';
% system(cmd)
%AAL
AAL_lab = importdata('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/aal_ROIs_cerebellum.txt');
for ii = 1:116 %over AAL ROIs
    cmd = ['fslmaths /projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/aal_8mm_try5.nii.gz -thr ' num2str(ii) ' -uthr ' num2str(ii) [' -bin /scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/Papers/Systematic_Variations_DCM/NatureCommunications_Revision/Revision_I/MEGSensors_FunctionalParcellation/AAL/AAL_parcels_nifti/' num2str(ii) '.nii.gz']];
    system(cmd)
end
%non-optimised but effective solution to rename the AAL ROIs
dirdir = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/Papers/Systematic_Variations_DCM/NatureCommunications_Revision/Revision_I/MEGSensors_FunctionalParcellation/AAL/AAL_parcels_nifti/';
for ii = 11:20 %over AAL ROIs
    dum = load_nii([dirdir '/' num2str(ii) '.nii.gz']);
    dum2 = dum.img;
    nii = make_nii(dum2,[8 8 8]); %making nifti image from 3D data matrix (and specifying the 8 mm of resolution)
    nii.img = dum2; %storing matrix within image structure
    nii.hdr.hist = dum.hdr.hist; %copying some information from maskk
    disp(['saving nifti image - AAL ROI ' num2str(ii)])
    save_nii(nii,['/aux/MINDLAB2021_MEG-TempSeqAges/26_09_2023/' AAL_lab{ii}(3:end) '.nii.gz']); %printing image
    disp(ii)
end
%AAL grey matter only
maskGM = load_nii('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/osl/std_masks/lb_combinedGM8mm_mask_thr.nii.gz');
dirdir = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/Papers/Systematic_Variations_DCM/NatureCommunications_Revision/Revision_I/MEGSensors_FunctionalParcellation/AAL/AAL_parcels_nifti/';
for ii = 1:116 %over AAL ROIs
    dum = load_nii([dirdir '/' num2str(ii) '.nii.gz']);
    dum2 = dum.img;
    dum2(maskGM.img~=1) = 0; %if there is no correspondance between AAL ROI and grey matter image of the brain, assigning 0
    nii = make_nii(dum2,[8 8 8]); %making nifti image from 3D data matrix (and specifying the 8 mm of resolution)
    nii.img = dum2; %storing matrix within image structure
    nii.hdr.hist = maskGM.hdr.hist; %copying some information from maskk
    disp(['saving nifti image - AAL ROI ' num2str(ii)])
    save_nii(nii,[dirdir '/' AAL_lab{ii} '_GM.nii.gz']); %printing image
    disp(ii)
end

% 2)check overlapping between AAL parcels and our parcellation 
%loading our parcellation
load('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/Beam_abs_1_sens_1_freq_broadband/DCM/ROIs_MainActivity/FinalROIs/ROIsDCM.mat')
%loading AAL parcels label
AAL_lab = importdata('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/aal_ROIs_cerebellum.txt');
MAT = zeros(90,6); %ROIs (AAL) x ROIs (our parcellation)
for ii = 1:90 %over AAL ROIs
    %getting LBPD indices of AAL ROIs
    S = [];
    S.input = 2; %1 = MNI coordinates; 2 = AAL ROIs; 3 = general image with non-zero values
    S.coordd = []; %coordinates in MNI space (x,y,z)
    S.AAL_ROIs = [ii]; %AAL ROIs numbers you want to use
    % S.image = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/Beam_abs_0_sens_1_freq_broadband/Contr_1_abs_0.nii.gz';
    S.image = [];
    %actual function
    idx_LBPD = From3DNifti_OrMNICoords_2_CoordMatrix_8mm_LBPD_D(S); %it provides only the indices
    dum = zeros(3559,1);
    dum(idx_LBPD) = 1; %getting indices as 0 and 1 in 3559 vector
    for pp = 1:size(ROIs_DCM,2) %over parcels in our parcellation
        MAT(ii,pp) = length(find(dum + ROIs_DCM(:,pp)==2)); %computing the number of overlapping voxels between AAL ROIs ii and our parcellation parcel pp
    end
end
[~,maxi] = max(MAT'); maxi = maxi'; %maxi shows the correspondance between parcellations
maxi(find(sum(MAT,2)==0)) = 7; %getting voxels that did not belong to any of our original parcel and assigning 7 to them
ROIs_to_AAL = cell(7,2);
for ii = 1:7 %over parcels in our parcellation
    ROIs_to_AAL{ii,1} = find(maxi==ii); %AAL indeces
    ROIs_to_AAL{ii,2} = AAL_lab(ROIs_to_AAL{ii,1}); %AAL labels
end
names2 = cell(7,1); names2(1:6) = names; names2(7) = {'No_Parcel'};
save('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/Papers/Systematic_Variations_DCM/NatureCommunications_Revision/Revision_I/MEGSensors_FunctionalParcellation/AAL/AAL_parcels_nifti/ROIs_to_AAL.mat','ROIs_to_AAL','names2')

%% SOURCE LEAKAGE CORRECTION - ORIGINAL PARCELLATION

cond = 5; %condition that you want; 1 = Old; 2-5 = NewTX (X = 1-4)
ROIs = [2 3 4 6]; %selected ROIs; 1 = MC, 2 = HITR, 3 = HITL, 4 = VMPFC, 5 = ACL, 6 = ACR
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

%% TIME-FREQUENCY ANALYSIS - FUNCTIONAL PARCELLATION

%here you provide a matrix in LBPD coordinates with 0s and 1s to index specific ROIs 

ROII{1} = 'MC'; ROII{2} = 'HITR'; ROII{3} = 'HITL'; ROII{4} = 'VMPFC'; ROII{5} = 'ACL'; ROII{6} = 'ACR';
%loading matrix in LBPD coordinates with 0s and 1s to index specific ROIs 
load('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/Beam_abs_1_sens_1_freq_broadband/DCM/ROIs_MainActivity/FinalROIs/ROIsDCM.mat');
for ii = 3%:6 %over ROIs
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
    S.outdir = ['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/Papers/Systematic_Variations_DCM/NatureCommunications_Revision/Revision_I/Time_Frequency_Analysis/LBPD_Parcels/' ROII{ii} '_SingleSubjects_SingleVoxels.mat'];
    %actual function
    jobid = job2cluster(@InducedResponses_Morlet_ROIs_LBPD_D,S);
end

%% STATISTICS AND PLOTTING - TF - ORIGINAL PARCELLATION
%LOADING INDEPENDENT SUBJECTS (COMPUTED INDEPENDENT VOXELS AND THEN AVERAGED IN THE ROIs)

condition = 1; % 1 = Old; 2 = NewT1; 3 = NewT2; 4 = NewT3; 5 = NewT5; 0 = Old - NewT1
contr = [1 2 3 4 5]; %conditions to be contrasted (NOW IT'S FROM 1 TO 5 BECAUSE I CREATED A LOOP FOR DOING ALL CONTRASTS), meaningul only if condition == 0
fff = 1:60; %frequencies to be plotted
ROIn = [1:6]; % 1 = MC; 2 = HITL; 3 = HITR; 4 = VMPFC; 5 = ACL; 6 = ACR with Morlet computed independently on each voxel and each trial and then the output was averaged (twice)
loadl = 1; %1 for loading; 0 for not loading, in case you already loaded the data and want to plot different conditions
exportl = 1; %1 for exporting figures

outdir = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/Papers/Systematic_Variations_DCM/NatureCommunications_Revision/Revision_I/Time_Frequency_Analysis/LBPD_Parcels/FinalImages';
ROII{1} = 'MC'; ROII{2} = 'HITL'; ROII{3} = 'HITR'; ROII{4} = 'VMPFC'; ROII{5} = 'ACL'; ROII{6} = 'ACR';
for ss = 1:length(ROIn)
    ROI = ROIn(ss);
    if loadl == 1
        if ROI == 1
            list = dir(['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/Papers/Systematic_Variations_DCM/NatureCommunications_Revision/Revision_I/Time_Frequency_Analysis/LBPD_Parcels/MC*.mat']);
        elseif ROI == 2
            list = dir(['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/Papers/Systematic_Variations_DCM/NatureCommunications_Revision/Revision_I/Time_Frequency_Analysis/LBPD_Parcels/HITL*.mat']);
        elseif ROI == 3
            list = dir(['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/Papers/Systematic_Variations_DCM/NatureCommunications_Revision/Revision_I/Time_Frequency_Analysis/LBPD_Parcels/HITR*.mat']);
        elseif ROI == 4
            list = dir(['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/Papers/Systematic_Variations_DCM/NatureCommunications_Revision/Revision_I/Time_Frequency_Analysis/LBPD_Parcels/VMPFC*.mat']);
        elseif ROI == 5
            list = dir(['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/Papers/Systematic_Variations_DCM/NatureCommunications_Revision/Revision_I/Time_Frequency_Analysis/LBPD_Parcels/ACL*.mat']);
        elseif ROI == 6
            list = dir(['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/Papers/Systematic_Variations_DCM/NatureCommunications_Revision/Revision_I/Time_Frequency_Analysis/LBPD_Parcels/ACR*.mat']);
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
        caxis([10 1400])
        %colormap with white for 0 values
        x = []; x.bottom = [0 0 0.5]; x.botmiddle = [0 0.5 1]; x.middle = [1 1 1]; x.topmiddle = [1 0 0]; x.top = [0.6 0 0]; %red - blue
        colormap(bluewhitered_PD(0,x))
        title([ROII{ROI} ' - cond ' num2str(condition)])
        if exportl == 1
            export_fig([outdir '/ROIs_Induced/ROIs_Induced_' ROII{ROI} '_Cond' num2str(condition) '.png'])
            export_fig([outdir '/ROIs_Induced/ROIs_Induced_' ROII{ROI} '_Cond' num2str(condition) '.eps'])
        end
    end
    if condition == 0
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
            %testing contrast results (correcting for multiple comparison) by performing Monte Carlo simulations
            P3 = zeros(size(T,1),size(T,2));
            P3(abs(T)>2) = 1; %threshold t-val = 2.6 corresponding to p-val < 0.01 (obtained by dividing 0.05 by the 4 comparisons employed here)
            thresh = 0;
            permut = 1000;
            threshMC = 0.001;
            perm_max = 1;
            t1 = f(fff); t2 = time;
            [ OUT ] = twoD_MCS_LBPD_D( P3, thresh, permut, threshMC, perm_max, t1 , t2 )
            %             PDn = cell2table(OUT); %table
            %             writetable(PDn,[outdir '/TF_cond' num2str(contr(1)) '_vs_cond' num2str(contr(cc+1)) '_' ROII{ROIn(ss)} '.xlsx'],'Sheet',1); %printing excel file
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
            if exportl == 1
                export_fig([outdir '/ROIs_Induced/TF_SignOnly_cond' num2str(contr(1)) '_vs_cond' num2str(contr(cc+1)) '_' ROII{ROIn(ss)} '.eps'])
                export_fig([outdir '/ROIs_Induced/TF_SignOnly_cond' num2str(contr(1)) '_vs_cond' num2str(contr(cc+1)) '_' ROII{ROIn(ss)} '.png'])
            end
        end            
    end
end

%%

%% *** ADDITIONAL FINE-TUNINGS OF PREVIOUS ANALYSES ***

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
