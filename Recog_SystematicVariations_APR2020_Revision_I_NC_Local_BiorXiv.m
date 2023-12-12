%%

%% SPATIOTEMPORAL BRAIN HIERARCHIES OF AUDITORY MEMORY RECOGNITION AND PREDICTIVE CODING - SYSTEMATIC VARIATION OF AUDITORY SEQUENCES (AUDITORY PATTERN RECOGNITION 2020) - LEONARDO BONETTI

%%

%% LOCAL COMPUTER - REVISION I - NATURE COMMUNICATIONS

%These codes were used locally to prepare plots with higher resolution compared to the solutions offered by the Aarhus cluster of computers

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

%%

%DO THAT WITH THE VOXELS

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
