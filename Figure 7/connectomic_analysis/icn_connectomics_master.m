clear all, close all, clc
root = 'D:\Dropbox (Brain Modulation Lab)\Shared Lab Folders\CRCNS\MOVEMENT DATA';
leaddbs = 'C:\code\leaddbs';
wjn_toolbox = 'C:\code\wjn_toolbox';
spm12 = 'C:\code\spm12';
icn_connectomics = 'C:\code\icn_connectomics';
addpath(spm12)
addpath(genpath(leaddbs))
addpath(wjn_toolbox)
addpath(icn_connectomics)


%% CREATE SPHERICAL ROIS FROM ELECTRODE LOCATIONS FOR ALL SUBJECTS
root = 'C:\tmp\connectomics_ROIs\';
cd(root)
csvfile='C:\code\icn\ECOG_vs_STN\UPDATE_NEWCV_PLOTS\df_all.csv';
T=readtable(csvfile);
connectomics_dir = fullfile('C:\tmp','connectomics_ROIv2');
mkdir(connectomics_dir)

for a =1:size(T,1)
    n=0;
    rois={};
    sub = ['sub-' sprintf( '%03d', T.sub(a) )];
    roiname = fullfile(connectomics_dir,[sus '_' T.sess_{a} '_ROI_'  T.ch{a} '.nii']);
    if strcmp(T.ch{a}(1:4),'ECOG')
        n=n+1;
        mni = [-abs(T.x(a)) T.y(a) T.z(a)];
        rois{n}=wjn_spherical_roi(roiname,mni,4);
    elseif strcmp(T.sess_{a},'right') && ismember(T.ch{a},{'STN_RIGHT_0','STN_RIGHT_1','STN_RIGHT_2'})
        n=n+1;
        mni = nanmean([-abs(T.x(b)) T.y(b) T.z(b);T.x(b+1),T.y(b+1),T.z(b+1)]);
        rois{n}=wjn_spherical_roi(roiname,mni,4);
    else
        n=n+1;
        mni = nanmean([-abs(T.x(b)) T.y(b) T.z(b);T.x(b+1),T.y(b+1),T.z(b+1)]);
        rois{n}=wjn_spherical_roi(roiname,mni,3);
    end
end


%% SEARCH INDIVIDUAL FIBERS
clear all, close all, clc
root = 'C:\tmp\connectomics_ROIs\';
cd(root)
csvfile='C:\code\icn\ECOG_vs_STN\UPDATE_NEWCV_PLOTS\df_all.csv';
T=readtable(csvfile);
connectomics_dir = fullfile('C:\tmp','connectomics_fibers');
mkdir(connectomics_dir)
connectome = load(fullfile(ea_getearoot,'connectomes','dMRI','PPMI_85 (Ewert 2017)','data.mat'));
mkdir(connectomics_dir)
for a =1:size(T,1)
    side = T.sess_{a};
    sub = ['sub-' sprintf( '%03d', T.sub(a) )];
    fibername = fullfile(connectomics_dir,[sub '_' T.sess_{a} '_' T.ch{a} '_dMRI_fibers.mat']);
    if strcmp(T.ch{a}(1:4),'ECOG')
        mni = [-abs(T.x(a)) T.y(a) T.z(a)];
        wjn_searchfibers(fibername,mni,connectome,5,1);
    elseif strcmp(T.sess_{a},'right') && ismember(T.ch{a},{'STN_RIGHT_0','STN_RIGHT_1','STN_RIGHT_2'})
        mni = nanmean([-abs(T.x(a)) T.y(a) T.z(b);T.x(a+1),T.y(a+1),T.z(a+1)]);
        wjn_searchfibers(fibername,mni,connectome,5,1);
    else
        mni = nanmean([-abs(T.x(a)) T.y(a) T.z(a);T.x(a+1),T.y(a+1),T.z(a+1)]);
        wjn_searchfibers(fibername,mni,connectome,5,1);
    end
end



%% Decoding performance interpolated to Surface
root = 'D:\Dropbox (Brain Modulation Lab)\Shared Lab Folders\CRCNS\MOVEMENT DATA';
csvfile='C:\code\icn\ECOG_vs_STN\UPDATE_NEWCV_PLOTS\df_all.csv'
[files,folders] = wjn_subdir(fullfile(root,'sub*ECOG*Z*.nii'));

T=readtable(csvfile);
istn = ci('STN',T.ch);
iecog = ci('ECOG',T.ch);

stn = [ T.r2_con(istn) -abs(T.x(istn)) T.y(istn) T.z(istn)];
figure
wjn_plot_surface('STN.surf.gii')
hold on
wjn_plot_surface('STN.surf.gii',stn)
hold on
caxis([0 .6])
caire = [-12.58, -13.41, -5.87];
plot3(caire(1),caire(2),caire(3),'linestyle','none','marker','x','MarkerEdgeColor','r','Markersize',40)
alpha .5
myprint('STN_overlay05')
alpha 1
myprint('STN_overlay1')
%

ecog = [ T.r2_con(iecog) -abs(T.x(iecog)) T.y(iecog) T.z(iecog)];
hk = [-37 -25 62];
ctx = load('CortexLowRes_15000V.mat');
ctx.Faces = ctx.Faces_lh;
ctx.Vertices = ctx.Vertices_lh;

figure
wjn_plot_surface(ctx)
hold on
wjn_plot_surface(ctx,ecog)
% plot3(ecog(:,2),ecog(:,3),ecog(:,4),'rx')
hold on
plot3(hk(1),hk(2),hk(3),'linestyle','none','marker','x','MarkerEdgeColor','r','Markersize',40)
caxis([0 .6])
alpha .5
myprint('ECOG_overlay05')
alpha 1
myprint('ECOG_overlay1')

%% Write out spatial interpolation nifti
root = 'D:\Dropbox (Brain Modulation Lab)\Shared Lab Folders\CRCNS\MOVEMENT DATA';
csvfile='C:\code\icn\ECOG_vs_STN\UPDATE_NEWCV_PLOTS\df_all.csv';
[files,folders] = wjn_subdir(fullfile(root,'sub*ECOG*Z*.nii'));

T=readtable(csvfile);
istn = ci('STN',T.ch);
iecog = ci('ECOG',T.ch);

stn = [ T.r2_con(istn) -abs(T.x(istn)) T.y(istn) T.z(istn)];

ecog = [ T.r2_con(iecog) -abs(T.x(iecog)) T.y(iecog) T.z(iecog)];
hk = [-37 -25 62];
wjn_heatmap('STN_XGB_performance.nii',stn(:,2:end),stn(:,1),'C:\code\spm12\canonical\single_subj_T1.nii')
wjn_heatmap('ECOG_XGB_performance.nii',ecog(:,2:end),ecog(:,1),'C:\code\spm12\canonical\single_subj_T1.nii')
%% Spatial Fingerprint analysis single channel LOOM

clear all, close all
root = 'C:\tmp\connectomics_ROIs\';
resultsfolder='C:\tmp\connectomics_ROIs\Fingerprints';
mkdir(resultsfolder)
cd(root)
csvfile='C:\code\icn\ECOG_vs_STN\UPDATE_NEWCV_PLOTS\df_all.csv';
T=readtable(csvfile);
n=0;nm=0;missing={};
T=T(ci('ECOG',T.ch),:);

for a = 1:size(T,1)
%     try
        [fmri_file,fdir] = ffind(fullfile(root,'PPMI 74_15 (Horn 2017)_Patients',['sub-' sprintf( '%03d', T.sub(a) )  '*' T.sess_{a} '*' T.ch{a} '*AvgR_Fz.nii']));
        nii = ea_load_nii(fullfile(fdir{1},fmri_file{1}));
        n=n+1;
        fmrimat(:,a) = nii.img(:);
        [dti_file,fdir] = ffind(fullfile(root,'PPMI_85 (Ewert 2017)',['ssub-' sprintf( '%03d', T.sub(a) )  '*' T.sess_{a} '*' T.ch{a} '*struc_seed.nii']));
        nii = ea_load_nii(fullfile(fdir{1},dti_file{1}));
        dtimat(:,n) = nii.img(:);
        regressor(n,1) = T.r2_con(a);
%     catch
%         nm=nm+1;
%         missing{nm} = ['sub-' sprintf( '%03d', T.sub(a) )  '*' T.sess_{a} '*' T.ch{a} '*.nii'];
%         warning(['sub-' sprintf( '%03d', T.sub(a) )  '*' T.sess_{a} '*' T.ch{a} '*.nii is missing.'])
%     end
    disp(a)
end
dtimat(dtimat==0)=nan;
rcv=[];
for leftout = 1:length(T.sub)
    remaining = setdiff(1:size(T,1),leftout);
    % Step 1 create R-Map without leftout 
    rfmri=corr(fmrimat(:,remaining)',regressor(remaining),'rows','pairwise');  
    rdti=corr(dtimat(:,remaining)',regressor(remaining),'rows','pairwise');  
    % Step 2 calculate spatial similarity of  individual fingerprints with R-Map
    rrdti=[];
    rrfmri=[];
    rrdti=corr(dtimat(:,remaining),rdti,'rows','pairwise');
    rrfmri=corr(fmrimat(:,remaining),rfmri,'rows','pairwise');
    % Step 3 calculate fingerprint similarity between leftout connectivity profile with R-Map from remaining contacts
    rldti=corr(dtimat(:,leftout),rdti,'rows','pairwise');
    rlfmri=corr(fmrimat(:,leftout),rfmri,'rows','pairwise');
    % OR train linear model with more than one fingerprint (e.g. DTI+fMRI or different brain parcellations)
    mdl = fitlm([rrdti rrfmri],regressor(remaining));
    r2 = mdl.predict([rldti rlfmri]);
    rcv = [rcv;regressor(leftout) rldti rlfmri r2];
    disp(['Currently left out: ' num2str(leftout)])
end

% Step 4 write out R-Maps for visualization of optimal fingerprint
or=corr(dtimat',regressor,'rows','pairwise');  
Rmap = nii;
Rmap.fname = fullfile(resultsfolder,'DTI_Rmap_R2con_ECOG.nii');
Rmap.img(:) = or(:);
spm_write_vol(Rmap,Rmap.img)

CV.DTI_nifti = Rmap;


or=corr(fmrimat',regressor,'rows','pairwise');  
Rmap = nii;
Rmap.fname = fullfile(resultsfolder,'fMRI_Rmap_R2con_ECOG.nii');
Rmap.img(:) = or(:);
spm_write_vol(Rmap,Rmap.img)

CV.fMRI_nifti = Rmap;
CV.regressor = rcv(:,1);
CV.fMRI = rcv(:,3);
CV.DTI = rcv(:,2);
CV.LM = rcv(:,4);
save(fullfile(resultsfolder,'CV'));

% Step 5 Evaluate whether connectivity fingerprint is predictive for decoding performance

figure
wjn_corr_plot(rcv(:,3),rcv(:,1))
title('Prediction network mapping - Leave one channel out')
xlabel('Spatial correlation')
ylabel('Decoding performance')


%% Spatial Fingerprint analysis whole subject LOOM

clear all, close all
root = 'C:\tmp\connectomics_ROIs\PPMI 74_15 (Horn 2017)_Patients';
cd(root)
csvfile='C:\code\icn\ECOG_vs_STN\UPDATE_NEWCV_PLOTS\df_all.csv';
T=readtable(csvfile);
n=0;nm=0;missing={};
T=T(ci('ECOG',T.ch),:);

for a = 1:length(T.sub)
    try
        nii_file = ffind(['sub-' sprintf( '%03d', T.sub(a) )  '*' T.sess_{a} '*' T.ch{a} '*AvgR_Fz.nii']);
        nii = ea_load_nii(nii_file{1});
        n=n+1;
        niimat(:,n) = nii.img(:);
        regressor(n,1) = T.r2_con(a);
    catch
        nm=nm+1;
        missing{nm} = ['sub-' sprintf( '%03d', T.sub(a) )  '*' T.sess_{a} '*' T.ch{a} '*.nii'];
        warning(['sub-' sprintf( '%03d', T.sub(a) )  '*' T.sess_{a} '*' T.ch{a} '*.nii is missing.'])
    end
    disp(a)
end
rcv=[];
usubs = unique(T.sub);
figure
for a = 1:length(usubs)
    leftout = find(T.sub==usubs(a));
    remaining = setdiff(1:size(T,1),leftout);
    rr=corr(niimat(:,remaining)',regressor(remaining),'rows','pairwise');
    rrr=[];
    for b = 1:length(remaining)
        rrr(:,b)=corr(niimat(:,remaining(b)),rr,'rows','pairwise');
        disp(b)
    end
    mdl = fitlm(rrr,regressor(remaining));
    rl=corr(niimat(:,leftout),rr,'rows','pairwise');
    r2 = mdl.predict(rl);
    rcv = [rcv;regressor(leftout) rl r2];
    disp(['Currently left out: ' num2str(leftout')])
    [~,best(a,1)] = max(regressor(leftout));
    [~,best(a,2)] = max(rl);
    [~,best(a,3)] = max(r2);
end

figure
wjn_corr_plot(rcv(:,2),rcv(:,1))
title('Prediction network mapping - Leave one subject out')
xlabel('Spatial correlation')
ylabel('Decoding performance')


%% Spatial Fingerprint analysis session LOOM

clear all, close all
root = 'C:\tmp\connectomics_ROIs\PPMI 74_15 (Horn 2017)_Patients';
cd(root)
csvfile='C:\code\icn\ECOG_vs_STN\UPDATE_NEWCV_PLOTS\df_all.csv';
T=readtable(csvfile);
n=0;nm=0;missing={};
T=T(ci('ECOG',T.ch),:);

for a = 1:length(T.sub)
    try
        nii_file = ffind(['sub-' sprintf( '%03d', T.sub(a) )  '*' T.sess_{a} '*' T.ch{a} '*AvgR_Fz.nii']);
        nii = ea_load_nii(nii_file{1});
        n=n+1;
        niimat(:,n) = nii.img(:);
        regressor(n,1) = T.r2_con(a);
    catch
        nm=nm+1;
        missing{nm} = ['sub-' sprintf( '%03d', T.sub(a) )  '*' T.sess_{a} '*' T.ch{a} '*.nii'];
        warning(['sub-' sprintf( '%03d', T.sub(a) )  '*' T.sess_{a} '*' T.ch{a} '*.nii is missing.'])
    end
    disp(a)
end
rcv=[];
subsess = unique(strcat(num2str(T.sub),'_',T.sess_));
figure
for a = 1:length(subsess)
    sub = str2double(strtok(subsess{a},'_'));
    [~,rem] = strtok(subsess{a},'_');
    leftout = find(T.sub==sub & strcmp(rem(2:end),T.sess_));
    remaining = setdiff(1:size(T,1),leftout);
    rr=corr(niimat(:,remaining)',regressor(remaining),'rows','pairwise');
    rrr=[];
    for b = 1:length(remaining)
        rrr(:,b)=corr(niimat(:,remaining(b)),rr,'rows','pairwise');
        disp(b)
    end
    mdl = fitlm(rrr,regressor(remaining));
    rl=corr(niimat(:,leftout),rr,'rows','pairwise');
    r2 = mdl.predict(rl);
    rcv = [rcv;regressor(leftout) rl r2];
    disp(['Currently left out: ' num2str(leftout')])
    [~,best(a,1)] = max(regressor(leftout));
    [~,best(a,2)] = max(rl);
    [~,best(a,3)] = max(r2);
end

figure
wjn_corr_plot(rcv(:,2),rcv(:,1))
title('Prediction network mapping - Leave one subject out')
xlabel('Spatial correlation')
ylabel('Decoding performance')

%% Clinical correlation [needs correction to include fisher transform and ECOG vs. STN separately]
clear all, close all
root = 'D:\Dropbox (Brain Modulation Lab)\Shared Lab Folders\CRCNS\MOVEMENT DATA\AvgR_Fz';
cd(root)
csvfile='C:\code\icn\ECOG_vs_STN\UPDATE_NEWCV_PLOTS\df_all.csv';
T=readtable(csvfile);
subs = unique(T.sub);
for a = 1:length(subs),i=T.sub==subs(a);perf(a) = nanmean([T.r2_con(i);T.r2_ips(i)]);U(a)=nanmean(T.UPDRS_total(i));end
figure
wjn_corr_plot(U',perf')
title('Prediction network mapping - Leave one subject out')
xlabel('Motor sign severity [UPDRS-III total]')
ylabel('XGB Decoding performance [R²]')
myprint('UPDRS-III correlation')

%% SPM analysis

% create factorial matrix
clear all, close all
root = 'D:\Dropbox (Brain Modulation Lab)\Shared Lab Folders\CRCNS\MOVEMENT DATA\AvgR_Fz';
cd(root)
csvfile='C:\code\icn\ECOG_vs_STN\UPDATE_NEWCV_PLOTS\df_all.csv';
T=readtable(csvfile);
T=T(ci('ECOG',T.ch),:);

decoding_performance = T.r2_con;

wjn_gaussianize(decoding_performance)

sub = [T.sub; T.sub];
unique_subjects = unique(sub);
contacts=[];
for a = 1:length(unique_subjects)
    contacts=[contacts;[1:size(find(T.sub==unique_subjects(a)))]'];
end
double_contacts = [contacts;contacts];
cond = [ones(200,1);ones(200,1).*2];

factor_matrix = [ones(400,1) sub cond double_contacts]

