function BOLD_func_conn_rat(varargin)
%Input: 4D sw* Raw Data Path, Template Path
%edited by Yi-Tien Li 2020/6/30
%
global Tha_Sig BA_Sig data data_d data_c temp
%% Setting
TR = 1;
BP_range = [0.01 0.1];
%% load data
CurrentDir = pwd;
cd(varargin{1});
fprintf('Loading fMRI data...\n');
fn = cellstr(ls('fMRI2mni.ni*'));
if strcmp(fn{1}(end-1:end),'gz')>0
        gunzip(fn{1});
        fn = cellstr(ls('fMRI2mni.nii'));
end
epi = MRIread(fn{1});
epi_r = reshape(epi.vol,size(epi.vol,1)*size(epi.vol,2)*size(epi.vol,3),size(epi.vol,4));
rp = cellstr(ls('rp*.txt'));
motion6 = load(rp{1});

if nargin <2
    mask = MRIread('E:\Scripts\TT2_Tohoku_mask.nii');
     csf = MRIread('H:\我的雲端硬碟\MTBI_rat_preprocess\rat_template\sigle_rat_template\wpriors_CSF.nii');
    label = MRIread('E:\Scripts\TT2_Tohoku_labels.nii');
else
    mask = MRIread(sprintf('%s\fmask.nii',varargin{2}));
    csf = MRIread(sprintf('%s\wpriors_CSF.nii',varargin{2}));
    label = MRIread(sprintf('%s\TT2_Tohoku_labels.nii',varargin{2}));
end

mask_r = reshape(mask.vol,size(mask.vol,1)*size(mask.vol,2)*size(mask.vol,3),1);
csf_r = reshape(csf.vol,size(csf.vol,1)*size(csf.vol,2)*size(csf.vol,3),1);
label_r = reshape(label.vol,size(label.vol,1)*size(label.vol,2)*size(label.vol,3),1);

%%%%%%%加條件 少運算
%% detrend & band-pass filter
filetest1 = cellstr(ls('d_*.nii'));
if isempty(filetest1{1})
fprintf('Processing fMRI data : Detrend + Band-pass Filter ...\n');
data_d = detrend(epi_r,1);
epi.vol = reshape(data_d,size(epi.vol,1),size(epi.vol,2),size(epi.vol,3),size(epi.vol,4));
fn_new = ['d_' fn{1}];
MRIwrite(epi,fn_new,'float');
else
temp_epi = MRIread(filetest1{1});
data_d = reshape(temp_epi.vol,size(epi.vol,1)*size(epi.vol,2)*size(epi.vol,3),size(epi.vol,4));
fn_new = ['d' fn{1}];
end

filetest2 = cellstr(ls('bd_*.nii'));
if isempty(filetest2{1})
data = bandpass(data_d',BP_range,1/TR)';
epi.vol = reshape(data,size(epi.vol,1),size(epi.vol,2),size(epi.vol,3),size(epi.vol,4));
fn_new = ['b' fn_new];
MRIwrite(epi,fn_new,'float');
else
temp_epi = MRIread(filetest2{1});
data = reshape(temp_epi.vol,size(epi.vol,1)*size(epi.vol,2)*size(epi.vol,3),size(epi.vol,4));
fn_new = ['b' fn_new];
end
%% regress out nuisance covariates
fprintf('Processing fMRI data : Regress out nuisance covariates ...\n');
cov = [motion6 mean(data_d(find(csf_r>0.2),:))'];
%cov = motion6;
beta = inv(cov'*cov)*cov'*data';
data_c = data - (cov*beta)';
epi.vol = reshape(data_c,size(epi.vol,1),size(epi.vol,2),size(epi.vol,3),size(epi.vol,4));
fn_new = ['p' fn_new];
MRIwrite(epi,fn_new,'float');

%% func conn
mask = MRIread('E:\Scripts\NewAtlas.nii');
mask_r = reshape(mask.vol,1,size(mask.vol,1)*size(mask.vol,2)*size(mask.vol,3));

for ii = 1:max(mask_r)
    ROI_Sig(ii,:) = nanmedian(data_c(mask_r==ii,:));   
end

CorrMat = corrcoef(ROI_Sig');
for kk = 1:size(CorrMat,1)
    CorrMat(kk,kk) = 0;
end
sCorrMat = squareform(CorrMat);

save zCorrMat_csf.mat ROI_Sig CorrMat sCorrMat
cd(CurrentDir);
end
