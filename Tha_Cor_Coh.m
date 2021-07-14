function Tha_Cor_Coh(varargin)
%Input: 4D sw* Raw Data Path, Template Path
%edited by Yi-Tien Li 2020/6/30
%
% global Tha_Sig BA_Sig data data_d data_c temp
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
    label = MRIread('E:\Scripts\NewAtlas.nii');
else
    mask = MRIread(sprintf('%s\fmask.nii',varargin{2}));
    csf = MRIread(sprintf('%s\wpriors_CSF.nii',varargin{2}));
    label = MRIread(sprintf('%s\NewAtlas.nii',varargin{2}));
end

mask_r = reshape(mask.vol,size(mask.vol,1)*size(mask.vol,2)*size(mask.vol,3),1);
csf_r = reshape(csf.vol,size(csf.vol,1)*size(csf.vol,2)*size(csf.vol,3),1);
label_r = reshape(label.vol,size(label.vol,1)*size(label.vol,2)*size(label.vol,3),1);

%%%%%%%加條件 少運算
%% detrend & band-pass filter
fprintf('Processing fMRI data : Detrend + Band-pass Filter ...\n');
data_d = detrend(epi_r,1);
epi.vol = reshape(data_d,size(epi.vol,1),size(epi.vol,2),size(epi.vol,3),size(epi.vol,4));
fn_new = ['d_' fn{1}];
MRIwrite(epi,fn_new,'float');

% data = bandpass(data_d',BP_range,1/TR)';
% epi.vol = reshape(data,size(epi.vol,1),size(epi.vol,2),size(epi.vol,3),size(epi.vol,4));
% fn_new = ['b' fn_new];
% MRIwrite(epi,fn_new,'float');

%% regress out nuisance covariates
fprintf('Processing fMRI data : Regress out nuisance covariates ...\n');
cov = [motion6 median(data_d(find(csf_r>0.2),:))'];
% cov = motion6;
beta = inv(cov'*cov)*cov'*data_d';
data_c = data_d - (cov*beta)';
epi.vol = reshape(data_c,size(epi.vol,1),size(epi.vol,2),size(epi.vol,3),size(epi.vol,4));
fn_new = ['p' fn_new];
MRIwrite(epi,fn_new,'float');

%% thalamocortical coherence
Tha_Sig(1,:) = nanmedian(data_c(find(label_r>=35),:));
% sk_list = [4 13 52 61];
h2 = waitbar(0,'Please wait...');
for ii = 1:34 %BA max(max(max(label_r)))
%     if ismember(ii,sk_list)
%     else
        temp = nanmedian(data_c(label_r==ii,:),1);
        BA_Sig(ii,:) = temp;
        [cxy, f] =mscohere(Tha_Sig,temp,[],[],[],1);
        TCcxy(ii,:) = cxy;
        TCf(ii,:) = f;
        waitbar(ii/102,h2,sprintf('Cal for No.%d, %0.2f%%\n',ii,100*ii/102));
%     end
end
close(h2);
save zTCcoh.mat TCcxy TCf Tha_Sig BA_Sig;
cd(CurrentDir);
end
