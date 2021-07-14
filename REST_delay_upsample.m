function REST_delay_upsample(data_dir,etco2_dir)
% Input Format 1: SPM8_Preprocess_Full(SMSFileDir)
% Input Format 2: SPM8_Preprocess_Full(SMSFileDir,TR,SliceNumbers)
% Default Setting: TR=0.1, SliceNumbers=24
%%
% global data_stc_lh cvr_lh LHFile
global temp latency
clear cvr_lh cvr_rh latency_lh latency_rh

CurrentDir = pwd;
cd(data_dir)
LHFile = cellstr(ls('4D*nii'));

for subj = 1:1%length(LHFile)
tic;
% cd(CurrentDir)
data_stc_lh=strcat(data_dir,'\',LHFile{subj});

nii = MRIread(data_stc_lh);
data_lh = reshape(nii.vol,size(nii.vol,1)*size(nii.vol,2)*size(nii.vol,3),size(nii.vol,4));
[vox t]=size(data_lh);

mask_da=MRIread('\\10.103.1.160\Study_VCI\fBrainMask_05_61x73x61.nii');
mask=mask_da.vol;
mask_r=reshape(mask,73*61*61,1);

load(fullfile(etco2_dir,'rest_reg'));
%% Load Subj Motion Log
% MotionFile = cellstr(ls('rp*.txt'));
% fclose('all');
% fid = fopen(MotionFile{1});
% File = textscan(fid,'%s');
% fclose('all');
% m = reshape(File{1},6,t);
% motionlog = cellfun(@str2double,m);

% data_lh=interpft(data_lh,size(nii.vol,4)*4,2);


% TR=0.1;
% timeVec = ((1:2+1)-1).*TR;
t=size(nii.vol,4)*4;
%%
%Generate block regressor
% fprintf('CVR analysis -- Generate Block Reg!!\n')
ind = find(mask_r>0);
data_mean=bandpass(mean(data_lh(ind,:)),[0.01 0.08],2);
data_lh=interpft(data_lh,size(nii.vol,4)*4,2);
data_mean=interpft(data_mean,size(nii.vol,4)*4,2);
% data_mean=rest_reg'; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%etCO2 regressor

% [r lags]=crosscorr(block,data_mean,round(t/20));
% [~,I]=max(r);
% lag=lags(I)
% % block_reg=zeros(1,t);
% if lag>0
%     block_reg=[zeros(1,abs(lag)) block];
% %     block1d_reg=[zeros(1,abs(lag)) block1d];
% %     block_reg(1,:)=block0_reg(1,1:t);
% %     block_reg(2,:)=block1d_reg(1,1:t);
% %     block_reg(3,:)=ones(1,t);
% else
%     block_reg=block(abs(lag)+1:end);
%     if length(block_reg)<720
%         block_reg(end:720)=0;
%     end
% %     block1d_reg=block1d(lag+1:end);
% %     block_reg(1,:)=block0_reg(1,1:t);
% %     block_reg(2,:)=block1d_reg(1,1:t);
% %     block_reg(3,:)=ones(1,t);
% end
%%
%GLM CVR calculation -lh
% fprintf('CVR analysis --GLM cross-reg!!\n')
%cvr_lh = mask_r';
latency_lh = mask_r';
% h = waitbar(0,'CVR calculation');
% for ii =1:size(data_lh,1)
%     waitbar(ii/size(data_lh,1));
%     if mask_r(ii)>0
% %     fprintf('Processing lh %02d%%...\r',round(100*ii/size(data_lh,1)));
%     temp=data_lh(ii,:);
%     [cor latencies]=crosscorr(block_reg,temp,round(t/40));
%     [~,II]=max(cor);
%     latency=latencies(II);
%     block_reg1=zeros(8,t);
% %     block_reg2=zeros(1,t);
% %     block_reg1=zeros(2,t);
%     if latency>0
%         block_reg2=[zeros(1,abs(latency)) block_reg];
%         block_reg2=block_reg2(1,1:t);
%         block_reg1(1,1:t)=block_reg2(1,1:t);
%     else
%         block_reg2=block_reg(abs(latency)+1:end);
%         if length(block_reg2)>t
%         block_reg1(1,1:t)=block_reg2(1,1:t);
%         else
%         block_reg1(1,1:length(block_reg2))=block_reg2(1,1:end);
%         block_reg1(1,length(block_reg2)+1:end)=0;
%         end
%     end
%     
%     block_reg1(2,1:t)=ones(1,t);
%     block_reg1(3:8,:)=motionlog;
%     
%     beta_lh=inv(block_reg1*block_reg1')*block_reg1*temp';
%     cvr_lh(1,ii)=100*beta_lh(1)./beta_lh(2).*max(block_reg1(1,:));
% %     cvr_lh(2,ii)=100*beta_lh(1)./beta_lh(2).*max(block_reg1(1,:));
% %     latency_lh(1,ii)=latency;
% %     latency_lh(2,ii)=latency;
%     else
%         cvr_lh(1,ii)=nan;
% %         latency_lh(1,ii)=nan;
%     end
% end
% close(h);


h = waitbar(0,'Latency calculation');
for ii =1:size(data_lh,1)
    waitbar(ii/size(data_lh,1));
    if mask_r(ii)>0
%     fprintf('Processing lh %02d%%...\r',round(100*ii/size(data_lh,1)));
    temp=data_lh(ii,:);
    [cor latencies]=crosscorr(data_mean,temp,round(t/40));
    [~,II]=max(cor);
    latency=latencies(II);
%     block_reg1=zeros(8,t);
%     block_reg2=zeros(1,t);
    block_reg1=zeros(2,t);
    if latency>0
        block_reg2=[median(data_mean).*ones(1,abs(latency)) data_mean];
        block_reg2=block_reg2(1,1:t);
        block_reg1(1,1:t)=block_reg2(1,1:t);
    else
        block_reg2=data_mean(abs(latency)+1:end);
        if length(block_reg2)>t
        block_reg1(1,1:t)=block_reg2(1,1:t);
        else
        block_reg1(1,1:length(block_reg2))=block_reg2(1,1:end);
        block_reg1(1,length(block_reg2)+1:end)=median(data_mean);
        end
    end
    
    block_reg1(2,1:t)=ones(1,t);
%     block_reg1(3:8,:)=motionlog;
    
    beta_lh=inv(block_reg1*block_reg1')*block_reg1*temp';
    cvr_lh(1,ii)=100*beta_lh(1)./beta_lh(2);
%     cvr_lh(2,ii)=100*beta_lh(1)./beta_lh(2).*max(block_reg1(1,:));
    latency_lh(1,ii)=latency;
%     latency_lh(2,ii)=latency;
    else
        cvr_lh(1,ii)=nan;
        latency_lh(1,ii)=nan;
    end
end
close(h);


% beta_lh=inv(block_reg*block_reg')*block_reg*data_lh';
% H_lh=sqrt(((beta_lh(1,:).^2)*sum(block_reg(1,:).^2))+((beta_lh(2,:).^2)*sum(block_reg(2,:).^2))).*sign(beta_lh(1,:));
% 
% 
% cvr_lh(1,:)=100*H_lh./beta_lh(3,:).*max(block_reg(1,:));
% cvr_lh(2,:)=100*H_lh./beta_lh(3,:).*max(block_reg(1,:));
% latency_lh(1,:)=beta_lh(2,:)./beta_lh(1,:);
% latency_lh(2,:)=beta_lh(2,:)./beta_lh(1,:);


%%
%GLM CVR calculation -rh
% for ii =1:size(data_rh,1)
%     fprintf('Processing rh %02d%%...\r',round(100*ii/size(data_rh,1)));
%     temp=data_rh(ii,:);
%     [cor latencies]=crosscorr(block_reg,temp,round(t/60));
%     [~,II]=max(cor);
%     latency=latencies(II);
% %     block_reg1=zeros(8,t);
%     block_reg1=zeros(2,t);
%     if latency>0
%         block_reg2=[zeros(1,abs(latency)) block_reg];
%     else
%         block_reg2=block_reg(abs(latency)+1:end);
%     end
%     block_reg1(1,:)=block_reg2(1,1:t);
%     block_reg1(2,:)=ones(1,t);
% %     block_reg1(3:8,:)=motionlog;
%     
%     beta_rh=inv(block_reg1*block_reg1')*block_reg1*temp';
%     cvr_rh(1,ii)=100*beta_rh(1)./beta_rh(2).*max(block_reg1(1,:));
%     cvr_rh(2,ii)=100*beta_rh(1)./beta_rh(2).*max(block_reg1(1,:));
%     latency_rh(1,ii)=latency;
%     latency_rh(2,ii)=latency;
% end
% % beta_rh=inv(block_reg*block_reg')*block_reg*data_rh';
% % H_rh=sqrt(((beta_rh(1,:).^2)*sum(block_reg(1,:).^2))+((beta_rh(2,:).^2)*sum(block_reg(2,:).^2))).*sign(beta_rh(1,:));
% % 
% % 
% % cvr_rh(1,:)=100*H_rh./beta_rh(3,:).*max(block_reg(1,:));
% % cvr_rh(2,:)=100*H_rh./beta_rh(3,:).*max(block_reg(1,:));
% % latency_rh(1,:)=beta_rh(2,:)./beta_rh(1,:);
% % latency_rh(2,:)=beta_rh(2,:)./beta_rh(1,:);

%%
%write stc
%   timeVec=((1:size(data_epi,2))-1).*dt;
%   cd(data_dir)
%   filename=['ruCVR_' LHFile{1}];
%   inverse_write_stc(cvr_lh',Vv_lh,timeVec(1).*1e3,mean(diff(timeVec)).*1e3,filename);
%   filename=['uLatency_' LHFile{1}];
%   inverse_write_stc(latency_lh',Vv_lh,timeVec(1).*1e3,mean(diff(timeVec)).*1e3,filename);
%   
%   filename=['ruCVR_' RHFile{1}];
%   inverse_write_stc(cvr_rh',Vv_rh,timeVec(1).*1e3,mean(diff(timeVec)).*1e3,filename);
%   filename=['uLatency_' RHFile{1}];
%   inverse_write_stc(latency_rh',Vv_rh,timeVec(1).*1e3,mean(diff(timeVec)).*1e3,filename);
%   cd(CurrentDir)
%   filename=['CVR_' LHFile{1}];
%   mask_da.vol = reshape(cvr_lh,size(nii.vol,1),size(nii.vol,2),size(nii.vol,3));
%   err = MRIwrite(mask_da,filename,'double');
%   inverse_write_stc(cvr_lh',Vv_lh,timeVec(1).*1e3,mean(diff(timeVec)).*1e3,filename);
  filename=['Latency_' LHFile{subj}];
  mask_da.vol = reshape(latency_lh,size(nii.vol,1),size(nii.vol,2),size(nii.vol,3));
  err = MRIwrite(mask_da,filename,'double');
 
  filename=['CVR_' LHFile{subj}];
  mask_da.vol = reshape(cvr_lh,size(nii.vol,1),size(nii.vol,2),size(nii.vol,3));
  err = MRIwrite(mask_da,filename,'double');
  fprintf('Finish %s_%s, Elapsed time = %0.2f mins\n',LHFile{subj}(1:4),LHFile{subj}(end-5:end-4),toc/60);
%   inverse_write_stc(latency_lh',Vv_lh,timeVec(1).*1e3,mean(diff(timeVec)).*1e3,filename);
end
cd(CurrentDir);
end
