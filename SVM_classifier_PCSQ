clear; close all;
Currentdir = pwd;
load('D:\mTBI_REST\PCSQ_1Y_2Y_SHHadd_new3.mat');
time = '2y';
cats = categorical(X_name);

X = eval(['X_' time]);
Y = cell(size(eval(['y_' time])));
Ytemp = eval(['y_' time]);
for ii = 1:length(Ytemp)
    if Ytemp(ii,1)<3
        Y{ii,1} = 'good';
    else
        Y{ii,1} = 'bad';
    end
end
fold = length(Ytemp); %leave-one-out
Mdl = fitcsvm(X,Y,'KFold',fold,'Standardize',true);
beta = [];
for ii = 1:fold
    beta = [beta Mdl.Trained{ii}.Beta];
%     temp = Mdl.ModelParameters.Generator.UseObsForIter(:,ii);
%     testset(ii) = find(temp==0);
%     Ynew(ii,1) = X(find(temp==0),:)*Mdl.Trained{ii}.Beta;
end
[fit,posterior] = kfoldPredict(Mdl);
[fpr,tpr,~,auc] = perfcurve(Y,posterior(:,1),Mdl.ClassNames(1));
cvtrainError = kfoldLoss(Mdl);
cvtrainAccuracy = 1-cvtrainError

figure;
plot([0 1],[0 1],'k--');hold on;
plot(fpr,tpr,'Linewidth',2); 
xlabel('False positive rate');
ylabel('True positive rate');
title(sprintf('ROC Curve, AUC = %0.4f, Accuracy = %0.4f',auc,cvtrainAccuracy));axis square;
if auc>cri
saveas(gcf,['ROC_' num2str(fold) 'fold_' time sprintf('_%d_%d',round(100*auc),round(100*cvtrainAccuracy)) '.png']);
end
% mse = kfoldLoss(Mdl)
figure;imagesc(beta);colormap(jet);
% figure;scatter(Y,fit);
% figure;scatter(0.5*(Y+fit),Y-fit);
% corrcoef(Y,fit)
% std(Y-fit)
m = mean(beta,2);
s = std(beta')'/2;
[~, I] = sort(m, 'ascend');
figure;hold on;
barh(linspace(1,length(m),length(m))',m(I),...
    'FaceColor',[0.6350 0.0780 0.1840],...
    'EdgeColor','r',...
    'LineWidth',1.5,...
    'BaseValue',0);
errorbar(m(I),linspace(1,length(m),length(m))',s(I),s(I),'k.','horizontal');
set(gca,'ytick',linspace(1,length(m),length(m))');
set(gca,'yticklabel',cats(I));
set(gca,'YDir','reverse');
if auc>cri
saveas(gcf,['Cmap_' num2str(fold) 'fold_' time sprintf('_%d_%d',round(100*auc),round(100*cvtrainAccuracy)) '.png']);
end
figure;
cm = confusionchart(confusionmat(Mdl.Y,fit,'Order',{'bad','good'}),{'bad outcome','good outcome'},'Fontsize',14);
if auc>cri
saveas(gcf,['CM_' num2str(fold) 'fold_' time sprintf('_%d_%d',round(100*auc),round(100*cvtrainAccuracy)) '.png']);
end
% figure;
% p(1) = sum(NPT_ind<=3);
% p(2) = sum(NPT_ind>3);
% p(3) = sum(ind{1}<65);
% p(4) = sum(ind{1}>=65);
% p(5) = sum(ind{2}<65);
% p(6) = sum(ind{2}>=65);
% plabel = {sprintf('demographic\n%d, [%0.2f%%]',p(1),100*p(1)/sum(p)),sprintf('neuropsychological\n%d, [%0.2f%%]',p(2),100*p(2)/sum(p)),sprintf('WM 1b activation\n%d, [%0.2f%%]',p(3),100*p(3)/sum(p)),sprintf('WM 1b deactivation\n%d, [%0.2f%%]',p(4),100*p(4)/sum(p)),sprintf('WM 2b activation\n%d, [%0.2f%%]',p(5),100*p(5)/sum(p)),sprintf('WM 2b deactivation\n%d, [%0.2f%%]',p(6),100*p(6)/sum(p))};
% pc = pie(p,[1 1 1 1 1 1],plabel);

% saveas(gcf,['pie_' num2str(fold) 'fold_' time sprintf('_%d_%d',round(100*auc),round(100*cvtrainAccuracy)) '.png']);
save(['Mdl_' num2str(fold) 'fold_' time sprintf('_%d_%d',round(100*auc),round(100*cvtrainAccuracy)) '.mat'],'Mdl');

