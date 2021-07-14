Currentdir = pwd;
load('E:\mTBI_rat_processed_data\TRN_behav.mat');
lab = {'mNSS','Beam Walking (s)','NOR (RI)','MWM path length (m)','MWM escape latency (s)'};
fn = 'TCcoh_diff';
param = {'mNSS','BW','NOR','MWM_path','MWM_time'};
time = {'1-week','2-week','3-week','5-week','7-week'};
group = [ones(8,1);2*ones(8,1);3*ones(8,1)];
for jj = 1:length(param)
for ii = 5:5%1:size(mNSS,2)
figure(ii)
x = eval(fn);
temp = eval(param{jj});
y = temp(:,ii);
alpha = 0.05;
x_min = min(x);
x_max = max(x);
n_pts = 10000;
N=length(x);
XX = [ones(size(x)) x];
beta = XX\y;
X = x_min:(x_max-x_min)/n_pts:x_max;
Y = ones(size(X))*beta(1) + beta(2)*X;
SE_y_cond_x = sum((y - beta(1)*ones(size(y))-beta(2)*x).^2)/(N-2);
SSX = (N-1)*var(x);
SE_Y = SE_y_cond_x*(ones(size(X))*(1/N + (mean(x)^2)/SSX) + (X.^2 - 2*mean(x)*X)/SSX);
Yoff = (2*finv(1-alpha,2,N-2)*SE_Y).^0.5;
% SE_b0 = SE_y_cond_x*sum(x.^2)/(N*SSX)
% sqrt(SE_b0)
top_int = Y + Yoff;
bot_int = Y - Yoff;
% subplot(1,5,ii);
h = scatterhist(x,y,'Group',group,'Color',[0.4 0.4 1;1 0.6 0;0.7 0 0],'Marker','...','MarkerSIze',[30 30 30]);hold on; %,'Kernel','on','Direction','out','Color','br');hold on;
% clr = get(h(1),'Color');
[r,p]=corrcoef(x,y);
plot(X,Y,'k','LineWidth',2);
plot(X,top_int,'green--','LineWidth',1.5);
plot(X,bot_int,'green--','LineWidth',1.5);
boxplot(h(2),x,group,'orientation','horizontal','label',{'sham','single','double'},'Color',[0.4 0.4 1;1 0.6 0;0.7 0 0]);
boxplot(h(3),y,group,'orientation','horizontal','label',{'sham','single','double'},'Color',[0.4 0.4 1;1 0.6 0;0.7 0 0]);
set(h(2:3),'XTickLabel','');
view(h(3),[270,90]);
legend('sham control','single CHI','double CHI');

% title(sprintf('%s: y = (%0.3f) x + (%0.3f), r = %0.3f, p = %0.3f',fn,beta(2),beta(1),r(1,2),p(1,2)),'Fontsize',14);
title(sprintf('R^2 = %0.3f, p = %0.3f',(r(1,2))^2,p(1,2)),'Fontsize',14);
switch jj
    case 1
        ylim([0 10]);
    case 2
        ylim([0 180]);
    case 3
        ylim([0 1]);
    case 4
        ylim([0 10]);
    case 5
        ylim([0 30]);   
end
ylabel([lab{jj} ' : ' time{ii}]);
xlabel('RD(TRN) : 1-week');
% xlabel('FC between Thalamus <---> TPN (Z)');
grid on;
axis(h(1),'square');
% axis(h(2:3),'auto');
cd('E:\mTBI_rat_processed_data\TC_behav_plot');
saveas(gcf,sprintf('Scatterhist_%sX%s%s.png',fn,param{jj},time{ii}));
end
end
cd(Currentdir);
