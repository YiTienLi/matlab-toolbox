close all;
Currentdir = pwd;
% cd('D:\mTBI_REST');
fn = 'WM2';
load('D:\mTBI_REST\mTBI_NST.mat');
load(sprintf('D:\\mTBI_REST\\%s_NST_result_WMseed_conn_blreg_all.mat',fn));
lab = {'PSQI','DHI','PCSQ','DS','CAL','WMI','BAI','BDI'};

for ii = 1:8
figure(ii)
x = t_dmn;
y = NST(:,ii);
alpha = 0.05/length(x);
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
h = scatterhist(x,y,'Group',group,'Color','br');hold on; %,'Kernel','on','Direction','out','Color','br');hold on;
% clr = get(h(1),'Color');
[r,p]=corrcoef(x,y);
plot(X,Y,'k','LineWidth',2);
plot(X,top_int,'green--','LineWidth',1.5);
plot(X,bot_int,'green--','LineWidth',1.5);
boxplot(h(2),x,group,'orientation','horizontal','label',{'',''},'color','br');
boxplot(h(3),y,group,'orientation','horizontal','label',{'',''},'color','br');
set(h(2:3),'XTickLabel','');
view(h(3),[270,90]);
legend('HC','mTBI');

title(sprintf('%s: y = %0.3f x + %0.3f, r = %0.3f, p = %0.3f',fn,beta(2),beta(1),r(1,2),p(1,2)),'Fontsize',14);
% ylim([0 10]);
ylabel(lab{ii});
xlabel('FC between Thalamus <---> DMN (Z)');
% xlabel('FC between Thalamus <---> TPN (Z)');
grid on;
axis(h(1),'square');
% axis(h(2:3),'auto');
%cd('H:\我的雲端硬碟\mTBI_Fig');
%saveas(gcf,sprintf('%s_T_DMNx%s.png',fn,lab{ii}));
end
cd(Currentdir);
