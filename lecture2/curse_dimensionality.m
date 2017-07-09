clear all;
close all;
addpath(genpath('data'));
%get beta
load('Data_20150102.mat');
data = Data_20150102.data;
data_m = data(:,1:end-3);
to_rem = find(sum(isnan(data_m),1)==size(data_m,1));
data_m(:,to_rem) = [];
[r_ c_] = find(isnan(data_m(1:end-50,:)));
data_m_nonan = data_m(max(r_)+1:end,:);
[r_ c_] = find(isnan(data_m_nonan));
data_m_nonan = data_m_nonan(1:min(r_)-1,:); X = data_m_nonan;
%regress payrolls on everything else
lm = fitlm(data_m_nonan(:,2:end),data_m_nonan(:,1));
%in-sample fit
fits = [ones(size(data_m_nonan,1),1) data_m_nonan(:,2:end)]*lm.Coefficients.Estimate;

dates = datenum('2014-10-01'); dates = dates:-31:0;
dates = dates(1:length(data_m_nonan)); dates = flip(dates); 
dates = datenum(datestr(dates,'yyyy-mm'));
cut_off = (dates(end));
Figure = figure;
hold on
p1 = plot(dates, fits, 'r-.','LineWidth',2);

%new data set
load('Data_20170106.mat');
data = Data_20170106.data;
data_m = data(:,1:end-3);
data_m(:,to_rem) = [];
[r_ c_] = find(isnan(data_m(1:end-50,:)));
data_m_nonan = data_m(max(r_)+1:end,:);
[r_ c_] = find(isnan(data_m_nonan));
data_m_nonan = data_m_nonan(1:min(r_)-1,:); XX = data_m_nonan;
dates = datenum('2016-10-01'); dates = dates:-31:0;
dates = dates(1:length(data_m_nonan)); dates = flip(dates); 
dates = datenum(datestr(dates,'yyyy-mm'));
%out-of-sample fit
fits = [ones(size(data_m_nonan,1),1) data_m_nonan(:,2:end)]*lm.Coefficients.Estimate;
p2 = plot(dates(end-24:end), fits(end-24:end), 'g-.','LineWidth',2);
p3 = plot(dates, data_m_nonan(:,1),'b','LineWidth',2);
yl_ = ylim();
line([cut_off cut_off], [yl_(1) yl_(2)]); 
ylabel('Change (in Ks) from previous period');
xlabel('date')
datetick('x','yyyy-mmm','keeplimits');
legend([p3 p1 p2],{'Actual Nonfarm Payrolls','In-sample fit','Out-of-sample fit'});
cd('plots');
print('cursefit','-dpdf');
cd ..


%% PCA
[T n] = size(X);
X_tilda = (X-repmat(mean(X,1),T,1))./repmat(std(X,1),T,1);
Omega = cov(X_tilda);
[evec eval] = eig(Omega); evec = fliplr(evec); eval = fliplr(eval);
FV = evec(:,1:4);%FV are the loadings
PCA_X = FV'*X_tilda'; 
%regressing our series on the principal components
lm_pca = fitlm(PCA_X',X_tilda(:,1));

%getting principal components for newest data set
[T n] = size(XX);
XX_tilda = (XX-repmat(mean(XX,1),T,1))./repmat(std(XX,1),T,1);
PCA_XX = FV'*XX_tilda'; %new principal components

%fitting to principal components to the new dataset
fit_pca = [ones(size(PCA_XX',1),1) PCA_XX']*lm_pca.Coefficients.Estimate;
%transforming fit into non-standardized version
fit_pca = fit_pca.*std(X(:,1)) + mean(X(:,1));

p4 = plot(dates(end-24:end), fit_pca(end-24:end),'y-.','LineWidth',2);
legend([p3 p1 p2 p4],{'Actual Nonfarm Payrolls','In-sample fit',...
    'Out-of-sample fit','Out-of-sample fit PCA'});
hold off
cd('plots');
print('cursefitPCA','-dpdf');
cd ..