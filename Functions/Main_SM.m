%{ 
The Main Script for the relationship between FC and SM
%}

close all; clear

cd('~/Documents/fMRI_Real/Figures')
addpath('/home/yuchen/Documents/MATLAB/cbrewer/cbrewer');
addpath('2019_03_03_BCT') % for network properties
addpath('powerlaws_full_v0.0.10-2012-01-17') % for power law dist function fitting

HCP_dic = '~/Documents/fMRI_Real/HCP/HCP_PTN1200/Analysis/';

Title_list = {'cov','Pcov','\Deltac','\Deltap','\Deltas'};
Name_list = {'cov','Pcov','$\Delta c$','$\Delta p$','$\Delta s$'};
Tag_list = {'cov','Pcov','ivcov','Preg','A'};
thres_list = [0.05 0.01 0.005 0.001 0.0005 0.0001];


%% HCP data loading

thres = 0.05;
ConvMatrix = load([HCP_dic 'ConvMatrix_IC100.mat']);
Mat_list = {ConvMatrix.cov_list, ConvMatrix.Pcov_list,ConvMatrix.ivCOV_list,ConvMatrix.Preg_list,ConvMatrix.A_list};
clear ConvMatrix
load([HCP_dic 'Timeseries_IC100.mat'],'SID_list')
load([HCP_dic 'PSDAllBootStat_IC100/p_list.mat'])
ID_list = cell2mat(ID_list);
SID_list = cellfun(@str2double,SID_list);
Tag = 'Version1/HCP';

%% Figure 4ab: Compute the network efficiency of the binary network 
Effi = []; % efficiency of the original matrix
Bidir_effi = []; % efficiency of reciprocal connections
Unidir_effi = []; % efficiency of uni-directional connections
Bidir_effi2 = []; % efficiency of upper triangle connections
Bidir_effi3 = []; % efficiency of lower triangle connections
Anti_Bidir_effi = []; % efficiency of network with significant anti-symmetric connections
for kk = 1:1003
    for i = 1:5 % across methods
        tmp = Mat_list{i}{kk,5};
        bi_tmp = double(p_list{kk,i,5}<thres);
        sig_tmp = tmp .* bi_tmp; sig_tmp(sig_tmp>0)=300; sig_tmp(sig_tmp<0)=-100;
        Anti_sym = double((sig_tmp + transpose(sig_tmp))==200);
        Bidir = bi_tmp.*transpose(bi_tmp);
        Unidir = bi_tmp - Bidir;
        tmp_triu = triu(bi_tmp); tmp_tril = tril(bi_tmp);
        Bidir_Net2 = tmp_triu + transpose(tmp_triu) - diag(diag(bi_tmp));
        Bidir_Net3 = tmp_tril + transpose(tmp_tril) - diag(diag(bi_tmp));
        [~, Effi(kk,i)] = charpath(distance_bin(bi_tmp));
        [~, Bidir_effi(kk,i)] = charpath(distance_bin(Bidir));
        [~, Unidir_effi(kk,i)] = charpath(distance_bin(Unidir));
        [~, Bidir_effi2(kk,i)] = charpath(distance_bin(Bidir_Net2));
        [~, Bidir_effi3(kk,i)] = charpath(distance_bin(Bidir_Net3));
        [~, Anti_Bidir_effi(kk,i)] = charpath(distance_bin(Anti_sym));
    end
end

% Efficiency Pcov vs Efficiency of Preg+SL
m = 2; n = 5;
f = figure;
set(gcf,'Position',[0 0 300 300])
[tmp, tmp_p] = corrcoef([Effi(:,m) Effi(:,n)]);
scatter(Effi(:,m),Effi(:,n),10,'filled')
xlabel(['Global efficiency of ' Title_list{m}]); ylabel(['Global efficiency of ' Title_list{n}])
title(['r=' num2str(tmp(1,2),2) '; p=' num2str(tmp_p(1,2),2)])
set(gca,'FontSize',12,'XColor','k','YColor','k','LineWidth',1.5);
saveas(f,[Tag '_Efficiency_BetweenCorr.eps'],'epsc')

f=figure;
set(gcf,'Position',[0 0 400 400],'Renderer','painters');
[tmp, tmp_p] = corrcoef([Effi(:,m) Effi(:,n)]);
scatterhist(Effi(:,m),Effi(:,n),'Direction','out','MarkerSize',10,'Marker','.');
h = lsline; h.Color = 'r'; h.LineWidth= 3; h.LineStyle = '--';
xlabel(['Global efficiency of ' Title_list{m}]); ylabel(['Global efficiency of ' Title_list{n}]);
title(['r=' num2str(tmp(1,2),2) '; p=' num2str(tmp_p(1,2),2)]);
set(gca,'FontSize',12,'XColor','k','YColor','k','LineWidth',1.5,'color','none')
box off
saveas(f,[Tag '_Efficiency_BetweenCorr_scatterhist.eps'],'epsc')

%% Read in SM spreadsheets available from HCP website
fileID = fopen('Analysis/behavioral_data.csv');
psycho = textscan(fileID, ['%d' repmat('%s',1,114) repmat('%f',1,397) repmat('%s',1,60) repmat('%f',1,10) '%*[^\n]'],'delimiter',',','Headerlines',1);
fclose(fileID);
fileID = fopen('Analysis/behavioral_data.csv');
psycho_name = textscan(fileID,[repmat('%s',1,582) '%*[^\n]'],1,'delimiter',',');
fclose(fileID);
% Select from ’PicSeq_Unadj‘ to 'SelfEff_Unadj'
select_idx = [1 116:512 573:582];
psycho = psycho(select_idx); 
psycho_name = psycho_name(select_idx);
% restricted data
fileID = fopen('RESTRICTED_yuchen_6_25_2019_23_47_20_tag.csv');
restricted_name = textscan(fileID,[repmat('%s',1,201) '%*[^\n]'],1,'delimiter',',');
fclose(fileID);
fileID = fopen('RESTRICTED_yuchen_6_25_2019_23_47_20_tag.csv');
restricted_tag = textscan(fileID,[repmat('%s',1,201) '%*[^\n]'],1,'delimiter',',','Headerlines',1);
fclose(fileID);
tmp = cellfun(@char,restricted_tag);
readin_pattern = strrep(tmp,'1','%f');
readin_pattern = strrep(readin_pattern,'0','%s');
fileID = fopen('RESTRICTED_yuchen_6_25_2019_23_47_20_tag.csv');
restricted = textscan(fileID,[readin_pattern '%*[^\n]'],'delimiter',',','Headerlines',2);
fclose(fileID);
tag_idx = [];% select the numeric columns
for i = 1:length(restricted_tag)
    tag_idx(i) = str2num(restricted_tag{i}{1});
end
tag_idx = find(tag_idx~=0);

psycho_sort = nan([size(SID_list,2),size(psycho,2)-1]); % exclude the SID column
psycho_name = psycho_name(2:end); % exclude the 'Subject'
for k = 1:size(psycho{1},1)
    SID = psycho{1}(k);
    ID = find(SID_list==SID);
    if ~isempty(ID)
        idx = find(ID_list==ID);
        for i = 1:size(psycho_sort,2)
            psycho_sort(idx,i) = psycho{i+1}(k);
        end
    end
end
% exclude FS measurements
noFS_idx = [];
psycho_name_noFS = {};
for i = 1:length(psycho_name)
    if strfind(char(psycho_name{i}),'FS')
        continue
    else
        noFS_idx = [noFS_idx i];
        psycho_name_noFS = [psycho_name_noFS cell2mat(psycho_name{i})];
    end
end
psycho_sort_noFS = psycho_sort(:,noFS_idx);
FS_idx = setdiff(1:size(psycho_sort,2), noFS_idx);
FS_corr = corrcoef(psycho_sort(:,FS_idx),'Rows','pairwise');


% deal with restricted dataset 
restricted_num = restricted(tag_idx);
restricted_name_num = restricted_name(tag_idx);
restricted_sort = nan([size(SID_list,2),size(restricted_num,2)-1]); % exclude the SID column
restricted_name_num = restricted_name_num(2:end);
tmp = [];
for i = 1:length(restricted_name_num)
    tmp = [tmp restricted_name_num{i}];
end
restricted_name_num = tmp;
for k = 1:size(restricted_num{1},1) % across subjects
    SID = restricted_num{1}(k);
    ID = find(SID_list==SID);
    if ~isempty(ID)
        idx = find(ID_list==ID);
        for i = 1:size(restricted_sort,2)
            restricted_sort(idx,i) = restricted_num{i+1}(k);
        end
    end
end

% Set 999 in Correction column to nan
for i = 1:length(restricted_name_num)
    if strfind(char(restricted_name_num{i}),'Correction')
        correction_idx = i;
    end
end
tmp = restricted_sort(:,correction_idx);tmp(tmp==999) = nan;
restricted_sort(:,correction_idx) = tmp;

all_sort = [psycho_sort restricted_sort];
all_name = [psycho_name restricted_name_num];
for i = 1:length(all_name)
    if strfind(char(all_name{i}),'Age_in')
        age_idx = i;
    end
    if strfind(char(all_name{i}),'Height')
        height_idx = i;
    end
    if strfind(char(all_name{i}),'FS_TotCort_GM_Vol')
        TotGM_idx = i;
    end
    if strfind(char(all_name{i}),'Strength')
        strength_idx = i;
    end
    if strfind(char(all_name{i}),'PMAT24_A_RTCR')
        fi_idx = i;
    end
    if strfind(char(all_name{i}),'BPSystolic')
        BPsys_idx = i;
    end
end
% correlation of the variables to be regressed out
corrcoef(all_sort(:,[age_idx height_idx strength_idx fi_idx TotGM_idx BPsys_idx]),'Row','pairwise')

combine_sort = [psycho_sort_noFS restricted_sort]; % sorted according to the position in ID_list
combine_name = [psycho_name_noFS restricted_name_num];
nan_N = sum(isnan(combine_sort));
small_nan_idx = find(nan_N<800);
combine_sort_select = combine_sort(:,small_nan_idx);
combine_name_select = combine_name(small_nan_idx);

% regress out age, height, TotalGMVolume, BPsys, from efficiency
X = [ones(size(Effi,1),1) all_sort(:,[age_idx height_idx TotGM_idx BPsys_idx])];
for j = 1:size(X,2)
    tmp = X(:,j); tmp(isnan(tmp)) = mean(tmp,'omitnan');
    X(:,j) = tmp;
end
Effi_res = Effi;
for i = 1:size(Effi,2)
    y = Effi(:,i); 
    b = X\y;
    Effi_res(:,i) = y - X*b;
end
[psycho_variant_corr, psycho_variant_corr_pvalue] = corrcoef([Effi Effi_res Bidir_effi Unidir_effi Bidir_effi2 Bidir_effi combine_sort_select],'Rows','pairwise');


% sig_psycho = and(psycho_variant_corr_pvalue(2,26:end)<0.05,psycho_variant_corr_pvalue(5,26:end)<0.05);
% sig_psycho = psycho_variant_corr_pvalue(14,26:end)<0.05;
sig_idx = find(psycho_variant_corr_pvalue(10,31:end)<0.005); % sig based on residual A-Net efficiency
sig_name = {};
for j = 1:length(sig_idx)
    tmp_name = combine_name_select{sig_idx(j)};
    sig_name{j} = strrep(tmp_name,'_',' ');
end

% Heatmap of correlation and p value with SM names on the xticklabels
f = figure;
set(gcf,'Position',[0 0 250 400]);
subplot(211);imagesc(psycho_variant_corr(6:10,sig_idx+30));colorbar;title('Correlation')
xticks(1:length(sig_idx)); xticklabels(sig_name); xtickangle(60)
a = get(gca,'XTickLabel'); set(gca,'XTickLabel',a,'fontsize',5)
yticks(1:5); yticklabels(Title_list);
ylabel('FC-Net efficiency')
subplot(212);imagesc(-log(psycho_variant_corr_pvalue(6:10,sig_idx+30)));colorbar; title('p value (-log10)')
xticks(1:length(sig_idx)); xticklabels(sig_name); xtickangle(60)
a = get(gca,'XTickLabel'); set(gca,'XTickLabel',a,'fontsize',5)
yticks(1:5); yticklabels(Title_list);
ylabel('FC-Net efficiency')
saveas(f,[Tag '_Effi_SM_corr_res.eps'],'epsc')

% Heatmap w/o xticklabels
f = figure;
set(gcf,'Position',[0 0 300 300]);
subplot(211);imagesc(psycho_variant_corr(6:10,sig_idx+30));colorbar;title('Correlation')
xlabel('Subject measurements');ylabel('FC-Net efficiency')
yticks(1:5); yticklabels(Title_list);
set(gca,'xticklabel',{[]});
subplot(212);imagesc(-log(psycho_variant_corr_pvalue(6:10,sig_idx+30)));colorbar; title('p value (-log10)')
yticks(1:5); yticklabels(Title_list);
xlabel('Subject measurements');ylabel('FC-Net efficiency')
set(gca,'xticklabel',{[]});
saveas(f,[Tag '_Effi_SM_corr_res_nolabel.eps'],'epsc')

% Scatterplots of SM v.s efficiency
k = sig_idx(1); 
f = figure;
set(gcf,'Position',[0 0 900 600])
for i = 1:5
    subplot(2,3,i)
    set(gca,'FontSize',12,'XColor','k','YColor','k','LineWidth',1.5)
    hold on
    scatter(Effi(:,i),combine_sort_select(:,k),10,'filled');
    h = lsline; h.Color = 'r'; h.LineWidth= 3; h.LineStyle = '--';
    hold off
    [tmp_corr,tmp_pvalue] = corrcoef([Effi(:,i) combine_sort_select(:,k)],'Rows','pairwise');
    title(Title_list{i})
    xlimit = xlim;
    text(0.5*(xlimit(1)+xlimit(2)),max(combine_sort_select(:,k))*0.9,['r=' num2str(tmp_corr(1,2),2) '; p=' num2str(tmp_pvalue(1,2),2)])
end
saveas(f,[Tag 'Effi_' combine_name_select{k} '_corr.eps'], 'epsc')



%% Figure S4: other network measurements

% degree distribution
degree_list = {}; % row: subjects; column: methods
for kk = 1:size(p_list,1) % ID_idx in ID_list
    for i = 1:5
        sigbi = double(p_list{kk,i,5}<thres);
        [~,~,tmp2] = degrees_dir(sigbi);
        degree_list{kk,i} = tmp2;    
    end
end

kk = 1;
% the survivor form of cdf (log-log)
g = figure;
set(gcf,'Position',[0 0 900 600])
for i = 1:5
    x = degree_list{kk,i}; x = x';
    fd_exp = fitdist(x,'Exponential');
    fd_gauss = fitdist(x,'Normal');
    [alpha,xmin,L] = plfit(x);
    n = length(x); q = unique(x);
    % get the survivor ecdf
    c = hist(x,q)'./n;
    c = [[q;q(end)+1] 1-[0; cumsum(c)]]; c(c(:,2)<10^-10,:) = [];
    cf = ((xmin:q(end))'.^-alpha)./(zeta(alpha) - sum((1:xmin-1).^-alpha));
    cf = [(xmin:q(end)+1)' 1-[0; cumsum(cf)]];
    cf(:,2) = cf(:,2) .* c(c(:,1)==xmin,2);
    exp_cdf = [q 1-cdf(fd_exp,q)];
    gauss_cdf = [q 1-cdf(fd_gauss,q)];
    
    subplot(2,3,i)
    [x_ecdf,x_uni] = ecdf(x);
    loglog(c(:,1),c(:,2),'bo','MarkerSize',8,'MarkerFaceColor',[1 1 1]); hold on;
    loglog(cf(:,1),cf(:,2),'k--','LineWidth',2); hold on
    loglog(exp_cdf(:,1),exp_cdf(:,2),'r--','LineWidth',2); hold on
    loglog(gauss_cdf(:,1),gauss_cdf(:,2),'g--','LineWidth',2); hold on
    ymax = ylim;plot([xmin xmin],[ymax(1) ymax(2)],'k--','Linewidth',1);
    ylabel('P(X \geq x)'); xlabel('x: node degree'); title(Title_list{i});
    legend({'1-cdf','Power Law','Exponential','Gaussian'},'Location','southwest');
    set(gca,'FontSize',12)
end
saveas(g,[Tag '_ID' num2str(kk) '_Degree_fit_log.eps'],'epsc')


% DTI degree distribution
load('~/Documents/fMRI_Real/human/GT/DSI_release2_2011.mat')
Cij_bi = double(CIJ_fbden_average~=0);
DTI_degree = degrees_und(Cij_bi);
[alpha,xmin,L] = plfit(DTI_degree);
fd_exp = fitdist(DTI_degree','Exponential');
fd_gauss = fitdist(DTI_degree','Normal');
x = DTI_degree(DTI_degree~=0)';
n = length(x);        
q = unique(x);
c = hist(x,q)'./n;
c = [[q;q(end)+1] 1-[0; cumsum(c)]]; c(c(:,2)<10^-10,:) = [];
cf = ((xmin:q(end))'.^-alpha)./(zeta(alpha) - sum((1:xmin-1).^-alpha));
cf = [(xmin:q(end)+1)' 1-[0; cumsum(cf)]];
cf(:,2) = cf(:,2) .* c(c(:,1)==xmin,2);

fig = figure;
set(gcf,'Position',[0 0 400 400]);
loglog(c(:,1),c(:,2),'bo','MarkerSize',8,'MarkerFaceColor',[1 1 1]); hold on;
loglog(cf(:,1),cf(:,2),'k--','LineWidth',2); hold on
loglog(q,1-cdf(fd_exp,q),'r--','LineWidth',2); hold on
loglog(q,1-cdf(fd_gauss,q),'g--','LineWidth',2); hold on
ymax = ylim;plot([xmin xmin],[ymax(1) ymax(2)],'k--','Linewidth',1);
legend({'1-cdf','Power Law','Exponential','Gaussian'},'Location','southwest');
ylabel('Pr(X \geq x)'); xlabel('x: node degreee');
set(gca,'FontSize',12);
saveas(fig,[Tag '_DTI_degree_dist_log.eps'],'epsc')


% other network measurements
cc = []; % clustering coefficient
trans = []; % transitivity
loc = []; % local efficiency
modul = []; % modularity
for kk = 1:1003
    for i = 1:5
        tmp = double(p_list{kk,i,5}<thres);
        cc(kk,i) = mean(clustering_coef_bd(tmp));
        trans(kk,i) = transitivity_bd(tmp);
        loc(kk,i) = efficiency_bin(tmp);
        [~,modul(kk,i)] = modularity_dir(tmp);
    end
end

cc(cc==Inf) = nan;
measure_corr = corrcoef([cc trans loc modul],'Rows','pairwise');

f = figure;
set(gcf,'Position',[0 0 900 650])
m = 2; n = 5; 
subplot(2,3,1); 
[corr_tmp,p_tmp] = corrcoef([Effi(:,m) Effi(:,n)],'Rows','pairwise');
scatter(Effi(:,m),Effi(:,n),10,'filled'); hold on;
xlimit = xlim;ylimit=ylim;text(0.47*(xlimit(1)+xlimit(2)),ylimit(2),['r=' num2str(corr_tmp(1,2),2) '; p=' num2str(p_tmp(1,2),2)]);
title('Global Efficiency'); xlabel(Title_list{m}); ylabel(Title_list{n});set(gca,'FontSize',12);
subplot(2,3,2); 
[corr_tmp,p_tmp] = corrcoef([cc(:,m) cc(:,n)],'Rows','pairwise');
scatter(cc(:,m),cc(:,n),10,'filled'); hold on;
xlimit = xlim;ylimit=ylim;text(0.45*(xlimit(1)+xlimit(2)),ylimit(2),['r=' num2str(corr_tmp(1,2),2) '; p=' num2str(p_tmp(1,2),2)]);
title('Clustering Coefficient'); xlabel(Title_list{m}); ylabel(Title_list{n});set(gca,'FontSize',12);
subplot(2,3,3); 
[corr_tmp,p_tmp] = corrcoef([trans(:,m) trans(:,n)],'Rows','pairwise');
scatter(trans(:,m),trans(:,n),10,'filled'); hold on;
xlimit = xlim;ylimit=ylim;text(0.45*(xlimit(1)+xlimit(2)),ylimit(2),['r=' num2str(corr_tmp(1,2),2) '; p=' num2str(p_tmp(1,2),2)]);
title('Transitivity'); xlabel(Title_list{m}); ylabel(Title_list{n}); set(gca,'FontSize',12);
subplot(2,3,4); 
[corr_tmp,p_tmp] = corrcoef([loc(:,m) loc(:,n)],'Rows','pairwise');
scatter(loc(:,m),loc(:,n),10,'filled'); hold on;
xlimit = xlim;ylimit=ylim;text(0.47*(xlimit(1)+xlimit(2)),ylimit(2),['r=' num2str(corr_tmp(1,2),2) '; p=' num2str(p_tmp(1,2),2)]);
title('Local Efficiency'); xlabel(Title_list{m}); ylabel(Title_list{n});set(gca,'FontSize',12)
subplot(2,3,5); 
[corr_tmp,p_tmp] = corrcoef([modul(:,m) modul(:,n)],'Rows','pairwise');
scatter(modul(:,m),modul(:,n),10,'filled'); hold on;
xlimit = xlim;ylimit=ylim;text(0.4*(xlimit(1)+xlimit(2)),ylimit(2),['r=' num2str(corr_tmp(1,2),2) '; p=' num2str(p_tmp(1,2),2)]);
title('Modularity'); xlabel(Title_list{m}); ylabel(Title_list{n});set(gca,'FontSize',12);
saveas(f,[Tag '_NM_BetweenCorr.eps'], 'epsc')
