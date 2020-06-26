% MAIN SCRIPT comparing SC and FC

close all
clear
% addpath('/home/yuchen/Documents/fMRI_Real/SupportFunctions');
% addpath('/home/yuchen/Documents/fMRI_Real/SupportFunctions/Violinplot-Matlab-master/');

%% load in time series data
dic = '/home/yuchen/Documents/fMRI_Real/human/subject1/';
cx = load([dic 'subject1_cx.mat']);
cx = cx.cx;

% loading GT matrix
GT = load('GT/human_GT_Hag.mat'); % Hagmann GT

% load in component selection file
fileID = fopen([dic 'IC_Annotation.txt']); % Load annotation file
Name = textscan(fileID, '%d%q%d','delimiter','-','HeaderLines',1); % 1st col: non-artefact ICs; 2nd col: non-artefact IC names; 3rd col: DMN or not
fclose(fileID);

name = {'cov';'Pcov';'ivcov';'Preg';'Preg+SL'};
Name_list = {'cov','Pcov','$\Delta c$','$\Delta p$','$\Delta s$'};


%% Compare bootstrap results and select the significant regions
thres = 0.001;

% Backward reconstruction + DiffCov 
[back,paras] = BackReconstr_human_v2(cx,'');
[cov,Pcov,ivcov,Preg,A,A2] = mri_model_human(back,dic,1.16,'');

% PSD bootstrap (fitting AR models separately)
load([dic 'BootConnectivity_PSDShuffle.mat'])
cov_boot_select = cov_boot;cov_select = cov(Name{1},Name{1});
Pcov_boot_select = Pcov_boot;Pcov_select = Pcov(Name{1},Name{1});
ivcov_boot_select = ivcov_boot; ivcov_select = ivcov(Name{1},Name{1});
Preg_boot_select = Preg_boot; Preg_select = Preg(Name{1},Name{1});
A_boot_select = A_boot; A_select = A(Name{1},Name{1});
A2_boot_select = A2_boot; A2_select = A2(Name{1},Name{1});
clear cov_boot Pcov_boot ivcov_boot Preg_boot A_boot
N = size(cov_select,1);

[cov_p,cov_f,~,~,cov_GT] = boot_stat(cov_boot_select,cov_select,thres,'cov',GT,1);savefig([dic 'Plots/CovBoot.fig']);
[Pcov_p,Pcov_f,~,~,Pcov_GT] = boot_stat(Pcov_boot_select,Pcov_select,thres,'Pcov',GT,1);savefig([dic 'Plots/PCovBoot.fig']);
[ivcov_p,ivcov_f,~,~,ivcov_GT] = boot_stat(ivcov_boot_select,ivcov_select,thres,'ivcov',GT,1);savefig(ivcov_f,[dic 'Plots/ivcovBoot.fig']);
[Preg_p,Preg_f,~,~,Preg_GT] = boot_stat(Preg_boot_select,Preg_select,thres,'Preg',GT,1);savefig(Preg_f,[dic 'Plots/PregBoot.fig']);
[A_p,A_f,~,~,A_GT] = boot_stat(A_boot_select,A_select,thres,'Preg+SL',GT,1);savefig(A_f,[dic 'Plots/ABoot.fig']);
[A2_p,~,~,~,A2_GT] = boot_stat(A2_boot_select,A2_select,thres,'A2',GT,1); 



% Create a table of thres and GTmean
thres_list = [0.05 0.01 0.005 0.001 0.0005 0.0001];
GTmean_list = zeros(6,length(thres_list)); % each row stands for cov, Pcov, ivcov, Preg, A and A2
AbsOrimean_list = zeros(6, length(thres_list));
mean_p_list = []; % each column is different across thresholds, each row is the significance of mean value difference between cov_GT and *_GT (ivcov_GT/Preg_GT/A_GT/A2_GT)
N_Sig_list = []; % sample size of significant connections  
for j = 1:length(thres_list)
    close all
    p = thres_list(j);
    [~,~,GTmean_list(1,j),AbsOrimean_list(1,j),cov_GT] = boot_stat(cov_boot_select,cov_select,p,'cov',GT,0);
    [~,~,GTmean_list(2,j),AbsOrimean_list(2,j),Pcov_GT] = boot_stat(Pcov_boot_select,Pcov_select,p,'Pcov',GT,0);
    [~,~,GTmean_list(3,j),AbsOrimean_list(3,j),ivcov_GT] = boot_stat(ivcov_boot_select,ivcov_select,p,'ivcov',GT,0);
    [~,~,GTmean_list(4,j),AbsOrimean_list(4,j),Preg_GT] = boot_stat(Preg_boot_select,Preg_select,p,'Preg',GT,0);
    [~,~,GTmean_list(5,j),AbsOrimean_list(5,j),A_GT] = boot_stat(A_boot_select,A_select,p,'Preg+SL',GT,0);
    [~,~,GTmean_list(6,j),AbsOrimean_list(6,j),A2_GT] = boot_stat(A2_boot_select,A2_select,p,'A2',GT,0);
    [mean_p_list(1,j),~] = ranksum(cov_GT,ivcov_GT);
    [mean_p_list(2,j),~] = ranksum(cov_GT,Preg_GT);
    [mean_p_list(3,j),~] = ranksum(cov_GT,A_GT);
    [mean_p_list(4,j),~] = ranksum(cov_GT,A2_GT);
    N_Sig_list(1,j) = length(cov_GT); N_Sig_list(2,j) = length(Pcov_GT); N_Sig_list(3,j) = length(ivcov_GT);
    N_Sig_list(4,j) = length(Preg_GT); N_Sig_list(5,j) = length(A_GT); N_Sig_list(6,j) = length(A2_GT);
    j
end
csvwrite([dic 'PSDbootstat_GTmean_Hag_all.csv'], GTmean_list)
csvwrite([dic 'PSDbootstat_GTmean_Hag_all_p.csv'], mean_p_list)
csvwrite([dic 'PSDbootstat_GTmean_Hag_all_N.csv'], N_Sig_list)




