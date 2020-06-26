% To simulate surrogate dataset for significance test

close all; clear
% addpath('/home/yuchen/Documents/fMRI_Real/mouse');    

% load in dataset and AR model
dic = '/home/yuchen/Documents/fMRI_Real/human/subject1/';
cx = load([dic 'subject1_cx.mat']);
cx = cx.cx;
load([dic 'ARmdl_list.mat'])
fileID = fopen('subject1/IC_Annotation.txt'); % Load annotation file
Name = textscan(fileID, '%d%q%d','delimiter','-','HeaderLines',1);
fclose(fileID);
cx_select = cx(:,Name{1});

% % backward reconstruction and FC estimation
% [back,paras] = BackReconstr_human_v2(cx_select,'');
% [cov,Pcov,ivcov,Preg,A,A2] = mri_model_human(back,dic,1.16,'');

% Bootstrap 1000 times and save
A_boot = zeros(size(cx_select,2),size(cx_select,2),1000);
cov_boot = A_boot;
Pcov_boot = A_boot;
ivcov_boot = A_boot;
Preg_boot = A_boot;
A2_boot = A_boot;
for j = 1:1000
	tic
	cx_PSDshuffle = cx_select;
	for i = 1:size(cx_select,2)
		u = idinput(size(cx_select,1),'rgs');
		cx_PSDshuffle(:,i) = sim(mld_list{i},u);
	end
    [back,paras]=BackReconstr_human_v2(cx_PSDshuffle,'_PSDShuffle');
    [cov_boot(:,:,j),Pcov_boot(:,:,j),ivcov_boot(:,:,j),Preg_boot(:,:,j),A_boot(:,:,j),A2_boot(:,:,j)] = mri_model_human(back,dic,1.16,'_PSDShuffle');
    j
    toc
end
save([dic 'BootConnectivity_PSDShuffle.mat'],'cov_boot','Pcov_boot','ivcov_boot','Preg_boot','A_boot','A2_boot');
