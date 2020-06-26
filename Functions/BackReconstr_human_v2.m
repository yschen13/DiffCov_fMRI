function [newV,paras_back] = BackReconstr_human_v2(cx,Vtag)
	%{
	INPUT: 
		cx: region-specific traces
		Vtag: version tag
    VERSION
        1.0: Y.C. 06/17/2019, restore from 'BackReconstr_HCP_v1.m'
	%}

	addpath('/home/yuchen/Documents/fMRI_Real/SupportFunctions/inexact_alm_rpca');
	addpath('/home/yuchen/Documents/fMRI_Real/SupportFunctions/inexact_alm_rpca/PROPACK');
	addpath('/home/yuchen/Documents/fMRI_Real/SupportFunctions/ivCOV_fMRI_paper');
	addpath('/home/yuchen/Documents/fMRI_Real/SupportFunctions');

	h = 1.16;

	%% laplace backward Balloon
	% parameters upgraded 03/25/2019
	eta = .5;
	kappa = .65;
	gamma = .41;
	tau = .98;
	alpha = .32;   
	Vo = .02;
	tau_s = 1.54; % inverse of kappa
	tau_f = 2.46; % inverse of gamma

	TE  = 0.030; % real human model: TE=30ms
	ep = 0.5;
	r0 = 167; % Specfic to 3T
	nu0 = 40.3/1.5*3; % Specific to 3T
	E0  = 0.8; 
	k1  = 4.3*nu0*E0*TE;
	k2  = ep*r0*E0*TE;
	k3  = 1 - ep;

	paras_back = struct('eta',eta,'kappa',kappa,'gamma',gamma,'tau',tau,...
		'alpha',alpha,'Vo',Vo,'tau_s',tau_s,'tau_f',tau_f,'TE',TE,'ep',ep,...
		'r0',r0,'nu0',nu0,'E0',E0,'k1',k1,'k2',k2,'k3',k3);

	z1 = ((k1+k2)*tau*(-eta) + 1*(k3-k2)*(-eta))*Vo; % no Vo in the paper appendix
	z0 = ((k1+k2)*(tau-1+alpha)/(alpha*tau)*(-eta) + (1/tau)*(k3-k2)*(-eta))*Vo; % no Vo in the paper appendix
	p4 = tau;
	p3 = 1+1/alpha+tau/tau_s;
	p2 = 1/(tau*alpha) + tau/tau_f+(1+1/alpha)/tau_s;
	p1 = (1+1/alpha)/tau_f+1/(tau*alpha*tau_s);
	p0 = 1/(tau*alpha*tau_f);


	%% newV: get rid of hemodynamic effects
	N = size(cx,2);
	%
	y = cx;
	% dy1 = diff(y);
	dy1 = -1/2*y(1:end-2,:) + 1/2*y(3:end,:);
	% dy1 = 1/280*y(1:end-8,:) - 4/105*y(2:end-7,:) + 1/5*y(3:end-6,:) - 4/5*y(4:end-5,:)  + 4/5*y(6:end-3,:) - 1/5*y(7:end-2,:) + 4/105*y(8:end-1,:) - 1/280*y(9:end,:);
	dy1 = dy1/h;
	% dy2 = diff(dy1);
	dy2 = (-1/2*dy1(1:end-2,:) + 1/2*dy1(3:end,:))*h;
	% dy2 = -1/560*y(1:end-8,:) - 8/315*y(2:end-7,:) - 1/5*y(3:end-6,:) + 8/5*y(4:end-5,:) - 205/72*y(5:end-4,:)  + 8/5*y(6:end-3,:) - 1/5*y(7:end-2,:) + 8/315*y(8:end-1,:) - 1/560*y(9:end,:);
	dy2 = dy2/(h*h);
	% dy3 = diff(dy2);
	dy3 = (-1/2*dy2(1:end-2,:) + 1/2*dy2(3:end,:))*h*h;
	% dy3 = -7/240*y(1:end-8,:) + 3/10*y(2:end-7,:) - 169/120*y(3:end-6,:) + 61/30*y(4:end-5,:)  - 61/30*y(6:end-3,:) + 169/120*y(7:end-2,:) - 3/10*y(8:end-1,:) + 7/240*y(9:end,:);
	dy3 = dy3/(h*h*h);
	% dy4 = diff(dy3);
	dy4 = (-1/2*dy3(1:end-2,:) + 1/2*dy3(3:end,:))*h*h*h;
	% dy4 = 7/240*y(1:end-8,:) - 2/5*y(2:end-7,:) + 169/60*y(3:end-6,:) - 122/15*y(4:end-5,:) + 91/8*y(5:end-4,:)  - 122/15*y(6:end-3,:) + 169/60*y(7:end-2,:) - 2/5*y(8:end-1,:) + 7/240*y(9:end,:);
	dy4 = dy4/(h*h*h*h);
	dy5 = (-1/2*dy4(1:end-2,:) + 1/2*dy4(3:end,:))*h*h*h*h;
	% dy5 = 1*y(2:end-7,:) - 6*y(3:end-6,:) + 15*y(4:end-5,:) - 20*y(5:end-4,:)  + 15*y(6:end-3,:) - 6*y(7:end-2,:) + 1*y(8:end-1,:);
	dy5 = dy5/(h*h*h*h*h);

	% dy3 = dy3(1:end-1,:);
	% dy2 = dy2(1:end-2,:);
	% dy1 = dy1(1:end-3,:);
	% y = y(5:end-4,:);
	dy1 = [mean(dy1);dy1 ; mean(dy1)];
	dy2 = [repmat(mean(dy2),2,1); dy2; repmat(mean(dy2),2,1)];
	dy3 = [repmat(mean(dy3),3,1); dy3; repmat(mean(dy3),3,1)];
	dy4 = [repmat(mean(dy4),4,1); dy4; repmat(mean(dy4),4,1)];
	dy5 = [repmat(mean(dy5),5,1); dy5; repmat(mean(dy5),5,1)];

	newV = p4*dy4+p3*dy3+p2*dy2+p1*dy1+p0*y;
end