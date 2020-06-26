function [data_deconv, event, HRF, adjust_global, PARA] = Deconv(data,h,SaveOn)
  % For blind deconvolution, Version 1: 02/20/19
  % Modified for real data, Version 2: 07/19/19
  % Input:
  %   data: hemodynamic traces 
  %   h: time interval (second)
  %   SaveOn: if save as *.mat file or not
  % Output:
  %   Deconv_*.mat

  addpath('/home/yuchen/Documents/Balloon/cxcx_56/spm12')
  addpath('/home/yuchen/Documents/Balloon/cxcx_56/rsHRF-master')

  N = size(data,2);
  L = size(data,1);

  data_norm = zscore(data);

  TR = h; %  second
  thr = 1; % (mean+) thr*standard deviation threshold to detect event.
  event_lag_max_seconds=9;    % the (estimated) maximum lagged time from neural event to BOLD event, in seconds.
  event_lag_max =round(event_lag_max_seconds/TR);

  [data_deconv,event,HRF,adjust_global,PARA] = wgr_deconv_canonhrf_par2(data_norm,thr,event_lag_max,TR);
 
  if SaveOn
    save([dic 'Deconv_' num2str(h) '.mat'],'data_norm','data_deconv','event','HRF','adjust_global','PARA')
  end
