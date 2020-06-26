function [BOLD_cov,BOLD_Pcov,BOLD_ivCOV,BOLD_Preg,A] = mri_human_model(cx, h)
%{
    Calculate diffCov, Ireg, +SL (fixed lambda) and plotting
    Applicable dataset: 
        Real human data B=3T, TR=1.16s, 100 IC 
    Keep the diagonal values
    Add Pcov (precision matrix of the covariance matrix) (05/22/19)
    INPUT: 
        cx: backreconstructed signal: T * N 
        h: sampled time interval
    OUTPUT:
        BOLD_cov: covariance matrix
        BOLD_Pcov: Partial covariance matrix
        BOLD_ivCOV: differential covariance matrix
        BOLD_Preg: Partial differential covariance matrix
        A: Preg+Sparse Regularization
%}

    addpath('inexact_alm_rpca');
    addpath('inexact_alm_rpca/PROPACK');

    N = size(cx,2);

    %% Deconv BOLD cov
    BOLD = cx;
    diff_cx = -1/2*cx(1:end-2,:) + 1/2*cx(3:end,:); % (z(t+1)-z(t-1))/2
    diff_cx = diff_cx/h; % (z(t+1)-z(t-1))/2h
    diff_cx = [mean(diff_cx);diff_cx ; mean(diff_cx)]; % substitute the first and last row with mean
    dBOLD = diff_cx;

    BOLD_cov = corrcoef(BOLD); % covariance matrix
    BOLD_Pcov = inv(BOLD_cov); % partial covariance matrix
    BOLD_Csample = corrcoef([BOLD dBOLD]);
    BOLD_ivCOV = BOLD_Csample(N+1:N+N,1:N); % differential covariance matrix
    BOLD_Preg = regC2(BOLD_Csample); % partial differential covariance matrix
    nodiagBOLD_Preg = BOLD_Preg - diag(diag(BOLD_Preg));
    lambda = 1/sqrt(N);
    [~,A,~] = inexact_alm_rpca(nodiagBOLD_Preg,lambda); % Preg+Sparse Regularization

end
