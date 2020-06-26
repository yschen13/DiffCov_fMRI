function [P,f,SigGTMean,SigAbsOriMean,Binary_GT] = boot_stat(boot,ori,thres,Tag,GT,PlotOn)
%{
    12/6/19: Transform p value to be smaller than 0.5
    INPUT:
        boot: bootstrap results, 30*30*1000
        ori: orignal connectivity matrix, 30*30
        thres: p value threshold
        Tag: which connectivity matrix
        GT: ground truth matrix
    OUTPUT:
        P: cdf of normal distribution with bootstrap-generated mean and std
        f: figure
        SigGTMean: the mean GT connectivity value from significant
        connections
        SigAbsMean: the mean absolute connection value from significant
        connections
        Binary_GT: the significant GT values
%}
    mu = mean(boot,3);
    sigma = std(boot,0,3);
    P = normcdf(ori,mu,sigma);
    P(P>0.5) = 1-P(P>0.5);
    Binary = double(P<thres);
%     Binary = double(Binary)-diag(diag(double(Binary)));
    
    Binary_GT = GT(Binary==1);
    SigGTMean = mean(Binary_GT,'omitnan');
    tmp = abs(ori);Binary_absori = tmp(Binary==1);
    SigAbsOriMean = mean(Binary_absori,'omitnan');
    
    if PlotOn
        f = figure;       
        addpath('/home/yuchen/Documents/MATLAB/cbrewer/cbrewer');
        customcolor = flipud(cbrewer('div','RdBu',20));
        subplot(221);
        heatmap(Binary);
        colormap(customcolor)
        title([Tag ' p<' num2str(thres)]);
        caxis([-1 1]);

        subplot(222);
        heatmap(ori.*Binary);
        colormap(customcolor);
        title([Tag ' Significant original connections']);
        caxis([-max(max(abs(ori.*Binary))) max(max(abs(ori.*Binary)))]);

        subplot(223);
        heatmap(abs(ori).*Binary);
        colormap(customcolor);
        title([Tag ' Absolute value of significant original connections']);
        caxis([-max(max(abs(ori).*Binary)) max(max(abs(ori).*Binary))]);

        subplot(224);
        heatmap(GT.*Binary);
        colormap(customcolor);
        title([Tag 'Significant Ground Truth Connectivity']);
        caxis([-max(max(abs(GT.*Binary))) max(max(abs(GT.*Binary)))]);
    else f = 0;
    end
end
