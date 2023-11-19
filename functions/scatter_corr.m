function md = scatter_corr(x,y, color, marker)
%SCATTER_CORR Plot correlation with scatter plot of the samples
%   
    scatter(x , y, 30, marker, "filled", "MarkerFaceColor",color);
    md = fitlm(x, y);
    p_val = md.Coefficients.pValue(2);
    coefs = md.Coefficients.Estimate;
     
    num_sample_plot = 1000;
    x_plot = linspace(min(x),max(x), num_sample_plot)';
    y_plot = [ones(num_sample_plot,1), x_plot]*coefs;
    plot(x_plot, y_plot, 'Color', color, 'LineWidth',2);
    

    % Bootstrap for ci
    num_bootstrap = 5000;
    Beta = zeros(num_bootstrap, 2);
    for i=1:num_bootstrap
        % Sample data
        bs_idx = datasample(1:length(x), length(x));
        % Fit line
        x_bs = [ones(length(x),1), x(bs_idx)];
        y_bs = y(bs_idx);
        beta = pinv(x_bs'*x_bs)*x_bs'*y_bs;
        Beta(i, :) = beta;
    end 

    y_pred_bs = [ones(length(x_plot),1),x_plot]*Beta';

    lb = prctile(y_pred_bs, 2.5, 2);
    hb = prctile(y_pred_bs,97.5, 2);
    plot(x_plot,  lb, 'Color', color);
    plot(x_plot,  hb,  'Color', color);
    fill([x_plot; flip(x_plot)]', [lb; flip(hb)]',color, 'FaceAlpha',0.1,'EdgeAlpha',0);

    text(x_plot(num_sample_plot/2),y_plot(num_sample_plot/2),['\leftarrow',' r = ',num2str(corr2(x,y)), ' p = ',num2str(p_val)])



end

