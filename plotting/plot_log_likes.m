function fig = plot_log_likes(config, results)
    % Plots log-likelihoods over MCMC samples
    
    [x1, x2] = get_mcmc_x_coords(config, results);

    fig = figure;
    hold on;

    if isfield(config, 'mask_lengths')
        y1 = results.anneal_test_log_p_hist;
        y2 = [y1(end); results.test_log_p_hist];
        y3 = results.anneal_train_log_p_hist;
        y4 = [y3(end); results.train_log_p_hist];
        
        plot(x1, y1, 'DisplayName', 'test anneal');
        plot(x2, y2, 'DisplayName', 'test after anneal');
        plot(x1, y3, 'DisplayName', 'train anneal');
        plot(x2, y4, 'DisplayName', 'train after anneal');
    else
        y1 = results.anneal_log_p_hist(:);
        y2 = [y1(end); results.log_p_hist];

        plot(x1, y1, 'DisplayName', 'anneal');
        plot(x2, y2, 'DisplayName', 'after anneal');
    end

    ylabel('log-likelihood');
    xlabel('MCMC samples');
    legend;
    hold off;
end
