function [x1, x2] = get_mcmc_x_coords(config, results)
    % Helper function to get x-axis coordinates for MCMC samples

    if isfield(config, 'mask_lengths')
        s1 = config.save_every_during_anneal;
        e1 = length(results.anneal_train_log_p_hist) * s1;
        x1 = s1:s1:e1;

        s2 = config.save_every_after_anneal;
        e2 = e1 + length(results.train_log_p_hist) * s2;
        x2 = e1:s2:e2;
    else
        s1 = config.save_every_during_anneal;
        e1 = length(results.anneal_log_p_hist) * s1;
        x1 = s1:s1:e1;

        s2 = config.save_every_after_anneal;
        e2 = e1 + length(results.log_p_hist) * s2;
        x2 = e1:s2:e2;
    end
end