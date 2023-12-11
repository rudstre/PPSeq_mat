function sampled_amplitude = rand_posterior(spike_count, seq_event_amplitude)
    % This function samples from the posterior distribution of a gamma
    % distributed parameter, given observed spike count data.
    
    % Extract the alpha (shape) and beta (rate) parameters from the prior
    alpha_prior = seq_event_amplitude.alpha;
    beta_prior = seq_event_amplitude.beta;

    % Calculate the posterior parameters
    alpha_post = alpha_prior + spike_count;
    beta_post = beta_prior + 1;

    % Sample from the posterior gamma distribution
    sampled_amplitude = gamrnd(alpha_post, 1/beta_post);
end