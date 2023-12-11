function [model, assignments, assignment_hist, log_p_hist, latent_event_hist, globals_hist] = annealed_gibbs(model, spikes, initial_assignments, num_anneals, samples_per_anneal, max_temperature, extra_split_merge_moves, split_merge_window, save_every, verbose)
    if nargin < 10
        verbose = false;
    end

    % Initialize storage
    assignment_hist = zeros(length(spikes), 0, 'int64');
    log_p_hist = [];
    latent_event_hist = {};
    globals_hist = {};

    % Return early if no samples
    if num_anneals == 0
        assignments = initial_assignments;
        return;
    end

    % Final amplitude for anneal
    model_priors = model.priors; 
    target_mean = mean(model_priors.seq_event_amplitude);
    target_var = var(model_priors.seq_event_amplitude);

    % Begin annealing
    temperatures = 10 .^ linspace(log10(max_temperature), 0, num_anneals);
    assignments = initial_assignments;

    for temp = temperatures
        % Print progress
        if verbose
            fprintf('TEMP: %f\n', temp);
        end

        % Anneal prior on sequence amplitude
        prior = get_priors(model);
        model = model.anneal_priors(specify_gamma(target_mean, target_var * temp));

        % Recompute probability of introducing a new cluster
        max_time = get_max_time(model);
        alpha = prior.seq_event_amplitude.alpha;
        beta = prior.seq_event_amplitude.beta;
        lambda = prior.seq_event_rate;
        log_prob = log(alpha) + log(lambda) + log(max_time) + alpha*(log(beta) - log(1+beta));
        model = set_new_cluster_log_prob(model, log_prob);

        % Draw Gibbs samples
        [model, assignments, assgn_hist, logp, latents, globals] = ...
            gibbs_sample(model, spikes, assignments, samples_per_anneal, extra_split_merge_moves, split_merge_window, save_every, verbose); % Assuming gibbs_sample is a custom function

        % Save samples
        assignment_hist = [assignment_hist, assgn_hist];
        log_p_hist = [log_p_hist; logp];
        latent_event_hist = [latent_event_hist; latents];
        globals_hist = [globals_hist; globals];
    end
end
