function [model, assignments, assignment_hist, log_p_hist, latent_event_hist, globals_hist] = gibbs_sample(model, spikes, initial_assignments, num_samples, extra_split_merge_moves, split_merge_window, save_every, verbose)
    if nargin < 8
        verbose = false;
    end

    % Initialize spike assignments and recompute model
    assignments = initial_assignments;
    model = recompute(model, spikes, assignments);

    % Initialize storage
    n_saved_samples = round(num_samples / save_every);
    log_p_hist = zeros(1, n_saved_samples);
    assignment_hist = zeros(length(spikes), n_saved_samples, 'int64');
    latent_event_hist = cell(1, n_saved_samples);
    globals_hist = cell(1, n_saved_samples);
    spike_order = 1:length(spikes);

    %% Main loop for Gibbs sampling
    for s = 1:num_samples
        % Update spike assignments in random order
        spike_order = spike_order(randperm(length(spike_order)));
        for i = spike_order
            model = remove_datapoint(model, spikes(i), assignments(i));
            [assignments(i),model] = gibbs_add_datapoint(model, spikes(i));
        end

        % Extra split-merge moves
        [assignments,model] = split_merge_sample(model, spikes, extra_split_merge_moves, assignments, split_merge_window);

        % Update latent events and globals
        model = gibbs_update_latents(model);
        model = gibbs_update_globals(model, spikes, assignments);

        % Store results
        if mod(s, save_every) == 0
            j = s / save_every;
            log_p_hist(j) = log_like(model, spikes);
            assignment_hist(:, j) = assignments;
            latent_event_hist{j} = event_list_summary(model);
            globals_hist{j} = model.globals;
            if verbose
                fprintf('%d-', s);
            end
        end
    end

    if verbose && (n_saved_samples > 0)
        fprintf('Done\n');
    end
end

%%
function model = gibbs_update_latents(model)
    log_probs = model.RW_buffer; % Assuming RW_buffer is pre-allocated

    for it = 1:length(model.sequence_events.evnts) % Loop over sequence events
        event = model.sequence_events.evnts(it);

        if event.spike_count == 0
            continue
        end

        % Sample sequence type
        log_probs(:) = event.seq_type_posterior;
        ind = sample_logprobs(log_probs); 

        % Convert linear index to matrix subscripts (row, column)
        [r, w] = ind2sub(size(log_probs), ind);
        event.sampled_type = r;
        event.sampled_warp = w;

        % Sample event time
        sigma2 = 1 / event.summed_precisions(r, w);
        mu = event.summed_potentials(r, w) * sigma2;
        event.sampled_timestamp = mu + sqrt(sigma2) * randn();

        % Sample event amplitude
        event.sampled_amplitude = rand_posterior(event.spike_count, model.priors.seq_event_amplitude);
        model.sequence_events.evnts(it) = event;
    end
end

function [z,model] = gibbs_add_datapoint(model, spike)
    % Create log-probability vector to sample assignments.
    %
    %  - We need to sample K + 2 possibilities. There are K existing clusters
    %    we could assign to. We could also form a new cluster (index K + 1),
    %    or assign the spike to the background (index K + 2).
    
    % Shape and rate parameters of gamma prior on latent event amplitude
    alpha = model.priors.seq_event_amplitude.alpha;
    beta = model.priors.seq_event_amplitude.beta;

    events = model.sequence_events;
    indices = events.indices;

    K = length(indices);

    % Create log-probability vector to sample assignments
    log_probs = zeros(1, K + 2);

    % Iterate over non-empty model events
    for it = 1:length(indices)
        k = indices(it);
        event = model.sequence_events.evnts(k);

        % Check if event is too far away to be considered
        too_far = abs(spike.timestamp - event.sampled_timestamp) > model.max_sequence_length;
        if too_far && (event.sampled_type > 0)
            log_probs(k) = -Inf;
        else
            % Compute probability of adding spike to cluster k
            log_probs(k) = log(event.spike_count + alpha) + model.log_posterior_predictive(event, spike);
        end
    end

    % New cluster probability
    log_probs(K + 1) = model.new_cluster_log_prob + log_posterior_predictive_singleton(model,spike);

    % Background probability
    log_probs(K + 2) = model.bkgd_log_prob + model.globals.bkgd_log_proportions(spike.neuron) - log(model.max_time);

    % Sample new assignment for spike x
    z = sample_logprobs(log_probs); % Assuming sample_logprobs is a custom function

    % Handle new sample
    if z == (K + 2)
        z = -1; % Background
    elseif z == (K + 1) || (K == 1 && isempty(model.sequence_events.indices))
        [z,model] = add_event(model, spike); % New cluster
    else
        k = model.sequence_events.indices(z); % Existing sequence event
        [z,model] = add_datapoint(model, spike, k);
    end
end

function model = gibbs_update_globals(model, spikes, assignments)
    N = num_neurons(model);
    R = num_sequence_types(model);

    priors = model.priors;
    globals = model.globals;

    %% RESAMPLE BACKGROUND SPIKE PARAMETERS
    num_bkgd_spikes = sum(assignments < 0);
    bkgd_counts = hist([spikes(assignments < 0).neuron], 1:N);

    % Sample proportions of background spikes across neurons and map to log space
    prior_bkg = priors.bkgd_proportions;
    globals.bkgd_log_proportions = log(rand(prior_bkg.posterior(bkgd_counts)));

    % Sample the rate of background spikes
    alpha_bg = priors.bkgd_amplitude.alpha;
    beta_bg = priors.bkgd_amplitude.beta / model.max_time;
    globals.bkgd_amplitude = gamrnd(num_bkgd_spikes + alpha_bg, 1 / (1 + beta_bg));

    %% RESAMPLE SEQUENCE TYPE PROPORTIONS
    seq_type_counts = zeros(1, R);
    indices = model.sequence_events.indices;
    for it = 1:length(indices)
        k = indices(it);
        event = model.sequence_events.evnts(k);
        seq_type_counts(event.sampled_type) = seq_type_counts(event.sampled_type) + 1;
    end

    globals.seq_type_log_proportions = log(rand(...
        posterior(priors.seq_type_proportions,seq_type_counts)));

    %% RESAMPLE NEURON RESPONSE PROFILES
    spk_count = zeros(N, R);
    spk_fm = zeros(N, R);
    spk_sm = zeros(N, R);

    for i = 1:length(spikes)
        % spike assignment to latent event
        x = spikes(i);
        k = assignments(i);

        % skip if spike is assigned to background
        if k < 0
            continue;
        end

        event = model.sequence_events.evnts(k);

        n = x.neuron;
        r = event.sampled_type;
        w = event.sampled_warp;
        offset = (x.timestamp - event.sampled_timestamp) / ...
            model.priors.warp_values(w);

        % compute sufficient statistics
        spk_count(n, r) = spk_count(n, r) + 1;
        spk_fm(n, r) = spk_fm(n, r) + offset;
        spk_sm(n, r) = spk_sm(n, r) + offset^2;
    end

    for r = 1:R
        neuron_response_proportions(:, r) = rand(...
            posterior(priors.neuron_response_proportions,spk_count(:, r)));

        for n = 1:N
            % Sample from the updated parameters of the normal-inverse-chi-squared distribution
            [mu, sigma2] = rand(posterior(...
                priors.neuron_response_profile,...
                spk_count(n, r), ...
                spk_fm(n, r), ...
                spk_sm(n, r)...
                ));

            % Store the sampled values in the corresponding entries of globals
            globals.neuron_response_offsets(n, r) = mu;
            globals.neuron_response_widths(n, r) = sigma2;
        end
    end

    globals.neuron_response_log_proportions = log(neuron_response_proportions);

    %% RECOMPUTE NEW CLUSTER PROBABILITIES
    alpha = priors.seq_event_amplitude.alpha;
    beta = priors.seq_event_amplitude.beta;
    lambda = priors.seq_event_rate;

    model.new_cluster_log_prob = log(alpha) + log(lambda) + log(model.max_time) + alpha * (log(beta) - log(1 + beta));
    model.bkgd_log_prob = log(globals.bkgd_amplitude) + log(model.max_time) + log(1 + beta);

    %% RECOMPUTE SUFFICIENT STATISTICS
    model = recompute(model, spikes, assignments);
end