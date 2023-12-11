classdef SeqModel
    properties
        max_time
        max_sequence_length
        priors
        globals
        sequence_events
        new_cluster_log_prob
        bkgd_log_prob
        R_buffer
        RW_buffer
        K_buffer
    end

    methods
        %% Constructor
        function obj = SeqModel(max_time, max_sequence_length, ...
                num_warp_values, max_warp, warp_variance, ...
                seq_event_rate, seq_type_proportions, ...
                seq_event_amplitude, neuron_response_proportions, ...
                neuron_response_profile, bkgd_amplitude, ...
                bkgd_proportions)
            if nargin == 5
                obj.max_time = max_time;
                obj.max_sequence_length = max_sequence_length;
                obj.priors = num_warp_values;
                obj.globals = max_warp;
                obj.sequence_events = warp_variance;

                % Initialize R_buffer, RW_buffer, K_buffer based on obj dimensions
                obj.R_buffer = zeros(obj.priors.seq_type_proportions.dim, 1);
                obj.RW_buffer = zeros(obj.priors.seq_type_proportions.dim, length(obj.priors.warp_values));
                obj.K_buffer = zeros(1, length(obj.priors.warp_values) + 2);

            else
                % Number of neurons (N), sequence types (R), and warp values (W)
                N = neuron_response_proportions.dim;
                R = seq_type_proportions.dim;
                W = num_warp_values;

                % Specify discrete distribution over warp values
                [warp_values, warp_log_proportions] = default_warps(num_warp_values, max_warp, warp_variance);

                % Create prior distributions
                priors = SeqPriors(seq_event_rate, seq_type_proportions, seq_event_amplitude, ...
                    neuron_response_proportions, neuron_response_profile, ...
                    bkgd_amplitude, bkgd_proportions, warp_values, warp_log_proportions);

                % Initialize global variables
                globals = sample_priors(priors);

                % Initialize sequence events
                sequence_events = SeqEventList(R, W);

                % Compute probabilities for new cluster and background spike
                alpha = seq_event_amplitude.alpha;
                beta = seq_event_amplitude.beta;
                lambda = seq_event_rate;
                bkgd_log_prob = log(globals.bkgd_amplitude) + log(max_time) + log(1 + beta);
                new_cluster_log_prob = log(alpha) + log(lambda) + log(max_time) + alpha * (log(beta) - log(1 + beta));

                % Assign properties
                obj.max_time = max_time;
                obj.max_sequence_length = max_sequence_length;
                obj.priors = priors;
                obj.globals = globals;
                obj.sequence_events = sequence_events;
                obj.new_cluster_log_prob = new_cluster_log_prob;
                obj.bkgd_log_prob = bkgd_log_prob;
                obj.R_buffer = zeros(1, R);
                obj.RW_buffer = zeros(R, W);
                obj.K_buffer = [];
            end
        end

        %% Core object methods
        function [spikes, assignments, events] = sample_model(obj, resample_latents, resample_globals)
            % Check if resampling of global parameters is required
            if resample_globals
                glbls = sample_priors(obj.priors);
            end

            events = EventSummaryInfo.empty; % Initialize empty EventSummaryInfo array

            if resample_latents
                % Sample the number of sequence events (K) from a Poisson distribution
                K = poissrnd(obj.priors.seq_event_rate * obj.max_time);

                % Initialize distributions for sequence types and warps
                seq_type_dist = makedist('Multinomial', 'probabilities', exp(obj.globals.seq_type_log_proportions));
                warp_dist = makedist('Multinomial', 'probabilities', exp(obj.priors.warp_log_proportions));

                for k = 1:K
                    % Sample event properties
                    tau = rand() * obj.max_time;
                    r = random(seq_type_dist);
                    omega = obj.priors.warp_values(random(warp_dist));
                    A = random(obj.priors.seq_event_amplitude);

                    % Add sampled event to the list
                    events(end+1) = EventSummaryInfo(k, tau, r, omega, A);
                end
            else
                events = obj.event_list_summary();
            end

            % Sample spikes and assignments
            S_bkgd = poissrnd(glbls.bkgd_amplitude * obj.max_time);
            probabilities = exp(glbls.bkgd_log_proportions);
            probabilities = probabilities / sum(probabilities); % Normalizing to get probabilities
            bkgd_dist = makedist('Multinomial', 'probabilities', probabilities);

            n_bkgd = random(bkgd_dist, [1, S_bkgd]);
            t_bkgd = rand(1, S_bkgd) * obj.max_time;

            for idx = 1:S_bkgd
                n = n_bkgd(idx);
                t = t_bkgd(idx);
                spikes = [spikes; Spike(n, t)]; % Assuming Spike takes two arguments: neuron and timestamp
                assignments = [assignments, -1];
            end

            % Compute neuron probabilities
            neuron_rel_amps = exp(glbls.neuron_response_log_proportions);

            % Sample sequence-evoked spikes
            for e = events
                % Num spikes evoked by latent event
                S = poissrnd(e.amplitude);

                % Neuron response proportions
                probabilities = neuron_rel_amps(:, e.seq_type);
                probabilities = probabilities / sum(probabilities); % Normalizing
                nrn_dist = makedist('Multinomial', 'probabilities', probabilities);

                for idx = 1:S
                    n = random(nrn_dist);

                    % Sample spike time
                    mu = glbls.neuron_response_offsets(n, e.seq_type);
                    sigma = sqrt(glbls.neuron_response_widths(n, e.seq_type));
                    t = e.seq_warp * (sigma * randn() + mu) + e.timestamp;

                    % Exclude spikes outside of time window
                    if (0 < t) && (t < model.max_time)
                        spikes = [spikes; Spike(n, t)];
                        assignments = [assignments, e.assignment_id];
                    end
                end
            end
        end

        function obj = anneal_priors(obj,gmdist)
            pr = obj.priors;
            pr_annealed = pr; pr_annealed.seq_event_amplitude = gmdist;
            obj.set_priors(pr_annealed);
        end

        % Recompute sufficient statistics for all sequence events
        function obj = recompute(obj, spikes, assignments)
            % Reset all events to zero spikes
            obj.sequence_events = obj.sequence_events.reset_events();
            obj.sequence_events.indices = [];

            % Add spikes back to their previously assigned event
            for i = 1:length(spikes)
                s = spikes(i);
                k = assignments(i);

                % Skip background spikes
                if k < 0
                    continue;
                end

                % Check that event k exists
                while k > length(obj.sequence_events.evnts)
                    num_sequence_types = numel(obj.priors.seq_type_proportions);
                    num_warp_values = length(obj.priors.warp_values);
                    obj.sequence_events.evnts(end + 1) = SeqEvent(num_sequence_types, num_warp_values);
                end

                % Add datapoint to event k
                [k,obj] = add_datapoint(obj, s, k, false);

                % Make sure that event k is marked as non-empty
                if ~ismember(k, obj.sequence_events.indices)
                    obj.sequence_events.indices(end + 1) = k;
                end
            end

            % Set the posterior for each event
            for k = obj.sequence_events.indices
                obj = set_posterior(obj, k);
            end
        end

        %% Add/Remove datapoint methods
        function [k,obj] = add_datapoint(obj, s, k, recompute_posterior)
            if nargin < 4
                recompute_posterior = true;
            end

            n = s.neuron;
            t = s.timestamp;
            event = obj.sequence_events.evnts(k);

            log_p_neuron = obj.globals.neuron_response_log_proportions;
            offsets = obj.globals.neuron_response_offsets;
            widths = obj.globals.neuron_response_widths;
            warps = obj.priors.warp_values;

            event.spike_count = event.spike_count + 1;

            for r = 1:num_sequence_types(obj)
                for w = 1:num_warp_values(obj)
                    v = 1 / (widths(n, r) * warps(w)^2);
                    m = (t - offsets(n, r) * warps(w)) * v;

                    event.summed_potentials(r, w) = event.summed_potentials(r, w) + m;
                    event.summed_precisions(r, w) = event.summed_precisions(r, w) + v;
                    event.summed_logZ(r, w) = event.summed_logZ(r, w) + (log_p_neuron(n, r) - gauss_info_logZ(m, v));
                end
            end

            if recompute_posterior
                for r = 1:num_sequence_types(obj)
                    for w = 1:num_warp_values(obj)
                        event.seq_type_posterior(r, w) = obj.globals.seq_type_log_proportions(r) + obj.priors.warp_log_proportions(w) + event.summed_logZ(r, w) + gauss_info_logZ(event.summed_potentials(r, w), event.summed_precisions(r, w));
                    end
                end

                event.seq_type_posterior = event.seq_type_posterior - logsumexp(event.seq_type_posterior);
            end
            obj.sequence_events.evnts(k) = event;
        end

        function obj = remove_datapoint(obj, s, k)
            if k == -1
                return;
            end

            n = s.neuron;
            t = s.timestamp;
            event = obj.sequence_events.evnts(k);

            log_p_neuron = obj.globals.neuron_response_log_proportions;
            offsets = obj.globals.neuron_response_offsets;
            widths = obj.globals.neuron_response_widths;
            warps = obj.priors.warp_values;

            if event.spike_count == 1
                obj.sequence_events = remove_event(obj.sequence_events, k);
                return;
            end

            event.spike_count = event.spike_count - 1;

            for r = 1:num_sequence_types(obj)
                for w = 1:num_warp_values(obj)
                    v = 1 / (widths(n, r) * warps(w)^2);
                    m = (t - offsets(n, r) * warps(w)) * v;
                    sv = event.summed_precisions(r, w);

                    event.summed_potentials(r, w) = event.summed_potentials(r, w) - m;
                    event.summed_precisions(r, w) = max(0, sv - v);
                    event.summed_logZ(r, w) = event.summed_logZ(r, w) - (log_p_neuron(n, r) - gauss_info_logZ(m, v)); % Assuming gauss_info_logZ is a custom function

                    event.seq_type_posterior(r, w) = obj.globals.seq_type_log_proportions(r) + obj.priors.warp_log_proportions(w) + event.summed_logZ(r, w) + gauss_info_logZ(event.summed_potentials(r, w), event.summed_precisions(r, w));
                end
            end

            event.seq_type_posterior = event.seq_type_posterior - logsumexp(event.seq_type_posterior); % Assuming logsumexp is a custom function

            obj.sequence_events.evnts(k) = event;
        end

        function [k,obj] = add_event(obj, s)
            [k,obj.sequence_events] = add_event(obj.sequence_events);
            [k,obj] = add_datapoint(obj, s, k);
        end

        %% Probability methods
        function log_prob = log_posterior_predictive(obj, event, spike)
            % Shorter names for convenience
            n = spike.neuron;
            log_p_neuron = obj.globals.neuron_response_log_proportions;
            offsets = obj.globals.neuron_response_offsets;
            widths = obj.globals.neuron_response_widths;
            warps = obj.priors.warp_values;

            % Grab size (R,W) matrix (already pre-allocated)
            log_prob = obj.RW_buffer;

            % Compute p({x1...xm} | r)
            for r = 1:obj.priors.seq_type_proportions.dim
                for w = 1:length(obj.priors.warp_values)
                    m = event.summed_potentials(r, w);
                    v = event.summed_precisions(r, w);
                    mu = m / v;
                    sigma2 = 1 / v;

                    % Calculating the log probability
                    if event.spike_count == 0
                        log_prob(r, w) = log_p_neuron(n, r);
                    else
                        log_prob(r, w) = (event.seq_type_posterior(r, w) + log_p_neuron(n, r) + log(normpdf(spike.timestamp,mu + offsets(n, r) * warps(w), sqrt(sigma2 + widths(n, r) * warps(w)))));
                    end
                end
            end
            log_prob = logsumexp(log_prob);
        end

        function log_prob = log_posterior_predictive_singleton(obj,spike)
            R = num_sequence_types(obj);
            W = num_warp_values(obj);
            n = spike.neuron;

            log_p_neuron = obj.globals.neuron_response_log_proportions;
            log_p_seqtype = obj.globals.seq_type_log_proportions;
            log_p_warps = obj.priors.warp_log_proportions;

            % Grab length-R array
            log_prob = obj.R_buffer;

            % Compute probability of singleton: p(x) = sum_r ( p(n | r) p(r) )
            for r = 1:R
                log_prob(r) = log_p_neuron(n, r) + log_p_seqtype(r);
            end

            % Combine probabilities in log space using logsumexp
            log_prob = logsumexp(log_prob);
        end

        function log_p = log_prior(obj)
            pr = obj.priors;
            glbls = obj.globals;

            R = num_sequence_types(obj); % num sequence types
            N = num_neurons(obj);        % num neurons

            % Rate of background spikes
            lp = logpdf(pr.bkgd_amplitude, glbls.bkgd_amplitude); % logpdf for gamma distribution

            % Background spike proportions across neurons using Dirichlet distribution
            for n = 1:N
                lp = lp + logpdf(pr.bkgd_proportions, exp(glbls.bkgd_log_proportions(n)));
            end

            % Add contribution of each sequence type
            for r = 1:R
                % Log prior on neuron amplitudes for each sequence type
                lp = lp + logpdf(pr.neuron_response_proportions, exp(glbls.neuron_response_log_proportions(:, r)));

                % Neuron offsets and widths for each sequence type
                for n = 1:N
                    lp = lp + logpdf(pr.neuron_response_profile, glbls.neuron_response_offsets(n, r), glbls.neuron_response_widths(n, r));
                end
            end

            log_p = lp; % Sum of all log probabilities
        end

        function lp = log_p_latents(obj)
            pr = obj.priors;
            glbls = obj.globals;
            events = obj.sequence_events;

            lp = 0.0;

            % For each latent event add to the log probability
            for event_idx = 1:length(events)
                event = events(event_idx);

                % Log prior on event amplitude (gamma distribution)
                lp = lp + logpdf_gamma(pr.seq_event_amplitude, event.sampled_amplitude);

                % Log probability of event type
                lp = lp + glbls.seq_type_log_proportions(event.sampled_type);

                % Log probability of event warp
                lp = lp + pr.warp_log_proportions(event.sampled_warp);
            end
        end

        % Starting translation of the 'log_like' function
        function ll = log_like(obj, spikes)
            glbls = obj.globals;
            seq_events = obj.sequence_events;
            t_max = obj.max_time;

            bkgd_amplitude = glbls.bkgd_amplitude;
            offsets = glbls.neuron_response_offsets;
            widths = glbls.neuron_response_widths;
            warps = obj.priors.warp_values;
            ll = 0.0;

            % Sum of Poisson Process intensity at all spikes
            for spike_idx = 1:length(spikes)
                x = spikes(spike_idx);
                g = bkgd_amplitude * exp(glbls.bkgd_log_proportions(x.neuron));
                for it = 1:length(seq_events.indices)
                    event_idx = seq_events.indices(it);
                    event = seq_events.evnts(event_idx);
                    w = warps(event.sampled_warp);
                    mu = event.sampled_timestamp + offsets(x.neuron, event.sampled_type) * w;
                    sigma2 = widths(x.neuron, event.sampled_type) * (w ^ 2);
                    g = g + event.sampled_amplitude * normpdf(mu, sqrt(sigma2), x.timestamp);
                end
                ll = ll + log(g);
            end

            % Penalty on integrated intensity function
            ll = ll - bkgd_amplitude * t_max;
            for event_idx = 1:length(seq_events.evnts)
                event = seq_events.evnts(event_idx);
                ll = ll - event.sampled_amplitude;
            end
        end

        % Starting translation of the 'event_log_like' function
        function ll = event_log_like(obj, spikes)
            % Compute likelihood of initial spike
            ll = log_posterior_predictive_singleton(obj,spikes(1));

            % Return early if one spike
            if length(spikes) == 1
                return;
            end

            % Create singleton cluster and account for remaining spikes
            [k,obj] = add_event(obj, spikes(1));
            event = obj.sequence_events.evnts(k);
            i = 2;
            while i <= length(spikes)
                ll = ll + log_posterior_predictive(obj, event, spikes(i));
                if i == length(spikes)
                    return;
                end
                [~,obj] = add_datapoint(obj, spikes(i), k);
                i = i + 1;
            end
        end

        %% Getters and setters for priors, globals, and log probabilities
        function max_time = get_max_time(obj)
            max_time = obj.max_time;
        end

        function bkgd_amplitude = get_bkgd_amplitude(obj)
            bkgd_amplitude = obj.globals.bkgd_amplitude;
        end

        function infos = event_list_summary(obj)
            % Generates a summary list of event information from Seqobj

            ev = obj.sequence_events;
            warp_vals = obj.priors.warp_values;
            infos = [];  % Initialize empty array for EventSummaryInfo

            for ind = 1:length(ev.indices)
                % Accessing event information by index
                event = ev.evnts(ev.indices(ind));

                % Create and append new EventSummaryInfo
                info.assignment_id = ev.indices(ind);
                info.timestamp = event.sampled_timestamp;
                info.seq_type = event.sampled_type;
                info.seq_warp = warp_vals(event.sampled_warp);
                info.amplitude = event.sampled_amplitude;

                infos = [infos, info];  % Append info to the list
            end
        end

        function priors = get_priors(obj)
            priors = obj.priors;
        end

        function obj = set_priors(obj, new_priors)
            obj.priors = new_priors;
        end

        function obj = set_posterior(obj, k)
            % Extract the specific event
            event = obj.sequence_events.evnts(k);

            % Get the number of sequence types and warp values from the model
            R = num_sequence_types(obj);
            W = num_warp_values(obj);

            % Update seq_type_posterior
            for r = 1:R
                for w = 1:W
                    event.seq_type_posterior(r, w) = obj.globals.seq_type_log_proportions(r) ...
                        + obj.priors.warp_log_proportions(w) ...
                        + event.summed_logZ(r, w) ...
                        + gauss_info_logZ(event.summed_potentials(r, w), ...
                        event.summed_precisions(r, w));
                end
            end

            % Normalize using logsumexp
            event.seq_type_posterior = event.seq_type_posterior - logsumexp(event.seq_type_posterior(:));

            % Update the event in the model
            obj.sequence_events.evnts(k) = event;
        end

        function obj = set_globals(obj, new_globals)
            obj.globals = new_globals;
        end

        function obj = set_new_cluster_log_prob(obj, prob)
            obj.new_cluster_log_prob = prob;
        end

        function obj = set_bkgd_log_prob(obj, prob)
            obj.bkgd_log_prob = prob;
        end

        function out = num_neurons(obj), out = obj.priors.bkgd_proportions.dim; end

        function out = num_sequence_types(obj), out = obj.priors.seq_type_proportions.dim; end

        function out = num_warp_values(obj), out = length(obj.priors.warp_values); end

        function out = num_sequence_events(obj), out = length(obj.sequence_events); end
        
    end
end

%% Helper functions

% Initialize warp probabilities
function [warp_values, warp_log_proportions] = default_warps(num_warp_values, tau_max, warp_variance)
    assert(num_warp_values >= 1, 'Number of warp values must be at least 1.');
    assert(tau_max >= 1, 'Max warp value must be >= 1 by convention.');
    assert(warp_variance >= 0, 'Warp variance must be non-negative.');

    if num_warp_values == 1
        warp_values = 1;
        warp_log_proportions = 0;
    else
        warp_values = logspace(-1, 1, num_warp_values) * tau_max;
        warp_log_proportions = -0.5 / warp_variance * log(warp_values).^2;
        warp_log_proportions = warp_log_proportions - logsumexp(warp_log_proportions);
    end
end

% Generates random samples based on the given priors
function glbls = sample_priors(pr)

    % Extract dimensions
    N = pr.neuron_response_proportions.dim;
    R = pr.seq_type_proportions.dim;

    % Initialize length-R probability vector over sequence types
    seq_type_proportions_samples = pr.seq_type_proportions.rand();
    seq_type_log_proportions = log(seq_type_proportions_samples);

    % Draw N x R matrix holding neuron response amplitudes
    neuron_response_log_proportions = log(pr.neuron_response_proportions.rand([N R]));

    % Draw N x R matrix holding neuron response offsets and widths
    neuron_response_offsets = zeros(N, R);
    neuron_response_widths = zeros(N, R);
    for n = 1:N
        for r = 1:R
            [mu, sigma2] = rand(pr.neuron_response_profile);
            neuron_response_offsets(n, r) = mu;
            neuron_response_widths(n, r) = sigma2;
        end
    end

    % Draw background rate parameter
    bkgd_amplitude = rand(pr.bkgd_amplitude);

    % Draw length-N probability vector holding normalized background firing rates
    bkgd_log_proportions = log(pr.bkgd_proportions.rand([1, N]));

    % Construct the SeqGlobals object
    glbls = SeqGlobals(seq_type_log_proportions, neuron_response_log_proportions, ...
        neuron_response_offsets, neuron_response_widths, ...
        bkgd_amplitude, bkgd_log_proportions);

end

% Utility function for log-sum-exp
function logsum = logsumexp(log_probs)
    max_log = max(log_probs);
    if max_log == -inf, logsum = -inf; return, end
    logsum = max_log + log(sum(exp(log_probs - max_log)));
end

function out = gauss_info_logZ(m, v), out = 0.5 * (log(2*pi) - log(v) + m * m / v); end