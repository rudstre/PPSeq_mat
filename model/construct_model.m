function model = construct_model(config, max_time, num_neurons)
    % Construct a SeqModel object from the configuration dictionary

    % Prior on sequence type proportions
    seq_type_proportions = SymmetricDirichlet(config.seq_type_conc_param, config.num_sequence_types);

    % Prior on expected number of spikes induced by sequence events
    seq_event_amplitude = specify_gamma(config.mean_event_amplitude, config.var_event_amplitude);

    % Prior on relative response amplitudes per neuron to each sequence type
    neuron_response_proportions = SymmetricDirichlet(config.neuron_response_conc_param, num_neurons);

    % Prior on the response offsets and widths for each neuron
    neuron_response_profile = NormalInvChisq(config.neuron_offset_pseudo_obs, 0.0, config.neuron_width_pseudo_obs, config.neuron_width_prior);

    % Prior on expected number of background spikes in a unit time interval
    bkgd_amplitude = specify_gamma(config.mean_bkgd_spike_rate, config.var_bkgd_spike_rate);

    % Prior on relative background firing rates across neurons
    bkgd_proportions = SymmetricDirichlet(config.bkgd_spikes_conc_param, num_neurons);

    % Construct the SeqModel
    model = SeqModel(max_time, config.max_sequence_length, config.num_warp_values, config.max_warp, config.warp_variance, config.seq_event_rate, seq_type_proportions, seq_event_amplitude, neuron_response_proportions, neuron_response_profile, bkgd_amplitude, bkgd_proportions);

    % If distributed processing is enabled, return a DistributedSeqModel
    if isfield(config, 'num_threads') && config.num_threads > 0
        model = DistributedSeqModel(model, config.num_threads);
    end
end

function result = specify_gamma(mean_val, var_val)
    % Utility function to specify a gamma distribution based on mean and variance
    alpha = (mean_val^2) / var_val;
    beta = mean_val / var_val;
    result = RateGamma(alpha, beta);
end