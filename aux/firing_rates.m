function F = firing_rates(globals, events, time_bins)
    num_neurons = length(globals.bkgd_log_proportions);
    bkgd_firing_rates = globals.bkgd_amplitude * exp(globals.bkgd_log_proportions);

    % Initialize the firing rate matrix
    F = repmat(bkgd_firing_rates, 1, length(time_bins));

    for e = events
        for n = 1:num_neurons
            mu = e.timestamp + globals.neuron_response_offsets(n, e.seq_type);
            sigma = sqrt(globals.neuron_response_widths(n, e.seq_type));
            F(n, :) = F(n, :) + e.amplitude * exp(globals.neuron_response_log_proportions(n, e.seq_type)) * normpdf(time_bins, mu, sigma);
        end
    end
end
