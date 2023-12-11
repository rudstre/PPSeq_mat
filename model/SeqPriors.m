classdef SeqPriors
    properties
        seq_event_rate
        seq_type_proportions
        seq_event_amplitude
        neuron_response_proportions
        neuron_response_profile
        bkgd_amplitude
        bkgd_proportions
        warp_values
        warp_log_proportions
    end
    
    methods
        function obj = SeqPriors(seq_event_rate, seq_type_proportions, seq_event_amplitude, neuron_response_proportions, neuron_response_profile, bkgd_amplitude, bkgd_proportions, warp_values, warp_log_proportions)
            obj.seq_event_rate = seq_event_rate;
            obj.seq_type_proportions = seq_type_proportions;
            obj.seq_event_amplitude = seq_event_amplitude;
            obj.neuron_response_proportions = neuron_response_proportions;
            obj.neuron_response_profile = neuron_response_profile;
            obj.bkgd_amplitude = bkgd_amplitude;
            obj.bkgd_proportions = bkgd_proportions;
            obj.warp_values = warp_values;
            obj.warp_log_proportions = warp_log_proportions;
        end
    end
end