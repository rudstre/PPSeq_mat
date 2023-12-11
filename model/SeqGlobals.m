classdef SeqGlobals
    properties
        seq_type_log_proportions
        neuron_response_log_proportions
        neuron_response_offsets
        neuron_response_widths
        bkgd_amplitude
        bkgd_log_proportions
    end
    
    methods
        function obj = SeqGlobals(seq_type_log_proportions, neuron_response_log_proportions, neuron_response_offsets, neuron_response_widths, bkgd_amplitude, bkgd_log_proportions)
            obj.seq_type_log_proportions = seq_type_log_proportions;
            obj.neuron_response_log_proportions = neuron_response_log_proportions;
            obj.neuron_response_offsets = neuron_response_offsets;
            obj.neuron_response_widths = neuron_response_widths;
            obj.bkgd_amplitude = bkgd_amplitude;
            obj.bkgd_log_proportions = bkgd_log_proportions;
        end
    end
end