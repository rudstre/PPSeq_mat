classdef SeqEvent
    properties
        spike_count            % Integer
        summed_potentials      % Matrix
        summed_precisions      % Matrix
        summed_logZ            % Matrix
        seq_type_posterior     % Matrix
        sampled_type           % Integer
        sampled_warp           % Integer
        sampled_timestamp      % Float
        sampled_amplitude      % Float
    end

    methods
        % Constructor for SeqEvent
        function obj = SeqEvent(spike_count, summed_potentials, summed_precisions, summed_logZ, seq_type_posterior, sampled_type, sampled_warp, sampled_timestamp, sampled_amplitude)

            if nargin == 2
                num_sequence_types = spike_count; 
                num_warp_values = summed_potentials;

                obj.spike_count = 0;
                obj.summed_potentials = zeros(num_sequence_types, num_warp_values);
                obj.summed_precisions = zeros(num_sequence_types, num_warp_values);
                obj.summed_logZ = zeros(num_sequence_types, num_warp_values);
                obj.seq_type_posterior = zeros(num_sequence_types, num_warp_values);
                obj.sampled_type = -1;
                obj.sampled_warp = -1;
                obj.sampled_timestamp = 0.0;
                obj.sampled_amplitude = 0.0;
            elseif nargin == 9
                obj.spike_count = spike_count;
                obj.summed_potentials = summed_potentials;
                obj.summed_precisions = summed_precisions;
                obj.summed_logZ = summed_logZ;
                obj.seq_type_posterior = seq_type_posterior;
                obj.sampled_type = sampled_type;
                obj.sampled_warp = sampled_warp;
                obj.sampled_timestamp = sampled_timestamp;
                obj.sampled_amplitude = sampled_amplitude;
            end

        end

        function obj = reset(obj)
            % Reset the SeqEvent object to its initial state
            obj.spike_count = 0;
            obj.summed_potentials(:) = 0;
            obj.summed_precisions(:) = 0;
            obj.summed_logZ(:) = 0;
            obj.seq_type_posterior(:) = 0;
        end
    end
end
