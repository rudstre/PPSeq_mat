classdef SeqEventList

    properties
        evnts    % Array of SeqEvent objects
        indices   % Array of indices for non-empty events
    end

    methods

        function obj = SeqEventList(num_sequence_types, num_warp_values)
            % Constructor to initialize the SeqEventList
            if isa(num_sequence_types,"SeqEvent")
                obj.evnts = num_sequence_types;
                obj.indices = num_warp_values;
            else
                obj.evnts = SeqEvent(num_sequence_types, num_warp_values);
                obj.indices = [];
            end
        end

        function event = get_event(obj, i)
            % Get event by index
            event = obj.evnts(i);
        end

        function len = get_length(obj)
            % Get the number of non-empty events
            len = length(obj.indices);
        end

        function [idx,obj] = add_event(obj)
            % Add a new event or find an empty event
            i = 1;
            for it = 1:length(obj.indices)
                j = obj.indices(it);
                if i ~= j
                    obj.indices = [obj.indices(1:i-1), i, obj.indices(i:end)];
                    idx = i;
                    return;
                end
                i = i + 1;
            end
            obj.indices(end + 1) = length(obj.indices) + 1;
            if length(obj.evnts) < length(obj.indices)
                [R,W] = size(obj.evnts(1).summed_potentials);
                obj.evnts(end + 1) = SeqEvent(R,W);
            end
            idx = obj.indices(end);
        end

        function obj = remove_event(obj, index)
            % Remove an event and reset its statistics
            obj.evnts(index) = obj.evnts(index).reset();
            obj.indices(obj.indices == index) = [];
        end

        function obj = recompute_distributed(obj, model, spikes, assignments)
            % Partition spikes and pass assignments to submodels
            [spk_partition, assgn_partition, ~] = partition_spikes(model, spikes, assignments); % Assuming partition_spikes is a custom function
            for i = 1:length(model.submodels)
                model.submodels(i) = recompute(model.submodels(i), spk_partition{i}, assgn_partition{i});
            end
        end

        function obj = reset_events(obj)
            for it = 1:length(obj.evnts)
                obj.evnts(it) = obj.evnts(it).reset;
            end
        end
    end
end