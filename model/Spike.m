classdef Spike
    properties
        neuron      % Integer
        timestamp   % Float
    end

    methods
        function obj = Spike(neuron, timestamp)
            obj.neuron = neuron;
            obj.timestamp = timestamp;
        end
    end
end