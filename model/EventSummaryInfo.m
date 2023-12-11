classdef EventSummaryInfo
    properties
        assignment_id  % Integer
        timestamp      % Float
        seq_type       % Integer
        seq_warp       % Float
        amplitude      % Float
    end
    
    methods
        function obj = EventSummaryInfo(assignment_id, timestamp, seq_type, seq_warp, amplitude)
            obj.assignment_id = assignment_id;
            obj.timestamp = timestamp;
            obj.seq_type = seq_type;
            obj.seq_warp = seq_warp;
            obj.amplitude = amplitude;
        end
    end
end