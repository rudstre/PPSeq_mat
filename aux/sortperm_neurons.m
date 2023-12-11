function sorted_indices = sortperm_neurons(globals, thres)
    % Sort neurons by preferred sequence type and offset

    if nargin < 2
        thres = 0.2;  % Default threshold
    end

    resp_props = exp(globals.neuron_response_log_proportions);
    offsets = globals.neuron_response_offsets;
    peak_resp = max(resp_props, [], 2);

    % Identify neurons with response above the threshold
    has_response = peak_resp > quantile(peak_resp, thres);

    % Get preferred sequence type and delay for each neuron
    [~, preferred_type] = max(resp_props, [], 2);
    preferred_delay = arrayfun(@(n, r) offsets(n, r), (1:size(offsets, 1))', preferred_type);

    % Combine into a matrix and sort based on the criteria
    Z = [has_response, preferred_type, preferred_delay];
    [~, sorted_indices] = sortrows(Z, [1, -2, 3]);
end
