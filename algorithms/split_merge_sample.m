function [assignments,model] = split_merge_sample(model, spikes, num_samples, assignments, split_merge_window)
    if num_samples == 0
        return;
    end

    % Spike indices not assigned to the background
    non_bkgd_idx = find(assignments ~= -1);

    if length(non_bkgd_idx) < 3
        return;
    end

    % Concentration parameter
    alpha = model.priors.seq_event_amplitude.alpha;

    for s = 1:num_samples

        % assgn_hist = hist(assignments,max(assignments)+2);
        % if
        % ~isequal([model.sequence_events.evnts.spike_count],assgn_hist(3:end))
        %     a = 5;
        % end

        % Sample random spike
        i = non_bkgd_idx(randi(length(non_bkgd_idx)));
        ti = spikes(i).timestamp;

        % Sample another random spike within 'split_merge_window' time units from spike i
        A = non_bkgd_idx(abs([spikes(non_bkgd_idx).timestamp] - ti) < split_merge_window & non_bkgd_idx' ~= i);

        if isempty(A)
            continue;
        end

        % Sample second spike and empty set A
        j = A(randi(length(A)));

        % Cluster assignment for spike i and spike j
        ci = assignments(i);
        cj = assignments(j);

        A = i; % Put spikes[i] in set A
        B = j; % Put spikes[j] in set B

        %% Propose split or merge
        if ci == cj
            % Propose split

            % Randomly assign other spikes in cluster to A or B
            for it = 1:length(non_bkgd_idx)
                k = non_bkgd_idx(it);
                if assignments(k) == ci && k ~= i && k ~= j
                    if rand() > 0.5
                        A(end + 1) = k;
                    else
                        B(end + 1) = k;
                    end
                end
            end
            C = [A, B];

            % Compute log acceptance probability
            accept_log_prob = event_log_like(model, spikes(A)) ...
                + event_log_like(model, spikes(B)) ...
                - event_log_like(model, spikes(C)) ...
                + gammaln(length(A) + alpha) ...
                + gammaln(length(B) + alpha) ...
                - gammaln(length(C) + alpha) ...
                - gammaln(alpha) + model.new_cluster_log_prob - (length(C) - 2) * log(0.5);

            % Accept move
            if rand() < exp(accept_log_prob)
                % Create new event for spike j and move it there
                for it = 1:length(B)
                    b = B(it);
                    model = remove_datapoint(model, spikes(b), assignments(b));
                    if b == j
                        [new_event,model] = add_event(model, spikes(b));
                        assignments(b) = new_event;
                    else
                        [assignments(b),model] = add_datapoint(model, spikes(b), new_event);
                    end
                end
            end
        else
            % Propose merge
            
            % Assign spikes to A or B based on current cluster assignment
            for it = 1:length(non_bkgd_idx)
                k = non_bkgd_idx(it);
                if k ~= i && k ~= j
                    if assignments(k) == ci
                        A(end + 1) = k;
                    elseif assignments(k) == cj
                        B(end + 1) = k;
                    end
                end
            end
            C = [A, B];

            % Compute log acceptance probability
            accept_log_prob = event_log_like(model, spikes(C)) - event_log_like(model, spikes(A)) ...
                - event_log_like(model, spikes(B)) + gammaln(length(C) + alpha) ...
                + gammaln(alpha) - gammaln(length(A) + alpha) - gammaln(length(B) + alpha) ...
                - model.new_cluster_log_prob + (length(C) - 2) * log(0.5);

            % Accept move
            if rand() < exp(accept_log_prob)
                % Move spikes in set B to cluster ci
                for it = 1:length(B)
                    b = B(it);
                    model = remove_datapoint(model, spikes(b), assignments(b));
                    [assignments(b),model] = add_datapoint(model, spikes(b), ci);
                end
            end
        end
    end
end

