function index = sample_logprobs(log_probs)
    % Converts log probabilities to probabilities using softmax
    probs = exp(log_probs - max(log_probs)); % for numerical stability
    probs = probs / sum(probs);

    % Sample from the distribution
    index = randsample(1:length(probs), 1, true, probs);
end