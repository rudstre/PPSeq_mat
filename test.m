% Songbird metadata
num_neurons = 75;
max_time = 22.2;

% Randomly permute neuron labels
% (This hides the sequences, to make things interesting)
p = randperm(num_neurons);

% Load spikes from a tab-delimited text file
[file,path] = uigetfile('*.txt');
data = readmatrix(fullfile(path,file), 'Delimiter', '\t');
neuron = data(:, 1);
time = data(:, 2);

% Apply the permutation to neuron labels
permuted_neurons = arrayfun(@(x) p(x), neuron);

% Combine permuted neuron labels and spike times
for i = 1:size(data,1)
    spikes(i) = Spike(permuted_neurons(i),time(i));
end

% Configuration structure
config = struct(...
    'num_sequence_types', 2, ...
    'seq_type_conc_param', 1.0, ...
    'seq_event_rate', 1.0, ...
    'mean_event_amplitude', 100.0, ...
    'var_event_amplitude', 1000.0, ...
    'neuron_response_conc_param', 0.1, ...
    'neuron_offset_pseudo_obs', 1.0, ...
    'neuron_width_pseudo_obs', 1.0, ...
    'neuron_width_prior', 0.5, ...
    'num_warp_values', 1, ...
    'max_warp', 1.0, ...
    'warp_variance', 1.0, ...
    'mean_bkgd_spike_rate', 30.0, ...
    'var_bkgd_spike_rate', 30.0, ...
    'bkgd_spikes_conc_param', 0.3, ...
    'max_sequence_length', Inf, ...
    'num_anneals', 10, ...
    'samples_per_anneal', 100, ...
    'max_temperature', 40.0, ...
    'save_every_during_anneal', 10, ...
    'samples_after_anneal', 2000, ...
    'save_every_after_anneal', 10, ...
    'split_merge_moves_during_anneal', 10, ...
    'split_merge_moves_after_anneal', 10, ...
    'split_merge_window', 1.0 ...
);

% Initialize all spikes to background process.
init_assignments = -1 * ones(length(spikes), 1);

% Construct model struct
model = construct_model(config, max_time, num_neurons);

% Run Gibbs sampling with an initial annealing period.
results = easy_sample(model, spikes, init_assignments, config);

% Grab the final MCMC sample
final_globals = results.globals_hist{end};
final_events = results.latent_event_hist{end};
final_assignments = results.assignment_hist(:, end);

% Helpful utility function that sorts the neurons to reveal sequences.
neuron_ordering = sortperm_neurons(final_globals);

% Plot model-annotated raster
fig = plot_raster_labeled(spikes, final_events, final_assignments, neuron_ordering);