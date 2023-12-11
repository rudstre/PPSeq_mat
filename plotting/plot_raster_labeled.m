function fig = plot_raster_labeled(spikes, events, spike_assignments, neuron_order, color_cycle, varargin)
    % Plots a model-labeled spike raster
    
    if nargin < 5
        color_cycle = {'#E41A1C', '#377EB8', '#4DAF4A', '#984EA3', '#FF7F00', '#FFFF33', '#A65628', '#F781BF'};
    end
    
    if all(ishex(color_cycle))
        color_cycle = [0,0,0; cell2mat(cellfun(@hex2rgb,color_cycle,'uni',false)')];
    end
    
    cmap = colormap(color_cycle);

    fig = figure;
    x = [spikes.timestamp];
    y = arrayfun(@(n) find(neuron_order == n, 1), [spikes.neuron]);
    c = repmat(1, length(spikes), 1); % Default color black

    typemap = containers.Map([events.assignment_id], [events.seq_type]);
    
    for i = 1:length(spikes)
        if spike_assignments(i) ~= -1
            k = typemap(spike_assignments(i));
            c(i) = 1 + mod(k - 1, length(color_cycle));
        end
    end

    scatter(x, y, 4, c, varargin{2:end}); % Additional arguments passed via varargin
    ylabel('neurons');
    xlabel('time (s)');
end
