function fig = plot_raster(spikes, varargin)
    % Plots an unlabeled spike raster
    
    fig = figure;
    x = [spikes.timestamp];
    y = [spikes.neuron];
    
    if nargin > 1 && varargin{1} == true
        y = y(randperm(length(y))); % Randomly permute neuron labels
    end
    
    scatter(x, y, 4, varargin{2:end}); % Additional arguments passed via varargin
    ylabel('neurons');
    xlabel('time (s)');
end
