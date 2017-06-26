
%% 1. Run watershed segmentation & classification on your image
% when you are done, save the data and load it back into MATLAB, or leave
% the GUI window open, so "h" is still accessible in memory
h = watershed_cells_gui 

%% 2. After running the segmentation & classification, gather the results 

% Parameters & results are stored here
h.UserData.params
h.UserData.results

% The "label matrix" defines the segmentation results (i.e. cell regions)
label_matrix = h.UserData.results.segmentation.label_matrix;
nRegions = h.UserData.results.segmentation.number;

% The classification function computes one value of f(R,G,B) for each cell
f_values = h.UserData.results.classification.f_values;

% The threshold determines if that cell is in state 1 or 2
f_threshold = h.UserData.params.classification.threshold.value;
in_state1 = f_values > f_threshold;
in_state2 = f_values <= f_threshold;

%% 3. Analyze the location of each cell, based on the classification

% Use the label matrix to compute the (x,y) centroid of each cell
stats = regionprops(label_matrix, 'Centroid');
centroids = cat(1, stats.Centroid);

% Divide the centroids based on the classification result
centroids_state1 = centroids(in_state1,:);
centroids_state2 = centroids(in_state2,:);

% bin the centroids along the X direction for cells in state 1
x_state1 = centroids_state1(:,1); % units: pixels from the top-left corner
x_state2 = centroids_state2(:,1); % units: pixels from the top-left corner
x_all = centroids(:,1); % units: pixels from the top-left corner
delta = 25; % bin spacing, units: pixels
max_x = max(x_all);
bin_edges = 0 : delta : (max_x+delta);
bin_centers = bin_edges(1:end-1) + delta/2;
numcells_state1 = histcounts(x_state1, bin_edges);
numcells_state2 = histcounts(x_state2, bin_edges);
numcells_total = histcounts(x_all, bin_edges);

% plot the binned centroids 
fig = figure('color', 'w');
ax = axes('parent', fig);
ax.NextPlot = 'add';
plot(ax, bin_centers, numcells_state1, 'm-o', 'displayname', 'State 1');
plot(ax, bin_centers, numcells_state2, 'c-d', 'displayname', 'State 2');
plot(ax, bin_centers, numcells_total, 'k--', 'displayname', 'Total');
ax.XLabel.String = 'X position (pixels)';
ax.YLabel.String = 'Frequency';
ax.XLim = bin_edges([1 end]);
ax.Title.String = 'Histogram of cell x-positions by state';
lg = legend('show', 'location', 'best');







