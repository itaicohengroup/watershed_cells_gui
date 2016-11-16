% WATERSHED_CELLS_GUI finds cells in an image using watershed segmentation
%
% Programmed by: Lena Bartell (lrb89@cornell.edu)
% Principle Investigator: Prof. Itai Cohen, Cornell University
%
% Quick-start guide:
%   To open the GUI, make sure the watershed_cells_gui folder is in your
%   MATLAB path and then run "watershed_cells_gui" from the command window.
%   Use the GUI that opens to:
%   (1) Select and import images via the "Select Images" panel
%   (2) Setup and run image segmentation via the "Segmentation" panel
%   (3) Setup and run region classification via the "Classification" panel
%   (4) Save the resulting data via the "Save Current Data" button
%   (5) Use the "Select Images" panel to select multiple image files and
%       then batch-process them with the current parameters using the
%       "Batch Process" button. This will loop through each listed image
%       file, load it, segment it, classify it, and save resulting data,
%       before moving on to the next image on the list.
%
% Usage:
%   h = watershed_cells_gui  
%       This command runs the GUI and returns the GUI figure handle h.
%       Within this handle h, you can access the current state of the GUI
%       data in the 'UserData' field. In particular, h.UserData is a struct
%       with the field 'params', which holds the current paramers from the
%       GUI and 'results', which holds the current analysis results. 
%       For example, run:
%           h.UserData.params.image.path
%       to return the list of image paths, or run:
%           h.UserData.results.segmentation.number
%       to return the total number of regions found during segmentation. 
%
% Output:
%   Click "Save Current Data" to prompt the user to select a folder. Then, 
%   inside this folder, the GUI will save:
%   (1) "<image name>_params.mat", a MATLAB file containing the structure
%       'params' with analysis parameters. This can be loaded into MATLAB
%       later for more detailed post-processing. 
%   (2) "<image name>_results.mat", a MATLAB file containing the structure
%       'resutls' with the analysis results. This can be loaded into MATLAB
%       later for more detailed post-processing.
%   (3) "<image name>_display.tif", a TIF image showing the original image
%       with segmented/classified regions outlined. This can be opened and
%       viewed to manually check the results.
%   Alternatively, running a batch process will automatically save the
%   above data in the same folder as the associated image file.
%
% See watershed_cells_gui README and documentation for more information.
%
% Requirements:
%   MATLAB R2014b or later
%   Image Processing Toolbox
%
% Tested with:
%   MATLAB R2016a & Image Processing Toolbox v9.4
%   MATLAB R2015b & Image Processing Toolbox v9.3
%   MATLAB R2014b & Image Processing Toolbox v9.1
%
% See also: imread adapthisteq medfilt2 imfilter fspecial multithresh
%           watershed bwconncomp bwperim regionprops labelmatrix

function varargout = watershed_cells_gui(varargin)
% WATERSHED_CELLS_GUI MATLAB code for watershed_cells_gui.fig
%      WATERSHED_CELLS_GUI, by itself, creates a new WATERSHED_CELLS_GUI or raises the existing
%      singleton*.
%
%      H = WATERSHED_CELLS_GUI returns the handle to a new WATERSHED_CELLS_GUI or the handle to
%      the existing singleton*.
%
%      WATERSHED_CELLS_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in WATERSHED_CELLS_GUI.M with the given input arguments.
%
%      WATERSHED_CELLS_GUI('Property','Value',...) creates a new WATERSHED_CELLS_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before watershed_cells_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to watershed_cells_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help watershed_cells_gui

% Last Modified by GUIDE v2.5 15-Nov-2016 23:40:09

%% ========== Begin initialization code - DO NOT EDIT ================== %%
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @watershed_cells_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @watershed_cells_gui_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
%% =========== End initialization code - DO NOT EDIT ABOVE ============= %%

%% ===================== Other functions =========================== %%

%% --------------------- General purpose

function objs = disable_gui(figure_handle)
% set objects in figure_handle as inactive (or off, if inactive is not
% valid) and return a list of disabled object handles
% but never disable the "cancel" button

objs = findobj(figure_handle, 'Enable', 'on', '-not', 'Tag', 'cancelbatch');

for ii = 1:length(objs)
    try
        set(objs(ii), 'Enable', 'inactive');
    catch
        set(objs(ii), 'Enable', 'off');
    end
end

function handles = reset_plots(handles, seg_im, seg_edges, state_im, state_edges, hist)
% set all plot data to be blank again

if seg_im
    % grayscale image
    handles.plots.grayim.CData = [];
end

% region outlines 
if seg_edges
    handles.plots.edges.CData = [];
    handles.numregions.String = sprintf('Number of regions: %d', []);
end

% raw image 
if state_im
    handles.plots.colorim.CData = [];
end

% region outlines by state image
if state_edges
    handles.plots.state1.CData = [];
    handles.plots.state2.CData = [];
    handles.numclass1.String = sprintf('State 1: %d', []);
    handles.numclass2.String = sprintf('State 2: %d', []);
end

% histogram of function value in each region
if hist
    handles.plots.hist.XData = 0;
    handles.plots.hist.YData = 0;
    handles.plots.hist.BarWidth = 0;
    handles.plots.thresh.XData = 0;
    handles.plots.thresh.YData = 0;
end

function handles = update_listbox(handles)

% udpate strings

if isempty(handles.output.UserData.params.image.path)
    % if we are out of images in the list... 
    % update the display
    handles.imagelistbox.String = sprintf('%s', '');

else
    parta = sprintf('%s\n', handles.output.UserData.params.image.path{1:end-1});
    partb = sprintf('%s', handles.output.UserData.params.image.path{end});
    handles.imagelistbox.String = sprintf('%s%s', parta, partb);
end

% update selection location to the clicked value 
handles.imagelistbox.Value = handles.output.UserData.params.image.index;

%% --------------------- Opening & Closing

function watershed_cells_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% --- Executes just before watershed_cells_gui is made visible.
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to watershed_cells_gui (see VARARGIN)
% UIWAIT makes watershed_cells_gui wait for user response (see UIRESUME)
% uiwait(handles.watershed_cells_gui);


% Default command line output is the GUI figure handle
handles.output = hObject;

% Store processing parameters and data in the figure handle's UserData, so
% it is accessible dynamically. Also store data to be displayed in the
% handles structure
[handles.output.UserData, handles.displaydata] = initialize_data();

% setup the listbox
handles = update_listbox(handles);

% Initialize the plots, create handles (with empty data) to all the thing
% swe will be plotting
handles.plots = initialize_display(...
    handles.segmentation_axes, handles.classify_axes, handles.hist_axes);

% Update gui data
guidata(hObject, handles)

function [data, displaydata] = initialize_data()

% data structure with parameters and results
data = struct(...
    'params', struct(... % default parameters struct
        'image', struct('path', {{}}, 'color', true, 'index', 1),...
        'segmentation', struct('equalization_cliplim', struct('value', 0.01, 'on', true),...
                               'background_size', struct('value', 19, 'on', true),...
                               'median_size', struct('value', 7, 'on', true),...
                               'gaussian_sigma', struct('value', 7, 'on', true),...
                               'minimum_area', struct('value', 35, 'on', true),...
                               'maximum_area', struct('value', 2000, 'on', true),...
                               'minimum_signal', struct('value', 0.2, 'on', true)...
                               ),...
        'classification', struct('f', struct('str', 'mean(R)-mean(G)', 'fun', @(R,G,B)mean(R)-mean(G), 'new', true),...
                                 'threshold', struct('value', 0, 'automatic', false, 'new', true)...
                                 )...
        ),...
    'results', struct(... % empty results struct
        'segmentation', struct('label_matrix', [],...
                               'number', []...
                               ),...
        'classification', struct('state1', struct('bw', [], 'number', []),...
                                 'state2', struct('bw', [], 'number', []),...
                                 'f_values', []...
                                 ),...
        'cancelbatch', false)...
    );

% empty display data struct
displaydata = struct(...
    'image', struct('color', [], 'gray', []),...
    'edges', struct('all', [], 'state1', [], 'state2', [])...
    );

function plots = initialize_display(segmentation_axes, classify_axes, hist_axes)
% initialize the axes & image data, return structure of plot/image handles

% hold on and link axes so they zoom/pan together
segmentation_axes.NextPlot = 'add';
classify_axes.NextPlot = 'add';
hist_axes.NextPlot = 'add';
linkaxes([segmentation_axes classify_axes]);

% grayscale image handles
plots.grayim = imshow([], 'parent', segmentation_axes);
plots.grayim.Tag = 'gray';

% region outlines image
plots.edges = imshow([], 'parent', segmentation_axes);
plots.edges.Tag = 'edges';

% raw image handles
plots.colorim = imshow([], 'parent', classify_axes);
plots.colorim.Tag = 'color';

% region outlines by state image
plots.state1 = imshow([], 'parent', classify_axes);
plots.state1.Tag = 'state1';
plots.state2 = imshow([], 'parent', classify_axes);
plots.state2.Tag = 'state2';

% histogram of function value in each region
hist_axes.AmbientLightColor = [0 0 0];
plots.hist = bar(hist_axes, 0, 0, 0, 'edgecolor', [1 1 1]*0.5, 'facecolor', [1 1 1]*0);
plots.hist.Tag = 'histogram';
plots.thresh = plot(hist_axes, 0, 0, 'color', [1 0 0]*0.75, 'linewidth', 1.5);
plots.thresh.Tag = 'thresh';
hist_axes.XLabel.String = 'f(R,G,B)';
hist_axes.YLabel.String = 'count';

function varargout = watershed_cells_gui_OutputFcn(hObject, eventdata, handles) 
% --- Outputs from this function are returned to the command line.
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function watershed_cells_gui_CloseRequestFcn(hObject, eventdata, handles)
% --- Executes when user attempts to close watershed_cells_gui.
% hObject    handle to watershed_cells_gui (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% ask the user if they really want to quit
button = questdlg('Are you sure you want to quit? Unsaved data will be lost.','Close Watershed Cells GUI','Quit','Cancel', 'Quit');

switch button
    case 'Quit'
        delete(hObject);
end

%% --------------------- Callbacks - List and Import Images

function paths = add_file(paths, fullpath)

% check if this file is already in the list
match = false;
for ii = 1:length(paths)
    if strcmp(fullpath, paths{ii});
        match = true;
        break
    end
end

% if not, add it to the list
if ~match
    paths{end+1} = fullpath;
end
        
function addbutton_Callback(hObject, eventdata, handles)
% --- Executes on button press in addbutton.
% hObject    handle to addbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Temporarily disable active gui components
objs = disable_gui(handles.watershed_cells_gui);
drawnow

% list of current file paths
paths = handles.output.UserData.params.image.path;

% get one or more file path(s)
[filename, pathname, success] = uigetfile({'*.tif;*.tiff'}, 'select image(s) to analyze', 'MultiSelect', 'on');
if success && ~isa(filename, 'cell')
    filename = {filename};
end

% if we got one or more files
if success
    
    % construct the full path to each file and add it to the list
    for ff = 1:length(filename)
        fullpath = fullfile(pathname, filename{ff});
        paths = add_file(paths, fullpath);
    end
    
    % Update the display & ui data & select the last path
    handles.output.UserData.params.image.index = max(1, length(paths));
    handles.output.UserData.params.image.path = paths;
    handles = update_listbox(handles);
    
    % Update handles structure
%     guidata(hObject, handles);
end

% Re-enable disabled gui components
set(objs, 'Enable', 'on');
drawnow

function removebutton_Callback(hObject, eventdata, handles)
% --- Executes on button press in removebutton.
% hObject    handle to removebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

paths = handles.output.UserData.params.image.path;
ix = handles.output.UserData.params.image.index;

% number of paths currently
maxix = length(paths);

% delete the selected path
if ix <= maxix
    paths(ix) = [];
    maxix = length(paths);
end

% move selection to the end of the list, if we are past the end
if ix > maxix && maxix > 0
    ix = maxix;
end

% reset path and index parameters
handles.output.UserData.params.image.path = paths;
handles.output.UserData.params.image.index = ix;

% update the listbox
handles = update_listbox(handles);

function imagelistbox_Callback(hObject, eventdata, handles)
% --- Executes on selection change in imagelistbox.
% hObject    handle to imagelistbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: contents = cellstr(get(hObject,'String')) returns imagelistbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from imagelistbox

% store the current selection index
handles.output.UserData.params.image.index = hObject.Value;

function [valid, iscolor] = check_image(rawim)
% check if valid color or grayscale image

if size(rawim, 3) == 3
    iscolor = true;
    valid = true;
    
elseif size(rawim, 3) == 1
    iscolor = false;
    valid = true;
    
else
    iscolor = false;
    valid = false;
end

function handles = show_images(handles)
% show color and grayscale images

% gray image
grayim = handles.displaydata.image.gray;
handles.plots.grayim.CData = grayim;
if ~isempty(grayim)
    handles.segmentation_axes.XLim = [0 size(grayim, 2)]+0.5;
    handles.segmentation_axes.YLim = [0 size(grayim, 1)]+0.5;
end

% color image
colorim = handles.displaydata.image.color;
handles.plots.colorim.CData = colorim;
if ~isempty(colorim)
    handles.classify_axes.XLim = [0 size(colorim, 2)]+0.5;
    handles.classify_axes.YLim = [0 size(colorim, 1)]+0.5;
end

function handles = loadimagebutton_Callback(hObject, eventdata, handles)
% --- Executes on button press in loadimagebutton.
% hObject    handle to loadimagebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Show user we are working on stuff and temporarily disable ui objects
handles.computingsegmentation.Visible = 'on';
handles.computingclassification.Visible = 'on';
objs = disable_gui(handles.watershed_cells_gui);
drawnow

% get the image patch 
if isempty(handles.output.UserData.params.image.path)
    % if we ask to load from an empty list, send a dummy/empty path to stop
    % the analysis
    impath = [];
    
    % also reset the results data and plots
    tmp = initialize_data();
    handles.output.UserData.results = tmp.results;
    handles = reset_plots(handles, 1, 1, 1, 1, 1);
else
    % get the currently-selected path
    impath = handles.output.UserData.params.image.path{...
        handles.output.UserData.params.image.index};
end

% check that the file exists
if exist(impath, 'file') 

    % import the image
    valid = false;
    try
        % import the image and convert to double [0 1] 
        rawim = im2double(imread(impath));

        % check the image for valid color channels
        [valid, iscolor] = check_image(rawim);

    catch err
        warning(err.message);
    end

    if valid

        % create grayscale image
        if iscolor
            colorim = rawim;
            grayim = rgb2gray(rawim);
        else
            colorim = cat(3, rawim, rawim, rawim);
            grayim = rawim;
        end

        % store the display info
        handles.output.UserData.params.image.color = iscolor;
        handles.displaydata.image.color = colorim;
        handles.displaydata.image.gray = grayim;
        
        % reset 'results' data
        tmp = initialize_data();
        handles.output.UserData.results = tmp.results;
        
        % update plots
        handles = reset_plots(handles, 1, 1, 1, 1, 1);
        handles = show_images(handles);
        
        % Update handles structure
        guidata(hObject, handles);

    else
        warning('Invalid image.')
    end
end

% Show user we are done working on stuff and enable ui objects
set(objs, 'Enable', 'on');
handles.computingsegmentation.Visible = 'off';
handles.computingclassification.Visible = 'off';
drawnow

%% --------------------- Callbacks - Segmentation

function sz_equalization_Callback(hObject, eventdata, handles)
% hObject    handle to sz_equalization (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of sz_equalization as text
%        str2double(get(hObject,'String')) returns contents of sz_equalization as a double

% get new value
val = str2double(hObject.String);

% check/correct value, positive number in the range [0 1]
val = min(abs(val), 1);

% reset string
hObject.String = num2str(val);

% save new value into parameters
handles.output.UserData.params.segmentation.equalization_cliplim.value = val;

function run_equalization_Callback(hObject, eventdata, handles)
% --- Executes on button press in run_equalization.
% hObject    handle to run_equalization (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of run_equalization

% save new value into parameters
handles.output.UserData.params.segmentation.equalization_cliplim.on = hObject.Value;

function sz_background_Callback(hObject, eventdata, handles)
% hObject    handle to sz_background (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of sz_background as text
%        str2double(get(hObject,'String')) returns contents of sz_background as a double

% get new value
val = str2double(hObject.String);

% check/correct value, positive odd integer
val = floor(abs(val)/2)*2+1;

% reset string
hObject.String = num2str(val);

% save new value into parameters
handles.output.UserData.params.segmentation.background_size.value = val;

function run_background_Callback(hObject, eventdata, handles)
% --- Executes on button press in run_background.
% hObject    handle to run_background (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of run_background

% save new value into parameters
handles.output.UserData.params.segmentation.background_size.on = hObject.Value;

function sz_median_Callback(hObject, eventdata, handles)
% hObject    handle to sz_median (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of sz_median as text
%        str2double(get(hObject,'String')) returns contents of sz_median as a double

% get new value
val = str2double(hObject.String);

% check/correct value, positive odd integer
val = floor(abs(val)/2)*2+1;

% reset string
hObject.String = num2str(val);

% save new value into parameters
handles.output.UserData.params.segmentation.median_size.value = val;

function run_median_Callback(hObject, eventdata, handles)
% --- Executes on button press in run_median.
% hObject    handle to run_median (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of run_median

% save new value into parameters
handles.output.UserData.params.segmentation.median_size.on = hObject.Value;

function sz_gaussian_Callback(hObject, eventdata, handles)
% hObject    handle to sz_gaussian (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of sz_gaussian as text
%        str2double(get(hObject,'String')) returns contents of sz_gaussian as a double

% get new value
val = str2double(hObject.String);

% check/correct value, positive 
val = abs(val);

% reset string
hObject.String = num2str(val);

% save new value into parameters
handles.output.UserData.params.segmentation.gaussian_sigma.value = val;

function run_gaussian_Callback(hObject, eventdata, handles)
% --- Executes on button press in run_gaussian.
% hObject    handle to run_gaussian (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of run_gaussian

% save new value into parameters
handles.output.UserData.params.segmentation.gaussian_sigma.on = hObject.Value;

function sz_minarea_Callback(hObject, eventdata, handles)
% hObject    handle to sz_minarea (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of sz_minarea as text
%        str2double(get(hObject,'String')) returns contents of sz_minarea as a double

% get new value
val = str2double(hObject.String);

% check/correct value, positive integer
val = round(abs(val));

% reset string
hObject.String = num2str(val);

% save new value into parameters
handles.output.UserData.params.segmentation.minimum_area.value = val;

function run_minarea_Callback(hObject, eventdata, handles)
% --- Executes on button press in run_minarea.
% hObject    handle to run_minarea (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of run_minarea

% save new value into parameters
handles.output.UserData.params.segmentation.minimum_area.on = hObject.Value;

function sz_maxarea_Callback(hObject, eventdata, handles)
% hObject    handle to sz_maxarea (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of sz_maxarea as text
%        str2double(get(hObject,'String')) returns contents of sz_maxarea as a double

% get new value
val = str2double(hObject.String);

% check/correct value, positive integer
val = round(abs(val));

% reset string
hObject.String = num2str(val);

% save new value into parameters
handles.output.UserData.params.segmentation.maximum_area.value = val;

function run_maxarea_Callback(hObject, eventdata, handles)
% --- Executes on button press in run_maxarea.
% hObject    handle to run_maxarea (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of run_maxarea

% save new value into parameters
handles.output.UserData.params.segmentation.maximum_area.on = hObject.Value;

function sz_minsignal_Callback(hObject, eventdata, handles)
% hObject    handle to sz_minsignal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of sz_minsignal as text
%        str2double(get(hObject,'String')) returns contents of sz_minsignal as a double

% get new value
val = str2double(hObject.String);

% check/correct value, positive number in the range [0 1]
val = min(abs(val), 1);

% reset string
hObject.String = num2str(val);

% save new value into parameters
handles.output.UserData.params.segmentation.minimum_signal.value = val;

function run_minsignal_Callback(hObject, eventdata, handles)
% --- Executes on button press in run_minsignal.
% hObject    handle to run_minsignal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of run_minsignal

% save new value into parameters
handles.output.UserData.params.segmentation.minimum_signal.on = hObject.Value;

function handles = segmentationbutton_Callback(hObject, eventdata, handles)
% --- Executes on button press in segmentationbutton.
% hObject    handle to segmentationbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Show user we are working on stuff and temporarily disable ui objects
handles.computingsegmentation.Visible = 'on';
objs = disable_gui(handles.watershed_cells_gui);
drawnow

% run & store segmentation
if ~isempty(handles.displaydata.image.color)
    
    [label_matrix, edges] = find_regions(handles.displaydata.image.color, ...
        handles.output.UserData.params.segmentation);
    handles.output.UserData.results.segmentation.number = max(label_matrix(:));
    handles.output.UserData.results.segmentation.label_matrix = label_matrix;
    handles.displaydata.edges.all = edges;

    % Reset classification calculations
    handles.output.UserData.params.classification.f.new = true;
    handles.output.UserData.params.classification.threshold.new = true;

    % Reset classification results and Display the segmentation results
    tmp = initialize_data();
    handles.output.UserData.results.classification = tmp.results.classification;
    handles.output.UserData.params.classification.threshold.new = true;
    handles.output.UserData.params.classification.f.new = true;
    handles = reset_plots(handles, 0, 0, 0, 1, 1);
    handles = show_segmentation(handles);
    
else
    warning('No valid image data. Load image before running segmentation.')
end

% Update gui data
guidata(hObject, handles)

% Show user we are done working on stuff and enable ui objects
set(objs, 'Enable', 'on');
handles.computingsegmentation.Visible = 'off';
drawnow

function handles = show_segmentation(handles)

% show segmentation with yellow outlines
edges = handles.displaydata.edges.all;
handles.plots.edges.CData = cat(3, edges, edges, edges*0);
handles.plots.edges.AlphaData = edges * 0.5;
handles.numregions.String = sprintf('Number of regions: %d', ...
    max(handles.output.UserData.results.segmentation.number));

%% --------------------- Callbacks - Classification

function classifyeq_Callback(hObject, eventdata, handles)
% hObject    handle to classifyeq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of classifyeq as text
%        str2double(get(hObject,'String')) returns contents of classifyeq as a double

% get new string
str = hObject.String;

% construct and check equation
try
    fun = eval(sprintf('@(R,G,B)(%s)', str));
    
    if ~isa(fun, 'function_handle')
        str = 'NaN';
        fun = @(R,G,B) NaN;
        warning('Invalid function.')
    end
    
catch
    str = 'NaN';
    fun = @(R,G,B) NaN;
    warning('Invalid function.')
end

% reset string
hObject.String = str;

% save new value into parameters
handles.output.UserData.params.classification.f.str = str;
handles.output.UserData.params.classification.f.fun = fun;
handles.output.UserData.params.classification.f.new = true;

function threshold_Callback(hObject, eventdata, handles)
% hObject    handle to threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of threshold as text
%        str2double(get(hObject,'String')) returns contents of threshold as a double

% get new value
val = str2double(hObject.String);

% reset string
hObject.String = num2str(val);

% save new value into parameters
handles.output.UserData.params.classification.threshold.value = val;
handles.output.UserData.params.classification.threshold.new = true;

function otsuthresh_Callback(hObject, eventdata, handles)
% --- Executes on button press in otsuthresh.
% hObject    handle to otsuthresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of otsuthresh

% save new value into parameters
handles.output.UserData.params.classification.threshold.automatic = hObject.Value;
handles.output.UserData.params.classification.threshold.new = true;

% Set the threshold box to inactive
if hObject.Value
    handles.threshold.Enable = 'inactive';
else
    handles.threshold.Enable = 'on';
end

function handles = classificationbutton_Callback(hObject, eventdata, handles)
% --- Executes on button press in classificationbutton.
% hObject    handle to classificationbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Show user we are working on stuff and temporarily disable ui objects
handles.computingclassification.Visible = 'on';
objs = disable_gui(handles.watershed_cells_gui);
drawnow

% Check - have we applied this function and/or threshold already? (if so,
% we can skip the calculation)
fnew = handles.output.UserData.params.classification.f.new;
tnew = handles.output.UserData.params.classification.threshold.new;

% apply function, if we haven't already
if fnew
    if ~isempty(handles.displaydata.image.color) && ...
            ~isempty(handles.output.UserData.results.segmentation.label_matrix)
        
        % apply function
        handles.output.UserData.results.classification.f_values = ...
            apply_function(handles.displaydata.image.color,...
                handles.output.UserData.results.segmentation.label_matrix, ...
                handles.output.UserData.params.classification.f.fun);

        % mark that we've now applied the function (it's no longer new)
        handles.output.UserData.params.classification.f.new = false;
        
    else
        warning('No valid image and/or segmentation data. Load image and run segmentation before classification.')
    end
end
    
% threshold function values
if fnew || tnew
    
    % apply and store threshold
    [bw1, bw2, n1, n2, thresh] = apply_threshold(...
        handles.output.UserData.results.segmentation.label_matrix,...
        handles.output.UserData.results.classification.f_values, ...
        handles.output.UserData.params.classification.threshold.value,...
        handles.output.UserData.params.classification.threshold.automatic);
        
    % store the results
    handles.output.UserData.results.classification.state1.bw = bw1;
    handles.output.UserData.results.classification.state2.bw = bw2;
    handles.output.UserData.results.classification.state1.number = n1;
    handles.output.UserData.results.classification.state2.number = n2; 
    handles.output.UserData.params.classification.threshold.value = thresh;
    
    % update the display data
    handles.threshold.String = num2str(thresh);
    handles.displaydata.edges.state1 = bwperim(bw1);
    handles.displaydata.edges.state2 = bwperim(bw2);
        
    % mark that we've now applied the threshold (it's no longer new)
    handles.output.UserData.params.classification.threshold.new = false;
end

% update classification display
if fnew || tnew
    handles = show_classification(handles);
end

% Update gui data
guidata(hObject, handles)

% Show user we are done working on stuff and enable ui objects
handles.computingclassification.Visible = 'off';
set(objs, 'Enable', 'on');
drawnow

function handles = show_classification(handles)

% display output histogram
values = handles.output.UserData.results.classification.f_values;
nbins = max(10, round(length(values)/100)*5);
[N, binedges] = histcounts(values, nbins);
handles.plots.hist.YData = N;
handles.plots.hist.XData = binedges(1:end-1) + diff(binedges)/2;
handles.plots.hist.BarWidth = 1;

% display threshold
handles.plots.thresh.XData = [1 1] * handles.output.UserData.params.classification.threshold.value;
handles.plots.thresh.YData = [0 1.1] * max(N);

% display new output edge images for state1 and state2
edges1 = handles.displaydata.edges.state1;
edges2 = handles.displaydata.edges.state2;
handles.plots.state1.CData = cat(3, edges1, edges1*0, edges1); % state 1 has magenta outlines
handles.plots.state1.AlphaData = edges1;
handles.plots.state2.CData = cat(3, edges2*0, edges2, edges2); % state 2 has cyan outlines
handles.plots.state2.AlphaData = edges2;

% display state counts
handles.numclass1.String = sprintf('State 1: %d', handles.output.UserData.results.classification.state1.number);
handles.numclass2.String = sprintf('State 2: %d', handles.output.UserData.results.classification.state2.number);

%% --------------------- Callbacks - Save Data

function handles = savedatabutton_Callback(hObject, eventdata, handles, varargin)
% --- Executes on button press in savedatabutton.
% hObject    handle to savedatabutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Show user we are working on stuff and temporarily disable ui objects
handles.computingsegmentation.Visible = 'on';
handles.computingclassification.Visible = 'on';
objs = disable_gui(handles.watershed_cells_gui);
drawnow

% If there is a 4th argument, it is the basename for saving. If this is not
% given, asking the user to select a save location
if nargin > 3
    % Use the given basename
    basename = varargin{1};
    
else
    % bring up a dialog box to pick the save name/location
    PathName = uigetdir('*.*', 'Select folder to save results');
    
    % construct the basename, using the image filename as part of the base
    if PathName
        [~,basename,~] = fileparts(handles.output.UserData.params.image.path{handles.output.UserData.params.image.index});
        if isempty(basename)
            basename = 'output';
        end
        basename = fullfile(PathName, basename);
    end
end

% save the data and classification display image
if basename
        
    % save params and results as mat files
    params = handles.output.UserData.params;
    results = handles.output.UserData.results;
    save([basename '_params.mat'], 'params', '-mat');
    save([basename '_results.mat'], 'results', '-mat');
    
    % save image showing the classification output
    im = handles.displaydata.image.color;
    if ~isempty(im)
        imR = im(:,:,1);
        imG = im(:,:,2);
        imB = im(:,:,3);
        imR(handles.displaydata.edges.all) = 1; % all outlines are yellow
        imG(handles.displaydata.edges.all) = 1;
        imB(handles.displaydata.edges.all) = 0;
        imR(handles.displaydata.edges.state1) = 1; % state 1 outlines are magenta
        imG(handles.displaydata.edges.state1) = 0;
        imB(handles.displaydata.edges.state1) = 1;
        imR(handles.displaydata.edges.state2) = 0; % state 2 outlines are cyan
        imG(handles.displaydata.edges.state2) = 1;
        imB(handles.displaydata.edges.state2) = 1;
        im = cat(3, imR, imG, imB);
        imwrite(im, [basename '_display.tif']) 
    end
    
    % Print confirmation to the screen
    fprintf('Data saved. Basename: %s\n', basename);
end

% Show user we are working on stuff and enable ui objects
handles.computingsegmentation.Visible = 'off';
handles.computingclassification.Visible = 'off';
set(objs, 'Enable', 'on');
drawnow

%% --------------------- Callbacks - Batch Process

function batchprocessbutton_Callback(hObject, eventdata, handles)
% --- Executes on button press in batchprocessbutton.
% hObject    handle to batchprocessbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Temporarily disable ui objects
objs = disable_gui(handles.watershed_cells_gui);
drawnow

% check if the user wants to proceed
button = questdlg('Batch process all images using the current parameters? Results will be saved in the same folder as the associated image.',...
    'Run Batch Process', 'OK', 'Cancel', 'Cancel');

% if we want to proceed
if strcmp(button, 'OK')
    
    % show user we are working
    handles.computingbatch.Visible = 'on';
    handles.cancelbatch.Visible = 'on';
    
    % process each image
    num_images = length(handles.output.UserData.params.image.path);
    for ii = 1:num_images
        
        % Select the image & update the listbox
        if handles.output.UserData.results.cancelbatch, break, end;
        handles.output.UserData.params.image.index = ii;
        handles = update_listbox(handles);
        impath = handles.output.UserData.params.image.path{ii};
        
        % Load the image
        if handles.output.UserData.results.cancelbatch, break, end;
        handles = loadimagebutton_Callback(handles.loadimagebutton, eventdata, handles);
        
        % Run the segmentation
        if handles.output.UserData.results.cancelbatch, break, end;
        handles = segmentationbutton_Callback(handles.segmentationbutton, eventdata, handles);
        
        % Run the classification
        if handles.output.UserData.results.cancelbatch, break, end;
        handles = classificationbutton_Callback(handles.classificationbutton, eventdata, handles);
        
        % Construct base filename for saving (same folder as the image
        % being processed)
        if handles.output.UserData.results.cancelbatch, break, end;
        [pathstring, filename, ~] = fileparts(impath); 
        basename = fullfile(pathstring, filename);
        
        % Save the data and Print a mini report to the command line
        if handles.output.UserData.results.cancelbatch, break, end;
        handles = savedatabutton_Callback(handles.savedatabutton, eventdata, handles, basename);
        fprintf('%s, Image %d of %d, Filename: %s, Regions: %d, State 1: %d, State 2: %d\n',...
            datestr(now), ii, num_images, impath, ...
            handles.output.UserData.results.segmentation.number,...
            handles.output.UserData.results.classification.state1.number,...
            handles.output.UserData.results.classification.state2.number)
        
    end
    
    % show user we are done working
    handles.computingbatch.Visible = 'off';
    handles.cancelbatch.Visible = 'off';
    handles.output.UserData.results.cancelbatch = false;
    
end

% Update gui data
guidata(hObject, handles)

% Enable ui objects 
set(objs, 'Enable', 'on');
drawnow

function cancelbatch_Callback(hObject, eventdata, handles)
% --- Executes on button press in cancelbatch.
% hObject    handle to cancelbatch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)cancelbatch


handles.output.UserData.results.cancelbatch = true;

% Update gui data
guidata(hObject, handles)














