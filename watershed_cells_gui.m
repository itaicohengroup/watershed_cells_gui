% WATERSHED_CELLS_GUI finds cells in an image using watershed segmentation
%
% Programmed by: Lena Bartell (lrb89@cornell.edu)
% Principle Investigator: Prof. Itai Cohen, Cornell University
%
% Quick-start guide:
%   To open the GUI, make sure the watershed_cells_gui folder is in your
%   MATLAB path and then run "watershed_cells_gui" from the command window.
%   Use the GUI that opens to:
%   (1) Read in an image via the "Import Images" panel
%   (2) Setup and run image segmentation via the "segmentation" panel
%   (3) Setup and run region classification via the "Classification" panel
%   (4) Save the resulting data via the "Save Data" button
%
% Usage:
%   h = watershed_cells_gui  
%       This command runs the GUI and returns the GUI figure handle h.
%       Within this handle h, you can access the current state of the GUI
%       data in the 'UserData' field. In particular, h.UserData is a struct
%       with the field 'params', which holds the current paramers from the
%       GUI and 'results', which holds the current analysis results. 
%       For example, run
%           h.UserData.params.image.path
%       to return the path to the image being analyzed, or run
%           h.UserData.results.segmentation.number
%       to return the total number of regions found during segmentation. 
%
% Output:
%   Click "Save Data" to prompt the user to select a folder. Then, inside
%   the folder, the GUI will save:
%   (1) "<image name>_params.mat", a MATLAB file containing the structure
%       'params' with analysis parameters. This can be loaded into MATLAB
%       later for more detailed post-processing. 
%   (2) "<image name>_results.mat", a MATLAB file containing the structure
%       'resutls' with the analysis results. This can be loaded into MATLAB
%       later for more detailed post-processing.
%   (3) "<image name>_display.tif", a TIF image showing the original image
%       with segmented/classified regions outlined. This can be opened and
%       viewed to manually check the results.
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

% Last Modified by GUIDE v2.5 14-Nov-2016 22:32:30

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

function push_data(hObject, handles)

handles.watershed_cells_gui.UserData = handles.data;
guidata(hObject, handles);

function objs = disable_gui(figure_handle)
% set objects in figure_handle as inactive (or off, if inactive is not
% valid) and return a list of disabled object handles

objs = findobj(figure_handle, 'Enable', 'on');
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
    handles.numclass1.Visible = 'off';
    handles.numclass2.Visible = 'off';
end

% histogram of function value in each region
if hist
    handles.plots.hist.XData = 0;
    handles.plots.hist.YData = 0;
    handles.plots.hist.BarWidth = 0;
    handles.plots.thresh.XData = 0;
    handles.plots.thresh.YData = 0;
end

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

% get default data structures
[handles.data, handles.displaydata] = initialize_data();

% Override default command line output
handles.output = hObject;

% Initialize the plots
handles.plots = initialize_display(...
    handles.segmentation_axes, handles.classify_axes, handles.hist_axes);

% Initialize the batch process
handles.batch = [];

% Update gui data
push_data(hObject, handles)

function [data, displaydata] = initialize_data()

% data structure with parameters and results
data = struct(...
    'params', struct(... % default parameters struct
        'image', struct('path', '', 'color', true),...
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
                                 )...
        )...
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

function watershed_cells_gui_DeleteFcn(hObject, eventdata, handles)
% --- Executes during object deletion, before destroying properties.
% hObject    handle to watershed_cells_gui (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes when user attempts to close watershed_cells_gui.
function watershed_cells_gui_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to watershed_cells_gui (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% ask the user if they really want to quit
button = questdlg('Are you sure you want to quit?','Close GUI','Yes','No', 'Yes');

switch button
    case 'Yes'
        if isa(handles.batch, 'handle')
            delete(handles.batch)
        end
        delete(hObject);
    otherwise
end

%% --------------------- Callbacks - Import Image

function browsebutton_Callback(hObject, eventdata, handles)
% --- Executes on button press in browsebutton.
% hObject    handle to browsebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Show user we are working on stuff and temporarily disable ui objects
handles.computingsegmentation.Visible = 'on';
handles.computingclassification.Visible = 'on';
objs = disable_gui(handles.watershed_cells_gui);
drawnow

% get the file path
[filename, pathname, success] = uigetfile({'*.tif;*.tiff'}, 'select image to analyze');
if success
    
    % construct the full path
    impath = [pathname filename];
    
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
            
            % store and display image info
            handles.data.params.image.path = impath;
            handles.data.params.image.color = iscolor;
            handles.displaydata.image.color = colorim;
            handles.displaydata.image.gray = grayim;
            handles.pathtoimage.String = impath;
            
            % reset results and update plots
            tmp = initialize_data();
            handles.data.results = tmp.results;
            handles.data.params.classification.threshold.new = true;
            handles.data.params.classification.f.new = true;
            handles = reset_plots(handles, 1, 1, 1, 1, 1);
            handles = show_images(handles);
            
        else
            warning('Invalid image. Skipping.')
        end
    end
end

% Reset classification calculations
handles.data.params.classification.f.new = true;
handles.data.params.classification.threshold.new = true;

% Update gui data
push_data(handles.output, handles)

% Show user we are working on stuff and enable ui objects
handles.computingsegmentation.Visible = 'off';
handles.computingclassification.Visible = 'off';
set(objs, 'Enable', 'on');
drawnow

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
handles.data.params.segmentation.equalization_cliplim.value = val;

% Update gui data
push_data(handles.output, handles)

function run_equalization_Callback(hObject, eventdata, handles)
% --- Executes on button press in run_equalization.
% hObject    handle to run_equalization (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of run_equalization

% save new value into parameters
handles.data.params.segmentation.equalization_cliplim.on = hObject.Value;

% Update gui data
push_data(handles.output, handles)

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
handles.data.params.segmentation.background_size.value = val;

% Update gui data
push_data(handles.output, handles)

function run_background_Callback(hObject, eventdata, handles)
% --- Executes on button press in run_background.
% hObject    handle to run_background (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of run_background

% save new value into parameters
handles.data.params.segmentation.background_size.on = hObject.Value;

% Update gui data
push_data(handles.output, handles)

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
handles.data.params.segmentation.median_size.value = val;

% Update gui data
push_data(handles.output, handles)

function run_median_Callback(hObject, eventdata, handles)
% --- Executes on button press in run_median.
% hObject    handle to run_median (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of run_median

% save new value into parameters
handles.data.params.segmentation.median_size.on = hObject.Value;

% Update gui data
push_data(handles.output, handles)

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
handles.data.params.segmentation.gaussian_sigma.value = val;

% Update gui data
push_data(handles.output, handles)

function run_gaussian_Callback(hObject, eventdata, handles)
% --- Executes on button press in run_gaussian.
% hObject    handle to run_gaussian (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of run_gaussian

% save new value into parameters
handles.data.params.segmentation.gaussian_sigma.on = hObject.Value;

% Update gui data
push_data(handles.output, handles)

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
handles.data.params.segmentation.minimum_area.value = val;

% Update gui data
push_data(handles.output, handles)

function run_minarea_Callback(hObject, eventdata, handles)
% --- Executes on button press in run_minarea.
% hObject    handle to run_minarea (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of run_minarea

% save new value into parameters
handles.data.params.segmentation.minimum_area.on = hObject.Value;

% Update gui data
push_data(handles.output, handles)

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
handles.data.params.segmentation.maximum_area.value = val;

% Update gui data
push_data(handles.output, handles)

function run_maxarea_Callback(hObject, eventdata, handles)
% --- Executes on button press in run_maxarea.
% hObject    handle to run_maxarea (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of run_maxarea

% save new value into parameters
handles.data.params.segmentation.maximum_area.on = hObject.Value;

% Update gui data
push_data(handles.output, handles)

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
handles.data.params.segmentation.minimum_signal.value = val;

% Update gui data
push_data(handles.output, handles)

function run_minsignal_Callback(hObject, eventdata, handles)
% --- Executes on button press in run_minsignal.
% hObject    handle to run_minsignal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of run_minsignal

% save new value into parameters
handles.data.params.segmentation.minimum_signal.on = hObject.Value;

% Update gui data
push_data(handles.output, handles)

function segmentationbutton_Callback(hObject, eventdata, handles)
% --- Executes on button press in segmentationbutton.
% hObject    handle to segmentationbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Show user we are working on stuff and temporarily disable ui objects
handles.computingsegmentation.Visible = 'on';
handles.computingclassification.Visible = 'on';
objs = disable_gui(handles.watershed_cells_gui);
drawnow

% run & store segmentation
if ~isempty(handles.displaydata.image.color)
    [label_matrix, edges] = ...
        find_regions(handles.displaydata.image.color, handles.data.params.segmentation);
    handles.data.results.segmentation.number = max(label_matrix(:));
    handles.data.results.segmentation.label_matrix = label_matrix;
    handles.displaydata.edges.all = edges;

    % Reset classification calculations
    handles.data.params.classification.f.new = true;
    handles.data.params.classification.threshold.new = true;

    % Reset classification results and Display the segmentation results
    tmp = initialize_data();
    handles.data.results.classification = tmp.results.classification;
    handles.data.params.classification.threshold.new = true;
    handles.data.params.classification.f.new = true;
    handles = reset_plots(handles, 0, 0, 0, 1, 1);
    handles = show_segmentation(handles);
    
else
    warning('No valid image data.')
end

% Update gui data
push_data(handles.output, handles)

% Show user we are done working on stuff and enable ui objects
handles.computingsegmentation.Visible = 'off';
handles.computingclassification.Visible = 'off';
set(objs, 'Enable', 'on');
drawnow

function handles = show_segmentation(handles)

% show segmentation with yellow outlines
edges = handles.displaydata.edges.all;
handles.plots.edges.CData = cat(3, edges, edges, edges*0);
handles.plots.edges.AlphaData = edges * 0.5;
handles.numregions.String = sprintf('Number of regions: %d', ...
    max(handles.data.results.segmentation.number));

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
handles.data.params.classification.f.str = str;
handles.data.params.classification.f.fun = fun;
handles.data.params.classification.f.new = true;

% Update gui data
push_data(handles.output, handles)

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
handles.data.params.classification.threshold.value = val;
handles.data.params.classification.threshold.new = true;

% Update gui data
push_data(handles.output, handles)

function otsuthresh_Callback(hObject, eventdata, handles)
% --- Executes on button press in otsuthresh.
% hObject    handle to otsuthresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of otsuthresh

% save new value into parameters
handles.data.params.classification.threshold.automatic = hObject.Value;
handles.data.params.classification.threshold.new = true;

% Set the threshold box to inactive
if hObject.Value
    handles.threshold.Enable = 'inactive';
else
    handles.threshold.Enable = 'on';
end

% Update gui data
push_data(handles.output, handles)

function classificationbutton_Callback(hObject, eventdata, handles)
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
fnew = handles.data.params.classification.f.new;
tnew = handles.data.params.classification.threshold.new;

% apply function, if we haven't already
if fnew
    if ~isempty(handles.displaydata.image.color) && ...
            ~isempty(handles.data.results.segmentation.label_matrix)
        
        % apply function
        handles.data.results.classification.f_values = ...
            apply_function(handles.displaydata.image.color,...
                handles.data.results.segmentation.label_matrix, ...
                handles.data.params.classification.f.fun);

        % mark that we've now applied the function (it's no longer new)
        handles.data.params.classification.f.new = false;
        
    else
        warning('No valid image and/or segmentation data.')
    end
end
    
% threshold function values
if fnew || tnew
    
    % apply and store threshold
    [bw1, bw2, n1, n2, thresh] = apply_threshold(...
        handles.data.results.segmentation.label_matrix,...
        handles.data.results.classification.f_values, ...
        handles.data.params.classification.threshold.value,...
        handles.data.params.classification.threshold.automatic);
        
    % store the results
    handles.data.results.classification.state1.bw = bw1;
    handles.data.results.classification.state2.bw = bw2;
    handles.data.results.classification.state1.number = n1;
    handles.data.results.classification.state2.number = n2; 
    handles.data.params.classification.threshold.value = thresh;
    
    % update the display data
    handles.threshold.String = num2str(thresh);
    handles.displaydata.edges.state1 = bwperim(bw1);
    handles.displaydata.edges.state2 = bwperim(bw2);
        
    % mark that we've now applied the threshold (it's no longer new)
    handles.data.params.classification.threshold.new = false;
end

% update classification display
if fnew || tnew
    handles = show_classification(handles);
end

% Update gui data
push_data(handles.output, handles)

% Show user we are done working on stuff and enable ui objects
handles.computingclassification.Visible = 'off';
set(objs, 'Enable', 'on');
drawnow

function handles = show_classification(handles)

% display output histogram
values = handles.data.results.classification.f_values;
nbins = max(10, round(length(values)/100)*5);
[N, binedges] = histcounts(values, nbins);
handles.plots.hist.YData = N;
handles.plots.hist.XData = binedges(1:end-1) + diff(binedges)/2;
handles.plots.hist.BarWidth = 1;

% display threshold
handles.plots.thresh.XData = [1 1] * handles.data.params.classification.threshold.value;
handles.plots.thresh.YData = [0 1.1] * max(N);

% display new output edge images for state1 and state2
edges1 = handles.displaydata.edges.state1;
edges2 = handles.displaydata.edges.state2;
handles.plots.state1.CData = cat(3, edges1, edges1*0, edges1); % state 1 has magenta outlines
handles.plots.state1.AlphaData = edges1;
handles.plots.state2.CData = cat(3, edges2*0, edges2, edges2); % state 2 has cyan outlines
handles.plots.state2.AlphaData = edges2;

% display state counts
handles.numclass1.String = sprintf('State 1: %d', handles.data.results.classification.state1.number);
handles.numclass1.Visible = 'on';
handles.numclass2.String = sprintf('State 2: %d', handles.data.results.classification.state2.number);
handles.numclass2.Visible = 'on';

%% --------------------- Callbacks - Save Data

function savedatabutton_Callback(hObject, eventdata, handles)
% --- Executes on button press in savedatabutton.
% hObject    handle to savedatabutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Show user we are working on stuff and temporarily disable ui objects
handles.computingsegmentation.Visible = 'on';
handles.computingclassification.Visible = 'on';
objs = disable_gui(handles.watershed_cells_gui);
drawnow

% bring up a dialog box to pick the save name/location
PathName = uigetdir('*.*', 'Select folder to save results');
if PathName
    
    % save results using the image filename as the base name
    [~,basename,~] = fileparts(handles.data.params.image.path);
    if isempty(basename)
        basename = 'output';
    end
    basename = fullfile(PathName, basename);
    
    % save params and results as mat files
    params = handles.data.params;
    results = handles.data.results;
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
end

% print some quick info to the screen
if ~isempty(handles.data.params.image.path) && ...
        ~isempty(handles.data.results.segmentation.number)
    fprintf('%s: %d regions (%d state 1, %d state 2)\n', ...
        handles.data.params.image.path,...
        handles.data.results.segmentation.number,...
        handles.data.results.classification.state1.number,...
        handles.data.results.classification.state2.number);
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


data.params = handles.watershed_cells_gui.UserData.params;
tmp = initialize_data();
data.results = tmp.results;


handles.batch = batch_process(data);

% Update gui data
push_data(handles.output, handles)

uiwait(handles.batch)

% Enable ui objects
set(findobj(objs), 'Enable', 'on');






