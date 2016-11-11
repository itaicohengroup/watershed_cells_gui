% WATERSHED_CELLS_GUI finds cells in an image using watershed segmentation
%
% Date created: October 26, 2016
% Programmed by: Lena Bartell (lrb89@cornell.edu)
% Principle Investigator: Prof. Itai Cohen, Cornell University
%
% Quick-start guide
%   Top open the GUI, make sure the watershed_Cells_gui folder is in your
%   MATLAB path and then run "watershed_cells_gui". Use the GUI to set the
%   parameters then click "Find Cells". This will run the analysis and
%   display the result. To skip a step in the analysis process, un-check
%   the associated box. If you change the cell outline alpha value, click
%   "Update Display" to show the change. To save the resulting data, click
%   "Save Data".
%
% Output
%   "Save Data" will prompt the user to save the structure "data" to a
%   ".mat" file. To view the data, load the .mat file into Matlab. The
%   stucture "data" contains:
%       data.raw_image  Input image (before any processing)
%       data.params     Parameters used during image analysis
%       data.cells      Structure with a Label Matrix "label_matrix" and
%                       associated Connected Components Structure ".CC"
%                       containing the segmentation output. For more
%                       information about label matrices and connected
%                       components structures, see Matlab documentation.
%
% See documentation PDF for more information
%
% Last updated: 
%   Version 1.0 October 28, 2016
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

% Last Modified by GUIDE v2.5 10-Nov-2016 16:50:21

%% ========== Begin initialization code - DO NOT EDIT ================== %%
gui_Singleton = 1;
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

%% ===================== Customized functions ========================== %%

function watershed_cells_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% --- Executes just before watershed_cells_gui is made visible.
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to watershed_cells_gui (see VARARGIN)
% UIWAIT makes watershed_cells_gui wait for user response (see UIRESUME)
% uiwait(handles.watershed_cells_gui);

% Choose default command line output for watershed_cells_gui
handles.output = hObject;

% setup default params structure
handles.UserData.params.values = struct(...
    'path_to_image', '',...
    'equalization_cliplim', 0.01,...
    'background_size', 19,...
    'median_size', 7,...
    'gaussian_sigma', 7,...
    'minimum_area', 35,...
    'maximum_area', 1000,...
    'minimum_signal', 0.2,...
    'edge_alpha', 0.5,...
    'f', 'mean(R)-mean(G)',...
    'threshold', 0);
handles.UserData.params.on = struct(...
    'path_to_image', [],...
    'equalization_cliplim', true,...
    'background_size', true,...
    'median_size', true,...
    'gaussian_sigma', true,...
    'minimum_area', true,...
    'maximum_area', true,...
    'minimum_signal', true);
handles.UserData.results = struct(...
    'num_regions', [],...
    'num_state1', [],...
    'num_state2', [],...
    'raw_image', [],...
    'gray_image', [], ...
    'label_matrix', [],...
    'edge_image', [],...
    'state1_image', [], ...
    'state2_image', []);

% populate GUI with default values
handles = update_params_display(handles);

% update axes display
handles = initialize_display(handles);

% Update handles structure
guidata(hObject, handles);

function handles = initialize_display(handles)
% initialize the axes image data

handles.segmentation_axes.NextPlot = 'add';
handles.classify_axes.NextPlot = 'add';
linkaxes([handles.segmentation_axes handles.classify_axes]);

% grayscale image handles
handles.UserData.h_grayim = imshow(handles.UserData.results.gray_image, 'parent', handles.segmentation_axes);
handles.UserData.h_grayim.Tag = 'gray';

% region outlines image
handles.UserData.h_edges = imshow(handles.UserData.results.edge_image, 'parent', handles.segmentation_axes);
handles.UserData.h_edges.Tag = 'edges';

% raw image handles
handles.UserData.h_rawim = imshow(handles.UserData.results.raw_image, 'parent', handles.classify_axes);
handles.UserData.h_rawim.Tag = 'raw';

% region outlines by state image
handles.UserData.h_state1 = imshow(handles.UserData.results.state1_image, 'parent', handles.classify_axes);
handles.UserData.h_state1.Tag = 'state1';
handles.UserData.h_state2 = imshow(handles.UserData.results.state2_image, 'parent', handles.classify_axes);
handles.UserData.h_state2.Tag = 'state2';

function handles = update_params_display(handles)
% display parameter values in the gui

params = handles.UserData.params.values;
params_on = handles.UserData.params.on;
results = handles.UserData.results;

% segmentation parameter string values
handles.pathtoimage.String = params.path_to_image;
handles.sz_equalization.String = num2str(params.equalization_cliplim);
handles.sz_background.String = num2str(params.background_size);
handles.sz_median.String = num2str(params.median_size);
handles.sz_gaussian.String = num2str(params.gaussian_sigma);
handles.sz_minarea.String = num2str(params.minimum_area);
handles.sz_maxarea.String = num2str(params.maximum_area);
handles.sz_minsignal.String = num2str(params.minimum_signal);
handles.sz_edgealpha.String = num2str(params.edge_alpha);

% segmentation parameters on/off values
handles.run_equalization.Value = params_on.equalization_cliplim;
handles.run_background.Value = params_on.background_size;
handles.run_median.Value = params_on.median_size;
handles.run_gaussian.Value = params_on.gaussian_sigma;
handles.run_minarea.Value = params_on.minimum_area;
handles.run_maxarea.Value = params_on.maximum_area;
handles.run_minsignal.Value = params_on.minimum_signal;

% classification parameters string values
handles.classifyeq.String = params.f;
handles.threshold.String = num2str(params.threshold);

% number of regions
handles.numregions.String = sprintf('Number of regions: %d', results.num_regions);
handles.numclasses.String = sprintf('State 1: %d, State 2: %d', ...
    results.num_state1, results.num_state2);

function handles = get_params(handles)
% get parameter values from the gui

% segmentation parameter string values
handles.UserData.params.values.path_to_image = handles.pathtoimage.String;
handles.UserData.params.values.equalization_cliplim = str2double(handles.sz_equalization.String);
handles.UserData.params.values.background_size = str2double(handles.sz_background.String);
handles.UserData.params.values.median_size = str2double(handles.sz_median.String);
handles.UserData.params.values.gaussian_sigma = str2double(handles.sz_gaussian.String);
handles.UserData.params.values.minimum_area = str2double(handles.sz_minarea.String);
handles.UserData.params.values.maximum_area = str2double(handles.sz_maxarea.String);
handles.UserData.params.values.edge_alpha = str2double(handles.sz_edgealpha.String);

% segmentation parameter on/off values
handles.UserData.params.on.equalization_cliplim = handles.run_equalization.Value;
handles.UserData.params.on.background_size = handles.run_background.Value;
handles.UserData.params.on.median_size = handles.run_median.Value;
handles.UserData.params.on.gaussian_size = handles.run_gaussian.Value;
handles.UserData.params.on.minimum_area = handles.run_minarea.Value;
handles.UserData.params.on.maximum_area = handles.run_maxarea.Value;
handles.UserData.params.on.minimum_signal = handles.run_minsignal.Value;

% classification parameters
handles.UserData.params.values.f = handles.classifyeq.String;
handles.UserData.params.values.threshold = str2double(handles.threshold.String);

% check/fix parameter values
params = handles.UserData.params.values;
params.equalization_cliplim = min(1, abs(params.equalization_cliplim)); % positive <= 1
params.background_size = floor(abs(params.background_size)/2)*2+1; % positive odd integer
params.median_size = floor(abs(params.median_size)/2)*2+1; % positive odd integer
params.gaussian_sigma = abs(params.gaussian_sigma); % positive
params.minimum_area = round(abs(params.minimum_area)); % positive integer
params.maximum_area = round(abs(params.maximum_area)); % positive integer
params.minimum_signal = min(1, abs(params.minimum_signal)); % positive <= 1
params.edge_alpha = min(1, abs(params.edge_alpha)); % positive <= 1
handles.UserData.params.values = params;

% update display with fixed values
handles = update_params_display(handles);

function browsebutton_Callback(hObject, eventdata, handles)
% --- Executes on button press in browsebutton.
% hObject    handle to browsebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% get the file path
[filename, pathname] = uigetfile({'*.tif;*.tiff'}, 'select image to analyze');
if filename
    handles.UserData.params.values.path_to_image = [pathname filename];

    % update parameter display
    handles = update_params_display(handles);

    % import image
    handles = import_image(handles);

    % Update handles structure
    guidata(hObject, handles);

    % Show the image
    updatedisplaybutton_Callback(hObject, eventdata, handles)
end

function handles = import_image(handles)
% import image from the path and store grayscale image

if exist(handles.UserData.params.values.path_to_image, 'file')    
    % raw image
    handles.UserData.results.raw_image = imread(handles.UserData.params.values.path_to_image);
    
    % grayscale image
    if size(handles.UserData.results.raw_image, 3) > 1 
        handles.UserData.results.gray_image = mat2gray(rgb2gray(handles.UserData.results.raw_image));
    else
        handles.UserData.results.gray_image = mat2gray(handles.UserData.results.raw_image);
    end
    
else
    warning('Image not found')
end

function findcellsbutton_Callback(hObject, eventdata, handles)
% --- Executes on button press in findcellsbutton.
% hObject    handle to findcellsbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% show user we are working on stuff
handles.computing.Visible = 'on';
drawnow

% get current parameters
handles = get_params(handles);
params = handles.UserData.params.values;
params_on = handles.UserData.params.on;

% read in the image
handles = import_image(handles);

% if we have an image, run the analysis to find cells and store result
if ~isempty(handles.UserData.results.gray_image)
    [label_matrix, edge_image] = find_cells(...
        handles.UserData.results.gray_image, params, params_on);
    handles.UserData.results.label_matrix = label_matrix;
    handles.UserData.results.edge_image = edge_image;
    handles.UserData.results.num_regions = max(label_matrix(:));
end

% automatically call the updatedisplay button
updatedisplaybutton_Callback(hObject, eventdata, handles);

% show user we are done computing stuff
handles.computing.Visible = 'off';
drawnow

% Update handles structure
guidata(hObject, handles);

function updatedisplaybutton_Callback(hObject, eventdata, handles)
% --- Executes on button press in updatedisplaybutton.
% hObject    handle to updatedisplaybutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% show user we are working on stuff
handles.computing.Visible = 'on';
drawnow

% update parameters
handles = get_params(handles);

% show grayscale image
if ~isempty(handles.UserData.results.gray_image) && ...
        ~isequal(handles.UserData.h_grayim.CData, handles.UserData.results.gray_image)
    handles.UserData.h_grayim.CData = handles.UserData.results.gray_image;
    handles.segmentation_axes.XTick = [];
    handles.segmentation_axes.YTick = [];
    handles.segmentation_axes.XLim = [0 size(handles.UserData.results.gray_image, 2)]+0.5;
    handles.segmentation_axes.YLim = [0 size(handles.UserData.results.gray_image, 1)]+0.5;
end

% show edge image
if ~isempty(handles.UserData.results.edge_image)
    
    % deal with different image classes 
    edgeim = handles.UserData.results.edge_image;
    switch class(handles.UserData.results.gray_image)
        case 'uint8'
            edgeim = uint8(edgeim*2^8);
        case 'double'
            edgeim = double(edgeim);
        otherwise
            edgeim = double(edgeim);
    end
    handles.UserData.h_edges.CData = edgeim;
    
    % set alpha
    handles.UserData.h_edges.AlphaData = ...
        handles.UserData.results.edge_image * ...
        handles.UserData.params.values.edge_alpha;
    
end


% show raw image in classify_axes
if ~isempty(handles.UserData.results.raw_image) && ...
        ~isequal(handles.UserData.h_rawim.CData, handles.UserData.results.raw_image)
    handles.UserData.h_rawim.CData = handles.UserData.results.raw_image;
    handles.classify_axes.XTick = [];
    handles.classify_axes.YTick = [];
    handles.classify_axes.XLim = [0 size(handles.UserData.results.raw_image, 2)]+0.5;
    handles.classify_axes.YLim = [0 size(handles.UserData.results.raw_image, 1)]+0.5;
end

% show state1 image in classify_axes
if ~isempty(handles.UserData.results.state1_image)
    
    % deal with different image classes 
    edgeim = handles.UserData.results.state1_image;
    switch class(handles.UserData.results.raw_image)
        case 'uint8'
            edgeim = uint8(edgeim*2^8);
        case 'double'
            edgeim = double(edgeim);
        otherwise
            edgeim = double(edgeim);
    end
    handles.UserData.h_state1.CData = cat(3, edgeim, edgeim, edgim*0);
    
    % set alpha
    handles.UserData.h_state1.AlphaData = ...
        handles.UserData.results.state1_image * ...
        handles.UserData.params.values.edge_alpha;
    
end

% show state2 image in classify_axes
if ~isempty(handles.UserData.results.state2_image)
    
    % deal with different image classes 
    edgeim = handles.UserData.results.state2_image;
    switch class(handles.UserData.results.raw_image)
        case 'uint8'
            edgeim = uint8(edgeim*2^8);
        case 'double'
            edgeim = double(edgeim);
        otherwise
            edgeim = double(edgeim);
    end
    handles.UserData.h_state2.CData = cat(3, edgeim*0, edgeim, edgim);
    
    % set alpha
    handles.UserData.h_state2.AlphaData = ...
        handles.UserData.results.state2_image * ...
        handles.UserData.params.values.edge_alpha;
    
end

% show user we are working on stuff
handles.computing.Visible = 'off';
drawnow

% Update handles structure
guidata(hObject, handles);

function savedatabutton_Callback(hObject, eventdata, handles)
% --- Executes on button press in savedatabutton.
% hObject    handle to savedatabutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% bring up a dialog box to pick the save name/location
[FileName,PathName] = uiputfile('*.*', 'Basename to save results', 'output');
if FileName
    basename = regexp(FileName, '\.\S\S\S$', 'split');
    basename = basename{1};

    % save parameters as matlab files
    params = handles.UserData.params;
    save([PathName basename '_params.mat'], 'params', '-mat');

    % print parameters & region totals to text file
    params = handles.UserData.params.values;
    params_on = handles.UserData.params.on;
    keys = fieldnames(params);
    f = fopen([PathName basename '_params.txt'], 'wt+');
    fprintf(f, 'NAME\tVALUE\tON\n');
    for ii = 1:length(keys)
        value = params.(keys{ii});
        try
            on = params_on.(keys{ii});
        catch
            on = 1;
        end
        if isnumeric(value)
            fprintf(f, '%s\t%f\t%d\n', keys{ii}, value, on);
        elseif ischar(value)
            fprintf(f, '%s\t%s\t%d\n', keys{ii}, value, on);
        end
    end
    fprintf(f, '%s\t%d\t\n', 'num_regions', handles.UserData.results.num_regions);
    fprintf(f, '%s\t%d\t\n', 'num_state1', handles.UserData.results.num_state1);
    fprintf(f, '%s\t%d\t\n', 'num_state2', handles.UserData.results.num_state2);
    fclose(f);

    % save label matrix (identifying all regions)
    label_matrix = handles.UserData.results.label_matrix;
    if ~isempty(label_matrix) && ~isempty(handles.UserData.results.num_regions)
        % matlab file
        save([PathName basename '_labelmatrix.mat'], 'label_matrix', '-mat');

        % tiff image
        imwrite(label_matrix, gray(handles.UserData.results.num_regions), [PathName basename '_labelmatrix.tif'], 'tif');
    end
    
    % save binary images identifying state1 and state2 regions
    disp('add stuff here!')
    
end

function classifybutton_Callback(hObject, eventdata, handles)
% --- Executes on button press in classifybutton.
% hObject    handle to classifybutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% show user we are working on stuff
handles.computing.Visible = 'on';
drawnow

% get current parameters & results
handles = get_params(handles);
params = handles.UserData.params.values;
results = handles.UserData.results;

% construct classification equation
classify_fun = eval(sprintf('@(R,G,B)(%s)', params.f));

% run classification
if isa(classify_fun, 'function_handle')
    
    
    % get each channel. if no color, set all channels as the raw intensity
    switch size(results.raw_image, 3)
        case 3
            imR = results.raw_image(:,:,1);
            imG = results.raw_image(:,:,1);
            imB = results.raw_image(:,:,1);
        case 1
            [imR, imG, imB] = deal(results.raw_image);
    end
    
    
    % gather linear pixel indices in each region
    stats = regionprops(results.label_matrix, 'PixelIdxList');
    
    
    
else
    warning('Invalid classification function.')
end



%% ===================== Unaltered functions =========================== %%

function varargout = watershed_cells_gui_OutputFcn(hObject, eventdata, handles) 
% --- Outputs from this function are returned to the command line.
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

%% --------------------- CreateFcn

function sz_median_CreateFcn(hObject, eventdata, handles)
% --- Executes during object creation, after setting all properties.
% hObject    handle to sz_median (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function sz_background_CreateFcn(hObject, eventdata, handles)
% --- Executes during object creation, after setting all properties.
% hObject    handle to sz_background (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function sz_gaussian_CreateFcn(hObject, eventdata, handles)
% --- Executes during object creation, after setting all properties.
% hObject    handle to sz_gaussian (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function sz_minarea_CreateFcn(hObject, eventdata, handles)
% --- Executes during object creation, after setting all properties.
% hObject    handle to sz_minarea (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function sz_maxarea_CreateFcn(hObject, eventdata, handles)
% --- Executes during object creation, after setting all properties.
% hObject    handle to sz_maxarea (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function sz_minsignal_CreateFcn(hObject, eventdata, handles)
% --- Executes during object creation, after setting all properties.
% hObject    handle to sz_minsignal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function pathtoimage_CreateFcn(hObject, eventdata, handles)
% --- Executes during object creation, after setting all properties.
% hObject    handle to pathtoimage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function sz_equalization_CreateFcn(hObject, eventdata, handles)
% --- Executes during object creation, after setting all properties.
% hObject    handle to sz_equalization (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function sz_edgealpha_CreateFcn(hObject, eventdata, handles)
% --- Executes during object creation, after setting all properties.
% hObject    handle to sz_edgealpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function classifyeq_CreateFcn(hObject, eventdata, handles)
% --- Executes during object creation, after setting all properties.
% hObject    handle to classifyeq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function threshold_CreateFcn(hObject, eventdata, handles)
% --- Executes during object creation, after setting all properties.
% hObject    handle to threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% --------------------- Callback

function run_background_Callback(hObject, eventdata, handles)
% --- Executes on button press in run_background.
% hObject    handle to run_background (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of run_background

function run_median_Callback(hObject, eventdata, handles)
% --- Executes on button press in run_median.
% hObject    handle to run_median (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of run_median

function run_gaussian_Callback(hObject, eventdata, handles)
% --- Executes on button press in run_gaussian.
% hObject    handle to run_gaussian (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of run_gaussian

function run_minarea_Callback(hObject, eventdata, handles)
% --- Executes on button press in run_minarea.
% hObject    handle to run_minarea (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of run_minarea

function run_maxarea_Callback(hObject, eventdata, handles)
% --- Executes on button press in run_maxarea.
% hObject    handle to run_maxarea (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of run_maxarea

function run_minsignal_Callback(hObject, eventdata, handles)
% --- Executes on button press in run_minsignal.
% hObject    handle to run_minsignal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of run_minsignal

function sz_background_Callback(hObject, eventdata, handles)
% hObject    handle to sz_background (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of sz_background as text
%        str2double(get(hObject,'String')) returns contents of sz_background as a double

function sz_median_Callback(hObject, eventdata, handles)
% hObject    handle to sz_median (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of sz_median as text
%        str2double(get(hObject,'String')) returns contents of sz_median as a double

function sz_gaussian_Callback(hObject, eventdata, handles)
% hObject    handle to sz_gaussian (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of sz_gaussian as text
%        str2double(get(hObject,'String')) returns contents of sz_gaussian as a double

function sz_minarea_Callback(hObject, eventdata, handles)
% hObject    handle to sz_minarea (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of sz_minarea as text
%        str2double(get(hObject,'String')) returns contents of sz_minarea as a double

function sz_maxarea_Callback(hObject, eventdata, handles)
% hObject    handle to sz_maxarea (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of sz_maxarea as text
%        str2double(get(hObject,'String')) returns contents of sz_maxarea as a double

function sz_minsignal_Callback(hObject, eventdata, handles)
% hObject    handle to sz_minsignal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of sz_minsignal as text
%        str2double(get(hObject,'String')) returns contents of sz_minsignal as a double

function pathtoimage_Callback(hObject, eventdata, handles)
% hObject    handle to pathtoimage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of pathtoimage as text
%        str2double(get(hObject,'String')) returns contents of pathtoimage as a double

function run_celloutlines_Callback(hObject, eventdata, handles)
% --- Executes on button press in run_celloutlines.
% hObject    handle to run_celloutlines (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of run_celloutlines

function run_coloroutlines_Callback(hObject, eventdata, handles)
% --- Executes on button press in run_coloroutlines.
% hObject    handle to run_coloroutlines (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of run_coloroutlines

function run_equalization_Callback(hObject, eventdata, handles)
% --- Executes on button press in run_equalization.
% hObject    handle to run_equalization (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of run_equalization

function sz_equalization_Callback(hObject, eventdata, handles)
% hObject    handle to sz_equalization (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of sz_equalization as text
%        str2double(get(hObject,'String')) returns contents of sz_equalization as a double

function sz_edgealpha_Callback(hObject, eventdata, handles)
% hObject    handle to sz_edgealpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of sz_edgealpha as text
%        str2double(get(hObject,'String')) returns contents of sz_edgealpha as a double

function classifyeq_Callback(hObject, eventdata, handles)
% hObject    handle to classifyeq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of classifyeq as text
%        str2double(get(hObject,'String')) returns contents of classifyeq as a double

function threshold_Callback(hObject, eventdata, handles)
% hObject    handle to threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of threshold as text
%        str2double(get(hObject,'String')) returns contents of threshold as a double


%% --------------------- ButtonDownFcn

function findcellsbutton_ButtonDownFcn(hObject, eventdata, handles)
% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over findcellsbutton.
% hObject    handle to findcellsbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function updatedisplaybutton_ButtonDownFcn(hObject, eventdata, handles)
% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over updatedisplaybutton.
% hObject    handle to updatedisplaybutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function classifybutton_ButtonDownFcn(hObject, eventdata, handles)
% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over classifybutton.
% hObject    handle to classifybutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
