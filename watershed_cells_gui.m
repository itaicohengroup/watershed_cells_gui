% WATERSHED_CELLS_GUI finds cells in an image using watershed segmentation
%
% Version: 1.0
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

% Last Modified by GUIDE v2.5 06-Nov-2016 20:18:36

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

function findcellsbutton_Callback(hObject, eventdata, handles)
% --- Executes on button press in findcellsbutton.
% hObject    handle to findcellsbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% show user we are working on stuff
handles.computing.Visible = 'on';
drawnow

% If the analysis step is un-checked (i.e. to be skipped), set the
% parameter to 0 (i.e. false, off). Otherwise, grab the number from the
% text box and check/fix it as appropriate
if handles.run_equalization.Value
    num = str2double(handles.sz_equalization.String);
    params.equalization = abs(num); % positive
else
    params.equalization = 0;
end
if handles.run_background.Value
    num = str2double(handles.sz_background.String);
    params.background = floor(abs(num)/2)*2+1; % positive odd integer
else
    params.background = 0;
end
if handles.run_median.Value
    num = str2double(handles.sz_median.String);
    params.median = floor(abs(num)/2)*2+1; % positive odd integer
else
    params.median = 0;
end
if handles.run_gaussian.Value
    num = str2double(handles.sz_gaussian.String);
    params.gaussian = abs(num); % positive
else
    params.gaussian = 0;
end
if handles.run_minarea.Value
    num = str2double(handles.sz_minarea.String);
    params.minarea = round(abs(num)); % positive integer
else
    params.minarea = 0;
end
if handles.run_maxarea.Value
    num = str2double(handles.sz_maxarea.String);
    params.maxarea = round(abs(num)); % positive integer
else
    params.maxarea = 0;
end
if handles.run_minsignal.Value
    num = str2double(handles.sz_minsignal.String);
    params.minsignal = abs(num); % positive 
else
    params.minsignal = 0;
end

% get input image file. if it can't be found, open a dialog box so the user
% can find it
params.path_to_image = handles.pathtoimage.String;
if ~exist(params.path_to_image, 'file')
    [filename, pathname] = uigetfile({'*.tiff';'*.tif'}, 'select image to analyze');
    if filename
        params.path_to_image = [pathname filename];
        handles.pathtoimage.String = params.path_to_image;
        % Update handles structure
        guidata(hObject, handles);
    end
end

if exist(params.path_to_image, 'file')
    % store parameters in handles
    handles.params = params;
    
    % read image from file and store it in handles
    handles.raw_image = imread(params.path_to_image);
    
    % run function to find cells using the provided parameters and store
    % results in handles
    [label_matrix, CC, edge_image] = find_cells(handles.raw_image, handles.params);
    handles.cells.label_matrix = label_matrix;
    handles.cells.CC = CC;
    handles.edge_image = edge_image;

    % automatically call the updatedisplay button
    updatedisplaybutton_Callback(hObject, eventdata, handles);
    
    % show user we are done computing stuff
    handles.computing.Visible = 'off';
    drawnow

    % Update handles structure
    guidata(hObject, handles);
    
else
    % show user we are done computing stuff
    warning('Image not found')
    handles.computing.Visible = 'off';
    drawnow

end

function updatedisplaybutton_Callback(hObject, eventdata, handles)
% --- Executes on button press in updatedisplaybutton.
% hObject    handle to updatedisplaybutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% show user we are working on stuff
handles.computing.Visible = 'on';
drawnow

% show raw image
cla(handles.result_axes)
if isfield(handles, 'raw_image')
    % show image
    imshow(handles.raw_image, 'parent', handles.result_axes)
    
    % update axes appearance
    handles.result_axes.XTick = [];
    handles.result_axes.YTick = [];
    handles.result_axes.XLim = [0 size(handles.raw_image,2)]+0.5;
    handles.result_axes.YLim = [0 size(handles.raw_image,1)]+0.5;
end

% show cell outlines
handles.result_axes.NextPlot = 'add';
if isfield(handles, 'edge_image')
    % deal with different image classes 
    edgeim = handles.edge_image;
    switch class(handles.raw_image)
        case 'uint8'
            edgeim = uint8(edgeim*2^8);
        case 'double'
            edgeim = double(edgeim);
        otherwise
            edgeim = double(edgeim);
    end
    
    % show the cell outlines
    hedge = imshow(edgeim, 'parent', handles.result_axes);
    edgealpha = abs(str2double(handles.sz_edgealpha.String));
    hedge.AlphaData = handles.edge_image * edgealpha;
    
    % display the number of cells found
    handles.numcells.String = sprintf('Number of cells: %d', handles.cells.CC.NumObjects);
end

% Update handles structure
guidata(hObject, handles);

% show user we are working on stuff
handles.computing.Visible = 'off';
drawnow

function savedata_Callback(hObject, eventdata, handles)
% --- Executes on button press in savedata.
% hObject    handle to savedata (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% gather data to be saved
data = struct();
if isfield(handles, 'raw_image')
    data.raw_image = handles.raw_image;
end
if isfield(handles, 'params')
    data.params = handles.params;
end
if isfield(handles, 'cells')
    data.cells = handles.cells;
end

% bring up a dialog box to save the data
uisave('data','output.mat')

function browsebutton_Callback(hObject, eventdata, handles)
% --- Executes on button press in browsebutton.
% hObject    handle to browsebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% get existing parameters, if any
if isfield(handles, 'params')
    params = handles.params;
else
    params = struct();
end

% get the file path
[filename, pathname] = uigetfile({'*.tiff';'*.tif'}, 'select image to analyze');
if filename
    params.path_to_image = [pathname filename];
    handles.pathtoimage.String = params.path_to_image;
end

% read image from file and store it in handles
if exist(params.path_to_image, 'file')    
    handles.raw_image = imread(params.path_to_image);
else
    warning('Image not found')
end

% store parameters in handles
handles.params = params;

% Update handles structure
guidata(hObject, handles);

% Show the image
updatedisplaybutton_Callback(hObject, eventdata, handles)

%% ===================== Unaltered functions =========================== %%
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

% Update handles structure
guidata(hObject, handles);

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
