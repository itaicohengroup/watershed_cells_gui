function varargout = batch_process(varargin)
% BATCH_PROCESS MATLAB code for batch_process.fig
%      BATCH_PROCESS, by itself, creates a new BATCH_PROCESS or raises the existing
%      singleton*.
%
%      H = BATCH_PROCESS returns the handle to a new BATCH_PROCESS or the handle to
%      the existing singleton*.
%
%      BATCH_PROCESS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BATCH_PROCESS.M with the given input arguments.
%
%      BATCH_PROCESS('Property','Value',...) creates a new BATCH_PROCESS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before batch_process_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to batch_process_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help batch_process

% Last Modified by GUIDE v2.5 14-Nov-2016 17:07:15

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @batch_process_OpeningFcn, ...
                   'gui_OutputFcn',  @batch_process_OutputFcn, ...
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

%% --------------------- Opening & Closing

function batch_process_OpeningFcn(hObject, eventdata, handles, varargin)
% --- Executes just before batch_process is made visible.
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to batch_process (see VARARGIN)
% 
% UIWAIT makes batch_process wait for user response (see UIRESUME)
% uiwait(handles.batch_process);

% Choose default command line output for batch_process
handles.output = hObject;

% Input argument should be the analysis parameters & (blank) results
data = varargin{1};
data.images.paths = {};
data.images.current_index = 1;
handles.output.UserData = data;

% update imagelistbox
handles = update_listbox(handles);

% Update handles structure
guidata(hObject, handles);

function handles = update_listbox(handles)

% udpate strings
handles.imagelistbox.String = ...
    sprintf('%s\n', handles.output.UserData.images.paths{:});

% update selection location
handles.imagelistbox.Value = handles.output.UserData.images.current_index;

function varargout = batch_process_OutputFcn(hObject, eventdata, handles) 
% --- Outputs from this function are returned to the command line.
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function imagelistbox_CreateFcn(hObject, eventdata, handles)
% --- Executes during object creation, after setting all properties.
% hObject    handle to imagelistbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% --------------------- Callbacks - Image selection

function imagelistbox_Callback(hObject, eventdata, handles)
% --- Executes on selection change in imagelistbox.
% hObject    handle to imagelistbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: contents = cellstr(get(hObject,'String')) returns imagelistbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from imagelistbox

handles.output.UserData.images.current_index = hObject.Value;

function addfilesbutton_Callback(hObject, eventdata, handles)
% --- Executes on button press in addfilesbutton.
% hObject    handle to addfilesbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
1;

images = handles.output.UserData.images.paths;

% get the file path
[filename, pathname, success] = uigetfile({'*.tif;*.tiff'}, 'select image(s) to analyze', 'MultiSelect', 'on');
if success
    
    % for each selected file
    if ~isa(filename, 'cell')
        filename = {filename};
    end
    for ff = 1:length(filename)
        
        % construct full path to file
        fullpath = fullfile(pathname, filename{ff});
    
        % check if this file is already in the list
        match = false;
        for ii = 1:length(images)
            if strcmp(fullpath, images{ii});
                match = true;
                break
            end
        end

        % add it to the list
        if ~match
            images{end+1} = fullpath;
        end
        
    end
    
    % Update the display & ui data
    handles.output.UserData.images.paths = images;
    handles = update_listbox(handles);
    
    % Update handles structure
    guidata(hObject, handles);
    
end

function removebutton_Callback(hObject, eventdata, handles)
% --- Executes on button press in removebutton.
% hObject    handle to removebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% number of paths currently
maxix = length(handles.output.UserData.images.paths);

% delete the selected path
ix = handles.output.UserData.images.current_index;
if ix <= maxix
    handles.output.UserData.images.paths(ix) = [];
    maxix = length(handles.output.UserData.images.paths);
end

% move selection to the end of the list, if we are past the end
if ix > maxix && maxix > 0
    handles.output.UserData.images.current_index = maxix;
end

% update the listbox
handles = update_listbox(handles);

% Update handles structure
guidata(hObject, handles);

function runbutton_Callback(hObject, eventdata, handles)
% --- Executes on button press in runbutton.
% hObject    handle to runbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
1;

% Update handles structure
guidata(hObject, handles);

