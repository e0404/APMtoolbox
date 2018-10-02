function varargout = covarianceInfluenceVisualizer(varargin)
% COVARIANCEINFLUENCEVISUALIZER MATLAB code for covarianceInfluenceVisualizer.fig
%      COVARIANCEINFLUENCEVISUALIZER, by itself, creates a new COVARIANCEINFLUENCEVISUALIZER or raises the existing
%      singleton*.
%
%      H = COVARIANCEINFLUENCEVISUALIZER returns the handle to a new COVARIANCEINFLUENCEVISUALIZER or the handle to
%      the existing singleton*.
%
%      COVARIANCEINFLUENCEVISUALIZER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in COVARIANCEINFLUENCEVISUALIZER.M with the given input arguments.
%
%      COVARIANCEINFLUENCEVISUALIZER('Property','Value',...) creates a new COVARIANCEINFLUENCEVISUALIZER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before covarianceInfluenceVisualizer_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to covarianceInfluenceVisualizer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help covarianceInfluenceVisualizer

% Last Modified by GUIDE v2.5 07-Jun-2018 16:38:59

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @covarianceInfluenceVisualizer_OpeningFcn, ...
                   'gui_OutputFcn',  @covarianceInfluenceVisualizer_OutputFcn, ...
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
% End initialization code - DO NOT EDIT

% --- Executes just before covarianceInfluenceVisualizer is made visible.
function covarianceInfluenceVisualizer_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to covarianceInfluenceVisualizer (see VARARGIN)

% Choose default command line output for covarianceInfluenceVisualizer
handles.output = hObject;

handles.covInfluence = evalin('base','covInfluence');

nSpots = size(handles.covInfluence,2);

set(handles.slider1,'SliderStep',[1 1]./nSpots);
set(handles.slider1,'Min',1);
set(handles.slider1,'Max',nSpots);
set(handles.slider1,'Value',1);


set(handles.slider2,'SliderStep',[1 1]./nSpots);
set(handles.slider2,'Min',1);
set(handles.slider2,'Max',nSpots);
set(handles.slider2,'Value',1);

% Update handles structure
guidata(hObject, handles);

% This sets up the initial plot - only do when we are invisible
% so window can get raised using covarianceInfluenceVisualizer.
if strcmp(get(hObject,'Visible'),'off')
    imagesc(squeeze(handles.covInfluence (:,1,:,1)));
end

% UIWAIT makes covarianceInfluenceVisualizer wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = covarianceInfluenceVisualizer_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --------------------------------------------------------------------
function FileMenu_Callback(hObject, eventdata, handles)
% hObject    handle to FileMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function OpenMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to OpenMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
file = uigetfile('*.fig');
if ~isequal(file, 0)
    open(file);
end

% --------------------------------------------------------------------
function PrintMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to PrintMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
printdlg(handles.figure1)

% --------------------------------------------------------------------
function CloseMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to CloseMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selection = questdlg(['Close ' get(handles.figure1,'Name') '?'],...
                     ['Close ' get(handles.figure1,'Name') '...'],...
                     'Yes','No','Yes');
if strcmp(selection,'No')
    return;
end

delete(handles.figure1)


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


hObject.Value = round(hObject.Value);
imagesc(squeeze(handles.covInfluence(:,hObject.Value,:,get(handles.slider2,'Value'))));


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider2_Callback(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

hObject.Value = round(hObject.Value);
imagesc(squeeze(handles.covInfluence(:,get(handles.slider1,'Value'),:,hObject.Value)));


% --- Executes during object creation, after setting all properties.
function slider2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
