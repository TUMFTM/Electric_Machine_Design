function varargout = GUI_LDS(varargin)
% figure1 MATLAB code for figure1.fig
%      figure1, by itself, creates a new figure1 or raises the existing
%      singleton*.
%
%      H = figure1 returns the handle to a new figure1 or the handle to
%      the existing singleton*.
%
%      figure1('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in figure1.M with the given input arguments.
%
%      figure1('Property','Value',...) creates a new figure1 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_LDS_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_LDS_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help figure1

% Last Modified by GUIDE v2.5 02-Jul-2019 12:51:14

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_LDS_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_LDS_OutputFcn, ...
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


% --- Executes just before figure1 is made visible.
function GUI_LDS_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to figure1 (see VARARGIN)

% Choose default command line output for figure1
handles.output = hObject;

handles.FZG_LDS.fz_m = str2double(handles.fz_m.String);
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes figure1 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_LDS_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles;

% --- Executes when user attempts to close figure1.
function GUI_LDS_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

button = questdlg('Wollen Sie wirklich beenden?', 'Beenden','Ja','Nein','Nein'); 
switch button 
    case 'Ja'
        uiresume(handles.GUI_LDS);
       % delete(handles.GUI_LDS);
    case 'Nein'
end


function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in BMWi3.
function BMWi3_Callback(hObject, eventdata, handles)
% hObject    handle to BMWi3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Set to default BMW i3 vehicle
handles.fz_m.String = num2str(1640); 
handles.cW.String = num2str(0.3);
handles.A.String = num2str(2.38);
handles.tyre_r.String = num2str(0.3498);
handles.battcap.String = num2str(22);
handles.aux.String = num2str(1500);
handles.GearRatio.String = num2str(9.7);
    
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in NissanLeaf.
function NissanLeaf_Callback(hObject, eventdata, handles)
% hObject    handle to NissanLeaf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Set to default Nissan Leaf vehicle [CHO13]
handles.fz_m.String = num2str(1505);        %1505 
handles.cW.String = num2str(0.3);
handles.A.String = num2str(2.27);           % eigentlich 2.27
handles.tyre_r.String = num2str(0.316);     %205/55R16
handles.battcap.String = num2str(30);       %kleine Version mit 24 kWh möglich, große mit 30 kWh
handles.aux.String = num2str(1300);         %Sebastian Jeschke: Grundlegende Untersuchungen von Elektrofahrzeugen im Bezug auf Energieeffizienz und EMV mit einer skalierbaren Power-HiL-Umgebung
handles.GearRatio.String = num2str(7.937);  %Wikipedia: 7.937

handles.FZG_LDS.fz_m = str2double(handles.fz_m.String);
% Aktualisieren
% handles_GUI_LDS.fz_m = str2double(handles.fz_m.String)*1000;
% handles_GUI_LDS.cW = str2double(handles.cW.String);
% handles_GUI_LDS.A = str2double(handles.A.String);
% handles_GUI_LDS.tyre_r = str2double(handles.tyre_r.String);
% handles_GUI_LDS.battcap = str2double(handles.battcap.String);
% handles_GUI_LDS.aux = str2double(handles.aux.String); 
    
% Update handles structure
guidata(hObject, handles);



function Massetest_Callback(hObject, eventdata, handles)
% hObject    handle to Massetest (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.massetest = str2double(get(hObject,'String'));



% --- Executes during object creation, after setting all properties.
function Massetest_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Massetest (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
delete(hObject);


% --- Executes on button press in Start_LDS.
function Start_LDS_Callback(~,~, handles)
% hObject    handle to Start_LDS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton_TeslaModel3.
function pushbutton_TeslaModel3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_TeslaModel3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Set to default Tesla Model 3
handles.fz_m.String = num2str(1610);        	 %Angabe Audi
handles.cW.String = num2str(0.228);               %Angabe Audi
handles.A.String = num2str(2.28);                %Angabe Audi
handles.tyre_r.String = num2str(0.3353);         %235/40R19 
handles.battcap.String = num2str(78.1);          %Angabe Audi
handles.aux.String = num2str(5000);         
handles.GearRatio.String = num2str(9.7);         %Angabe Audi

handles.FZG_LDS.fz_m = str2double(handles.fz_m.String);

guidata(hObject, handles);


% --- Executes on button press in VW_eGolf.
function VW_eGolf_Callback(hObject, eventdata, handles)
% hObject    handle to VW_eGolf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Set to default Tesla Model 3
handles.fz_m.String = num2str(1615);        	 
handles.cW.String = num2str(0.27);              
handles.A.String = num2str(2.19);               
handles.tyre_r.String = num2str(0.6319);              %205/55 R16 H
handles.battcap.String = num2str(35.8);         % in kWh      
handles.aux.String = num2str(5000);         
handles.GearRatio.String = num2str(9.76);         

handles.FZG_LDS.fz_m = str2double(handles.fz_m.String);

guidata(hObject, handles);
