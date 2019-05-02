% -------------------------------------------------------------------------
% TU Munich - Institute of Automotive Technology
% -------------------------------------------------------------------------
% Model for the design of a permanent magnet excited synchronous machine and
% subsequent efficiency map calculation
% -------------------------------------------------------------------------
% Autor:    Svenja Kalt (kalt@ftm.mw.tum.de)
%           Jonathan Erhard
%           Prof. Markus Lienkamp
% -------------------------------------------------------------------------

function varargout = GUI_Kennfeld_PMSM(varargin)
% GUI_Kennfeld_PMSM MATLAB code for GUI_Kennfeld_PMSM.fig
%      GUI_Kennfeld_PMSM, by itself, creates a new GUI_Kennfeld_PMSM or raises the existing
%      singleton*.
%
%      H = GUI_Kennfeld_PMSM returns the handle to a new GUI_Kennfeld_PMSM or the handle to
%      the existing singleton*.
%
%      GUI_Kennfeld_PMSM('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_Kennfeld_PMSM.M with the given input arguments.
%
%      GUI_Kennfeld_PMSM('Property','Value',...) creates a new GUI_Kennfeld_PMSM or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI_Kennfeld_PMSM before GUI_Kennfeld_PMSM_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_Kennfeld_PMSM_OpeningFcn via varargin.
%
%      *See GUI_Kennfeld_PMSM Options on GUIDE's Tools menu.  Choose "GUI_Kennfeld_PMSM allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_Kennfeld_PMSM

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_Kennfeld_PMSM_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_Kennfeld_PMSM_OutputFcn, ...
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

% --- Executes just before GUI_Kennfeld_PMSM is made visible.
function GUI_Kennfeld_PMSM_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_Kennfeld_PMSM (see VARARGIN)

% Choose default command line output for GUI_Kennfeld_PMSM
handles.output = hObject;

% Read data
handles.Verluste.Stromwaermeverluste = handles.checkbox_Stromwaermeverluste.Value;
handles.Verluste.Eisenverluste = handles.checkbox_Eisenverluste.Value;
handles.Verluste.mechanische_Verluste = handles.checkbox_mechanische_Verluste.Value;
handles.Verluste.Zusatzverluste = handles.checkbox_Zusatzverluste.Value;
handles.Optionen.tics_M = str2double(handles.edit_tics_M.String);
handles.Optionen.tics_omega = str2double(handles.edit_tics_omega.String);
handles.Optionen.n_max_String = handles.popupmenu_n_max.String{handles.popupmenu_n_max.Value};
handles.Optionen.n_max_Value = handles.popupmenu_n_max.Value;
handles.Optionen.n_max = str2double(handles.edit_n_max.String);
handles.Optionen.Elektroband_Value = handles.popupmenu_Elektroband.Value;
handles.Optionen.Elektroband_String = handles.popupmenu_Elektroband.String{handles.popupmenu_Elektroband.Value};
handles.Regelgroessen.u_max = str2double(handles.edit_u_max.String);
handles.Regelgroessen.i_max = str2double(handles.edit_i_max.String);

% Save default data
handles.default.Regelgroessen = handles.Regelgroessen;
handles.default.Verluste = handles.Verluste;
handles.default.Optionen = handles.Optionen;

% Center GUI_Kennfeld_PMSM
movegui([844,100]);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI_Kennfeld_PMSM wait for user response (see UIRESUME)
uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = GUI_Kennfeld_PMSM_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles;

% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: delete(hObject) closes the figure

button = questdlg('Do you really want to quit?', 'Quit','Yes','No','No'); 
switch button 
    case 'Yes'
        uiresume(handles.figure1);
        delete(handles.figure1);
    case 'No'
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOSSES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on button press in checkbox_Stromwaermeverluste.
function checkbox_Stromwaermeverluste_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_Stromwaermeverluste (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of checkbox_Stromwaermeverluste

handles.Verluste.Stromwaermeverluste = get(hObject,'Value');

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in checkbox_Eisenverluste.
function checkbox_Eisenverluste_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_Eisenverluste (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of checkbox_Eisenverluste

handles.Verluste.Eisenverluste = get(hObject,'Value');

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in checkbox_mechanische_Verluste.
function checkbox_mechanische_Verluste_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_mechanische_Verluste (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of checkbox_mechanische_Verluste

handles.Verluste.mechanische_Verluste = get(hObject,'Value');

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in checkbox_Zusatzverluste.
function checkbox_Zusatzverluste_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_Zusatzverluste (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of checkbox_Zusatzverluste

handles.Verluste.Zusatzverluste = get(hObject,'Value');

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in pushbutton_Verluste.
function pushbutton_Verluste_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Verluste (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

uiwait(msgbox({'The losses that should be calculated,' 'can be selected here.'},'Verluste','help','modal'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OPTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function edit_tics_M_Callback(hObject, eventdata, handles)
% hObject    handle to edit_tics_M (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of edit_tics_M as text
%        str2double(get(hObject,'String')) returns contents of edit_tics_M as a double

handles.Optionen.tics_M = str2double(get(hObject,'String'));

% Update handles structure
guidata(hObject, handles);

function edit_tics_omega_Callback(hObject, eventdata, handles)
% hObject    handle to edit_tics_omega (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of edit_tics_omega as text
%        str2double(get(hObject,'String')) returns contents of edit_tics_omega as a double

handles.Optionen.tics_omega = str2double(get(hObject,'String'));

% Update handles structure
guidata(hObject, handles);

% --- Executes on selection change in popupmenu_n_max.
function popupmenu_n_max_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_n_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_n_max contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_n_max

contents = cellstr(get(hObject,'String'));
handles.Optionen.n_max_String = contents{get(hObject,'Value')};
handles.Optionen.n_max_Value = handles.popupmenu_n_max.Value;

if(strcmp(handles.Optionen.n_max_String,'max. Drehzahl in U/min'))
    set(handles.edit_n_max,'String', num2str(10000));
else
    set(handles.edit_n_max,'String', num2str(3.0));
end

handles.Optionen.n_max = str2double(handles.edit_n_max.String);

% Update handles structure
guidata(hObject, handles);

function edit_n_max_Callback(hObject, eventdata, handles)
% hObject    handle to edit_i_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of edit_i_max as text
%        str2double(get(hObject,'String')) returns contents of edit_i_max as a double

handles.Optionen.n_max = str2double(get(hObject,'String'));

% Update handles structure
guidata(hObject, handles);

% --- Executes on selection change in popupmenu_Elektroband.
function popupmenu_Elektroband_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_Elektroband (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_Elektroband contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_Elektroband

contents = cellstr(get(hObject,'String'));
handles.Optionen.Elektroband = contents{get(hObject,'Value')};

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in pushbutton_Optionen.
function pushbutton_Optionen_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Optionen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

uiwait(msgbox({'Here, the resolution and the max. speed of the map as well as the' 'the electrical steel strip is determined for the loss calculation.' 'ATTENTION: the resolution influences the calculation timE immensely!'},'Optionen','help','modal'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONTROLLED VARIABLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function edit_u_max_Callback(hObject, eventdata, handles)
% hObject    handle to edit_i_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of edit_i_max as text
%        str2double(get(hObject,'String')) returns contents of edit_i_max as a double

handles.Regelgroessen.u_max = str2double(get(hObject,'String'));

% Update handles structure
guidata(hObject, handles);

function edit_i_max_Callback(hObject, eventdata, handles)
% hObject    handle to edit_u_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of edit_u_max as text
%        str2double(get(hObject,'String')) returns contents of edit_u_max as a double

handles.Regelgroessen.i_max = str2double(get(hObject,'String'));

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in pushbutton_Regelgroessen.
function pushbutton_Regelgroessen_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Regelgroessen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

uiwait(msgbox({'ATTENTION: There is no thermal check.' 'The user is therefore required to design the machine using only' 'realistic data.'},'Regelgroessen','help','modal'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on selection change in popupmenu_Auswahl_Plot.
function popupmenu_Auswahl_Plot_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_Auswahl_Plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_Auswahl_Plot contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_Auswahl_Plot

contents = cellstr(get(hObject,'String'));
handles.Plot.Auswahl_Plot = contents{get(hObject,'Value')};

if(~strcmp(handles.Plot.Auswahl_Plot,'bitte auswaehlen ...'))
    var = handles.popupmenu_Auswahl_Plot.String;
    for i = 1:length(var)
        if(strcmp(var{i},'bitte auswaehlen ...'))
            handles.popupmenu_Auswahl_Plot.String(i,:) = [];
            handles.popupmenu_Auswahl_Plot.Value = handles.popupmenu_Auswahl_Plot.Value - 1;
            break;
        end
    end
end

% Loading the results
listing = dir('Ergebnisse');
listing = listing([listing.isdir] & contains({listing.name},'_data'));
[~,idx] =min(now*ones(size(listing)) - [listing.datenum]');
listing = listing(idx);
var = replace(listing.name,'_data','');
load(['Ergebnisse/',var,'_data/Maschinendaten_',var,'.mat']);
handles.folder = ['Ergebnisse/',var,'_data'];

% Parameter re-storage for easier use
% Rated sizes
prim = Maschinendaten.Bemessungsgroessen.Primaerparameter;
sek = Maschinendaten.Bemessungsgroessen.Sekundaerparameter;
% Efficiency diagram
kenn = Maschinendaten.Kennfeld;

% Plot
if(strcmp(handles.Plot.Auswahl_Plot,'Wirkungsgrad gesamt'))
    plot_Kennfeld;
    title(handles.axes_plot,'Wirkungsgrad gesamt $\eta_{V,ges}$ in \%','interpreter','latex','FontSize', 20);
    Plot_Label_Vektor_eta = [0.7,0.75,0.8,0.84,0.86,0.88,0.9,0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.98,0.99,0.995];
    [C,h]=contourf(handles.axes_plot,kenn.omega_k_mesh(1,:)/prim.p/(2*pi)*60,kenn.M_max_mesh(:,1),kenn.eta_ges_mesh,Plot_Label_Vektor_eta);
    clabel(C,h,Plot_Label_Vektor_eta);
elseif(strcmp(handles.Plot.Auswahl_Plot,'Verluste gesamt'))
    plot_Kennfeld;
    title(handles.axes_plot,'Verluste gesamt $P_{V,ges}$ in W','interpreter','latex','FontSize', 20);
    contourf(handles.axes_plot,kenn.omega_k_mesh(1,:)/prim.p/(2*pi)*60,kenn.M_max_mesh(:,1),kenn.Verluste.P_vges_mesh);
elseif(strcmp(handles.Plot.Auswahl_Plot,'Wirkungsgrad Stromwaermeverluste'))
    plot_Kennfeld;
    title(handles.axes_plot,'Wirkungsgrad Stromwaermeverluste $\eta_{V,sw}$ in \%','interpreter','latex','FontSize', 20);
    Plot_Label_Vektor_eta = [0.7,0.75,0.8,0.84,0.86,0.88,0.9,0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.98,0.99,0.995];
    [C,h]=contourf(handles.axes_plot,kenn.omega_k_mesh(1,:)/prim.p/(2*pi)*60,kenn.M_max_mesh(:,1),kenn.eta_vw_mesh,Plot_Label_Vektor_eta);
    clabel(C,h,Plot_Label_Vektor_eta);
elseif(strcmp(handles.Plot.Auswahl_Plot,'Stromwaermeverluste'))
    plot_Kennfeld;
    title(handles.axes_plot,'Stromwaermeverluste $P_{V,sw}$ in W','interpreter','latex','FontSize', 20);
    contourf(handles.axes_plot,kenn.omega_k_mesh(1,:)/prim.p/(2*pi)*60,kenn.M_max_mesh(:,1),kenn.Verluste.P_vw_mesh);
elseif(strcmp(handles.Plot.Auswahl_Plot,'Wirkungsgrad Eisenverluste'))
    plot_Kennfeld;
    title(handles.axes_plot,'Wirkungsgrad Eisenverluste $\eta_{V,fe}$ in \%','interpreter','latex','FontSize', 20);
    Plot_Label_Vektor_eta = [0.7,0.75,0.8,0.84,0.86,0.88,0.9,0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.98,0.99,0.995];
    [C,h]=contourf(handles.axes_plot,kenn.omega_k_mesh(1,:)/prim.p/(2*pi)*60,kenn.M_max_mesh(:,1),kenn.eta_fe_mesh,Plot_Label_Vektor_eta);
    clabel(C,h,Plot_Label_Vektor_eta);
elseif(strcmp(handles.Plot.Auswahl_Plot,'Eisenverluste'))
    plot_Kennfeld;
    title(handles.axes_plot,'Eisenverluste $P_{V,fe}$ in W','interpreter','latex','FontSize', 20);
    contourf(handles.axes_plot,kenn.omega_k_mesh(1,:)/prim.p/(2*pi)*60,kenn.M_max_mesh(:,1),kenn.Verluste.P_fe_mesh);
elseif(strcmp(handles.Plot.Auswahl_Plot,'Wirkungsgrad mechanische Verluste'))
    plot_Kennfeld;
    title(handles.axes_plot,'Wirkungsgrad mechanische Verluste $\eta_{V,mech}$ in \%','interpreter','latex','FontSize', 20);
    Plot_Label_Vektor_eta = [0.7,0.75,0.8,0.84,0.86,0.88,0.9,0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.98,0.99,0.995];
    [C,h]=contourf(handles.axes_plot,kenn.omega_k_mesh(1,:)/prim.p/(2*pi)*60,kenn.M_max_mesh(:,1),kenn.eta_vme_mesh,Plot_Label_Vektor_eta);
    clabel(C,h,Plot_Label_Vektor_eta);
elseif(strcmp(handles.Plot.Auswahl_Plot,'mechanische Verluste'))
    plot_Kennfeld;
    title(handles.axes_plot,'mechanische Verluste $P_{V,mech}$ in W','interpreter','latex','FontSize', 20);
    contourf(handles.axes_plot,kenn.omega_k_mesh(1,:)/prim.p/(2*pi)*60,kenn.M_max_mesh(:,1),kenn.Verluste.P_vme_mesh);
elseif(strcmp(handles.Plot.Auswahl_Plot,'Wirkungsgrad Zusatzverluste'))
    plot_Kennfeld;
    title(handles.axes_plot,'Wirkungsgrad Zusatzverluste $\eta_{V,zus}$ in \%','interpreter','latex','FontSize', 20);
    Plot_Label_Vektor_eta = [0.7,0.75,0.8,0.84,0.86,0.88,0.9,0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.98,0.99,0.995];
    [C,h]=contourf(handles.axes_plot,kenn.omega_k_mesh(1,:)/prim.p/(2*pi)*60,kenn.M_max_mesh(:,1),kenn.eta_zus_mesh,Plot_Label_Vektor_eta);
    clabel(C,h,Plot_Label_Vektor_eta);
elseif(strcmp(handles.Plot.Auswahl_Plot,'Zusatzverluste'))
    plot_Kennfeld;
    title(handles.axes_plot,'Zusatzverluste $P_{V,zus}$ in W','interpreter','latex','FontSize', 20);
    contourf(handles.axes_plot,kenn.omega_k_mesh(1,:)/prim.p/(2*pi)*60,kenn.M_max_mesh(:,1),kenn.Verluste.P_zus_mesh);
elseif(strcmp(handles.Plot.Auswahl_Plot,'i_d'))
    plot_Kennfeld;
    title(handles.axes_plot,'$i_d$ in A','interpreter','latex','FontSize', 20);
    contourf(handles.axes_plot,kenn.omega_k_mesh(1,:)/prim.p/(2*pi)*60,kenn.M_max_mesh(:,1),kenn.i_d_mesh);
elseif(strcmp(handles.Plot.Auswahl_Plot,'i_q'))
    plot_Kennfeld;
    title(handles.axes_plot,'$i_q$ in A','interpreter','latex','FontSize', 20);
    contourf(handles.axes_plot,kenn.omega_k_mesh(1,:)/prim.p/(2*pi)*60,kenn.M_max_mesh(:,1),kenn.i_q_mesh);
elseif(strcmp(handles.Plot.Auswahl_Plot,'u_d'))
    plot_Kennfeld;
    title(handles.axes_plot,'$u_d$ in V','interpreter','latex','FontSize', 20);
    contourf(handles.axes_plot,kenn.omega_k_mesh(1,:)/prim.p/(2*pi)*60,kenn.M_max_mesh(:,1),kenn.u_d_mesh);
elseif(strcmp(handles.Plot.Auswahl_Plot,'u_q'))
    plot_Kennfeld;
    title(handles.axes_plot,'$u_q$ in V','interpreter','latex','FontSize', 20);
    contourf(handles.axes_plot,kenn.omega_k_mesh(1,:)/prim.p/(2*pi)*60,kenn.M_max_mesh(:,1),kenn.u_q_mesh);
end

handles.pushbutton_save.Enable = 'on';

% Update handles structure
guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% START UND RESET
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on button press in pushbutton_reset.
function pushbutton_reset_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Reset to default
handles.checkbox_Stromwaermeverluste.Value = handles.default.Verluste.Stromwaermeverluste;
handles.checkbox_Eisenverluste.Value = handles.default.Verluste.Eisenverluste;
handles.checkbox_mechanische_Verluste.Value = handles.default.Verluste.mechanische_Verluste;
handles.checkbox_Zusatzverluste.Value = handles.default.Verluste.Zusatzverluste;
handles.edit_tics_M.String = num2str(handles.default.Optionen.tics_M);
handles.edit_tics_omega.String = num2str(handles.default.Optionen.tics_omega);
handles.popupmenu_n_max.Value = handles.default.Optionen.n_max_Value;
handles.edit_n_max.String = num2str(handles.default.Optionen.n_max);
handles.popupmenu_Elektroband.Value = handles.default.Optionen.Elektroband_Value;
handles.edit_u_max.String = num2str(handles.default.Regelgroessen.u_max);
handles.edit_i_max.String = num2str(handles.default.Regelgroessen.i_max);

% Read data
handles.Verluste.Stromwaermeverluste = handles.checkbox_Stromwaermeverluste.Value;
handles.Verluste.Eisenverluste = handles.checkbox_Eisenverluste.Value;
handles.Verluste.mechanische_Verluste = handles.checkbox_mechanische_Verluste.Value;
handles.Verluste.Zusatzverluste = handles.checkbox_Zusatzverluste.Value;
handles.Optionen.tics_M = str2double(handles.edit_tics_M.String);
handles.Optionen.tics_omega = str2double(handles.edit_tics_omega.String);
handles.Optionen.n_max_String = handles.popupmenu_n_max.String{handles.popupmenu_n_max.Value};
handles.Optionen.n_max_Value = handles.popupmenu_n_max.Value;
handles.Optionen.n_max = str2double(handles.edit_n_max.String);
handles.Optionen.Elektroband_Value = handles.popupmenu_Elektroband.Value;
handles.Optionen.Elektroband_String = handles.popupmenu_Elektroband.String{handles.popupmenu_Elektroband.Value};
handles.Regelgroessen.u_max = str2double(handles.edit_u_max.String);
handles.Regelgroessen.i_max = str2double(handles.edit_i_max.String);

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in pushbutton_start.
function pushbutton_start_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

uiresume(handles.figure1);
set(findall(handles.uipanel1,'-property','Enable'),'Enable','off')
set(handles.pushbutton_start,'Enable','off')
handles.pushbutton_start.String = 'Bitte warten ...';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SAVE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on button press in pushbutton_save.
function pushbutton_save_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

var = figure;
copyobj(handles.axes_plot, var);
saveas(var, [handles.folder,'/',handles.Plot.Auswahl_Plot], 'fig');
close(var);
