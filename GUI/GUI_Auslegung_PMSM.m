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

function varargout = GUI_Auslegung_PMSM(varargin)
% GUI_Auslegung_PMSM MATLAB code for GUI_Auslegung_PMSM.fig
%      GUI_Auslegung_PMSM, by itself, creates a new GUI_Auslegung_PMSM or raises the existing
%      singleton*.
%
%      H = GUI_Auslegung_PMSM returns the handle to a new GUI_Auslegung_PMSM or the handle to
%      the existing singleton*.
%
%      GUI_Auslegung_PMSM('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_Auslegung_PMSM.M with the given input arguments.
%
%      GUI_Auslegung_PMSM('Property','Value',...) creates a new GUI_Auslegung_PMSM or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI_Auslegung_PMSM before GUI_Auslegung_PMSM_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_Auslegung_PMSM_OpeningFcn via varargin.
%
%      *See GUI_Auslegung_PMSM Options on GUIDE's Tools menu.  Choose "GUI_Auslegung_PMSM allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_Auslegung_PMSM

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_Auslegung_PMSM_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_Auslegung_PMSM_OutputFcn, ...
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

% --- Executes just before GUI_Auslegung_PMSM is made visible.
function GUI_Auslegung_PMSM_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_Auslegung_PMSM (see VARARGIN)

% Choose default command line output for GUI_Auslegung_PMSM
handles.output = hObject;

% Center GUI_Auslegung_PMSM
movegui([100,100]);

% Read data
% Primary Parameters
handles.Primaerparameter.P_N = str2double(handles.edit_P_N.String);
handles.Primaerparameter.n_N = str2double(handles.edit_n_N.String);
handles.Primaerparameter.U_N = str2double(handles.edit_U_N.String);
handles.f_N_p = handles.popupmenu_f_N_p.String{handles.popupmenu_f_N_p.Value};
if(strcmp(handles.f_N_p,'Polpaarzahl'))
    handles.Primaerparameter.p = str2double(handles.edit_f_N_p.String);
    handles.p = handles.Primaerparameter.p;
    %handles.f_N = ceil((handles.p * handles.Primaerparameter.n_N)/60);
else
    handles.Primaerparameter.f_N = str2double(handles.edit_f_N_p.String);
    handles.p = ceil((handles.Primaerparameter.f_N * 60) / handles.Primaerparameter.n_N);
end
handles.Primaerparameter.cos_phi_N = str2double(handles.edit_cos_phi_N.String);
handles.Primaerparameter.m = str2double(handles.edit_m.String);
% Secondary Parameters
handles.Sekundaerparameter.Material_Stator = handles.popupmenu_Material_Stator.String{handles.popupmenu_Material_Stator.Value};
handles.Sekundaerparameter.Schaltung = handles.popupmenu_Schaltung.String;
handles.Sekundaerparameter.Kuehlung_Stator = handles.popupmenu_Kuehlung_Stator.String{handles.popupmenu_Kuehlung_Stator.Value};
handles.Sekundaerparameter.Magnetanordnung = handles.popupmenu_Magnetanordnung.String{handles.popupmenu_Magnetanordnung.Value};
% Approx. Values
lambda = str2double(handles.edit_lambda.String);
if(handles.p > 1)
    if(lambda > 2.5 || lambda < 0.5)
        set(handles.edit_lambda,'BackgroundColor',[1 0.4 0.4]);
    else
        set(handles.edit_lambda,'BackgroundColor',[0.4 1 0.4]);
    end
else
    if(lambda > 4 || lambda < 1.0)
        set(handles.edit_lambda,'BackgroundColor',[1 0.4 0.4]);
    else
        set(handles.edit_lambda,'BackgroundColor',[0.4 1 0.4]);
    end
end
handles.Richtwerte.lambda = lambda;
B_delta = str2double(handles.edit_B_delta.String);
if(B_delta > 1.05 || B_delta < 0.75)
    set(handles.edit_B_delta,'BackgroundColor',[1 0.4 0.4]);
else
    set(handles.edit_B_delta,'BackgroundColor',[0.4 1 0.4]);
end
handles.Richtwerte.B_delta = B_delta;
B_r_max = str2double(handles.edit_B_r_max.String);
if(B_r_max > 1.5 || B_r_max < 1.0)
    set(handles.edit_B_r_max,'BackgroundColor',[1 0.4 0.4]);
else
    set(handles.edit_B_r_max,'BackgroundColor',[0.4 1 0.4]);
end
handles.Richtwerte.B_r_max = B_r_max;
B_z_max = str2double(handles.edit_B_z_max.String);
if(B_z_max > 2.0 || B_z_max < 1.6)
    set(handles.edit_B_z_max,'BackgroundColor',[1 0.4 0.4]);
else
    set(handles.edit_B_z_max,'BackgroundColor',[0.4 1 0.4]);
end
handles.Richtwerte.B_z_max = B_z_max;
phi_n = str2double(handles.edit_phi_n.String);
if(phi_n > 0.5 || phi_n < 0.3)
    set(handles.edit_phi_n,'BackgroundColor',[1 0.4 0.4]);
else
    set(handles.edit_phi_n,'BackgroundColor',[0.4 1 0.4]);
end
handles.Richtwerte.phi_n = phi_n;
xi_p = str2double(handles.edit_xi_p.String);
if(xi_p > 1.0 || xi_p < 0.9)
    set(handles.edit_xi_p,'BackgroundColor',[1 0.4 0.4]);
else
    set(handles.edit_xi_p,'BackgroundColor',[0.4 1 0.4]);
end
handles.Richtwerte.xi_p = xi_p;
S = str2double(handles.edit_S.String);
if(strcmp(handles.Sekundaerparameter.Kuehlung_Stator,'Luft'))
    if(S > 7.0 || S < 3.0)
        set(handles.edit_S,'BackgroundColor',[1 0.4 0.4]);
    else
        set(handles.edit_S,'BackgroundColor',[0.4 1 0.4]);
    end
else
    if(S > 18.0 || S < 7.0)
        set(handles.edit_S,'BackgroundColor',[1 0.4 0.4]);
    else
        set(handles.edit_S,'BackgroundColor',[0.4 1 0.4]);
    end
end
handles.Richtwerte.S = S;
l_v = str2double(handles.edit_l_v.String);
if(l_v > 0.01 || l_v < 0.006)
    set(handles.edit_l_v,'BackgroundColor',[1 0.4 0.4]);
else
    set(handles.edit_l_v,'BackgroundColor',[0.4 1 0.4]);
end
handles.Richtwerte.l_v = l_v;
handles.Richtwerte.Wicklung = handles.popupmenu_Wicklung.String{handles.popupmenu_Wicklung.Value};
tau_n_min = str2double(handles.edit_tau_n_min.String);
if(tau_n_min > 0.07 || tau_n_min < 0.007)
    set(handles.edit_tau_n_min,'BackgroundColor',[1 0.4 0.4]);
else
    set(handles.edit_tau_n_min,'BackgroundColor',[0.4 1 0.4]);
end
handles.Richtwerte.tau_n_min = tau_n_min;
phi_Fe = str2double(handles.edit_phi_Fe.String);
if(phi_Fe > 0.98 || phi_Fe < 0.92)
    set(handles.edit_phi_Fe,'BackgroundColor',[1 0.4 0.4]);
else
    set(handles.edit_phi_Fe,'BackgroundColor',[0.4 1 0.4]);
end
handles.Richtwerte.phi_Fe = phi_Fe;
theta = str2double(handles.edit_theta.String);
if(theta > 120 || theta < 20)
    set(handles.edit_theta,'BackgroundColor',[1 0.4 0.4]);
else
    set(handles.edit_theta,'BackgroundColor',[0.4 1 0.4]);
end
handles.Richtwerte.theta = theta;
handles.Richtwerte.Material_PM = handles.popupmenu_Material_PM.String{handles.popupmenu_Material_PM.Value};
B_r_PM = str2double(handles.edit_B_r_PM.String);
if(strcmp(handles.Richtwerte.Material_PM,'Neodym-Eisen-Bor / s., a.'))
    if(B_r_PM > 1.52 || B_r_PM < 1.05)
        set(handles.edit_B_r_PM,'BackgroundColor',[1 0.4 0.4]);
    else
        set(handles.edit_B_r_PM,'BackgroundColor',[0.4 1 0.4]);
    end
elseif(strcmp(handles.Richtwerte.Material_PM,'Neodym-Eisen-Bor / k., i.'))
    if(B_r_PM > 0.74 || B_r_PM < 0.52)
        set(handles.edit_B_r_PM,'BackgroundColor',[1 0.4 0.4]);
    else
        set(handles.edit_B_r_PM,'BackgroundColor',[0.4 1 0.4]);
    end
elseif(strcmp(handles.Richtwerte.Material_PM,'Neodym-Eisen-Bor / k., a.'))
    if(B_r_PM > 1.00 || B_r_PM < 0.85)
        set(handles.edit_B_r_PM,'BackgroundColor',[1 0.4 0.4]);
    else
        set(handles.edit_B_r_PM,'BackgroundColor',[0.4 1 0.4]);
    end
elseif(strcmp(handles.Richtwerte.Material_PM,'Hartferrite / s., i.'))
    if(B_r_PM > 0.23 || B_r_PM < 0.19)
        set(handles.edit_B_r_PM,'BackgroundColor',[1 0.4 0.4]);
    else
        set(handles.edit_B_r_PM,'BackgroundColor',[0.4 1 0.4]);
    end
elseif(strcmp(handles.Richtwerte.Material_PM,'Hartferrite / s., a.'))
    if(B_r_PM > 0.45 || B_r_PM < 0.38)
        set(handles.edit_B_r_PM,'BackgroundColor',[1 0.4 0.4]);
    else
        set(handles.edit_B_r_PM,'BackgroundColor',[0.4 1 0.4]);
    end
elseif(strcmp(handles.Richtwerte.Material_PM,'Hartferrite / k., i.'))
    if(B_r_PM > 0.15 || B_r_PM < 0.08)
        set(handles.edit_B_r_PM,'BackgroundColor',[1 0.4 0.4]);
    else
        set(handles.edit_B_r_PM,'BackgroundColor',[0.4 1 0.4]);
    end
elseif(strcmp(handles.Richtwerte.Material_PM,'Hartferrite / k., a.'))
    if(B_r_PM > 0.29 || B_r_PM < 0.22)
        set(handles.edit_B_r_PM,'BackgroundColor',[1 0.4 0.4]);
    else
        set(handles.edit_B_r_PM,'BackgroundColor',[0.4 1 0.4]);
    end
elseif(strcmp(handles.Richtwerte.Material_PM,'AlNiCo / g., i.'))
    if(B_r_PM > 0.62 || B_r_PM < 0.58)
        set(handles.edit_B_r_PM,'BackgroundColor',[1 0.4 0.4]);
    else
        set(handles.edit_B_r_PM,'BackgroundColor',[0.4 1 0.4]);
    end
elseif(strcmp(handles.Richtwerte.Material_PM,'AlNiCo / g., a.'))
    if(B_r_PM > 1.35 || B_r_PM < 0.80)
        set(handles.edit_B_r_PM,'BackgroundColor',[1 0.4 0.4]);
    else
        set(handles.edit_B_r_PM,'BackgroundColor',[0.4 1 0.4]);
    end
elseif(strcmp(handles.Richtwerte.Material_PM,'SmCo5 / s., a.'))
    if(B_r_PM > 1.05 || B_r_PM < 0.85)
        set(handles.edit_B_r_PM,'BackgroundColor',[1 0.4 0.4]);
    else
        set(handles.edit_B_r_PM,'BackgroundColor',[0.4 1 0.4]);
    end
elseif(strcmp(handles.Richtwerte.Material_PM,'Sm2Co17 / s., a.'))
    if(B_r_PM > 1.15 || B_r_PM < 0.98)
        set(handles.edit_B_r_PM,'BackgroundColor',[1 0.4 0.4]);
    else
        set(handles.edit_B_r_PM,'BackgroundColor',[0.4 1 0.4]);
    end
end
handles.Richtwerte.B_r_PM = B_r_PM;
mu_mr = str2double(handles.edit_mu_mr.String);
if(strcmp(handles.Richtwerte.Material_PM,'Neodym-Eisen-Bor / s., a.'))
    if(mu_mr ~= 1.05)
        set(handles.edit_mu_mr,'BackgroundColor',[1 0.4 0.4]);
    else
        set(handles.edit_mu_mr,'BackgroundColor',[0.4 1 0.4]);
    end
elseif(strcmp(handles.Richtwerte.Material_PM,'Neodym-Eisen-Bor / k., i.'))
    if(mu_mr ~= 1.20)
        set(handles.edit_mu_mr,'BackgroundColor',[1 0.4 0.4]);
    else
        set(handles.edit_mu_mr,'BackgroundColor',[0.4 1 0.4]);
    end
elseif(strcmp(handles.Richtwerte.Material_PM,'Neodym-Eisen-Bor / k., a.'))
    if(mu_mr > 1.80 || mu_mr < 1.20)
        set(handles.edit_mu_mr,'BackgroundColor',[1 0.4 0.4]);
    else
        set(handles.edit_mu_mr,'BackgroundColor',[0.4 1 0.4]);
    end
elseif(strcmp(handles.Richtwerte.Material_PM,'Hartferrite / s., i.'))
    if(mu_mr ~= 1.25)
        set(handles.edit_mu_mr,'BackgroundColor',[1 0.4 0.4]);
    else
        set(handles.edit_mu_mr,'BackgroundColor',[0.4 1 0.4]);
    end
elseif(strcmp(handles.Richtwerte.Material_PM,'Hartferrite / s., a.'))
    if(mu_mr > 1.10 || mu_mr < 1.05)
        set(handles.edit_mu_mr,'BackgroundColor',[1 0.4 0.4]);
    else
        set(handles.edit_mu_mr,'BackgroundColor',[0.4 1 0.4]);
    end
elseif(strcmp(handles.Richtwerte.Material_PM,'Hartferrite / k., i.'))
    if(mu_mr ~= 1.30)
        set(handles.edit_mu_mr,'BackgroundColor',[1 0.4 0.4]);
    else
        set(handles.edit_mu_mr,'BackgroundColor',[0.4 1 0.4]);
    end
elseif(strcmp(handles.Richtwerte.Material_PM,'Hartferrite / k., a.'))
    if(mu_mr ~= 1.10)
        set(handles.edit_mu_mr,'BackgroundColor',[1 0.4 0.4]);
    else
        set(handles.edit_mu_mr,'BackgroundColor',[0.4 1 0.4]);
    end
elseif(strcmp(handles.Richtwerte.Material_PM,'AlNiCo / g., i.'))
    if(mu_mr > 6.0 || mu_mr < 4.0)
        set(handles.edit_mu_mr,'BackgroundColor',[1 0.4 0.4]);
    else
        set(handles.edit_mu_mr,'BackgroundColor',[0.4 1 0.4]);
    end
elseif(strcmp(handles.Richtwerte.Material_PM,'AlNiCo / g., a.'))
    if(mu_mr > 5.0 || mu_mr < 3.0)
        set(handles.edit_mu_mr,'BackgroundColor',[1 0.4 0.4]);
    else
        set(handles.edit_mu_mr,'BackgroundColor',[0.4 1 0.4]);
    end
elseif(strcmp(handles.Richtwerte.Material_PM,'SmCo5 / s., a.'))
    if(mu_mr ~= 1.08)
        set(handles.edit_mu_mr,'BackgroundColor',[1 0.4 0.4]);
    else
        set(handles.edit_mu_mr,'BackgroundColor',[0.4 1 0.4]);
    end
elseif(strcmp(handles.Richtwerte.Material_PM,'Sm2Co17 / s., a.'))
    if(mu_mr ~= 1.08)
        set(handles.edit_mu_mr,'BackgroundColor',[1 0.4 0.4]);
    else
        set(handles.edit_mu_mr,'BackgroundColor',[0.4 1 0.4]);
    end
end
handles.Richtwerte.mu_mr = mu_mr;
% Options
handles.Optionen.Nutraumbilanz = handles.popupmenu_Nutraumbilanz.String{handles.popupmenu_Nutraumbilanz.Value};
handles.Optionen.Animate_Nut = handles.popupmenu_Animate_Nut.String{handles.popupmenu_Animate_Nut.Value};

% Save default data
handles.default.Richtwerte = handles.Richtwerte;
handles.default.Richtwerte.Wicklung_Value = handles.popupmenu_Wicklung.Value;
handles.default.Richtwerte.Material_PM_Value = handles.popupmenu_Material_PM.Value;

% Set TooltipString
set(handles.edit_lambda, 'TooltipString', ['p=1: zwischen 1.0 und 4.0' 10 'p>1: zwischen 0.5 und 2.5'])
set(handles.edit_B_delta, 'TooltipString', ['zwischen 0.75 und 1.05'])
set(handles.edit_B_r_max, 'TooltipString', ['zwischen 1.0 und 1.5'])
set(handles.edit_B_z_max, 'TooltipString', ['zwischen 1.6 und 2.0'])
set(handles.edit_phi_n, 'TooltipString', ['zwischen 0.3 und 0.5'])
set(handles.edit_xi_p, 'TooltipString', ['zwischen 0.9 und 1.0'])
set(handles.edit_S, 'TooltipString', ['Luftkuehlung: zwischen 3.0 und 7.0' 10 'Wasserkuehlung: zwischen 7.0 und 18.0'])
set(handles.edit_l_v, 'TooltipString', ['zwischen 0.006 und 0.01'])
set(handles.popupmenu_Wicklung, 'TooltipString', ['Wicklung A1: Ganzlochwicklung, Einschicht, ungesehnt, gezont' 10 ...
    'Wicklung A2: Ganzlochwicklung, Zweischicht, ungesehnt, gezont' 10 ...
    'Wicklung B: Ganzlochwicklung, Zweischicht, gesehnt, gezont' 10 ...
    'Wicklung C: Bruchlochwicklung, Zweischicht, gesehnt, gezont,' 10 ...
    'Zahnspulenwicklung (q<1) wird dann umgesetzt, wenn der Durchmesser q>1 nicht zulaesst'])
set(handles.edit_tau_n_min, 'TooltipString', ['zwischen 0.007 und 0.07'])
set(handles.edit_phi_Fe, 'TooltipString', ['zwischen 0.92 und 0.98'])
set(handles.edit_theta, 'TooltipString', ['zwischen 20 und 120'])
set(handles.edit_B_r_PM, 'TooltipString', ['siehe ? Button'])
set(handles.edit_mu_mr, 'TooltipString', ['siehe ? Button'])
set(handles.popupmenu_Nutraumbilanz, 'TooltipString', ['JA: genaue Ausgestaltung der Nutmaße' 10 'NEIN: grobe Abschaetzung der Nutmaße'])
set(handles.popupmenu_Animate_Nut, 'TooltipString', ['JA: Nutgenerierung animieren' 10 'NEIN: Nutgenerierung nicht animieren' 10 'nur verfuegbar wenn Nutraumbilanz JA gewaehlt ist'])

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI_Auslegung_PMSM wait for user response (see UIRESUME)
uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = GUI_Auslegung_PMSM_OutputFcn(hObject, eventdata, handles) 
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
% PRIMARY DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function edit_P_N_Callback(hObject, eventdata, handles)
% hObject    handle to edit_n_N (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of edit_n_N as text
%        str2double(get(hObject,'String')) returns contents of edit_n_N as a double

handles.Primaerparameter.P_N = str2double(get(hObject,'String'));

% Update handles structure
guidata(hObject, handles);

function edit_n_N_Callback(hObject, eventdata, handles)
% hObject    handle to edit_P_N (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of edit_P_N as text
%        str2double(get(hObject,'String')) returns contents of edit_P_N as a double

handles.Primaerparameter.n_N = str2double(get(hObject,'String'));

% Update handles structure
guidata(hObject, handles);

function edit_U_N_Callback(hObject, eventdata, handles)
% hObject    handle to edit_P_N (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of edit_P_N as text
%        str2double(get(hObject,'String')) returns contents of edit_P_N as a double

handles.Primaerparameter.U_N = str2double(get(hObject,'String'));

% Update handles structure
guidata(hObject, handles);

% --- Executes on selection change in popupmenu_f_N_p.
function popupmenu_f_N_p_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_f_N_p (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_f_N_p contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_f_N_p

var = handles.f_N_p;
contents = cellstr(get(hObject,'String'));
handles.f_N_p = contents{get(hObject,'Value')};

if(strcmp(handles.f_N_p,'Polpaarzahl'))
    if(~strcmp(handles.f_N_p,var))
        set(handles.edit_f_N_p,'String', num2str(4.0));
        handles.Primaerparameter.p = str2double(handles.edit_f_N_p.String);
        handles.Primaerparameter = rmfield(handles.Primaerparameter,'f_N');
    end
else
    if(~strcmp(handles.f_N_p,var))
        set(handles.edit_f_N_p,'String', num2str(500.0));
        handles.Primaerparameter.f_N = str2double(handles.edit__p.String);
        handles.Primaerparameter = rmfield(handles.Primaerparameter,'p');
        handles.p = ceil((handles.Primaerparameter.f_N * 60) / handles.Primaerparameter.n_N);
    end
end

lambda = str2double(handles.edit_lambda.String);
if(strcmp(handles.f_N_p,'Polpaarzahl'))
    if(handles.Primaerparameter.p > 1)
        if(lambda > 2.5 || lambda < 0.5)
            set(handles.edit_lambda,'BackgroundColor',[1 0.4 0.4]);
        else
            set(handles.edit_lambda,'BackgroundColor',[0.4 1 0.4]);
        end
    else
        if(lambda > 4 || lambda < 1.0)
            set(handles.edit_lambda,'BackgroundColor',[1 0.4 0.4]);
        else
            set(handles.edit_lambda,'BackgroundColor',[0.4 1 0.4]);
        end
    end
else
    if(handles.p > 1)
        if(lambda > 2.5 || lambda < 0.5)
            set(handles.edit_lambda,'BackgroundColor',[1 0.4 0.4]);
        else
            set(handles.edit_lambda,'BackgroundColor',[0.4 1 0.4]);
        end
    else
        if(lambda > 4 || lambda < 1.0)
            set(handles.edit_lambda,'BackgroundColor',[1 0.4 0.4]);
        else
            set(handles.edit_lambda,'BackgroundColor',[0.4 1 0.4]);
        end
    end
end

% Update handles structure
guidata(hObject, handles);

function edit_f_N_p_Callback(hObject, eventdata, handles)
% hObject    handle to edit_n_N (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of edit_n_N as text
%        str2double(get(hObject,'String')) returns contents of edit_n_N as a double

if(strcmp(handles.f_N_p,'Polpaarzahl'))
    handles.Primaerparameter.p = str2double(get(hObject,'String'));
else
    handles.Primaerparameter.f_N = str2double(get(hObject,'String'));
    handles.p = ceil((handles.Primaerparameter.f_N * 60) / handles.Primaerparameter.n_N);
end

lambda = str2double(handles.edit_lambda.String);
if(strcmp(handles.f_N_p,'Polpaarzahl'))
    if(handles.Primaerparameter.p > 1)
        if(lambda > 2.5 || lambda < 0.5)
            set(handles.edit_lambda,'BackgroundColor',[1 0.4 0.4]);
        else
            set(handles.edit_lambda,'BackgroundColor',[0.4 1 0.4]);
        end
    else
        if(lambda > 4 || lambda < 1.0)
            set(handles.edit_lambda,'BackgroundColor',[1 0.4 0.4]);
        else
            set(handles.edit_lambda,'BackgroundColor',[0.4 1 0.4]);
        end
    end
else
    if(handles.p > 1)
        if(lambda > 2.5 || lambda < 0.5)
            set(handles.edit_lambda,'BackgroundColor',[1 0.4 0.4]);
        else
            set(handles.edit_lambda,'BackgroundColor',[0.4 1 0.4]);
        end
    else
        if(lambda > 4 || lambda < 1.0)
            set(handles.edit_lambda,'BackgroundColor',[1 0.4 0.4]);
        else
            set(handles.edit_lambda,'BackgroundColor',[0.4 1 0.4]);
        end
    end
end

% Update handles structure
guidata(hObject, handles);

function edit_cos_phi_N_Callback(hObject, eventdata, handles)
% hObject    handle to edit_cos_phi_N (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of edit_cos_phi_N as text
%        str2double(get(hObject,'String')) returns contents of edit_cos_phi_N as a double

handles.Primaerparameter.cos_phi_N = str2double(get(hObject,'String'));

% Update handles structure
guidata(hObject, handles);

function edit_m_Callback(hObject, eventdata, handles)
% hObject    handle to edit_m (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of edit_m as text
%        str2double(get(hObject,'String')) returns contents of edit_m as a double

handles.Primaerparameter.m = str2double(get(hObject,'String'));

% Update handles structure
guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SECONDARY DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on selection change in popupmenu_Material_Stator.
function popupmenu_Material_Stator_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_Material_Stator (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_Material_Stator contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_Material_Stator

contents = cellstr(get(hObject,'String'));
handles.Sekundaerparameter.Material_Stator = contents{get(hObject,'Value')};

% Update handles structure
guidata(hObject, handles);

% --- Executes on selection change in popupmenu_Schaltung.
function popupmenu_Schaltung_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_Schaltung (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_Schaltung contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_Schaltung

contents = cellstr(get(hObject,'String'));
handles.Sekundaerparameter.Schaltung = contents{get(hObject,'Value')};

% Update handles structure
guidata(hObject, handles);

% --- Executes on selection change in popupmenu_Kuehlung_Stator.
function popupmenu_Kuehlung_Stator_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_Kuehlung_Stator (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_Kuehlung_Stator contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_Kuehlung_Stator

contents = cellstr(get(hObject,'String'));
handles.Sekundaerparameter.Kuehlung_Stator = contents{get(hObject,'Value')};

S = str2double(handles.edit_S.String);
if(strcmp(handles.Sekundaerparameter.Kuehlung_Stator,'Luft'))
    if(S > 7.0 || S < 3.0)
        set(handles.edit_S,'BackgroundColor',[1 0.4 0.4]);
    else
        set(handles.edit_S,'BackgroundColor',[0.4 1 0.4]);
    end
else
    if(S > 18.0 || S < 7.0)
        set(handles.edit_S,'BackgroundColor',[1 0.4 0.4]);
    else
        set(handles.edit_S,'BackgroundColor',[0.4 1 0.4]);
    end
end

% Update handles structure
guidata(hObject, handles);

% --- Executes on selection change in popupmenu_Magnetanordnung.
function popupmenu_Magnetanordnung_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_Magnetanordnung (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_Magnetanordnung contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_Magnetanordnung

contents = cellstr(get(hObject,'String'));
handles.Sekundaerparameter.Magnetanordnung = contents{get(hObject,'Value')};

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in pushbutton_Magnetanordnung_Rotor.
function pushbutton_Magnetanordnung_Rotor_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Magnetanordnung_Rotor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

figHandles = findall(groot, 'Type', 'figure');
if(any(contains({figHandles.Name},handles.popupmenu_Magnetanordnung.String)))
    var = {figHandles.Name};
    var = var(contains({figHandles.Name},handles.popupmenu_Magnetanordnung.String));
    close(var{1})
end

f = figure('Name',handles.Sekundaerparameter.Magnetanordnung,'NumberTitle','off');
set(f,'MenuBar','none');
set(f,'ToolBar','none');
set(f,'Units','pixels');
axes(f)
imshow([handles.Sekundaerparameter.Magnetanordnung, '.jpg']);
set(f,'Position',[481 390 336 336]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% APPROX. VALUES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function edit_lambda_Callback(hObject, eventdata, handles)
% hObject    handle to edit_lambda (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_lambda as text
%        str2double(get(hObject,'String')) returns contents of edit_lambda as a double

if(strcmp(handles.f_N_p,'Polpaarzahl'))
    p = handles.Primaerparameter.p;
else
    p = handles.p;
end

lambda = str2double(get(hObject,'String'));
if(p > 1)
    if(lambda > 2.5 | lambda < 0.5)
        set(hObject,'BackgroundColor',[1 0.4 0.4]);
    else
        set(hObject,'BackgroundColor',[0.4 1 0.4]);
    end
else
    if(lambda > 4 | lambda < 1.0)
        set(hObject,'BackgroundColor',[1 0.4 0.4]);
    else
        set(hObject,'BackgroundColor',[0.4 1 0.4]);
    end
end

handles.Richtwerte.lambda = lambda;

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in pushbutton_lambda.
function pushbutton_lambda_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_lambda (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

uiwait(msgbox({'The value for p>1 should be between 0.5 and 2.5' 'and für p=1 between 1.0 and 4.0.'},'Relative Ankerlaenge','help','modal'))

function edit_B_delta_Callback(hObject, eventdata, handles)
% hObject    handle to edit_B_delta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_B_delta as text
%        str2double(get(hObject,'String')) returns contents of edit_B_delta as a double

B_delta = str2double(get(hObject,'String'));
if(B_delta > 1.05 | B_delta < 0.75)
    set(hObject,'BackgroundColor',[1 0.4 0.4]);
else
    set(hObject,'BackgroundColor',[0.4 1 0.4]);
end

handles.Richtwerte.B_delta = B_delta;

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in pushbutton_B_delta.
function pushbutton_B_delta_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_B_delta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

uiwait(msgbox({'The value should be between 0.75 and 1.05 .'},'Luftspaltinduktion','help','modal'))

function edit_B_r_max_Callback(hObject, eventdata, handles)
% hObject    handle to edit_B_r_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_B_r_max as text
%        str2double(get(hObject,'String')) returns contents of edit_B_r_max as a double

B_r_max = str2double(get(hObject,'String'));
if(B_r_max > 1.5 | B_r_max < 1.0)
    set(hObject,'BackgroundColor',[1 0.4 0.4]);
else
    set(hObject,'BackgroundColor',[0.4 1 0.4]);
end

handles.Richtwerte.B_r_max = B_r_max;

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in pushbutton_B_r_max.
function pushbutton_B_r_max_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_B_r_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

uiwait(msgbox({'The value should be between 1.0 and 1.5.'},'max. Rueckeninduktion','help','modal'))

function edit_B_z_max_Callback(hObject, eventdata, handles)
% hObject    handle to edit_B_z_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_B_z_max as text
%        str2double(get(hObject,'String')) returns contents of edit_B_z_max as a double

B_z_max = str2double(get(hObject,'String'));
if(B_z_max > 2.0 | B_z_max < 1.6)
    set(hObject,'BackgroundColor',[1 0.4 0.4]);
else
    set(hObject,'BackgroundColor',[0.4 1 0.4]);
end

handles.Richtwerte.B_z_max = B_z_max;

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in pushbutton_B_z_max.
function pushbutton_B_z_max_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_B_z_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

uiwait(msgbox({'The value should be between 1.6 and 2.0.'},'max. Zahninduktion','help','modal'))

function edit_phi_n_Callback(hObject, eventdata, handles)
% hObject    handle to edit_phi_n (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_phi_n as text
%        str2double(get(hObject,'String')) returns contents of edit_phi_n as a double

phi_n = str2double(get(hObject,'String'));
if(phi_n > 0.5 | phi_n < 0.3)
    set(hObject,'BackgroundColor',[1 0.4 0.4]);
else
    set(hObject,'BackgroundColor',[0.4 1 0.4]);
end

handles.Richtwerte.phi_n = phi_n;

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in pushbutton_phi_n.
function pushbutton_phi_n_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_phi_n (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

uiwait(msgbox({'The value should be between 0.3 and 0.5.'},'Nutfuellfaktor','help','modal'))

function edit_xi_p_Callback(hObject, eventdata, handles)
% hObject    handle to edit_xi_p (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_xi_p as text
%        str2double(get(hObject,'String')) returns contents of edit_xi_p as a double

xi_p = str2double(get(hObject,'String'));
if(xi_p > 1.0 | xi_p < 0.9)
    set(hObject,'BackgroundColor',[1 0.4 0.4]);
else
    set(hObject,'BackgroundColor',[0.4 1 0.4]);
end

handles.Richtwerte.xi_p = xi_p;

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in pushbutton_xi_p.
function pushbutton_xi_p_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_xi_p (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

uiwait(msgbox({'The value should be between 0.9 and 1.0.'},'Wicklungsfaktor Grundwelle','help','modal'))

function edit_S_Callback(hObject, eventdata, handles)
% hObject    handle to edit_S (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_S as text
%        str2double(get(hObject,'String')) returns contents of edit_S as a double

Kuehlung_Stator = handles.Sekundaerparameter.Kuehlung_Stator;
S = str2double(get(hObject,'String'));
if(strcmp(Kuehlung_Stator,'Luft'))
    if(S > 7.0 | S < 3.0)
        set(hObject,'BackgroundColor',[1 0.4 0.4]);
    else
        set(hObject,'BackgroundColor',[0.4 1 0.4]);
    end
else
    if(S > 18.0 | S < 7.0)
        set(hObject,'BackgroundColor',[1 0.4 0.4]);
    else
        set(hObject,'BackgroundColor',[0.4 1 0.4]);
    end
end

handles.Richtwerte.S = S;

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in pushbutton_S.
function pushbutton_S_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_S (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

uiwait(msgbox({'The value should be between 3.0 und 7.0 for air cooling' 'and for water cooling between 7.0 und 18.0.'},'Stromdichte','help','modal'))

function edit_l_v_Callback(hObject, eventdata, handles)
% hObject    handle to edit_l_v (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_l_v as text
%        str2double(get(hObject,'String')) returns contents of edit_l_v as a double

l_v = str2double(get(hObject,'String'));
if(l_v > 0.01 | l_v < 0.006)
    set(hObject,'BackgroundColor',[1 0.4 0.4]);
else
    set(hObject,'BackgroundColor',[0.4 1 0.4]);
end

handles.Richtwerte.l_v = l_v;

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in pushbutton_l_v.
function pushbutton_l_v_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_l_v (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

uiwait(msgbox({'The value should be between 0.006 and 0.01.'},'Kanalbreite der Ventilationskanaele','help','modal'))

% --- Executes on selection change in popupmenu_Wicklung.
function popupmenu_Wicklung_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_Wicklung (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_Wicklung contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_Wicklung

contents = cellstr(get(hObject,'String'));
handles.Richtwerte.Wicklung = contents{get(hObject,'Value')};

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in pushbutton_Wicklung.
function pushbutton_Wicklung_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Wicklung (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

uiwait(msgbox({'Winding A1: integral slot winding, single-layer, non-chorded, zoning' ...
    'Winding A2: integral slot winding, double-layer, non-chorded, zoning' ...
    'Winding B: integral slot winding, double-layer, chorded, zoning' ...
    'Winding C: fractional slot winding, double-layer, chorded, zoning,' ...
    '          tooth coil winding (q<1) is implemented if diameter coil q>1 not possible'},'Wicklung','help','modal'))

function edit_tau_n_min_Callback(hObject, eventdata, handles)
% hObject    handle to edit_tau_n_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_tau_n_min as text
%        str2double(get(hObject,'String')) returns contents of edit_tau_n_min as a double

tau_n_min = str2double(get(hObject,'String'));
if(tau_n_min > 0.07 | tau_n_min < 0.007)
    set(hObject,'BackgroundColor',[1 0.4 0.4]);
else
    set(hObject,'BackgroundColor',[0.4 1 0.4]);
end

handles.Richtwerte.tau_n_min = tau_n_min;

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in pushbutton_tau_n_min.
function pushbutton_tau_n_min_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_tau_n_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

uiwait(msgbox({'The value should be between 0.007 and 0.07.'},'min. Nutteilung','help','modal'))

function edit_phi_Fe_Callback(hObject, eventdata, handles)
% hObject    handle to edit_phi_Fe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_phi_Fe as text
%        str2double(get(hObject,'String')) returns contents of edit_phi_Fe as a double

phi_Fe = str2double(get(hObject,'String'));
if(phi_Fe > 0.98 | phi_Fe < 0.92)
    set(hObject,'BackgroundColor',[1 0.4 0.4]);
else
    set(hObject,'BackgroundColor',[0.4 1 0.4]);
end

handles.Richtwerte.phi_Fe = phi_Fe;

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in pushbutton_phi_Fe.
function pushbutton_phi_Fe_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_phi_Fe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

uiwait(msgbox({'The value should be between 0.92 and 0.98.'},'Eisenfuellfaktor','help','modal'))

function edit_theta_Callback(hObject, eventdata, handles)
% hObject    handle to edit_theta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_theta as text
%        str2double(get(hObject,'String')) returns contents of edit_theta as a double

theta = str2double(get(hObject,'String'));
if(theta > 120 | theta < 20)
    set(hObject,'BackgroundColor',[1 0.4 0.4]);
else
    set(hObject,'BackgroundColor',[0.4 1 0.4]);
end

handles.Richtwerte.theta = theta;

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in pushbutton_theta.
function pushbutton_theta_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_theta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

uiwait(msgbox({'temperature for the calculation of the specific' 'conductivity of the selected conductor material'},'Temperatur Wicklung','help','modal'))

% --- Executes on selection change in popupmenu_Material_PM.
function popupmenu_Material_PM_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_Material_PM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_Material_PM contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_Material_PM

contents = cellstr(get(hObject,'String'));
handles.Richtwerte.Material_PM = contents{get(hObject,'Value')};

B_r_PM = str2double(handles.edit_B_r_PM.String);
if(strcmp(handles.Richtwerte.Material_PM,'Neodym-Eisen-Bor / s., a.'))
    if(B_r_PM > 1.52 || B_r_PM < 1.05)
        set(handles.edit_B_r_PM,'BackgroundColor',[1 0.4 0.4]);
    else
        set(handles.edit_B_r_PM,'BackgroundColor',[0.4 1 0.4]);
    end
elseif(strcmp(handles.Richtwerte.Material_PM,'Neodym-Eisen-Bor / k., i.'))
    if(B_r_PM > 0.74 || B_r_PM < 0.52)
        set(handles.edit_B_r_PM,'BackgroundColor',[1 0.4 0.4]);
    else
        set(handles.edit_B_r_PM,'BackgroundColor',[0.4 1 0.4]);
    end
elseif(strcmp(handles.Richtwerte.Material_PM,'Neodym-Eisen-Bor / k., a.'))
    if(B_r_PM > 1.00 || B_r_PM < 0.85)
        set(handles.edit_B_r_PM,'BackgroundColor',[1 0.4 0.4]);
    else
        set(handles.edit_B_r_PM,'BackgroundColor',[0.4 1 0.4]);
    end
elseif(strcmp(handles.Richtwerte.Material_PM,'Hartferrite / s., i.'))
    if(B_r_PM > 0.23 || B_r_PM < 0.19)
        set(handles.edit_B_r_PM,'BackgroundColor',[1 0.4 0.4]);
    else
        set(handles.edit_B_r_PM,'BackgroundColor',[0.4 1 0.4]);
    end
elseif(strcmp(handles.Richtwerte.Material_PM,'Hartferrite / s., a.'))
    if(B_r_PM > 0.45 || B_r_PM < 0.38)
        set(handles.edit_B_r_PM,'BackgroundColor',[1 0.4 0.4]);
    else
        set(handles.edit_B_r_PM,'BackgroundColor',[0.4 1 0.4]);
    end
elseif(strcmp(handles.Richtwerte.Material_PM,'Hartferrite / k., i.'))
    if(B_r_PM > 0.15 || B_r_PM < 0.08)
        set(handles.edit_B_r_PM,'BackgroundColor',[1 0.4 0.4]);
    else
        set(handles.edit_B_r_PM,'BackgroundColor',[0.4 1 0.4]);
    end
elseif(strcmp(handles.Richtwerte.Material_PM,'Hartferrite / k., a.'))
    if(B_r_PM > 0.29 || B_r_PM < 0.22)
        set(handles.edit_B_r_PM,'BackgroundColor',[1 0.4 0.4]);
    else
        set(handles.edit_B_r_PM,'BackgroundColor',[0.4 1 0.4]);
    end
elseif(strcmp(handles.Richtwerte.Material_PM,'AlNiCo / g., i.'))
    if(B_r_PM > 0.62 || B_r_PM < 0.58)
        set(handles.edit_B_r_PM,'BackgroundColor',[1 0.4 0.4]);
    else
        set(handles.edit_B_r_PM,'BackgroundColor',[0.4 1 0.4]);
    end
elseif(strcmp(handles.Richtwerte.Material_PM,'AlNiCo / g., a.'))
    if(B_r_PM > 1.35 || B_r_PM < 0.80)
        set(handles.edit_B_r_PM,'BackgroundColor',[1 0.4 0.4]);
    else
        set(handles.edit_B_r_PM,'BackgroundColor',[0.4 1 0.4]);
    end
elseif(strcmp(handles.Richtwerte.Material_PM,'SmCo5 / s., a.'))
    if(B_r_PM > 1.05 || B_r_PM < 0.85)
        set(handles.edit_B_r_PM,'BackgroundColor',[1 0.4 0.4]);
    else
        set(handles.edit_B_r_PM,'BackgroundColor',[0.4 1 0.4]);
    end
elseif(strcmp(handles.Richtwerte.Material_PM,'Sm2Co17 / s., a.'))
    if(B_r_PM > 1.15 || B_r_PM < 0.98)
        set(handles.edit_B_r_PM,'BackgroundColor',[1 0.4 0.4]);
    else
        set(handles.edit_B_r_PM,'BackgroundColor',[0.4 1 0.4]);
    end
end

mu_mr = str2double(handles.edit_mu_mr.String);
if(strcmp(handles.Richtwerte.Material_PM,'Neodym-Eisen-Bor / s., a.'))
    if(mu_mr ~= 1.05)
        set(handles.edit_mu_mr,'BackgroundColor',[1 0.4 0.4]);
    else
        set(handles.edit_mu_mr,'BackgroundColor',[0.4 1 0.4]);
    end
elseif(strcmp(handles.Richtwerte.Material_PM,'Neodym-Eisen-Bor / k., i.'))
    if(mu_mr ~= 1.20)
        set(handles.edit_mu_mr,'BackgroundColor',[1 0.4 0.4]);
    else
        set(handles.edit_mu_mr,'BackgroundColor',[0.4 1 0.4]);
    end
elseif(strcmp(handles.Richtwerte.Material_PM,'Neodym-Eisen-Bor / k., a.'))
    if(mu_mr > 1.80 || mu_mr < 1.20)
        set(handles.edit_mu_mr,'BackgroundColor',[1 0.4 0.4]);
    else
        set(handles.edit_mu_mr,'BackgroundColor',[0.4 1 0.4]);
    end
elseif(strcmp(handles.Richtwerte.Material_PM,'Hartferrite / s., i.'))
    if(mu_mr ~= 1.25)
        set(handles.edit_mu_mr,'BackgroundColor',[1 0.4 0.4]);
    else
        set(handles.edit_mu_mr,'BackgroundColor',[0.4 1 0.4]);
    end
elseif(strcmp(handles.Richtwerte.Material_PM,'Hartferrite / s., a.'))
    if(mu_mr > 1.10 || mu_mr < 1.05)
        set(handles.edit_mu_mr,'BackgroundColor',[1 0.4 0.4]);
    else
        set(handles.edit_mu_mr,'BackgroundColor',[0.4 1 0.4]);
    end
elseif(strcmp(handles.Richtwerte.Material_PM,'Hartferrite / k., i.'))
    if(mu_mr ~= 1.30)
        set(handles.edit_mu_mr,'BackgroundColor',[1 0.4 0.4]);
    else
        set(handles.edit_mu_mr,'BackgroundColor',[0.4 1 0.4]);
    end
elseif(strcmp(handles.Richtwerte.Material_PM,'Hartferrite / k., a.'))
    if(mu_mr ~= 1.10)
        set(handles.edit_mu_mr,'BackgroundColor',[1 0.4 0.4]);
    else
        set(handles.edit_mu_mr,'BackgroundColor',[0.4 1 0.4]);
    end
elseif(strcmp(handles.Richtwerte.Material_PM,'AlNiCo / g., i.'))
    if(mu_mr > 6.0 || mu_mr < 4.0)
        set(handles.edit_mu_mr,'BackgroundColor',[1 0.4 0.4]);
    else
        set(handles.edit_mu_mr,'BackgroundColor',[0.4 1 0.4]);
    end
elseif(strcmp(handles.Richtwerte.Material_PM,'AlNiCo / g., a.'))
    if(mu_mr > 5.0 || mu_mr < 3.0)
        set(handles.edit_mu_mr,'BackgroundColor',[1 0.4 0.4]);
    else
        set(handles.edit_mu_mr,'BackgroundColor',[0.4 1 0.4]);
    end
elseif(strcmp(handles.Richtwerte.Material_PM,'SmCo5 / s., a.'))
    if(mu_mr ~= 1.08)
        set(handles.edit_mu_mr,'BackgroundColor',[1 0.4 0.4]);
    else
        set(handles.edit_mu_mr,'BackgroundColor',[0.4 1 0.4]);
    end
elseif(strcmp(handles.Richtwerte.Material_PM,'Sm2Co17 / s., a.'))
    if(mu_mr ~= 1.08)
        set(handles.edit_mu_mr,'BackgroundColor',[1 0.4 0.4]);
    else
        set(handles.edit_mu_mr,'BackgroundColor',[0.4 1 0.4]);
    end
end

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in pushbutton_Material_PM.
function pushbutton_Material_PM_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Material_PM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

uiwait(msgbox({'Abbreviations:' 'g. = molded' 'k. = plastic-bound' 's. = sintered' 'a. = anisotrop' 'i. = isotrop'},'Material PM','help','modal'))

function edit_B_r_PM_Callback(hObject, eventdata, handles)
% hObject    handle to edit_B_r_PM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of edit_B_r_PM as text
%        str2double(get(hObject,'String')) returns contents of edit_B_r_PM as a double

B_r_PM = str2double(get(hObject,'String'));
if(strcmp(handles.Richtwerte.Material_PM,'Neodym-Eisen-Bor / s., a.'))
    if(B_r_PM > 1.52 || B_r_PM < 1.05)
        set(handles.edit_B_r_PM,'BackgroundColor',[1 0.4 0.4]);
    else
        set(handles.edit_B_r_PM,'BackgroundColor',[0.4 1 0.4]);
    end
elseif(strcmp(handles.Richtwerte.Material_PM,'Neodym-Eisen-Bor / k., i.'))
    if(B_r_PM > 0.74 || B_r_PM < 0.52)
        set(handles.edit_B_r_PM,'BackgroundColor',[1 0.4 0.4]);
    else
        set(handles.edit_B_r_PM,'BackgroundColor',[0.4 1 0.4]);
    end
elseif(strcmp(handles.Richtwerte.Material_PM,'Neodym-Eisen-Bor / k., a.'))
    if(B_r_PM > 1.00 || B_r_PM < 0.85)
        set(handles.edit_B_r_PM,'BackgroundColor',[1 0.4 0.4]);
    else
        set(handles.edit_B_r_PM,'BackgroundColor',[0.4 1 0.4]);
    end
elseif(strcmp(handles.Richtwerte.Material_PM,'Hartferrite / s., i.'))
    if(B_r_PM > 0.23 || B_r_PM < 0.19)
        set(handles.edit_B_r_PM,'BackgroundColor',[1 0.4 0.4]);
    else
        set(handles.edit_B_r_PM,'BackgroundColor',[0.4 1 0.4]);
    end
elseif(strcmp(handles.Richtwerte.Material_PM,'Hartferrite / s., a.'))
    if(B_r_PM > 0.45 || B_r_PM < 0.38)
        set(handles.edit_B_r_PM,'BackgroundColor',[1 0.4 0.4]);
    else
        set(handles.edit_B_r_PM,'BackgroundColor',[0.4 1 0.4]);
    end
elseif(strcmp(handles.Richtwerte.Material_PM,'Hartferrite / k., i.'))
    if(B_r_PM > 0.15 || B_r_PM < 0.08)
        set(handles.edit_B_r_PM,'BackgroundColor',[1 0.4 0.4]);
    else
        set(handles.edit_B_r_PM,'BackgroundColor',[0.4 1 0.4]);
    end
elseif(strcmp(handles.Richtwerte.Material_PM,'Hartferrite / k., a.'))
    if(B_r_PM > 0.29 || B_r_PM < 0.22)
        set(handles.edit_B_r_PM,'BackgroundColor',[1 0.4 0.4]);
    else
        set(handles.edit_B_r_PM,'BackgroundColor',[0.4 1 0.4]);
    end
elseif(strcmp(handles.Richtwerte.Material_PM,'AlNiCo / g., i.'))
    if(B_r_PM > 0.62 || B_r_PM < 0.58)
        set(handles.edit_B_r_PM,'BackgroundColor',[1 0.4 0.4]);
    else
        set(handles.edit_B_r_PM,'BackgroundColor',[0.4 1 0.4]);
    end
elseif(strcmp(handles.Richtwerte.Material_PM,'AlNiCo / g., a.'))
    if(B_r_PM > 1.35 || B_r_PM < 0.80)
        set(handles.edit_B_r_PM,'BackgroundColor',[1 0.4 0.4]);
    else
        set(handles.edit_B_r_PM,'BackgroundColor',[0.4 1 0.4]);
    end
elseif(strcmp(handles.Richtwerte.Material_PM,'SmCo5 / s., a.'))
    if(B_r_PM > 1.05 || B_r_PM < 0.85)
        set(handles.edit_B_r_PM,'BackgroundColor',[1 0.4 0.4]);
    else
        set(handles.edit_B_r_PM,'BackgroundColor',[0.4 1 0.4]);
    end
elseif(strcmp(handles.Richtwerte.Material_PM,'Sm2Co17 / s., a.'))
    if(B_r_PM > 1.15 || B_r_PM < 0.98)
        set(handles.edit_B_r_PM,'BackgroundColor',[1 0.4 0.4]);
    else
        set(handles.edit_B_r_PM,'BackgroundColor',[0.4 1 0.4]);
    end
end

handles.Richtwerte.B_r_PM = B_r_PM;

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in pushbutton_B_r_PM.
function pushbutton_B_r_PM_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_B_r_PM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

uiwait(msgbox({'for Neodym-Eisen-Bor / s., a. between 1.05 and 1.52'...
    'for Neodym-Eisen-Bor / k., i. between 0.52 and 0.74'...
    'for Neodym-Eisen-Bor / k., a. between 0.85 and 1.00'...
    'for Hartferrite / s., i. between 0.19 and 0.23'...
    'for Hartferrite / s., a. between 0.38 and 0.45'...
    'for Hartferrite / k., i. between 0.08 and 0.15'...
    'for Hartferrite / k., a. between 0.22 and 0.29'...
    'for AlNiCo / g., i. between 0.58 and 0.62'...
    'for AlNiCo / g., a. between 0.80 and 1.35'...
    'for SmCo5 / s., a. between 0.85 and 1.05'...
    'for Sm2Co17 / s., a. between 0.98 and 1.15'},'Remanenzinduktion PM','help','modal'))

function edit_mu_mr_Callback(hObject, eventdata, handles)
% hObject    handle to edit_mu_mr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of edit_mu_mr as text
%        str2double(get(hObject,'String')) returns contents of edit_mu_mr as a double

mu_mr = str2double(get(hObject,'String'));
if(strcmp(handles.Richtwerte.Material_PM,'Neodym-Eisen-Bor / s., a.'))
    if(mu_mr ~= 1.05)
        set(handles.edit_mu_mr,'BackgroundColor',[1 0.4 0.4]);
    else
        set(handles.edit_mu_mr,'BackgroundColor',[0.4 1 0.4]);
    end
elseif(strcmp(handles.Richtwerte.Material_PM,'Neodym-Eisen-Bor / k., i.'))
    if(mu_mr ~= 1.20)
        set(handles.edit_mu_mr,'BackgroundColor',[1 0.4 0.4]);
    else
        set(handles.edit_mu_mr,'BackgroundColor',[0.4 1 0.4]);
    end
elseif(strcmp(handles.Richtwerte.Material_PM,'Neodym-Eisen-Bor / k., a.'))
    if(mu_mr > 1.80 || mu_mr < 1.20)
        set(handles.edit_mu_mr,'BackgroundColor',[1 0.4 0.4]);
    else
        set(handles.edit_mu_mr,'BackgroundColor',[0.4 1 0.4]);
    end
elseif(strcmp(handles.Richtwerte.Material_PM,'Hartferrite / s., i.'))
    if(mu_mr ~= 1.25)
        set(handles.edit_mu_mr,'BackgroundColor',[1 0.4 0.4]);
    else
        set(handles.edit_mu_mr,'BackgroundColor',[0.4 1 0.4]);
    end
elseif(strcmp(handles.Richtwerte.Material_PM,'Hartferrite / s., a.'))
    if(mu_mr > 1.10 || mu_mr < 1.05)
        set(handles.edit_mu_mr,'BackgroundColor',[1 0.4 0.4]);
    else
        set(handles.edit_mu_mr,'BackgroundColor',[0.4 1 0.4]);
    end
elseif(strcmp(handles.Richtwerte.Material_PM,'Hartferrite / k., i.'))
    if(mu_mr ~= 1.30)
        set(handles.edit_mu_mr,'BackgroundColor',[1 0.4 0.4]);
    else
        set(handles.edit_mu_mr,'BackgroundColor',[0.4 1 0.4]);
    end
elseif(strcmp(handles.Richtwerte.Material_PM,'Hartferrite / k., a.'))
    if(mu_mr ~= 1.10)
        set(handles.edit_mu_mr,'BackgroundColor',[1 0.4 0.4]);
    else
        set(handles.edit_mu_mr,'BackgroundColor',[0.4 1 0.4]);
    end
elseif(strcmp(handles.Richtwerte.Material_PM,'AlNiCo / g., i.'))
    if(mu_mr > 6.0 || mu_mr < 4.0)
        set(handles.edit_mu_mr,'BackgroundColor',[1 0.4 0.4]);
    else
        set(handles.edit_mu_mr,'BackgroundColor',[0.4 1 0.4]);
    end
elseif(strcmp(handles.Richtwerte.Material_PM,'AlNiCo / g., a.'))
    if(mu_mr > 5.0 || mu_mr < 3.0)
        set(handles.edit_mu_mr,'BackgroundColor',[1 0.4 0.4]);
    else
        set(handles.edit_mu_mr,'BackgroundColor',[0.4 1 0.4]);
    end
elseif(strcmp(handles.Richtwerte.Material_PM,'SmCo5 / s., a.'))
    if(mu_mr ~= 1.08)
        set(handles.edit_mu_mr,'BackgroundColor',[1 0.4 0.4]);
    else
        set(handles.edit_mu_mr,'BackgroundColor',[0.4 1 0.4]);
    end
elseif(strcmp(handles.Richtwerte.Material_PM,'Sm2Co17 / s., a.'))
    if(mu_mr ~= 1.08)
        set(handles.edit_mu_mr,'BackgroundColor',[1 0.4 0.4]);
    else
        set(handles.edit_mu_mr,'BackgroundColor',[0.4 1 0.4]);
    end
end

handles.Richtwerte.mu_mr = mu_mr;

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in pushbutton_mu_mr.
function pushbutton_mu_mr_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_mu_mr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

uiwait(msgbox({'for Neodym-Eisen-Bor / s., a. 1.05'...
    'for Neodym-Eisen-Bor / k., i. 1.20'...
    'for Neodym-Eisen-Bor / k., a. between 1.20 and 1.80'...
    'for Hartferrite / s., i. 1.25'...
    'for Hartferrite / s., a. between 1.05 and 1.10'...
    'for Hartferrite / k., i. 1.30'...
    'for Hartferrite / k., a. 1.10'...
    'for AlNiCo / g., i. between 4.0 and 6.0'...
    'for AlNiCo / g., a. between 3.0 and 5.0'...
    'for SmCo5 / s., a. 1.08'...
    'for Sm2Co17 / s., a. 1.08'},'Remanenzinduktion PM','help','modal'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OPTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on selection change in popupmenu_Nutraumbilanz.
function popupmenu_Nutraumbilanz_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_Nutraumbilanz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_Nutraumbilanz contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_Nutraumbilanz

contents = cellstr(get(hObject,'String'));
handles.Optionen.Nutraumbilanz = contents{get(hObject,'Value')};

if(strcmp(handles.Optionen.Nutraumbilanz,'JA'))
    set(handles.popupmenu_Animate_Nut,'Enable','on')
else
    set(handles.popupmenu_Animate_Nut,'Enable','off')
end

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in pushbutton_Nutraumbilanz.
function pushbutton_Nutraumbilanz_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Nutraumbilanz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

uiwait(msgbox({'JA: exact design of the slot dimensions' 'NEIN: rough estimation of the groove dimensions'},'Nutraumbilanz','help','modal'))

% --- Executes on selection change in popupmenu_Animate_Nut.
function popupmenu_Animate_Nut_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_Animate_Nut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_Animate_Nut contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_Animate_Nut

contents = cellstr(get(hObject,'String'));
handles.Optionen.Animate_Nut = contents{get(hObject,'Value')};

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in pushbutton_Animate_Nut.
function pushbutton_Animate_Nut_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Animate_Nut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

uiwait(msgbox({'JA: Animate slot generation' 'NEIN: Do not animate groove generation' 'only available if groove space balance YES is selected'},'Animation Nutgenerierung','help','modal'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% START AND RESET
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on button press in pushbutton_reset.
function pushbutton_reset_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Reset to default
handles.edit_lambda.String = num2str(handles.default.Richtwerte.lambda);
handles.edit_B_delta.String = num2str(handles.default.Richtwerte.B_delta);
handles.edit_B_r_max.String = num2str(handles.default.Richtwerte.B_r_max);
handles.edit_B_z_max.String = num2str(handles.default.Richtwerte.B_z_max);
handles.edit_phi_n.String = num2str(handles.default.Richtwerte.phi_n);
handles.edit_xi_p.String = num2str(handles.default.Richtwerte.xi_p);
handles.edit_S.String = num2str(handles.default.Richtwerte.S);
handles.edit_l_v.String = num2str(handles.default.Richtwerte.l_v);
handles.popupmenu_Wicklung.Value = handles.default.Richtwerte.Wicklung_Value;
handles.edit_tau_n_min.String = num2str(handles.default.Richtwerte.tau_n_min);
handles.edit_phi_Fe.String = num2str(handles.default.Richtwerte.phi_Fe);
handles.edit_theta.String = num2str(handles.default.Richtwerte.theta);
handles.popupmenu_Material_PM.Value = handles.default.Richtwerte.Material_PM_Value;
handles.edit_B_r_PM.String = num2str(handles.default.Richtwerte.B_r_PM);
handles.edit_mu_mr.String = num2str(handles.default.Richtwerte.mu_mr);

% Read data
% Primary Parameters
handles.Primaerparameter.P_N = str2double(handles.edit_P_N.String);
handles.Primaerparameter.n_N = str2double(handles.edit_n_N.String);
handles.Primaerparameter.U_N = str2double(handles.edit_U_N.String);
handles.f_N_p = handles.popupmenu_f_N_p.String{handles.popupmenu_f_N_p.Value};
if(strcmp(handles.f_N_p,'Polpaarzahl'))
    handles.Primaerparameter.p = str2double(handles.edit_f_N_p.String);
    handles.p = handles.Primaerparameter.p;
else
    handles.Primaerparameter.f_N = str2double(handles.edit_f_N_p.String);
    handles.p = ceil((handles.Primaerparameter.f_N * 60) / handles.Primaerparameter.n_N);
end
handles.Primaerparameter.cos_phi_N = str2double(handles.edit_cos_phi_N.String);
handles.Primaerparameter.m = str2double(handles.edit_m.String);
% Secondary Parameters
handles.Sekundaerparameter.Material_Stator = handles.popupmenu_Material_Stator.String{handles.popupmenu_Material_Stator.Value};
handles.Sekundaerparameter.Schaltung = handles.popupmenu_Schaltung.String;
handles.Sekundaerparameter.Kuehlung_Stator = handles.popupmenu_Kuehlung_Stator.String{handles.popupmenu_Kuehlung_Stator.Value};
handles.Sekundaerparameter.Magnetanordnung = handles.popupmenu_Magnetanordnung.String{handles.popupmenu_Magnetanordnung.Value};
% Approx. Values
lambda = str2double(handles.edit_lambda.String);
if(handles.p > 1)
    if(lambda > 2.5 || lambda < 0.5)
        set(handles.edit_lambda,'BackgroundColor',[1 0.4 0.4]);
    else
        set(handles.edit_lambda,'BackgroundColor',[0.4 1 0.4]);
    end
else
    if(lambda > 4 || lambda < 1.0)
        set(handles.edit_lambda,'BackgroundColor',[1 0.4 0.4]);
    else
        set(handles.edit_lambda,'BackgroundColor',[0.4 1 0.4]);
    end
end
handles.Richtwerte.lambda = lambda;
B_delta = str2double(handles.edit_B_delta.String);
if(B_delta > 1.05 || B_delta < 0.75)
    set(handles.edit_B_delta,'BackgroundColor',[1 0.4 0.4]);
else
    set(handles.edit_B_delta,'BackgroundColor',[0.4 1 0.4]);
end
handles.Richtwerte.B_delta = B_delta;
B_r_max = str2double(handles.edit_B_r_max.String);
if(B_r_max > 1.5 || B_r_max < 1.0)
    set(handles.edit_B_r_max,'BackgroundColor',[1 0.4 0.4]);
else
    set(handles.edit_B_r_max,'BackgroundColor',[0.4 1 0.4]);
end
handles.Richtwerte.B_r_max = B_r_max;
B_z_max = str2double(handles.edit_B_z_max.String);
if(B_z_max > 2.0 || B_z_max < 1.6)
    set(handles.edit_B_z_max,'BackgroundColor',[1 0.4 0.4]);
else
    set(handles.edit_B_z_max,'BackgroundColor',[0.4 1 0.4]);
end
handles.Richtwerte.B_z_max = B_z_max;
phi_n = str2double(handles.edit_phi_n.String);
if(phi_n > 0.5 || phi_n < 0.3)
    set(handles.edit_phi_n,'BackgroundColor',[1 0.4 0.4]);
else
    set(handles.edit_phi_n,'BackgroundColor',[0.4 1 0.4]);
end
handles.Richtwerte.phi_n = phi_n;
xi_p = str2double(handles.edit_xi_p.String);
if(xi_p > 1.0 || xi_p < 0.9)
    set(handles.edit_xi_p,'BackgroundColor',[1 0.4 0.4]);
else
    set(handles.edit_xi_p,'BackgroundColor',[0.4 1 0.4]);
end
handles.Richtwerte.xi_p = xi_p;
S = str2double(handles.edit_S.String);
if(strcmp(handles.Sekundaerparameter.Kuehlung_Stator,'Luft'))
    if(S > 7.0 || S < 3.0)
        set(handles.edit_S,'BackgroundColor',[1 0.4 0.4]);
    else
        set(handles.edit_S,'BackgroundColor',[0.4 1 0.4]);
    end
else
    if(S > 18.0 || S < 7.0)
        set(handles.edit_S,'BackgroundColor',[1 0.4 0.4]);
    else
        set(handles.edit_S,'BackgroundColor',[0.4 1 0.4]);
    end
end
handles.Richtwerte.S = S;
l_v = str2double(handles.edit_l_v.String);
if(l_v > 0.01 || l_v < 0.006)
    set(handles.edit_l_v,'BackgroundColor',[1 0.4 0.4]);
else
    set(handles.edit_l_v,'BackgroundColor',[0.4 1 0.4]);
end
handles.Richtwerte.l_v = l_v;
handles.Richtwerte.Wicklung = handles.popupmenu_Wicklung.String{handles.popupmenu_Wicklung.Value};
tau_n_min = str2double(handles.edit_tau_n_min.String);
if(tau_n_min > 0.07 || tau_n_min < 0.007)
    set(handles.edit_tau_n_min,'BackgroundColor',[1 0.4 0.4]);
else
    set(handles.edit_tau_n_min,'BackgroundColor',[0.4 1 0.4]);
end
handles.Richtwerte.tau_n_min = tau_n_min;
phi_Fe = str2double(handles.edit_phi_Fe.String);
if(phi_Fe > 0.98 || phi_Fe < 0.92)
    set(handles.edit_phi_Fe,'BackgroundColor',[1 0.4 0.4]);
else
    set(handles.edit_phi_Fe,'BackgroundColor',[0.4 1 0.4]);
end
handles.Richtwerte.phi_Fe = phi_Fe;
theta = str2double(handles.edit_theta.String);
if(theta > 120 || theta < 20)
    set(handles.edit_theta,'BackgroundColor',[1 0.4 0.4]);
else
    set(handles.edit_theta,'BackgroundColor',[0.4 1 0.4]);
end
handles.Richtwerte.theta = theta;
handles.Richtwerte.Material_PM = handles.popupmenu_Material_PM.String{handles.popupmenu_Material_PM.Value};
B_r_PM = str2double(handles.edit_B_r_PM.String);
if(strcmp(handles.Richtwerte.Material_PM,'Neodym-Eisen-Bor / s., a.'))
    if(B_r_PM > 1.52 || B_r_PM < 1.05)
        set(handles.edit_B_r_PM,'BackgroundColor',[1 0.4 0.4]);
    else
        set(handles.edit_B_r_PM,'BackgroundColor',[0.4 1 0.4]);
    end
elseif(strcmp(handles.Richtwerte.Material_PM,'Neodym-Eisen-Bor / k., i.'))
    if(B_r_PM > 0.74 || B_r_PM < 0.52)
        set(handles.edit_B_r_PM,'BackgroundColor',[1 0.4 0.4]);
    else
        set(handles.edit_B_r_PM,'BackgroundColor',[0.4 1 0.4]);
    end
elseif(strcmp(handles.Richtwerte.Material_PM,'Neodym-Eisen-Bor / k., a.'))
    if(B_r_PM > 1.00 || B_r_PM < 0.85)
        set(handles.edit_B_r_PM,'BackgroundColor',[1 0.4 0.4]);
    else
        set(handles.edit_B_r_PM,'BackgroundColor',[0.4 1 0.4]);
    end
elseif(strcmp(handles.Richtwerte.Material_PM,'Hartferrite / s., i.'))
    if(B_r_PM > 0.23 || B_r_PM < 0.19)
        set(handles.edit_B_r_PM,'BackgroundColor',[1 0.4 0.4]);
    else
        set(handles.edit_B_r_PM,'BackgroundColor',[0.4 1 0.4]);
    end
elseif(strcmp(handles.Richtwerte.Material_PM,'Hartferrite / s., a.'))
    if(B_r_PM > 0.45 || B_r_PM < 0.38)
        set(handles.edit_B_r_PM,'BackgroundColor',[1 0.4 0.4]);
    else
        set(handles.edit_B_r_PM,'BackgroundColor',[0.4 1 0.4]);
    end
elseif(strcmp(handles.Richtwerte.Material_PM,'Hartferrite / k., i.'))
    if(B_r_PM > 0.15 || B_r_PM < 0.08)
        set(handles.edit_B_r_PM,'BackgroundColor',[1 0.4 0.4]);
    else
        set(handles.edit_B_r_PM,'BackgroundColor',[0.4 1 0.4]);
    end
elseif(strcmp(handles.Richtwerte.Material_PM,'Hartferrite / k., a.'))
    if(B_r_PM > 0.29 || B_r_PM < 0.22)
        set(handles.edit_B_r_PM,'BackgroundColor',[1 0.4 0.4]);
    else
        set(handles.edit_B_r_PM,'BackgroundColor',[0.4 1 0.4]);
    end
elseif(strcmp(handles.Richtwerte.Material_PM,'AlNiCo / g., i.'))
    if(B_r_PM > 0.62 || B_r_PM < 0.58)
        set(handles.edit_B_r_PM,'BackgroundColor',[1 0.4 0.4]);
    else
        set(handles.edit_B_r_PM,'BackgroundColor',[0.4 1 0.4]);
    end
elseif(strcmp(handles.Richtwerte.Material_PM,'AlNiCo / g., a.'))
    if(B_r_PM > 1.35 || B_r_PM < 0.80)
        set(handles.edit_B_r_PM,'BackgroundColor',[1 0.4 0.4]);
    else
        set(handles.edit_B_r_PM,'BackgroundColor',[0.4 1 0.4]);
    end
elseif(strcmp(handles.Richtwerte.Material_PM,'SmCo5 / s., a.'))
    if(B_r_PM > 1.05 || B_r_PM < 0.85)
        set(handles.edit_B_r_PM,'BackgroundColor',[1 0.4 0.4]);
    else
        set(handles.edit_B_r_PM,'BackgroundColor',[0.4 1 0.4]);
    end
elseif(strcmp(handles.Richtwerte.Material_PM,'Sm2Co17 / s., a.'))
    if(B_r_PM > 1.15 || B_r_PM < 0.98)
        set(handles.edit_B_r_PM,'BackgroundColor',[1 0.4 0.4]);
    else
        set(handles.edit_B_r_PM,'BackgroundColor',[0.4 1 0.4]);
    end
end
mu_mr = str2double(handles.edit_mu_mr.String);
if(strcmp(handles.Richtwerte.Material_PM,'Neodym-Eisen-Bor / s., a.'))
    if(mu_mr ~= 1.05)
        set(handles.edit_mu_mr,'BackgroundColor',[1 0.4 0.4]);
    else
        set(handles.edit_mu_mr,'BackgroundColor',[0.4 1 0.4]);
    end
elseif(strcmp(handles.Richtwerte.Material_PM,'Neodym-Eisen-Bor / k., i.'))
    if(mu_mr ~= 1.20)
        set(handles.edit_mu_mr,'BackgroundColor',[1 0.4 0.4]);
    else
        set(handles.edit_mu_mr,'BackgroundColor',[0.4 1 0.4]);
    end
elseif(strcmp(handles.Richtwerte.Material_PM,'Neodym-Eisen-Bor / k., a.'))
    if(mu_mr > 1.80 || mu_mr < 1.20)
        set(handles.edit_mu_mr,'BackgroundColor',[1 0.4 0.4]);
    else
        set(handles.edit_mu_mr,'BackgroundColor',[0.4 1 0.4]);
    end
elseif(strcmp(handles.Richtwerte.Material_PM,'Hartferrite / s., i.'))
    if(mu_mr ~= 1.25)
        set(handles.edit_mu_mr,'BackgroundColor',[1 0.4 0.4]);
    else
        set(handles.edit_mu_mr,'BackgroundColor',[0.4 1 0.4]);
    end
elseif(strcmp(handles.Richtwerte.Material_PM,'Hartferrite / s., a.'))
    if(mu_mr > 1.10 || mu_mr < 1.05)
        set(handles.edit_mu_mr,'BackgroundColor',[1 0.4 0.4]);
    else
        set(handles.edit_mu_mr,'BackgroundColor',[0.4 1 0.4]);
    end
elseif(strcmp(handles.Richtwerte.Material_PM,'Hartferrite / k., i.'))
    if(mu_mr ~= 1.30)
        set(handles.edit_mu_mr,'BackgroundColor',[1 0.4 0.4]);
    else
        set(handles.edit_mu_mr,'BackgroundColor',[0.4 1 0.4]);
    end
elseif(strcmp(handles.Richtwerte.Material_PM,'Hartferrite / k., a.'))
    if(mu_mr ~= 1.10)
        set(handles.edit_mu_mr,'BackgroundColor',[1 0.4 0.4]);
    else
        set(handles.edit_mu_mr,'BackgroundColor',[0.4 1 0.4]);
    end
elseif(strcmp(handles.Richtwerte.Material_PM,'AlNiCo / g., i.'))
    if(mu_mr > 6.0 || mu_mr < 4.0)
        set(handles.edit_mu_mr,'BackgroundColor',[1 0.4 0.4]);
    else
        set(handles.edit_mu_mr,'BackgroundColor',[0.4 1 0.4]);
    end
elseif(strcmp(handles.Richtwerte.Material_PM,'AlNiCo / g., a.'))
    if(mu_mr > 5.0 || mu_mr < 3.0)
        set(handles.edit_mu_mr,'BackgroundColor',[1 0.4 0.4]);
    else
        set(handles.edit_mu_mr,'BackgroundColor',[0.4 1 0.4]);
    end
elseif(strcmp(handles.Richtwerte.Material_PM,'SmCo5 / s., a.'))
    if(mu_mr ~= 1.08)
        set(handles.edit_mu_mr,'BackgroundColor',[1 0.4 0.4]);
    else
        set(handles.edit_mu_mr,'BackgroundColor',[0.4 1 0.4]);
    end
elseif(strcmp(handles.Richtwerte.Material_PM,'Sm2Co17 / s., a.'))
    if(mu_mr ~= 1.08)
        set(handles.edit_mu_mr,'BackgroundColor',[1 0.4 0.4]);
    else
        set(handles.edit_mu_mr,'BackgroundColor',[0.4 1 0.4]);
    end
end
% Options
handles.Optionen.Nutraumbilanz = handles.popupmenu_Nutraumbilanz.String{handles.popupmenu_Nutraumbilanz.Value};
handles.Optionen.Animate_Nut = handles.popupmenu_Animate_Nut.String{handles.popupmenu_Animate_Nut.Value};

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in pushbutton_start.
function pushbutton_start_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

figHandles = findall(groot, 'Type', 'figure');
if(any(contains({figHandles.Name},handles.popupmenu_Magnetanordnung.String)))
    var = {figHandles.Name};
    var = var(contains({figHandles.Name},handles.popupmenu_Magnetanordnung.String));
    close(var{1})
end

uiresume(handles.figure1);
set(findall(handles.uipanel1,'-property','Enable'),'Enable','off')
set(findall(handles.uipanel4,'-property','Enable'),'Enable','off')
set(findall(handles.uipanel5,'-property','Enable'),'Enable','off')
set(handles.pushbutton_start,'Enable','off')
if(strcmp(handles.Optionen.Nutraumbilanz,'NEIN'))
    set(handles.text_Plot_Nut,'Visible','on')
end
handles.pushbutton_start.String = 'Bitte warten ...';
