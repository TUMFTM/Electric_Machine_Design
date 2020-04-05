% -------------------------------------------------------------------------
% TU Muenchen - Lehrstuhl fuer Fahrzeugtechnik (FTM)
% -------------------------------------------------------------------------
% Modell fuer den Entwurf und die Analyse einer PMSM oder ASM (MEAPA)
% -------------------------------------------------------------------------
% Autor: Svenja Kalt (kalt@ftm.mw.tum.de)
%        Jonathan Erhard
% -------------------------------------------------------------------------

function varargout = GUI_Entwurf(varargin)
% GUI_ENTWURF MATLAB code for GUI_Entwurf.fig
%      GUI_ENTWURF, by itself, creates a new GUI_ENTWURF or raises the existing
%      singleton*.
%
%      H = GUI_ENTWURF returns the handle to a new GUI_ENTWURF or the handle to
%      the existing singleton*.
%
%      GUI_ENTWURF('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_ENTWURF.M with the given input arguments.
%
%      GUI_ENTWURF('Property','Value',...) creates a new GUI_ENTWURF or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_Entwurf_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_Entwurf_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_Entwurf

% Last Modified by GUIDE v2.5 23-Mar-2019 15:56:07

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_Entwurf_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_Entwurf_OutputFcn, ...
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

% --- Executes just before GUI_Entwurf is made visible.
function GUI_Entwurf_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_Entwurf (see VARARGIN)

% Font size optimization for Windows
decreaseFontSizesIfReq(handles);

% Choose default command line output for GUI_Entwurf
handles.output = hObject;

% Center GUI_Entwurf
set(handles.figure1,'Position',[100 100 740 885]);

% Load Logos
FTM_Logo = imread('FTM_Logo.png');
axes(handles.axes_FTM_Logo);
imshow(FTM_Logo);

TUM_Logo = imread('TUM_Logo.jpeg');
axes(handles.axes_TUM_Logo);
imshow(TUM_Logo);

handles = newLoad(handles);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI_Entwurf wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = GUI_Entwurf_OutputFcn(hObject, eventdata, handles) 
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

button = questdlg('Wollen Sie den Entwurf und die Analyse beenden? Alle nicht gespeicherten Aenderungen gehen verloren.', 'Beenden','Ja','Nein','Nein'); 
switch button 
    case 'Ja'
        figHandles = findall(groot, 'Type', 'figure');
        if(any(contains({figHandles.Name},'GUI_Analyse')))
            var = {figHandles.Name};
            var = var(contains({figHandles.Name},'GUI_Analyse'));
            close(var{1},'force')
        end
        delete(handles.figure1);
    case 'Nein'
end

% #########################################################################
% #   A) TOOLBAR                                                          #
% #########################################################################

% --------------------------------------------------------------------
function toolbar_NewFile_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to toolbar_NewFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

button = questdlg('Wollen Sie eine neue Datei anlegen? Alle nicht gespeicherten Aenderungen gehen verloren.', 'Neue Datei','Ja','Nein','Nein'); 
switch button 
    case 'Ja'
        if(handles.opt.Locked)
            if(isfield(handles,'Entwurf'))
                handles = rmfield(handles,'Entwurf');
            end
        end
        
        if(isfield(handles,'rated'))
            handles = rmfield(handles,'rated');
        end
        if(isfield(handles,'richt'))
            handles = rmfield(handles,'richt');
        end
        if(isfield(handles,'opt'))
            handles = rmfield(handles,'opt');
        end
        
        handles.popupmenu_Content.Value = 1;
        handles.opt.Content = handles.popupmenu_Content.String{handles.popupmenu_Content.Value};
        handles = newLoad(handles);
        
    case 'Nein'
end

% Update handles structure
guidata(hObject, handles);

% --------------------------------------------------------------------
function toolbar_OpenFile_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to toolbar_OpenFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[file,path] = uigetfile('*.mat',[pwd '/3_Ergebnisse/']);

if(isequal(file,0))
else
    % Loeschen der bisherigen Daten
    if(handles.opt.Locked)
        if(isfield(handles,'Entwurf'))
            handles = rmfield(handles,'Entwurf');
        end
    end
    if(isfield(handles,'rated'))
        handles = rmfield(handles,'rated');
    end
    if(isfield(handles,'richt'))
        handles = rmfield(handles,'richt');
    end
    if(isfield(handles,'opt'))
        handles = rmfield(handles,'opt');
    end
    
    % Laden
    load([path file]);

    % Speichern der geladenen Daten
    handles.rated = Entwurf.Bemessungswerte;
    handles.richt = Entwurf.Richtwerte;
    handles.opt = Entwurf.Optionen;
    if(Entwurf.Optionen.Locked)
        handles.Entwurf = Entwurf;
    end
    
    % Maschinentyp
    idx = find(strcmp(handles.popupmenu_Maschinentyp.String,handles.opt.Maschinentyp));
    handles.popupmenu_Maschinentyp.Value = idx;
    handles.opt.Maschinentyp = handles.popupmenu_Maschinentyp.String{handles.popupmenu_Maschinentyp.Value};
    
    % Load default values
    handles = loadDefault(handles);
    
    % Set default data
    handles = setDefaultOptionen(handles);
    
    % Set from file
    handles = setFromFile(handles);
    
    % Read data
    handles = readOptionen(handles);
    handles = readBemessungswerte(handles);

    % Set GUI
    handles = toggleMaschinentyp(handles);
    handles = toggleContent(handles);

    % Read data
    handles = readRichtwerte(handles);

    handles = toggleLocked(handles);
    handles.text_ID.String = ['ID: Entwurf_' handles.opt.file_id];
end

% Update handles structure
guidata(hObject, handles);

% --------------------------------------------------------------------
function toolbar_SaveFile_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to toolbar_SaveFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if(isfolder(['3_Ergebnisse/' handles.opt.folder_id]))
    button = questdlg('Speichern als ... ?', 'Datei speichern','Neue Datei','Ueberschreiben','Neue Datei'); 
    switch button 
        case 'Neue Datei'
            handles.opt.unique_id = datestr(now,'yyyymmdd_HHMMSS');
            handles.opt.folder_id = [handles.opt.Maschinentyp '_data_' handles.opt.unique_id];
            handles.opt.file_id = [handles.opt.Maschinentyp,'_',handles.opt.unique_id];
            handles.text_ID.String = ['ID: Entwurf_' handles.opt.file_id];
            mkdir('3_Ergebnisse',handles.opt.folder_id);
            mkdir(['3_Ergebnisse/',handles.opt.folder_id],'1_Entwurf');
            mkdir(['3_Ergebnisse/',handles.opt.folder_id],'2_Analyse');
        case 'Ueberschreiben'  
    end
else
    mkdir('3_Ergebnisse',handles.opt.folder_id);
    mkdir(['3_Ergebnisse/',handles.opt.folder_id],'1_Entwurf');
    mkdir(['3_Ergebnisse/',handles.opt.folder_id],'2_Analyse');
end

if(handles.opt.Locked)
    handles.opt.Saved = 1;
    handles.Entwurf.Optionen.Saved = handles.opt.Saved;
    handles.Entwurf.Optionen.unique_id = handles.opt.unique_id;
    handles.Entwurf.Optionen.folder_id = handles.opt.folder_id;
    handles.Entwurf.Optionen.file_id = handles.opt.file_id;
    Entwurf = handles.Entwurf;
else
    Entwurf.Bemessungswerte = handles.rated;
    Entwurf.Richtwerte = handles.richt;
    Entwurf.Optionen = handles.opt;
end

try
    % Save struct Entwurf
    save(['3_Ergebnisse/',handles.opt.folder_id,'/1_Entwurf/Entwurf_',handles.opt.file_id,'.mat'],'Entwurf');
    % Save DXF file
    if(handles.opt.Locked)
        saveDXF(handles,handles.Entwurf.Geometrie)
    end
    % Save Excel file
    saveExcel(handles, Entwurf)
    
    uiwait(msgbox({'Speichern erfolgreich.'},'Datei speichern','help','modal'))
catch
    uiwait(msgbox({'Speichern fehlgeschlagen.'},'Datei speichern','help','error'))
end

% Update handles structure
guidata(hObject, handles);

% --------------------------------------------------------------------
function toolbar_UnlockFile_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to toolbar_UnlockFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if(handles.opt.Locked)
    button = questdlg('Wollen Sie die Datei entsperren? Es wird eine neue Datei mit den bisherigen Eingabedaten erzeugt. Alle Ergebnisse (inkl. Analyse) werden dabei zurueckgesetzt.', 'Datei entsperren','Ja','Nein','Nein'); 
    switch button 
        case 'Ja'
            figHandles = findall(groot, 'Type', 'figure');
            if(any(contains({figHandles.Name},'GUI_Analyse')))
                var = {figHandles.Name};
                var = var(contains({figHandles.Name},'GUI_Analyse'));
                close(var{1},'force')
            end
            
            if(handles.opt.Locked)
                if(isfield(handles,'Entwurf'))
                    handles = rmfield(handles,'Entwurf');
                end
            end
            
            handles.opt.Locked = 0;
            handles.opt.Saved = 0;
            
            % Maschinentyp
            handles.opt.Maschinentyp = handles.popupmenu_Maschinentyp.String{handles.popupmenu_Maschinentyp.Value};

            % Read data
            handles = readOptionen(handles);
            handles = readBemessungswerte(handles);
            handles = readRichtwerte(handles);

            % Set GUI
            handles = toggleMaschinentyp(handles);
            handles.popupmenu_Content.Value = 1;
            handles.opt.Content = handles.popupmenu_Content.String{handles.popupmenu_Content.Value};
            handles = toggleContent(handles);            

            handles = toggleLocked(handles);
            handles.opt.unique_id = datestr(now,'yyyymmdd_HHMMSS');
            handles.opt.folder_id = [handles.opt.Maschinentyp '_data_' handles.opt.unique_id];
            handles.opt.file_id = [handles.opt.Maschinentyp,'_',handles.opt.unique_id];
            handles.text_ID.String = ['ID: Entwurf_' handles.opt.file_id];
            
        case 'Nein'
    end
end

% Update handles structure
guidata(hObject, handles);

% #########################################################################
% #   B) ENTWURF                                                          #
% #########################################################################

% --- Executes on selection change in popupmenu_Maschinentyp.
function popupmenu_Maschinentyp_Callback(hObject, eventdata, handles)

handles = newLoad(handles);

figHandles = findall(groot, 'Type', 'figure');
if(any(contains({figHandles.Name},handles.popupmenu_Maschinenausfuehrung.String)))
    var = {figHandles.Name};
    var = var(contains({figHandles.Name},handles.popupmenu_Maschinenausfuehrung.String));
    close(var{1})
end

% Update handles structure
guidata(hObject, handles);

% --- Executes on selection change in popupmenu_Content.
function popupmenu_Content_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_Content (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get popupmenu contents
contents = cellstr(get(hObject,'String'));
handles.opt.Content = contents{get(hObject,'Value')};

handles = toggleContent(handles);

% Update handles structure
guidata(hObject, handles);

% #########################################################################
% #   C) EINGABE                                                          #
% #########################################################################

% --- Executes on button press in pushbutton_Entwurf.
function pushbutton_Entwurf_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Entwurf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

figHandles = findall(groot, 'Type', 'figure');
if(any(contains({figHandles.Name},handles.popupmenu_Maschinenausfuehrung.String)))
    var = {figHandles.Name};
    var = var(contains({figHandles.Name},handles.popupmenu_Maschinenausfuehrung.String));
    close(var{1})
end

handles.pushbutton_Content2.Enable = 'off';
handles.popupmenu_Content.Enable = 'off';
handles.popupmenu_Content.Value = 2;
handles.opt.Content = handles.popupmenu_Content.String{handles.popupmenu_Content.Value};
handles = toggleContent(handles);

set(findall(handles.uipanel_Eingabe,'-property','Enable'),'Enable','off')
set(handles.popupmenu_Maschinentyp,'Enable','off');
handles.text_Wait.String = [{'Bitte warten ...'}; {'Die Berechnungen werden ausgefuehrt.'}];

pause(0.1)

if(strcmp(handles.opt.Maschinentyp,'ASM'))
    [handles.Entwurf] = Entwurf_ASM(handles);
    %[handles.Entwurf] = Copy_of_Entwurf_ASM(handles);
elseif(strcmp(handles.opt.Maschinentyp,'PMSM'))
    [handles.Entwurf] = Entwurf_PMSM(handles);
else
    error('Ungueltige Eingabe bei Variable "Entwurf.Optionen.Maschinentyp"');
end

% Berechnung erfolgreich
handles.opt.Locked = 1;
handles.text_Wait.String = [{'Success ...'}; {'Die Berechnungen waren erfolgreich.'}];
handles.pushbutton_Content2.Enable = 'on';
handles.popupmenu_Content.Enable = 'on';
pause(1.0)

handles = toggleLocked(handles);

handles = setErgebnisse(handles);

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in pushbutton_Content1.
function pushbutton_Content1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Content1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Content
handles.popupmenu_Content.Value = 2;
handles.opt.Content = handles.popupmenu_Content.String{handles.popupmenu_Content.Value};

handles = toggleContent(handles);

% Update handles structure
guidata(hObject, handles);

% #########################################################################
% #   C.1) BEMESSUNGSWERTE                                                #
% #########################################################################

function edit_P_N_Callback(hObject, eventdata, handles)
% hObject    handle to edit_P_N (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Read data
handles = readBemessungswerte(handles);
handles = readRichtwerte(handles);

% Update handles structure
guidata(hObject, handles);

function edit_n_N_Callback(hObject, eventdata, handles)
% hObject    handle to edit_n_N (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Read data
handles = readBemessungswerte(handles);
handles = readRichtwerte(handles);

% Update handles structure
guidata(hObject, handles);

function edit_U_N_Callback(hObject, eventdata, handles)
% hObject    handle to edit_U_N (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Read data
handles = readBemessungswerte(handles);
handles = readRichtwerte(handles);

% Update handles structure
guidata(hObject, handles);

function edit_p_Callback(hObject, eventdata, handles)
% hObject    handle to edit_p (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Read data
handles = readBemessungswerte(handles);
handles = readRichtwerte(handles);

% Update handles structure
guidata(hObject, handles);

function edit_f_N_Callback(hObject, eventdata, handles)
% hObject    handle to edit_f_N (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.rated.f_N = str2double(handles.edit_f_N.String);
handles.rated.p = round((handles.rated.f_N * 60) / handles.rated.n_N);
handles.rated.f_N = (handles.rated.p * handles.rated.n_N) / 60;
handles.edit_p.String = num2str(handles.rated.p);
handles.edit_f_N.String = num2str(handles.rated.f_N);

% Read data
handles = readBemessungswerte(handles);
handles = readRichtwerte(handles);

% Update handles structure
guidata(hObject, handles);

function edit_cos_phi_N_Callback(hObject, eventdata, handles)
% hObject    handle to edit_cos_phi_N (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Read data
handles = readBemessungswerte(handles);
handles = readRichtwerte(handles);

% Update handles structure
guidata(hObject, handles);

function edit_m_Callback(hObject, eventdata, handles)
% hObject    handle to edit_m (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Read data
handles = readBemessungswerte(handles);
handles = readRichtwerte(handles);

% Update handles structure
guidata(hObject, handles);

% #########################################################################
% #   C.2) RICHTWERTE                                                     #
% #########################################################################

function edit_lambda_Callback(hObject, eventdata, handles)
% hObject    handle to edit_lambda (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Read data
handles = readRichtwerte(handles);

% Update handles structure
guidata(hObject, handles);

function edit_l_v_Callback(hObject, eventdata, handles)
% hObject    handle to edit_l_v (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Read data
handles = readRichtwerte(handles);

% Update handles structure
guidata(hObject, handles);

% --- Executes on selection change in popupmenu_B_A.
function popupmenu_B_A_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_B_A (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Luftspaltinduktion bzw. Strombelag
handles.B_A = handles.popupmenu_B_A.String{handles.popupmenu_B_A.Value};
if(strcmp(handles.opt.Maschinentyp,'ASM'))
    if(strcmp(handles.B_A,'Mittelwert der Luftspaltinduktion in T'))
        handles.edit_B_A.String = num2str(handles.default_richt.B_m(1,3));
    elseif(strcmp(handles.B_A,'Strombelag in A/mm'))
        handles.edit_B_A.String = num2str(handles.default_richt.A(1,3));
    else
        error('Ungueltige Eingabe bei Variable "handles.B_A"');
    end
elseif(strcmp(handles.opt.Maschinentyp,'PMSM'))
    if(strcmp(handles.B_A,'Amplitude der Luftspaltinduktion in T'))
        handles.edit_B_A.String = num2str(handles.default_richt.B_p(1,3));
    elseif(strcmp(handles.B_A,'Strombelag in A/mm'))
        handles.edit_B_A.String = num2str(handles.default_richt.A(1,3));
    else
        error('Ungueltige Eingabe bei Variable "handles.B_A"');
    end
else
    error('Ungueltige Eingabe bei Variable "handles.opt.Maschinentyp"');
end

% Read data
handles = readRichtwerte(handles);

% Update handles structure
guidata(hObject, handles);

function edit_B_A_Callback(hObject, eventdata, handles)
% hObject    handle to edit_B_A (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Read data
handles = readRichtwerte(handles);

% Update handles structure
guidata(hObject, handles);

function edit_S_1_Callback(hObject, eventdata, handles)
% hObject    handle to edit_S_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Read data
handles = readRichtwerte(handles);

% Update handles structure
guidata(hObject, handles);

function edit_B_1r_max_Callback(hObject, eventdata, handles)
% hObject    handle to edit_B_1r_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Read data
handles = readRichtwerte(handles);

% Update handles structure
guidata(hObject, handles);

function edit_B_1z_max_Callback(hObject, eventdata, handles)
% hObject    handle to edit_B_1z_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Read data
handles = readRichtwerte(handles);

% Update handles structure
guidata(hObject, handles);

function edit_phi_1n_Callback(hObject, eventdata, handles)
% hObject    handle to edit_phi_1n (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Read data
handles = readRichtwerte(handles);

% Update handles structure
guidata(hObject, handles);

function edit_xi_1p_Callback(hObject, eventdata, handles)
% hObject    handle to edit_xi_1p (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Read data
handles = readRichtwerte(handles);

% Update handles structure
guidata(hObject, handles);

function edit_tau_1n_min_Callback(hObject, eventdata, handles)
% hObject    handle to edit_tau_1n_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Read data
handles = readRichtwerte(handles);

% Update handles structure
guidata(hObject, handles);

function edit_phi_1Fe_Callback(hObject, eventdata, handles)
% hObject    handle to edit_phi_1Fe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Read data
handles = readRichtwerte(handles);

% Update handles structure
guidata(hObject, handles);

function edit_S_2_Callback(hObject, eventdata, handles)
% hObject    handle to edit_S_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Read data
handles = readRichtwerte(handles);

% Update handles structure
guidata(hObject, handles);

function edit_S_2s_Callback(hObject, eventdata, handles)
% hObject    handle to edit_S_2s (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Read data
handles = readRichtwerte(handles);

% Update handles structure
guidata(hObject, handles);

function edit_S_2r_Callback(hObject, eventdata, handles)
% hObject    handle to edit_S_2r (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Read data
handles = readRichtwerte(handles);

% Update handles structure
guidata(hObject, handles);

function edit_B_2r_max_Callback(hObject, eventdata, handles)
% hObject    handle to edit_B_2r_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Read data
handles = readRichtwerte(handles);

% Update handles structure
guidata(hObject, handles);

function edit_B_2z_max_Callback(hObject, eventdata, handles)
% hObject    handle to edit_B_2z_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Read data
handles = readRichtwerte(handles);

% Update handles structure
guidata(hObject, handles);

function edit_phi_2n_Callback(hObject, eventdata, handles)
% hObject    handle to edit_phi_2n (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Read data
handles = readRichtwerte(handles);

% Update handles structure
guidata(hObject, handles);

function edit_xi_2p_Callback(hObject, eventdata, handles)
% hObject    handle to edit_xi_2p (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Read data
handles = readRichtwerte(handles);

% Update handles structure
guidata(hObject, handles);

function edit_tau_2n_min_Callback(hObject, eventdata, handles)
% hObject    handle to edit_tau_2n_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Read data
handles = readRichtwerte(handles);

% Update handles structure
guidata(hObject, handles);

function edit_phi_2Fe_Callback(hObject, eventdata, handles)
% hObject    handle to edit_phi_2Fe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Read data
handles = readRichtwerte(handles);

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in pushbutton_Reset.
function pushbutton_Reset_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Reset data
handles = setDefaultRichtwerte(handles);

% Read data
handles = readRichtwerte(handles);

% Update handles structure
guidata(hObject, handles);

% #########################################################################
% #   C.3) OPTIONEN MASCHINE                                              #
% #########################################################################

% --- Executes on selection change in popupmenu_Maschinenausfuehrung.
function popupmenu_Maschinenausfuehrung_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_Maschinenausfuehrung (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

figHandles = findall(groot, 'Type', 'figure');
if(any(contains({figHandles.Name},handles.popupmenu_Maschinenausfuehrung.String)))
    var = {figHandles.Name};
    var = var(contains({figHandles.Name},handles.popupmenu_Maschinenausfuehrung.String));
    close(var{1})
end

% Read data
handles = readOptionen(handles);
handles = readRichtwerte(handles);

handles = toggleMaschinenausfuehrung(handles);

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in pushbutton_Maschinenausfuehrung.
function pushbutton_Maschinenausfuehrung_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Maschinenausfuehrung (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

figHandles = findall(groot, 'Type', 'figure');
if(any(contains({figHandles.Name},handles.popupmenu_Maschinenausfuehrung.String)))
    var = {figHandles.Name};
    var = var(contains({figHandles.Name},handles.popupmenu_Maschinenausfuehrung.String));
    close(var{1})
end

f = figure('Name',handles.opt.Maschinenausfuehrung,'NumberTitle','off');
set(f,'Visible','Off')
set(f,'MenuBar','none');
set(f,'ToolBar','none');
set(f,'Units','pixels');
axes(f)
imshow([handles.opt.Maschinenausfuehrung, '.jpg']);
set(f,'Position',[840 680 200 200]);
set(f,'Visible','On')

% --- Executes on selection change in popupmenu_Schaltung.
function popupmenu_Schaltung_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_Schaltung (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Read data
handles = readOptionen(handles);
handles = readRichtwerte(handles);

% Update handles structure
guidata(hObject, handles);

% --- Executes on selection change in popupmenu_Spulenform_Stator.
function popupmenu_Spulenform_Stator_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_Spulenform_Stator (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Read data
handles = readOptionen(handles);
handles = readRichtwerte(handles);

% Update handles structure
guidata(hObject, handles);

% --- Executes on selection change in popupmenu_Spulenform_Rotor.
function popupmenu_Spulenform_Rotor_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_Spulenform_Rotor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Read data
handles = readOptionen(handles);
handles = readRichtwerte(handles);

% Update handles structure
guidata(hObject, handles);

% --- Executes on selection change in popupmenu_Nutform_Stator.
function popupmenu_Nutform_Stator_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_Nutform_Stator (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Read data
handles = readOptionen(handles);
handles = readRichtwerte(handles);

% Update handles structure
guidata(hObject, handles);

% --- Executes on selection change in popupmenu_Nutform_Rotor.
function popupmenu_Nutform_Rotor_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_Nutform_Rotor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Read data
handles = readOptionen(handles);
handles = readRichtwerte(handles);

% Update handles structure
guidata(hObject, handles);

% --- Executes on selection change in popupmenu_Kuehlungsart.
function popupmenu_Kuehlungsart_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_Kuehlungsart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Read data
handles = readOptionen(handles);
handles = readRichtwerte(handles);

% Update handles structure
guidata(hObject, handles);

% --- Executes on selection change in popupmenu_Stator_Eisenmaterial.
function popupmenu_Stator_Eisenmaterial_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_Stator_Eisenmaterial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Read data
handles = readOptionen(handles);
handles = readRichtwerte(handles);

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in pushbutton_Plot_BH_Stator.
function pushbutton_Plot_BH_Stator_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Plot_BH_Stator (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

plotBH(handles.opt.Stator_Eisenmaterial);

% --- Executes on button press in pushbutton_Plot_Eisenverlust_Stator.
function pushbutton_Plot_Eisenverlust_Stator_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Plot_Eisenverlust_Stator (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

plotVerluste(handles.opt.Stator_Eisenmaterial);

% --- Executes on selection change in popupmenu_Stator_Leitermaterial.
function popupmenu_Stator_Leitermaterial_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_Stator_Leitermaterial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Read data
handles = readOptionen(handles);
handles = readRichtwerte(handles);

% Update handles structure
guidata(hObject, handles);

function edit_theta_1_Callback(hObject, eventdata, handles)
% hObject    handle to edit_theta_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Read data
handles = readOptionen(handles);
handles = readRichtwerte(handles);

% Update handles structure
guidata(hObject, handles);

% --- Executes on selection change in popupmenu_Rotor_Eisenmaterial.
function popupmenu_Rotor_Eisenmaterial_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_Rotor_Eisenmaterial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Read data
handles = readOptionen(handles);
handles = readRichtwerte(handles);

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in pushbutton_Plot_BH_Rotor.
function pushbutton_Plot_BH_Rotor_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Plot_BH_Rotor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

plotBH(handles.opt.Rotor_Eisenmaterial);

% --- Executes on button press in pushbutton_Plot_Eisenverlust_Rotor.
function pushbutton_Plot_Eisenverlust_Rotor_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Plot_Eisenverlust_Rotor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

plotVerluste(handles.opt.Rotor_Eisenmaterial);

% --- Executes on selection change in popupmenu_Rotor_Leitermaterial.
function popupmenu_Rotor_Leitermaterial_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_Rotor_Leitermaterial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Read data
handles = readOptionen(handles);
handles = readRichtwerte(handles);

% Update handles structure
guidata(hObject, handles);

function edit_theta_2_Callback(hObject, eventdata, handles)
% hObject    handle to edit_theta_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Read data
handles = readOptionen(handles);
handles = readRichtwerte(handles);

% Update handles structure
guidata(hObject, handles);

% --- Executes on selection change in popupmenu_Rotor_Magnetmaterial.
function popupmenu_Rotor_Magnetmaterial_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_Rotor_Magnetmaterial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Read data
handles = readOptionen(handles);
handles = readRichtwerte(handles);

% Update handles structure
guidata(hObject, handles);

% #########################################################################
% #   C.4) OPTIONEN BERECHNUNG                                            #
% #########################################################################

% --- Executes on selection change in popupmenu_Mode_Wicklung.
function popupmenu_Mode_Wicklung_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_Mode_Wicklung (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Read data
handles = readOptionen(handles);
handles = readRichtwerte(handles);

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in pushbutton_Mode_Wicklung.
function pushbutton_Mode_Wicklung_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Mode_Wicklung (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% to-do

% Update handles structure
guidata(hObject, handles);

% --- Executes on selection change in popupmenu_Wicklungstyp.
function popupmenu_Wicklungstyp_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_Wicklungstyp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Read data
handles = readOptionen(handles);
handles = readRichtwerte(handles);

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in pushbutton_Wicklungstyp.
function pushbutton_Wicklungstyp_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Wicklungstyp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

uiwait(msgbox({'Wicklung A: Ganzlochwicklung, Einschicht, ungesehnt, gezont' ...
    'Wicklung B: Ganzlochwicklung, Zweischicht, gesehnt, gezont' ...
    'Wicklung C: Bruchlochwicklung, Zweischicht, gesehnt, gezont,' ...
    '          Zahnspulenwicklung (q<1) wird dann umgesetzt, wenn der Durchmesser q>1 nicht zulaesst'},'Wicklung','help','modal'))

% #########################################################################
% #   D) ERGEBNISSE                                                       #
% #########################################################################

function edit_Psi_PM_PMSM_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Psi_PM_PMSM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Verketteter Fluss der PM
handles.Entwurf.EMAG.Psi_PM = str2double(handles.edit_Psi_PM_PMSM.String);

handles.opt.Saved = 0;

if(handles.Entwurf.EMAG.Psi_PM~=0 && handles.Entwurf.EMAG.L_d~=0 && handles.Entwurf.EMAG.L_q~=0)
    handles.pushbutton_Analyse.Enable = 'on';
else
    handles.pushbutton_Analyse.Enable = 'off';
end

% Update handles structure
guidata(hObject, handles);

function edit_L_d_PMSM_Callback(hObject, eventdata, handles)
% hObject    handle to edit_L_d_PMSM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Induktivitaet d-Achse
handles.Entwurf.EMAG.L_d = str2double(handles.edit_L_d_PMSM.String)*1e-3;

handles.opt.Saved = 0;

if(handles.Entwurf.EMAG.Psi_PM~=0 && handles.Entwurf.EMAG.L_d~=0 && handles.Entwurf.EMAG.L_q~=0)
    handles.pushbutton_Analyse.Enable = 'on';
else
    handles.pushbutton_Analyse.Enable = 'off';
end

% Update handles structure
guidata(hObject, handles);

function edit_L_q_PMSM_Callback(hObject, eventdata, handles)
% hObject    handle to edit_L_q_PMSM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Induktivitaet q-Achse
handles.Entwurf.EMAG.L_q = str2double(handles.edit_L_q_PMSM.String)*1e-3;

handles.opt.Saved = 0;

if(handles.Entwurf.EMAG.Psi_PM~=0 && handles.Entwurf.EMAG.L_d~=0 && handles.Entwurf.EMAG.L_q~=0)
    handles.pushbutton_Analyse.Enable = 'on';
else
    handles.pushbutton_Analyse.Enable = 'off';
end

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in pushbutton_FEM.
function pushbutton_FEM_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_FEM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if(ispc)
    % Start FEMAG
    % Create Script for FEMAG
    % to-do
else
    opts = struct('WindowStyle','modal','Interpreter','tex');
    errordlg('Platform not supported for FEM Calculations', 'Platform Error', opts)
end

% Update handles structure
guidata(hObject, handles);

% --- Executes on selection change in popupmenu_Ergebnis_Plot.
function popupmenu_Ergebnis_Plot_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_Ergebnis_Plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get popupmenu contents
contents = cellstr(get(hObject,'String'));
handles.Ergebnis_Plot = contents{get(hObject,'Value')};

if(strcmp(handles.Ergebnis_Plot, 'Maschine'))
    cla(handles.axes_Ergebnis_Plot)
    plotMaschine(handles,handles.Entwurf.Geometrie)
    zoom(handles.axes_Ergebnis_Plot,'reset')
    zoom(handles.axes_Ergebnis_Plot,'off')
    pan(handles.axes_Ergebnis_Plot,'off')
elseif(strcmp(handles.Ergebnis_Plot, 'Statornut'))
    cla(handles.axes_Ergebnis_Plot)
    plotNut_1(handles,handles.Entwurf.Geometrie.Nut_1);
    zoom(handles.axes_Ergebnis_Plot,'reset')
    zoom(handles.axes_Ergebnis_Plot,'off')
    pan(handles.axes_Ergebnis_Plot,'off')
elseif(strcmp(handles.Ergebnis_Plot, 'Rotornut'))
    cla(handles.axes_Ergebnis_Plot)
    plotNut_2(handles,handles.Entwurf.Geometrie.Nut_2);
    zoom(handles.axes_Ergebnis_Plot,'reset')
    zoom(handles.axes_Ergebnis_Plot,'off')
    pan(handles.axes_Ergebnis_Plot,'off')
else
    error('Ungueltige Eingabe bei Variable "handles.Ergebnis_Plot"')
end

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in pushbutton_Magnetanpassung.
function pushbutton_Magnetanpassung_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Magnetanpassung (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.uipanel_Magnetanpassung.Visible = 'on';

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in pushbutton_Magnetanpassung_Plot.
function pushbutton_Magnetanpassung_Plot_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Magnetanpassung_Plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.opt.Saved = 0;

if(strcmp(handles.opt.Maschinenausfuehrung,'IPMSM (eingelassen)'))
    handles.Entwurf.Geometrie.b_PM = str2double(handles.edit_b_PM.String);
    handles.Entwurf.Geometrie.h_PM = str2double(handles.edit_h_PM.String);
elseif(strcmp(handles.opt.Maschinenausfuehrung,'IPMSM (tangential)'))
    handles.Entwurf.Geometrie.b_PM = str2double(handles.edit_b_PM.String);
    handles.Entwurf.Geometrie.h_PM = str2double(handles.edit_h_PM.String);
    handles.Entwurf.Geometrie.Abstand_PM_Rotoroberflaeche = str2double(handles.edit_Abstand_PM_Rotoroberflaeche.String);
elseif(strcmp(handles.opt.Maschinenausfuehrung,'IPMSM (V-Form)'))
    handles.Entwurf.Geometrie.b_PM = str2double(handles.edit_b_PM.String);
    handles.Entwurf.Geometrie.h_PM = str2double(handles.edit_h_PM.String);
    handles.Entwurf.Geometrie.Abstand_PM_Rotoroberflaeche = str2double(handles.edit_Abstand_PM_Rotoroberflaeche.String);
    handles.Entwurf.Geometrie.alpha_PM = str2double(handles.edit_alpha_PM.String);
    handles.Entwurf.Geometrie.Abstand_PM_unten = str2double(handles.edit_Abstand_PM_unten.String);
else
    error('Ungueltige Eingabe bei Variable "opt.Maschinenausfuehrung"')
end

handles.Entwurf.EMAG.Psi_PM = 0;
handles.Entwurf.EMAG.L_d = 0;
handles.Entwurf.EMAG.L_q = 0;

handles = setErgebnisse(handles);

handles.Entwurf.Geometrie = Update_x_y_Koordinaten(handles.Entwurf.Bemessungswerte, handles.Entwurf.Geometrie, handles.Entwurf.Wicklung, handles.Entwurf.Optionen);

object = handles.popupmenu_Ergebnis_Plot;
GUI_Entwurf('popupmenu_Ergebnis_Plot_Callback', object, eventdata, handles)

handles.uipanel_Magnetanpassung.Visible = 'off';

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in pushbutton_Analyse.
function pushbutton_Analyse_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Analyse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if(handles.opt.Saved)
else
    if(isfolder(['3_Ergebnisse/' handles.opt.folder_id]))
        button = questdlg('Speichern als ... ?', 'Datei speichern','Neue Datei','Ueberschreiben','Neue Datei'); 
        switch button 
            case 'Neue Datei'
                handles.opt.unique_id = datestr(now,'yyyymmdd_HHMMSS');
                handles.opt.folder_id = [handles.opt.Maschinentyp '_data_' handles.opt.unique_id];
                handles.opt.file_id = [handles.opt.Maschinentyp,'_',handles.opt.unique_id];
                handles.text_ID.String = ['ID: Entwurf_' handles.opt.file_id];
                mkdir('3_Ergebnisse',handles.opt.folder_id);
                mkdir(['3_Ergebnisse/',handles.opt.folder_id],'1_Entwurf');
                mkdir(['3_Ergebnisse/',handles.opt.folder_id],'2_Analyse');
            case 'Ueberschreiben'  
        end
    else
        mkdir('3_Ergebnisse',handles.opt.folder_id);
        mkdir(['3_Ergebnisse/',handles.opt.folder_id],'1_Entwurf');
        mkdir(['3_Ergebnisse/',handles.opt.folder_id],'2_Analyse');
    end

    if(handles.opt.Locked)
        handles.opt.Saved = 1;
        handles.Entwurf.Optionen.Saved = handles.opt.Saved;
        handles.Entwurf.Optionen.unique_id = handles.opt.unique_id;
        handles.Entwurf.Optionen.folder_id = handles.opt.folder_id;
        handles.Entwurf.Optionen.file_id = handles.opt.file_id;
        Entwurf = handles.Entwurf;
    else
        Entwurf.Bemessungswerte = handles.rated;
        Entwurf.Richtwerte = handles.richt;
        Entwurf.Optionen = handles.opt;
    end

    try
        % Save struct Entwurf
        save(['3_Ergebnisse/',handles.opt.folder_id,'/1_Entwurf/Entwurf_',handles.opt.file_id,'.mat'],'Entwurf');
        % Save DXF file
        if(handles.opt.Locked)
            saveDXF(handles,handles.Entwurf.Geometrie)
        end
        % Save Excel file
        saveExcel(handles, Entwurf)

        uiwait(msgbox({'Speichern erfolgreich.'},'Datei speichern','help','modal'))
    catch
        uiwait(msgbox({'Speichern fehlgeschlagen.'},'Datei speichern','help','error'))
    end
end

GUI_Analyse(handles.Entwurf);

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in pushbutton_Content2.
function pushbutton_Content2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Content2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Content
handles.popupmenu_Content.Value = 1;
handles.opt.Content = handles.popupmenu_Content.String{handles.popupmenu_Content.Value};

handles = toggleContent(handles);

% Update handles structure
guidata(hObject, handles);

% #########################################################################
% #   E) AUSWAHL WICKLUNG                                                 #
% #########################################################################

% --- Executes on selection change in popupmenu_Auswahl_Wicklung.
function popupmenu_Auswahl_Wicklung_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_Auswahl_Wicklung (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if(strcmp(handles.opt.Mode_Wicklung,'Manuell'))
    load(['3_Ergebnisse/1_Misc/Entwurf_',handles.opt.file_id,'_temp.mat']);
    
    set(handles.edit_Wicklungslayout, 'Enable', 'inactive');
    
    if(isfield(wick,'y_1v'))
        wick = rmfield(wick,'y_1v');
    end
    if(isfield(wick,'y_1'))
        wick = rmfield(wick,'y_1');
    end
    if(isfield(wick,'Sehnung_1'))
        wick = rmfield(wick,'Sehnung_1');
    end
    
    var1 = wick.table_1pos(handles.popupmenu_Auswahl_Wicklung.Value-1,:);
    var1 = table2struct(var1);
    if(isfield(wick,'Matrix_lay'))
        wick = rmfield(wick,'Matrix_lay');
    end
    var1 = [fieldnames(wick)' fieldnames(var1)'; struct2cell(wick)' struct2cell(var1)'];
    var1 = struct(var1{:});
    
    handles.edit_y_1v.String = num2str(var1.y_1v);
    handles.edit_y_1v.Enable = 'off';

    var1 = TingleyAlg(rated, var1);
    wick.Matrix = var1.Matrix;
    wick.Matrix_lay = var1.Matrix_lay;
    
    % Update edit_Wicklungslayout
    var = var1.N_1/rated.m;
    name_3 = [{'A'},{'C'},{'B'}; {'a'},{'c'},{'b'}];
    for j = 1:rated.m
        vec = var1.Matrix(j,:);
        for i = 1:length(vec)
            if(sign(vec(i)<1))
                M1(j,i) = name_3(2,j);
            else
                M1(j,i) = name_3(1,j);
            end
        end
    end
    for k = 1:var1.n_lay
        M2 = var1.Matrix(:,var*(k-1)+1:var*k);
        M2 = M2(:);
        M3 = M1(:,var*(k-1)+1:var*k);
        M3 = M3(:);
        [~, indices] = sort(M2,1,'ComparisonMethod','abs');
        M4(k,:) = M3(indices)';
    end
    M4 = M4(:)';
    
    n    =  var1.n_lay;   % Step width
    x    = {'|'};   % To insert
    len    = numel(M4);
    M4 = [M4; repmat({'x'}, 1, (var1.N_1 * var1.n_lay))];
    M4(2, n:n:len) = x;
    M4 = M4(~strcmp(M4,'x')).';
    M4(end) = [];
    
    string_Wicklungslayout = char(M4)';
    handles.edit_Wicklungslayout.String = string_Wicklungslayout;
    
    cla(handles.axes_Auswahl_Wicklung,'reset')
    plotWicklungslayout(handles,var1,rated)
    
    save(['3_Ergebnisse/1_Misc/Entwurf_',handles.opt.file_id,'_temp.mat'],'wick','rated');

end
% Update handles structure
guidata(hObject, handles);

% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over edit_Wicklungslayout.
function edit_Wicklungslayout_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to edit_Wicklungslayout (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(hObject, 'Enable', 'On');

handles.edit_y_1v.String = '-';
handles.edit_y_1v.Enable = 'on';
handles.popupmenu_Auswahl_Wicklung.Value = 1;

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in pushbutton_Wicklungslayout.
function pushbutton_Wicklungslayout_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Wicklungslayout (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if(~strcmp(handles.edit_Wicklungslayout.Enable, 'inactive'))

    load(['3_Ergebnisse/1_Misc/Entwurf_',handles.opt.file_id,'_temp.mat']);

    % Update wick.Matrix
    string_Wicklungslayout = handles.edit_Wicklungslayout.String;
    M1 = strsplit(string_Wicklungslayout,'|');
    M2 = char(M1);
    M2 = M2(:)';

    N_1 = length(M1);
    n_lay = cellfun('length',M1);
    n_lay = n_lay(1);

    % Consistency Check
    lett = unique(M2);
    for i = 1:length(lett)
        count_lett(i) = count(M2,lett(i));
    end
    if(~any(N_1==wick.table_1pos.N_1 & n_lay==wick.table_1pos.n_lay) || any(cellfun('length',M1)~=n_lay) || length(unique(count_lett))>1)
        uiwait(msgbox({'Nicht symmetrische Wicklung eingegeben'},'Wicklungslayout','help','error'))
        return
    end

    name_3 = [{'A'},{'C'},{'B'};{'a'},{'c'},{'b'}];
    M4 = num2cell(zeros(1,rated.m));
    for k = 1:n_lay
        M3 = M2(:,N_1*(k-1)+1:N_1*k);
        M3 = cellstr(M3')';
        for j = 1:rated.m
            M4(j) = {[M4{j}, find(contains(M3,name_3(1,j)))]};
            M4(j) = {[M4{j}, find(contains(M3,name_3(2,j))) * (-1)]};
            M5(k,j) = {find(contains(M3,name_3(1,j)))};
            M5(k,j) = {[M5{k,j}, find(contains(M3,name_3(2,j))) * (-1)]};
        end
    end
    M4 = M4';
    M5 = M5';
    M6 = cell2mat(M4);
    M6(:,1) = [];

    var1 = wick.table_1pos(wick.table_1pos.N_1==N_1 & wick.table_1pos.n_lay==n_lay,:);
    var1 = var1(1,1:9);
    wick.table_1sel = var1;
    var1 = table2struct(var1);
    if(isfield(wick,'Matrix_lay'))
        wick = rmfield(wick,'Matrix_lay');
    end
    var1 = [fieldnames(wick)' fieldnames(var1)'; struct2cell(wick)' struct2cell(var1)'];
    var1 = struct(var1{:});

    var1.Matrix = M6;
    var1.Matrix_lay = M5;
    wick.Matrix = var1.Matrix;
    wick.Matrix_lay = var1.Matrix_lay;

    if(~isnan(str2double(handles.edit_y_1v.String)))
        if(~mod(var1.y_1Durchmesser - str2double(handles.edit_y_1v.String),1))
            wick.y_1v = str2double(handles.edit_y_1v.String);
            wick.y_1 = var1.y_1Durchmesser - wick.y_1v;
            wick.Sehnung_1 = wick.y_1 ./ var1.y_1Durchmesser;
        else
            uiwait(msgbox({'Bitte gueltige Wicklungsschrittverkuerzung bzw -verlaengerung angeben.'},'Auswahl Wicklung','help','error'))
            return
        end
    else
        uiwait(msgbox({'Bitte gueltige Wicklungsschrittverkuerzung bzw -verlaengerung angeben.'},'Auswahl Wicklung','help','error'))
        return
    end

    cla(handles.axes_Auswahl_Wicklung,'reset')
    plotWicklungslayout(handles,var1,rated)

    save(['3_Ergebnisse/1_Misc/Entwurf_',handles.opt.file_id,'_temp.mat'],'wick','rated');

end

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in pushbutton_Auswahl_Wicklung.
function pushbutton_Auswahl_Wicklung_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Auswahl_Wicklung (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if(strcmp(handles.opt.Mode_Wicklung,'Manuell'))
    if(strcmp(handles.edit_Wicklungslayout.String,'-'))
        uiwait(msgbox({'Bitte Wicklung auswaehlen.'},'Auswahl Wicklung','help','error'))
    else
        cla(handles.axes_Auswahl_Wicklung,'reset')
        handles.edit_Wicklungslayout.String = '-';
        uiresume(gcbf);
    end
else
    if(strcmp(handles.popupmenu_Auswahl_Wicklung.String{handles.popupmenu_Auswahl_Wicklung.Value},'-'))
        uiwait(msgbox({'Bitte Wicklung auswaehlen.'},'Auswahl Wicklung','help','error'))
    else
        uiresume(gcbf);
    end
end

% Update handles structure
        guidata(hObject, handles);

% #########################################################################
% #   F) HILFSFUNKTIONEN                                                  #
% #########################################################################

% Background Color Richtwerte
function var = bgCl(var,limit)
    var2 = str2double(var.String);
    if(var2 > limit(1) || var2 < limit(2))
        set(var,'BackgroundColor',[1 0.8 0.8]);
    else
        set(var,'BackgroundColor',[0.8 1 0.8]);
    end

% Toggle Locked
function handles = toggleLocked(handles)
    handles.uipanel_Magnetanpassung.Visible = 'off';
    handles.pushbutton_Content2.Enable = 'on';
    handles.popupmenu_Content.Enable = 'on';
    if(handles.opt.Locked)
        set(findall(handles.uipanel_Eingabe,'-property','Enable'),'Enable','off')
        handles.pushbutton_Content1.Enable = 'on';
        handles.uipanel_Hauptabmessungen.Visible = 'on';
        handles.popupmenu_Ergebnis_Plot.Visible = 'on';
        handles.uipanel_Ergebnis_Plot.Visible = 'on';
        if(strcmp(handles.opt.Maschinentyp,'ASM'))
            handles.uipanel_Kenndaten_PMSM.Visible = 'off';
            handles.uipanel_Kenndaten_ASM.Visible = 'on';
        elseif(strcmp(handles.opt.Maschinentyp,'PMSM'))
            handles.uipanel_Kenndaten_PMSM.Visible = 'on';
            handles.uipanel_Kenndaten_ASM.Visible = 'off';
        else
            error('Ungueltige Eingabe bei Variable "handles.opt.Maschinentyp"');
        end
        if(strcmp(handles.opt.Maschinenausfuehrung,'IPMSM (eingelassen)') || strcmp(handles.opt.Maschinenausfuehrung,'IPMSM (tangential)') || strcmp(handles.opt.Maschinenausfuehrung,'IPMSM (V-Form)'))
            if(handles.Entwurf.EMAG.Psi_PM~=0 && handles.Entwurf.EMAG.L_d~=0 && handles.Entwurf.EMAG.L_q~=0)
                handles.pushbutton_Analyse.Enable = 'on';
            else
                handles.pushbutton_Analyse.Enable = 'off';
            end
        else
            handles.pushbutton_Analyse.Enable = 'on';
        end
        
        set(handles.popupmenu_Maschinentyp,'Enable','off');

        handles.text_Wait.Visible = 'off';
        handles.text_Wait.String = [{'Success ...'}; {'Die Berechnungen waren erfolgreich.'}];
        
        handles = setErgebnisse(handles);
        cla(handles.axes_Ergebnis_Plot)
        plotMaschine(handles,handles.Entwurf.Geometrie)
        zoom(handles.axes_Ergebnis_Plot,'reset')
        zoom(handles.axes_Ergebnis_Plot,'off')
        pan(handles.axes_Ergebnis_Plot,'off')
    else
        set(findall(handles.uipanel_Eingabe,'-property','Enable'),'Enable','on')
        handles.uipanel_Hauptabmessungen.Visible = 'off';
        handles.popupmenu_Ergebnis_Plot.Visible = 'off';
        handles.uipanel_Ergebnis_Plot.Visible = 'off';
        handles.uipanel_Kenndaten_PMSM.Visible = 'off';
        handles.uipanel_Kenndaten_ASM.Visible = 'off';
        handles.pushbutton_Analyse.Enable = 'off';
        
        set(handles.popupmenu_Maschinentyp,'Enable','on');
        
        handles.text_Wait.Visible = 'on';
        handles.text_Wait.String = 'Noch keine Berechnungen ausgefuehrt.';
        
        handles = toggleMaschinentyp(handles);
        
        if(strcmp(handles.opt.Mode_Wicklung,'Klassisch'))
            set(handles.popupmenu_Wicklungstyp,'Enable','on')
            set(handles.pushbutton_Wicklungstyp,'Enable','on')
        else
            set(handles.popupmenu_Wicklungstyp,'Enable','off')
            set(handles.pushbutton_Wicklungstyp,'Enable','off')
        end
    end

% Toggle Maschinentyp
function handles = toggleMaschinentyp(handles)
    if(strcmp(handles.opt.Maschinentyp,'ASM'))
        set(handles.edit_cos_phi_N,'Enable','off')

        set(handles.popupmenu_B_A,'Value',1)
        set(handles.popupmenu_B_A,'String',[{'Mittelwert der Luftspaltinduktion in T'}; {'Strombelag in A/mm'}])
        
        handles.edit_S_2s.Enable = 'on';
        handles.edit_S_2r.Enable = 'on';
        handles.edit_S_2.Enable = 'on';
        handles.edit_B_2r_max.Enable = 'on';
        handles.edit_B_2z_max.Enable = 'on';
        handles.edit_phi_2n.Enable = 'on';
        handles.edit_xi_2p.Enable = 'on';
        handles.edit_tau_2n_min.Enable = 'on';
        handles.edit_phi_2Fe.Enable = 'on';
        
        set(handles.pushbutton_Maschinenausfuehrung,'Enable','off')
        set(handles.popupmenu_Spulenform_Rotor,'Enable','on')
        set(handles.popupmenu_Nutform_Rotor,'Enable','on')

        set(handles.uipanel_Optionen_Maschine_Rotor_Leitermaterial,'Visible','on')
        set(handles.uipanel_Optionen_Maschine_Rotor_Magnetmaterial,'Visible','off')
        
        handles.popupmenu_Ergebnis_Plot.Value = 1;
        handles.popupmenu_Ergebnis_Plot.String = [{'Maschine'}; {'Statornut'}; {'Rotornut'}];
    elseif(strcmp(handles.opt.Maschinentyp,'PMSM'))
        set(handles.edit_cos_phi_N,'Enable','on')

        set(handles.popupmenu_B_A,'Value',1)
        set(handles.popupmenu_B_A,'String',[{'Amplitude der Luftspaltinduktion in T'}; {'Strombelag in A/mm'}])

        handles.edit_S_2s.Enable = 'off';
        handles.edit_S_2r.Enable = 'off';
        handles.edit_S_2.Enable = 'off';
        handles.edit_B_2r_max.Enable = 'on';
        handles.edit_B_2z_max.Enable = 'off';
        handles.edit_phi_2n.Enable = 'off';
        handles.edit_xi_2p.Enable = 'off';
        handles.edit_tau_2n_min.Enable = 'off';
        handles.edit_phi_2Fe.Enable = 'on';

        set(handles.pushbutton_Maschinenausfuehrung,'Enable','on')
        set(handles.popupmenu_Spulenform_Rotor,'Enable','off')
        set(handles.popupmenu_Nutform_Rotor,'Enable','off')

        set(handles.uipanel_Optionen_Maschine_Rotor_Leitermaterial,'Visible','off')
        set(handles.uipanel_Optionen_Maschine_Rotor_Magnetmaterial,'Visible','on')

        handles.popupmenu_Ergebnis_Plot.Value = 1;
        handles.popupmenu_Ergebnis_Plot.String = [{'Maschine'}; {'Statornut'}];
    else
        error('Ungueltige Eingabe bei Variable "handles.opt.Maschinentyp"');
    end
    
    handles = toggleMaschinenausfuehrung(handles);

% Toggle Maschinenausfuehrung
function handles = toggleMaschinenausfuehrung(handles)
    handles.pushbutton_Magnetanpassung.Visible = 'off';
    if(strcmp(handles.opt.Maschinenausfuehrung,'Kaefiglaeufer'))
        handles.text_S_2.String = 'Stab- bzw. Ringstromdichte (Rotor) in A/mm^2';
        handles.edit_S_2s.Visible = 'on';
        handles.edit_S_2r.Visible = 'on';
        handles.edit_S_2.Visible = 'off';
        handles.edit_S_2s.Enable = 'on';
        handles.edit_S_2r.Enable = 'on';
        handles.edit_S_2.Enable = 'off';
        handles.edit_xi_2p.Enable = 'off';
        handles.edit_phi_2n.Enable = 'off';
    elseif(strcmp(handles.opt.Maschinenausfuehrung,'Schleifringlaeufer'))
        handles.text_S_2.String = 'Stromdichte (Rotor) in A/mm^2';
        handles.edit_S_2s.Visible = 'off';
        handles.edit_S_2r.Visible = 'off';
        handles.edit_S_2.Visible = 'on';
        handles.edit_S_2s.Enable = 'off';
        handles.edit_S_2r.Enable = 'off';
        handles.edit_S_2.Enable = 'on';
        handles.edit_xi_2p.Enable = 'on';
        handles.edit_phi_2n.Enable = 'on';
    elseif(strcmp(handles.opt.Maschinenausfuehrung,'IPMSM (eingelassen)') || strcmp(handles.opt.Maschinenausfuehrung,'IPMSM (tangential)') || strcmp(handles.opt.Maschinenausfuehrung,'IPMSM (V-Form)'))
        handles.pushbutton_Magnetanpassung.Visible = 'on';
    end

% Toggle Content
function handles = toggleContent(handles)
    if(strcmp(handles.opt.Content,'Eingabe'))
        set(handles.uipanel_Eingabe,'Visible','on')
        set(handles.uipanel_Ergebnisse,'Visible','off')
        set(handles.uipanel_Auswahl_Wicklung,'Visible','off')
        set(handles.uipanel_Ergebnisse,'Position',[742,0,742,826])
    elseif(strcmp(handles.opt.Content,'Ergebnisse'))
        handles = toggleLocked(handles);
        set(handles.uipanel_Ergebnisse,'Position',[0,0,742,826])
        set(handles.uipanel_Eingabe,'Visible','off')
        set(handles.uipanel_Ergebnisse,'Visible','on')
        set(handles.uipanel_Auswahl_Wicklung,'Visible','off')
    else
        error('Ungueltige Eingabe bei Variable "handles.opt.Content"');
    end

% New Load
function handles = newLoad(handles)
    handles.opt.Locked = 0;
    handles.opt.Saved = 0;
    
    % Maschinentyp
    handles.opt.Maschinentyp = handles.popupmenu_Maschinentyp.String{handles.popupmenu_Maschinentyp.Value};

    % Load default values
    handles = loadDefault(handles);

    % Set default data
    handles = setDefaultOptionen(handles);

    % Read data
    handles = readOptionen(handles);
    handles = readBemessungswerte(handles);

    % Set GUI
    handles = toggleMaschinentyp(handles);
    handles = toggleContent(handles);

    % Set default data
    handles = setDefaultRichtwerte(handles);

    % Read data
    handles = readRichtwerte(handles);

    handles = toggleLocked(handles);
    handles.opt.unique_id = datestr(now,'yyyymmdd_HHMMSS');
    handles.opt.folder_id = [handles.opt.Maschinentyp '_data_' handles.opt.unique_id];
    handles.opt.file_id = [handles.opt.Maschinentyp,'_',handles.opt.unique_id];
    handles.text_ID.String = ['ID: Entwurf_' handles.opt.file_id];
    
% Load default
function handles = loadDefault(handles)
    if(strcmp(handles.opt.Maschinentyp,'ASM'))
        Richtwerte_ASM;
        Optionen_ASM;
        Beschraenkungen_ASM;
    elseif(strcmp(handles.opt.Maschinentyp,'PMSM'))
        Richtwerte_PMSM;
        Optionen_PMSM;
        Beschraenkungen_PMSM;
    else
        error('Ungueltige Eingabe bei Variable "handles.Maschinentyp"');
    end
    if(isfield(handles,'default_richt'))
        handles = rmfield(handles,'default_richt');
    end
    if(isfield(handles,'default_opt'))
        handles = rmfield(handles,'default_opt');
    end
    if(isfield(handles,'limit'))
        handles = rmfield(handles,'limit');
    end
    handles.default_richt = richt;
    handles.default_opt = opt;
    handles.limit = limit;

% Load material
function data = loadMaterial(var,typ)
    
    if(strcmp(var.String,'-'))
        data.String = var.String;
        return
    end
    switch typ
        case 'Elektroblech'
            filepath = '5_Materialien/1_Elektroblech/';
        case 'Leiter'
            filepath = '5_Materialien/2_Leiter/';
        case 'Magnet'
            filepath = '5_Materialien/3_Magnet/';
        otherwise
            error('Undefined material type')
    end
    
    load([filepath var.String '.mat']);
    data.String = var.String;

% Read Optionen
function handles = readOptionen(handles)
    % Maschinentyp
    handles.opt.Maschinentyp = handles.popupmenu_Maschinentyp.String{handles.popupmenu_Maschinentyp.Value};

    % Content
    handles.opt.Content = handles.popupmenu_Content.String{handles.popupmenu_Content.Value};
    
    % Maschinenausfuehrung
    handles.opt.Maschinenausfuehrung = handles.popupmenu_Maschinenausfuehrung.String{handles.popupmenu_Maschinenausfuehrung.Value};
    
    % Schaltung
    handles.opt.Schaltung = handles.popupmenu_Schaltung.String{handles.popupmenu_Schaltung.Value};
    
    % Spulenform Stator
    handles.opt.Spulenform_Stator = handles.popupmenu_Spulenform_Stator.String{handles.popupmenu_Spulenform_Stator.Value};
    
    % Spulenform Rotor
    handles.opt.Spulenform_Rotor = handles.popupmenu_Spulenform_Rotor.String{handles.popupmenu_Spulenform_Rotor.Value};
    
    % Nutform Stator
    handles.opt.Nutform_Stator = handles.popupmenu_Nutform_Stator.String{handles.popupmenu_Nutform_Stator.Value};
    
    % Nutform Rotor
    handles.opt.Nutform_Rotor = handles.popupmenu_Nutform_Rotor.String{handles.popupmenu_Nutform_Rotor.Value};
    
    % Kuehlungsart
    handles.opt.Kuehlungsart = handles.popupmenu_Kuehlungsart.String{handles.popupmenu_Kuehlungsart.Value};
    
    % Eisenmaterial Stator
    handles.opt.Stator_Eisenmaterial.String = handles.popupmenu_Stator_Eisenmaterial.String{handles.popupmenu_Stator_Eisenmaterial.Value};
    handles.opt.Stator_Eisenmaterial = loadMaterial(handles.opt.Stator_Eisenmaterial,'Elektroblech');
    
    % Leitermaterial Stator
    handles.opt.Stator_Leitermaterial.String = handles.popupmenu_Stator_Leitermaterial.String{handles.popupmenu_Stator_Leitermaterial.Value};
    handles.opt.Stator_Leitermaterial = loadMaterial(handles.opt.Stator_Leitermaterial,'Leiter');
    
    % Temperatur Leiter Stator
    handles.opt.theta_1 = str2double(handles.edit_theta_1.String);
    
    % Eisenmaterial Rotor
    handles.opt.Rotor_Eisenmaterial.String = handles.popupmenu_Rotor_Eisenmaterial.String{handles.popupmenu_Rotor_Eisenmaterial.Value};
    handles.opt.Rotor_Eisenmaterial = loadMaterial(handles.opt.Rotor_Eisenmaterial,'Elektroblech');
    
    % Leitermaterial Rotor
    handles.opt.Rotor_Leitermaterial.String = handles.popupmenu_Rotor_Leitermaterial.String{handles.popupmenu_Rotor_Leitermaterial.Value};
    handles.opt.Rotor_Leitermaterial = loadMaterial(handles.opt.Rotor_Leitermaterial,'Leiter');
    
    % Temperatur Leiter Rotor
    handles.opt.theta_2 = str2double(handles.edit_theta_2.String);
    
    % Magnetmaterial Rotor
    handles.opt.Rotor_Magnetmaterial.String = handles.popupmenu_Rotor_Magnetmaterial.String{handles.popupmenu_Rotor_Magnetmaterial.Value};
    handles.opt.Rotor_Magnetmaterial = loadMaterial(handles.opt.Rotor_Magnetmaterial,'Magnet');
    
    % Remanenzinduktion PM
    if(strcmp(handles.opt.Maschinentyp,'PMSM'))
        handles.edit_B_r_PM.String = num2str(handles.opt.Rotor_Magnetmaterial.B_r);
    end
    
    % Modus Wicklungsauslegung
    handles.opt.Mode_Wicklung = handles.popupmenu_Mode_Wicklung.String{handles.popupmenu_Mode_Wicklung.Value};
    
    % Wicklungstyp
    if(strcmp(handles.opt.Mode_Wicklung,'Klassisch'))
        set(handles.popupmenu_Wicklungstyp,'Enable','on')
        set(handles.pushbutton_Wicklungstyp,'Enable','on')
        handles.popupmenu_Wicklungstyp.String = handles.default_opt.Wicklungstyp;
        handles.opt.Wicklungstyp = handles.popupmenu_Wicklungstyp.String{handles.popupmenu_Wicklungstyp.Value};
    else
        set(handles.popupmenu_Wicklungstyp,'Enable','off')
        set(handles.pushbutton_Wicklungstyp,'Enable','off')
        handles.popupmenu_Wicklungstyp.Value = 1;
        handles.popupmenu_Wicklungstyp.String = {'-'};
        handles.opt.Wicklungstyp = handles.popupmenu_Wicklungstyp.String{handles.popupmenu_Wicklungstyp.Value};
    end
    
% Read Bemessungswerte
function handles = readBemessungswerte(handles)
    % Nenndrehzahl
    handles.rated.n_N = str2double(handles.edit_n_N.String);

    % Nennspannung
    handles.rated.U_N = str2double(handles.edit_U_N.String);

    % Polpaarzahl
    handles.rated.p = str2double(handles.edit_p.String);
    handles.rated.f_N = (handles.rated.p * handles.rated.n_N) / 60;
    handles.edit_f_N.String = num2str(handles.rated.f_N);

    % Nennleistung
    handles.rated.P_N = str2double(handles.edit_P_N.String)*1e3;
    
    % Nennfrequenz
    handles.rated.f_N = str2double(handles.edit_f_N.String);
    
    % Leistungsfaktor
    handles.rated.cos_phi_N = str2double(handles.edit_cos_phi_N.String);
    % Strangzahl
    handles.rated.m = str2double(handles.edit_m.String);

% Read Richtwerte
function handles = readRichtwerte(handles)
    if(strcmp(handles.opt.Maschinentyp,'ASM'))
        % Relative Ankerlaenge
        if(handles.rated.p > 1)
            bgCl(handles.edit_lambda,handles.default_richt.lambda(2,:));
        else
            bgCl(handles.edit_lambda,handles.default_richt.lambda(1,:));
        end
        handles.richt.lambda = str2double(handles.edit_lambda.String);

        % Kanalbreite der Ventilationskanaele
        bgCl(handles.edit_l_v,handles.default_richt.l_v(1,:));
        handles.richt.l_v = str2double(handles.edit_l_v.String);

        % Luftspaltinduktion bzw. Strombelag
        handles.B_A = handles.popupmenu_B_A.String{handles.popupmenu_B_A.Value};
        if(strcmp(handles.B_A,'Mittelwert der Luftspaltinduktion in T'))
            bgCl(handles.edit_B_A,handles.default_richt.B_m(1,:));
            handles.richt.B_m = str2double(handles.edit_B_A.String);
        elseif(strcmp(handles.B_A,'Strombelag in A/mm'))
            bgCl(handles.edit_B_A,handles.default_richt.A(1,:));
            handles.richt.A = str2double(handles.edit_B_A.String);
        else
            error('Ungueltige Eingabe bei Variable "handles.B_A"');
        end

        % Stromdichte (Stator)
        bgCl(handles.edit_S_1,handles.default_richt.S_1(1,:));
        handles.richt.S_1 = str2double(handles.edit_S_1.String);

        % max. Rueckeninduktion (Stator)
        bgCl(handles.edit_B_1r_max,handles.default_richt.B_1r_max(1,:));
        handles.richt.B_1r_max = str2double(handles.edit_B_1r_max.String);
        
        % max. Zahninduktion (Stator)
        bgCl(handles.edit_B_1z_max,handles.default_richt.B_1z_max(1,:));
        handles.richt.B_1z_max = str2double(handles.edit_B_1z_max.String);
        
        % Nutfuellfaktor (Stator)
        var = strcmp(handles.opt.Spulenform_Stator,handles.default_richt.phi_1n_opt);
        bgCl(handles.edit_phi_1n,handles.default_richt.phi_1n(var~=0,:));
        handles.richt.phi_1n = str2double(handles.edit_phi_1n.String);
        
        % Wicklungsfaktor Grundwelle (Stator)
        bgCl(handles.edit_xi_1p,handles.default_richt.xi_1p(1,:));
        handles.richt.xi_1p = str2double(handles.edit_xi_1p.String);
        
        % min. Nutteilung (Stator)
        bgCl(handles.edit_tau_1n_min,handles.default_richt.tau_1n_min(1,:));
        handles.richt.tau_1n_min = str2double(handles.edit_tau_1n_min.String);
        
        % Eisenfuellfaktor (Stator)
        bgCl(handles.edit_phi_1Fe,handles.default_richt.phi_1Fe(1,:));
        handles.richt.phi_1Fe = str2double(handles.edit_phi_1Fe.String);
        
        % Stabstromdichte (Rotor)
        var = strcmp(handles.opt.Rotor_Leitermaterial.String,handles.default_richt.S_2s_opt);
        bgCl(handles.edit_S_2s,handles.default_richt.S_2s(var~=0,:));
        handles.richt.S_2s = str2double(handles.edit_S_2s.String);
        % Ringstromdichte (Rotor)
        var = strcmp(handles.opt.Rotor_Leitermaterial.String,handles.default_richt.S_2r_opt);
        bgCl(handles.edit_S_2r,handles.default_richt.S_2r(var~=0,:));
        handles.richt.S_2r = str2double(handles.edit_S_2r.String);

        % Stromdichte (Rotor)
        var = strcmp(handles.opt.Rotor_Leitermaterial.String,handles.default_richt.S_2_opt);
        bgCl(handles.edit_S_2,handles.default_richt.S_2(var~=0,:));
        handles.richt.S_2 = str2double(handles.edit_S_2.String);

        % max. Rueckeninduktion (Rotor)
        bgCl(handles.edit_B_2r_max,handles.default_richt.B_2r_max(1,:));
        handles.richt.B_2r_max = str2double(handles.edit_B_2r_max.String);
        
        % max. Zahninduktion (Rotor)
        bgCl(handles.edit_B_2z_max,handles.default_richt.B_2z_max(1,:));
        handles.richt.B_2z_max = str2double(handles.edit_B_2z_max.String);
        
        % Nutfuellfaktor (Rotor)
        var = strcmp(handles.opt.Spulenform_Rotor,handles.default_richt.phi_2n_opt);
        bgCl(handles.edit_phi_2n,handles.default_richt.phi_2n(var~=0,:));
        handles.richt.phi_2n = str2double(handles.edit_phi_2n.String);
        
        % Wicklungsfaktor Grundwelle (Rotor)
        bgCl(handles.edit_xi_2p,handles.default_richt.xi_2p(1,:));
        handles.richt.xi_2p = str2double(handles.edit_xi_2p.String);
        
        % min. Nutteilung (Rotor)
        bgCl(handles.edit_tau_2n_min,handles.default_richt.tau_2n_min(1,:));
        handles.richt.tau_2n_min = str2double(handles.edit_tau_2n_min.String);
        
        % Eisenfuellfaktor (Rotor)
        bgCl(handles.edit_phi_2Fe,handles.default_richt.phi_2Fe(1,:));
        handles.richt.phi_2Fe = str2double(handles.edit_phi_2Fe.String);
        
    elseif(strcmp(handles.opt.Maschinentyp,'PMSM'))
        % Relative Ankerlaenge
        if(handles.rated.p > 1)
            bgCl(handles.edit_lambda,handles.default_richt.lambda(2,:));
        else
            bgCl(handles.edit_lambda,handles.default_richt.lambda(1,:));
        end
        handles.richt.lambda = str2double(handles.edit_lambda.String);

        % Kanalbreite der Ventilationskanaele
        bgCl(handles.edit_l_v,handles.default_richt.l_v(1,:));
        handles.richt.l_v = str2double(handles.edit_l_v.String);

        % Luftspaltinduktion bzw. Strombelag
        handles.B_A = handles.popupmenu_B_A.String{handles.popupmenu_B_A.Value};
        if(strcmp(handles.B_A,'Amplitude der Luftspaltinduktion in T'))
            bgCl(handles.edit_B_A,handles.default_richt.B_p(1,:));
            handles.richt.B_p = str2double(handles.edit_B_A.String);
        elseif(strcmp(handles.B_A,'Strombelag in A/mm'))
            var = strcmp(handles.opt.Kuehlungsart,handles.default_richt.A_opt);
            bgCl(handles.edit_B_A,handles.default_richt.A(var~=0,:));
            handles.richt.A = str2double(handles.edit_B_A.String);
        else
            error('Ungueltige Eingabe bei Variable "handles.B_A"');
        end

        % Stromdichte (Stator)
        var = strcmp(handles.opt.Kuehlungsart,handles.default_richt.S_1_opt);
        bgCl(handles.edit_S_1,handles.default_richt.S_1(var~=0,:));
        handles.richt.S_1 = str2double(handles.edit_S_1.String);

        % max. Rueckeninduktion (Stator)
        bgCl(handles.edit_B_1r_max,handles.default_richt.B_1r_max(1,:));
        handles.richt.B_1r_max = str2double(handles.edit_B_1r_max.String);
        
        % max. Zahninduktion (Stator)
        bgCl(handles.edit_B_1z_max,handles.default_richt.B_1z_max(1,:));
        handles.richt.B_1z_max = str2double(handles.edit_B_1z_max.String);
        
        % Nutfuellfaktor (Stator)
        var = strcmp(handles.opt.Spulenform_Stator,handles.default_richt.phi_1n_opt);
        bgCl(handles.edit_phi_1n,handles.default_richt.phi_1n(var~=0,:));
        handles.richt.phi_1n = str2double(handles.edit_phi_1n.String);
        
        % Wicklungsfaktor Grundwelle (Stator)
        bgCl(handles.edit_xi_1p,handles.default_richt.xi_1p(1,:));
        handles.richt.xi_1p = str2double(handles.edit_xi_1p.String);
        
        % min. Nutteilung (Stator)
        bgCl(handles.edit_tau_1n_min,handles.default_richt.tau_1n_min(1,:));
        handles.richt.tau_1n_min = str2double(handles.edit_tau_1n_min.String);
        
        % Eisenfuellfaktor (Stator)
        bgCl(handles.edit_phi_1Fe,handles.default_richt.phi_1Fe(1,:));
        handles.richt.phi_1Fe = str2double(handles.edit_phi_1Fe.String);
        
        % max. Rueckeninduktion (Rotor)
        bgCl(handles.edit_B_2r_max,handles.default_richt.B_2r_max(1,:));
        handles.richt.B_2r_max = str2double(handles.edit_B_2r_max.String);
        
        % Eisenfuellfaktor (Rotor)
        bgCl(handles.edit_phi_2Fe,handles.default_richt.phi_2Fe(1,:));
        handles.richt.phi_2Fe = str2double(handles.edit_phi_2Fe.String);
        
    else
        error('Ungueltige Eingabe bei Variable "handles.opt.Maschinentyp"');
    end
    
    % Save and Restore Richtwerte
    handles = MiscRichtwerte(handles);
    
    % Set Tooltip String
    handles = setTooltipString(handles);
    
    % Remove fields
    handles = removeFields(handles);

% Set default Richtwerte   
function handles = setDefaultRichtwerte(handles)
    if(strcmp(handles.opt.Maschinentyp,'ASM'))
        % Relative Ankerlaenge
        if(handles.rated.p > 1)
            handles.edit_lambda.String = num2str(handles.default_richt.lambda(2,3));
        else
            handles.edit_lambda.String = num2str(handles.default_richt.lambda(1,3));
        end
        
        % Kanalbreite der Ventilationskanaele
        handles.edit_l_v.String = num2str(handles.default_richt.l_v(1,3));
        
        % Luftspaltinduktion bzw. Strombelag
        handles.B_A = handles.popupmenu_B_A.String{handles.popupmenu_B_A.Value};
        if(strcmp(handles.B_A,'Mittelwert der Luftspaltinduktion in T'))
            handles.edit_B_A.String = num2str(handles.default_richt.B_m(1,3));
        elseif(strcmp(handles.B_A,'Strombelag in A/mm'))
            handles.edit_B_A.String = num2str(handles.default_richt.A(1,3));
        else
            error('Ungueltige Eingabe bei Variable "handles.B_A"');
        end
        
        % Stromdichte (Stator)
        handles.edit_S_1.String = num2str(handles.default_richt.S_1(1,3));
        
        % max. Rueckeninduktion (Stator)
        handles.edit_B_1r_max.String = num2str(handles.default_richt.B_1r_max(1,3));
        
        % max. Zahninduktion (Stator)
        handles.edit_B_1z_max.String = num2str(handles.default_richt.B_1z_max(1,3));
        
        % Nutfuellfaktor (Stator)
        var = strcmp(handles.opt.Spulenform_Stator,handles.default_richt.phi_1n_opt);
        handles.edit_phi_1n.String = num2str(handles.default_richt.phi_1n(var~=0,3));
        
        % Wicklungsfaktor Grundwelle (Stator)
        handles.edit_xi_1p.String = num2str(handles.default_richt.xi_1p(1,3));
        
        % min. Nutteilung (Stator)
        handles.edit_tau_1n_min.String = num2str(handles.default_richt.tau_1n_min(1,3));
        
        % Eisenfuellfaktor (Stator)
        handles.edit_phi_1Fe.String = num2str(handles.default_richt.phi_1Fe(1,3));
        
        % Stabstromdichte (Rotor)
        var = strcmp(handles.opt.Rotor_Leitermaterial.String,handles.default_richt.S_2s_opt);
        handles.edit_S_2s.String = num2str(handles.default_richt.S_2r(var~=0,3));
        % Ringstromdichte (Rotor)
        var = strcmp(handles.opt.Rotor_Leitermaterial.String,handles.default_richt.S_2r_opt);
        handles.edit_S_2r.String = num2str(handles.default_richt.S_2s(var~=0,3));

        % Stromdichte (Rotor)
        var = strcmp(handles.opt.Rotor_Leitermaterial.String,handles.default_richt.S_2_opt);
        handles.edit_S_2.String = num2str(handles.default_richt.S_2(var~=0,3));
        
        % max. Rueckeninduktion (Rotor)
        handles.edit_B_2r_max.String = num2str(handles.default_richt.B_2r_max(1,3));
        
        % max. Zahninduktion (Rotor)
        handles.edit_B_2z_max.String = num2str(handles.default_richt.B_2z_max(1,3));
        
        % Nutfuellfaktor (Rotor)
        var = strcmp(handles.opt.Spulenform_Rotor,handles.default_richt.phi_2n_opt);
        handles.edit_phi_2n.String = num2str(handles.default_richt.phi_2n(var~=0,3));
        
        % Wicklungsfaktor Grundwelle (Rotor)
        handles.edit_xi_2p.String = num2str(handles.default_richt.xi_2p(1,3));
        
        % min. Nutteilung (Rotor)
        handles.edit_tau_2n_min.String = num2str(handles.default_richt.tau_2n_min(1,3));
        
        % Eisenfuellfaktor (Rotor)
        handles.edit_phi_2Fe.String = num2str(handles.default_richt.phi_2Fe(1,3));
        
    elseif(strcmp(handles.opt.Maschinentyp,'PMSM'))
        % Relative Ankerlaenge
        if(handles.rated.p > 1)
            handles.edit_lambda.String = num2str(handles.default_richt.lambda(2,3));
        else
            handles.edit_lambda.String = num2str(handles.default_richt.lambda(1,3));
        end
        
        % Kanalbreite der Ventilationskanaele
        handles.edit_l_v.String = num2str(handles.default_richt.l_v(1,3));
        
        % Luftspaltinduktion bzw. Strombelag
        handles.B_A = handles.popupmenu_B_A.String{handles.popupmenu_B_A.Value};
        if(strcmp(handles.B_A,'Amplitude der Luftspaltinduktion in T'))
            handles.edit_B_A.String = num2str(handles.default_richt.B_p(1,3));
        elseif(strcmp(handles.B_A,'Strombelag in A/mm'))
            var = strcmp(handles.opt.Kuehlungsart,handles.default_richt.A_opt);
            handles.edit_B_A.String = num2str(handles.default_richt.A(var~=0,3));
        else
            error('Ungueltige Eingabe bei Variable "handles.B_A"');
        end
        
        % Stromdichte (Stator)
        var = strcmp(handles.opt.Kuehlungsart,handles.default_richt.S_1_opt);
        handles.edit_S_1.String = num2str(handles.default_richt.S_1(var~=0,3));
        
        % max. Rueckeninduktion (Stator)
        handles.edit_B_1r_max.String = num2str(handles.default_richt.B_1r_max(1,3));
        
        % max. Zahninduktion (Stator)
        handles.edit_B_1z_max.String = num2str(handles.default_richt.B_1z_max(1,3));
        
        % Nutfuellfaktor (Stator)
        var = strcmp(handles.opt.Spulenform_Stator,handles.default_richt.phi_1n_opt);
        handles.edit_phi_1n.String = num2str(handles.default_richt.phi_1n(var~=0,3));
        
        % Wicklungsfaktor Grundwelle (Stator)
        handles.edit_xi_1p.String = num2str(handles.default_richt.xi_1p(1,3));
        
        % min. Nutteilung (Stator)
        handles.edit_tau_1n_min.String = num2str(handles.default_richt.tau_1n_min(1,3));
        
        % Eisenfuellfaktor (Stator)
        handles.edit_phi_1Fe.String = num2str(handles.default_richt.phi_1Fe(1,3));
        
        % max. Rueckeninduktion (Rotor)
        handles.edit_B_2r_max.String = num2str(handles.default_richt.B_2r_max(1,3));
        
        % Eisenfuellfaktor (Rotor)
        handles.edit_phi_2Fe.String = num2str(handles.default_richt.phi_2Fe(1,3));
        
    else
        error('Ungueltige Eingabe bei Variable "handles.opt.Maschinentyp"');
    end

% Set default Optionen
function handles = setDefaultOptionen(handles)
    % Maschinenausfuehrung
    handles.popupmenu_Maschinenausfuehrung.String = handles.default_opt.Maschinenausfuehrung;
    handles.popupmenu_Maschinenausfuehrung.Value = 1;
    
    % Schaltung
    handles.popupmenu_Schaltung.String = handles.default_opt.Schaltung;
    handles.popupmenu_Schaltung.Value = 1;
    
    % Spulenform Stator
    handles.popupmenu_Spulenform_Stator.String = handles.default_opt.Spulenform_Stator;
    handles.popupmenu_Spulenform_Stator.Value = 1;
    
    % Spulenform Rotor
    handles.popupmenu_Spulenform_Rotor.String = handles.default_opt.Spulenform_Rotor;
    handles.popupmenu_Spulenform_Rotor.Value = 1;
    
    % Nutform Stator
    handles.popupmenu_Nutform_Stator.String = handles.default_opt.Nutform_Stator;
    handles.popupmenu_Nutform_Stator.Value = 1;
    
    % Nutform Rotor
    handles.popupmenu_Nutform_Rotor.String = handles.default_opt.Nutform_Rotor;
    handles.popupmenu_Nutform_Rotor.Value = 1;
    
    % Kuehlungsart
    handles.popupmenu_Kuehlungsart.String = handles.default_opt.Kuehlungsart;
    handles.popupmenu_Kuehlungsart.Value = 1;
    
    % Eisenmaterial Stator
    handles.popupmenu_Stator_Eisenmaterial.String = handles.default_opt.Stator_Eisenmaterial.String;
    handles.popupmenu_Stator_Eisenmaterial.Value = 1;
    
    % Leitermaterial Stator
    handles.popupmenu_Stator_Leitermaterial.String = handles.default_opt.Stator_Leitermaterial.String;
    handles.popupmenu_Stator_Leitermaterial.Value = 1;
    
    % Eisenmaterial Rotor
    handles.popupmenu_Rotor_Eisenmaterial.String = handles.default_opt.Rotor_Eisenmaterial.String;
    handles.popupmenu_Rotor_Eisenmaterial.Value = 1;
    
    % Leitermaterial Rotor
    handles.popupmenu_Rotor_Leitermaterial.String = handles.default_opt.Rotor_Leitermaterial.String;
    handles.popupmenu_Rotor_Leitermaterial.Value = 1;
    
    % Magnetmaterial Rotor
    handles.popupmenu_Rotor_Magnetmaterial.String = handles.default_opt.Rotor_Magnetmaterial.String;
    handles.popupmenu_Rotor_Magnetmaterial.Value = 1;
    
    % Modus Wicklungsauslegung
    handles.popupmenu_Mode_Wicklung.String = handles.default_opt.Mode_Wicklung;
    handles.popupmenu_Mode_Wicklung.Value = 1;
    
    % Wicklungstyp
    handles.popupmenu_Wicklungstyp.String = handles.default_opt.Wicklungstyp;
    handles.popupmenu_Wicklungstyp.Value = 1;

% Set Tooltip String
function handles = setTooltipString(handles)
    set(handles.edit_n_N, 'TooltipString', ['Beschraenkung zwischen ' num2str(handles.limit.n_N(1,2)) ' und ' num2str(handles.limit.n_N(1,1))])
    set(handles.edit_U_N, 'TooltipString', ['Beschraenkung zwischen ' num2str(handles.limit.U_N(1,2)) ' und ' num2str(handles.limit.U_N(1,1))])
    set(handles.edit_f_N, 'TooltipString', ['Beschraenkung zwischen ' num2str(handles.limit.f_N(1,2)) ' und ' num2str(handles.limit.f_N(1,1))])
    set(handles.edit_m, 'TooltipString', ['Beschraenkung zwischen ' num2str(handles.limit.m(1,2)) ' und ' num2str(handles.limit.m(1,1))])
    
    set(handles.popupmenu_Wicklungstyp, 'TooltipString', ['Wicklung A: Ganzlochwicklung, Einschicht, ungesehnt, gezont' 10 ...
    'Wicklung B: Ganzlochwicklung, Zweischicht, gesehnt, gezont' 10 ...
    'Wicklung C: Bruchlochwicklung, Zweischicht, gesehnt, gezont,' 10 ...
    'Zahnspulenwicklung (q<1) wird dann umgesetzt, wenn der Durchmesser q>1 nicht zulaesst'])
    
    if(handles.rated.p > 1)
        set(handles.edit_lambda, 'TooltipString', ['Richtwertbereich zwischen ' num2str(handles.default_richt.lambda(2,2)) ' und ' num2str(handles.default_richt.lambda(2,1))])
    else
        set(handles.edit_lambda, 'TooltipString', ['Richtwertbereich zwischen ' num2str(handles.default_richt.lambda(1,2)) ' und ' num2str(handles.default_richt.lambda(1,1))])
    end
    set(handles.edit_l_v, 'TooltipString', ['Richtwertbereich zwischen ' num2str(handles.default_richt.l_v(1,2)) ' und ' num2str(num2str(handles.default_richt.l_v(1,1)))])
    set(handles.edit_B_1r_max, 'TooltipString', ['Richtwertbereich zwischen ' num2str(handles.default_richt.B_1r_max(1,2)) ' und ' num2str(handles.default_richt.B_1r_max(1,1))])
    set(handles.edit_B_1z_max, 'TooltipString', ['Richtwertbereich zwischen ' num2str(handles.default_richt.B_1z_max(1,2)) ' und ' num2str(handles.default_richt.B_1z_max(1,1))])
    var = strcmp(handles.opt.Spulenform_Stator,handles.default_richt.phi_1n_opt);
    set(handles.edit_phi_1n, 'TooltipString', ['Richtwertbereich zwischen ' num2str(handles.default_richt.phi_1n(var~=0,2)) ' und ' num2str(handles.default_richt.phi_1n(var~=0,1))])
    set(handles.edit_xi_1p, 'TooltipString', ['Richtwertbereich zwischen ' num2str(handles.default_richt.xi_1p(1,2)) ' und ' num2str(handles.default_richt.xi_1p(1,1))])
    set(handles.edit_tau_1n_min, 'TooltipString', ['Richtwertbereich zwischen ' num2str(handles.default_richt.tau_1n_min(1,2)) ' und ' num2str(handles.default_richt.tau_1n_min(1,1))])
    set(handles.edit_phi_1Fe, 'TooltipString', ['Richtwertbereich zwischen ' num2str(handles.default_richt.phi_1Fe(1,2)) ' und ' num2str(handles.default_richt.phi_1Fe(1,1))])
    set(handles.edit_B_2r_max, 'TooltipString', ['Richtwertbereich zwischen ' num2str(handles.default_richt.B_2r_max(1,2)) ' und ' num2str(handles.default_richt.B_2r_max(1,1))])
    set(handles.edit_phi_2Fe, 'TooltipString', ['Richtwertbereich zwischen ' num2str(handles.default_richt.phi_2Fe(1,2)) ' und ' num2str(handles.default_richt.phi_2Fe(1,1))])
    
    if(strcmp(handles.opt.Maschinentyp,'ASM'))
        idx = find(handles.rated.p==handles.limit.p);
        if(isempty(idx))
            set(handles.edit_P_N, 'TooltipString', ['Beschraenkung zwischen ' num2str(NaN) ' und ' num2str(NaN)])
        else
            set(handles.edit_P_N, 'TooltipString', ['Beschraenkung zwischen ' num2str(handles.limit.P_N(idx,2)*1e-3) ' und ' num2str(handles.limit.P_N(idx,1)*1e-3)])
        end
        set(handles.edit_p, 'TooltipString', ['Beschraenkung auf die Polpaarzahlen ' num2str(handles.limit.p)])
        if(strcmp(handles.B_A,'Mittelwert der Luftspaltinduktion in T'))
            set(handles.edit_B_A, 'TooltipString', ['Richtwertbereich zwischen ' num2str(handles.default_richt.B_m(1,2)) ' und ' num2str(handles.default_richt.B_m(1,1))])
        elseif(strcmp(handles.B_A,'Strombelag in A/mm'))
            set(handles.edit_B_A, 'TooltipString', ['Richtwertbereich zwischen ' num2str(handles.default_richt.A(1,2)) ' und ' num2str(handles.default_richt.A(1,1))])
        else
            error('Ungueltige Eingabe bei Variable "handles.B_A"');
        end
        set(handles.edit_S_1, 'TooltipString', ['Richtwertbereich zwischen ' num2str(handles.default_richt.S_1(1,2)) ' und ' num2str(handles.default_richt.S_1(1,1))])
        var = strcmp(handles.opt.Rotor_Leitermaterial.String,handles.default_richt.S_2s_opt);
        set(handles.edit_S_2s, 'TooltipString', ['Richtwertbereich zwischen ' num2str(handles.default_richt.S_2s(var~=0,2)) ' und ' num2str(handles.default_richt.S_2s(var~=0,1))])
        var = strcmp(handles.opt.Rotor_Leitermaterial.String,handles.default_richt.S_2r_opt);
        set(handles.edit_S_2r, 'TooltipString', ['Richtwertbereich zwischen ' num2str(handles.default_richt.S_2r(var~=0,2)) ' und ' num2str(handles.default_richt.S_2r(var~=0,1))])
        var = strcmp(handles.opt.Rotor_Leitermaterial.String,handles.default_richt.S_2_opt);
        set(handles.edit_S_2, 'TooltipString', ['Richtwertbereich zwischen ' num2str(handles.default_richt.S_2(var~=0,2)) ' und ' num2str(handles.default_richt.S_2(var~=0,1))])
        set(handles.edit_B_2z_max, 'TooltipString', ['Richtwertbereich zwischen ' num2str(handles.default_richt.B_2z_max(1,2)) ' und ' num2str(handles.default_richt.B_2z_max(1,1))])
        var = strcmp(handles.opt.Spulenform_Stator,handles.default_richt.phi_2n_opt);
        set(handles.edit_phi_2n, 'TooltipString', ['Richtwertbereich zwischen ' num2str(handles.default_richt.phi_2n(var~=0,2)) ' und ' num2str(handles.default_richt.phi_2n(var~=0,1))])
        set(handles.edit_xi_2p, 'TooltipString', ['Richtwertbereich zwischen ' num2str(handles.default_richt.xi_2p(1,2)) ' und ' num2str(handles.default_richt.xi_2p(1,1))])
        set(handles.edit_tau_2n_min, 'TooltipString', ['Richtwertbereich zwischen ' num2str(handles.default_richt.tau_2n_min(1,2)) ' und ' num2str(handles.default_richt.tau_2n_min(1,1))])
    elseif(strcmp(handles.opt.Maschinentyp,'PMSM'))
        set(handles.edit_P_N, 'TooltipString', ['Beschraenkung zwischen ' num2str(handles.limit.P_N(1,2)*1e-3) ' und ' num2str(handles.limit.P_N(1,1)*1e-3)])
        set(handles.edit_n_N, 'TooltipString', ['Beschraenkung zwischen ' num2str(handles.limit.p(1,2)) ' und ' num2str(handles.limit.p(1,1))])
        set(handles.edit_cos_phi_N, 'TooltipString', ['Beschraenkung zwischen ' num2str(handles.limit.cos_phi_N(1,2)) ' und ' num2str(handles.limit.cos_phi_N(1,1))])
        if(strcmp(handles.B_A,'Amplitude der Luftspaltinduktion in T'))
            set(handles.edit_B_A, 'TooltipString', ['Richtwertbereich zwischen ' num2str(handles.default_richt.B_p(1,2)) ' und ' num2str(handles.default_richt.B_p(1,1))])
        elseif(strcmp(handles.B_A,'Strombelag in A/mm'))
            var = strcmp(handles.opt.Kuehlungsart,handles.default_richt.A_opt);
            set(handles.edit_B_A, 'TooltipString', ['Richtwertbereich zwischen ' num2str(handles.default_richt.A(var~=0,2)) ' und ' num2str(handles.default_richt.A(var~=0,1))])
        else
            error('Ungueltige Eingabe bei Variable "handles.B_A"');
        end
        var = strcmp(handles.opt.Kuehlungsart,handles.default_richt.S_1_opt);
        set(handles.edit_S_1, 'TooltipString', ['Richtwertbereich zwischen ' num2str(handles.default_richt.S_1(var~=0,2)) ' und ' num2str(handles.default_richt.S_1(var~=0,1))])
    else
        error('Ungueltige Eingabe bei Variable "handles.opt.Maschinentyp"');
    end

% Remove fields
function handles = removeFields(handles)
    if(strcmp(handles.opt.Maschinentyp,'ASM'))
        if(isfield(handles.rated, 'cos_phi_N'))
            handles.rated = rmfield(handles.rated, 'cos_phi_N');
        end
        if(strcmp(handles.B_A,'Mittelwert der Luftspaltinduktion in T'))
            if(isfield(handles.richt, 'B_p'))
                handles.richt = rmfield(handles.richt, 'B_p');
            end
            if(isfield(handles.richt, 'A'))
                handles.richt = rmfield(handles.richt, 'A');
            end
        elseif(strcmp(handles.B_A,'Strombelag in A/mm'))
            if(isfield(handles.richt, 'B_m'))
                handles.richt = rmfield(handles.richt, 'B_m');
            end
            if(isfield(handles.richt, 'B_p'))
                handles.richt = rmfield(handles.richt, 'B_p');
            end
        else
            error('Ungueltige Eingabe bei Variable "handles.B_A"');
        end
        if(isfield(handles.opt, 'Rotor_Magnetmaterial'))
            handles.opt = rmfield(handles.opt, 'Rotor_Magnetmaterial');
        end
    elseif(strcmp(handles.opt.Maschinentyp,'PMSM'))
        if(strcmp(handles.B_A,'Amplitude der Luftspaltinduktion in T'))
            if(isfield(handles.richt, 'B_m'))
                handles.richt = rmfield(handles.richt, 'B_m');
            end
            if(isfield(handles.richt, 'A'))
                handles.richt = rmfield(handles.richt, 'A');
            end
        elseif(strcmp(handles.B_A,'Strombelag in A/mm'))
            if(isfield(handles.richt, 'B_m'))
                handles.richt = rmfield(handles.richt, 'B_m');
            end
            if(isfield(handles.richt, 'B_p'))
                handles.richt = rmfield(handles.richt, 'B_p');
            end
        else
            error('Ungueltige Eingabe bei Variable "handles.B_A"');
        end
        if(isfield(handles.richt, 'S_2'))
            handles.richt = rmfield(handles.richt, 'S_2');
        end
        if(isfield(handles.richt, 'S_2s'))
            handles.richt = rmfield(handles.richt, 'S_2s');
        end
        if(isfield(handles.richt, 'S_2r'))
            handles.richt = rmfield(handles.richt, 'S_2r');
        end
        if(isfield(handles.richt, 'phi_2n'))
            handles.richt = rmfield(handles.richt, 'phi_2n');
        end
        if(isfield(handles.richt, 'xi_2p'))
            handles.richt = rmfield(handles.richt, 'xi_2p');
        end
        if(isfield(handles.richt, 'tau_2n_min'))
            handles.richt = rmfield(handles.richt, 'tau_2n_min');
        end
        if(isfield(handles.opt, 'Spulenform_Rotor'))
            handles.opt = rmfield(handles.opt, 'Spulenform_Rotor');
        end
        if(isfield(handles.opt, 'Nutform_Rotor'))
            handles.opt = rmfield(handles.opt, 'Nutform_Rotor');
        end
        if(isfield(handles.opt, 'Rotor_Leitermaterial'))
            handles.opt = rmfield(handles.opt, 'Rotor_Leitermaterial');
        end
        if(isfield(handles.opt, 'theta_2'))
            handles.opt = rmfield(handles.opt, 'theta_2');
        end
    else
        error('Ungueltige Eingabe bei Variable "handles.opt.Maschinentyp"');
    end

% Save and Restore Richtwerte
function handles = MiscRichtwerte(handles)
    % Relative Ankerlaenge
    var = str2double(handles.edit_lambda.String);
    handles.edit_lambda.String = num2str(0);
    handles.edit_lambda.String = num2str(var);

    % Kanalbreite der Ventilationskanaele
    var = str2double(handles.edit_l_v.String);
    handles.edit_l_v.String = num2str(0);
    handles.edit_l_v.String = num2str(var);

    % Luftspaltinduktion bzw. Strombelag
    var = str2double(handles.edit_B_A.String);
    handles.edit_B_A.String = num2str(0);
    handles.edit_B_A.String = num2str(var);

    % Stromdichte (Stator)
    var = str2double(handles.edit_S_1.String);
    handles.edit_S_1.String = num2str(0);
    handles.edit_S_1.String = num2str(var);

    % max. Rueckeninduktion (Stator)
    var = str2double(handles.edit_B_1r_max.String);
    handles.edit_B_1r_max.String = num2str(0);
    handles.edit_B_1r_max.String = num2str(var);

    % max. Zahninduktion (Stator)
    var = str2double(handles.edit_B_1z_max.String);
    handles.edit_B_1z_max.String = num2str(0);
    handles.edit_B_1z_max.String = num2str(var);

    % Nutfuellfaktor (Stator)
    var = str2double(handles.edit_phi_1n.String);
    handles.edit_phi_1n.String = num2str(0);
    handles.edit_phi_1n.String = num2str(var);

    % Wicklungsfaktor Grundwelle (Stator)
    var = str2double(handles.edit_xi_1p.String);
    handles.edit_xi_1p.String = num2str(0);
    handles.edit_xi_1p.String = num2str(var);

    % min. Nutteilung (Stator)
    var = str2double(handles.edit_tau_1n_min.String);
    handles.edit_tau_1n_min.String = num2str(0);
    handles.edit_tau_1n_min.String = num2str(var);

    % Eisenfuellfaktor (Stator)
    var = str2double(handles.edit_phi_1Fe.String);
    handles.edit_phi_1Fe.String = num2str(0);
    handles.edit_phi_1Fe.String = num2str(var);
    
    % Stabstromdichte (Rotor)
    var = str2double(handles.edit_S_2s.String);
    handles.edit_S_2s.String = num2str(0);
    handles.edit_S_2s.String = num2str(var);
    % Ringstromdichte (Rotor)
    var = str2double(handles.edit_S_2r.String);
    handles.edit_S_2r.String = num2str(0);
    handles.edit_S_2r.String = num2str(var);

    % Stromdichte (Rotor)
    var = str2double(handles.edit_S_2.String);
    handles.edit_S_2.String = num2str(0);
    handles.edit_S_2.String = num2str(var);

    % max. Rueckeninduktion (Rotor)
    var = str2double(handles.edit_B_2r_max.String);
    handles.edit_B_2r_max.String = num2str(0);
    handles.edit_B_2r_max.String = num2str(var);

    % max. Zahninduktion (Rotor)
    var = str2double(handles.edit_B_2z_max.String);
    handles.edit_B_2z_max.String = num2str(0);
    handles.edit_B_2z_max.String = num2str(var);

    % Nutfuellfaktor (Rotor)
    var = str2double(handles.edit_phi_2n.String);
    handles.edit_phi_2n.String = num2str(0);
    handles.edit_phi_2n.String = num2str(var);

    % Wicklungsfaktor Grundwelle (Rotor)
    var = str2double(handles.edit_xi_2p.String);
    handles.edit_xi_2p.String = num2str(0);
    handles.edit_xi_2p.String = num2str(var);

    % min. Nutteilung (Rotor)
    var = str2double(handles.edit_tau_2n_min.String);
    handles.edit_tau_2n_min.String = num2str(0);
    handles.edit_tau_2n_min.String = num2str(var);

    % Eisenfuellfaktor (Rotor)
    var = str2double(handles.edit_phi_2Fe.String);
    handles.edit_phi_2Fe.String = num2str(0);
    handles.edit_phi_2Fe.String = num2str(var);

% Set from file
function handles = setFromFile(handles)
    % Nennleistung
    handles.edit_P_N.String = num2str(handles.rated.P_N/1000);

    % Nenndrehzahl
    handles.edit_n_N.String = num2str(handles.rated.n_N);

    % Nennspannung
    handles.edit_U_N.String = num2str(handles.rated.U_N);

    % Polpaarzahl
    handles.edit_p.String = num2str(handles.rated.p);

    % Nennfrequenz
    handles.edit_f_N.String = num2str(handles.rated.f_N);

    % Strangzahl
    handles.edit_m.String = num2str(handles.rated.m);
    
    % Relative Ankerlaenge
    handles.edit_lambda.String = num2str(handles.richt.lambda);

    % Kanalbreite der Ventilationskanaele
    handles.edit_l_v.String = num2str(handles.richt.l_v);

    % Luftspaltinduktion bzw. Strombelag
    if(isfield(handles.richt,'B_m'))
        idx = find(strcmp(handles.popupmenu_B_A.String,'Mittelwert der Luftspaltinduktion in T'));
        handles.popupmenu_B_A.Value = idx;
        handles.edit_B_A.String = num2str(handles.richt.B_m);
    elseif(isfield(handles.richt,'B_p'))
        idx = find(strcmp(handles.popupmenu_B_A.String,'Amplitude der Luftspaltinduktion in T'));
        handles.popupmenu_B_A.Value = idx;
        handles.edit_B_A.String = num2str(handles.richt.B_p);
    elseif(isfield(handles.richt,'A'))
        idx = find(strcmp(handles.popupmenu_B_A.String,'Strombelag in A/mm'));
        handles.popupmenu_B_A.Value = idx;
        handles.edit_B_A.String = num2str(handles.richt.A);
    else
        error('Kein Wert fuer die Luftspaltinduktion "B_m" bzw. "B_p" oder den Strombelag "A" gefunden.')
    end

    % Stromdichte (Stator)
    handles.edit_S_1.String = num2str(handles.richt.S_1);

    % max. Rueckeninduktion (Stator)
    handles.edit_B_1r_max.String = num2str(handles.richt.B_1r_max);

    % max. Zahninduktion (Stator)
    handles.edit_B_1z_max.String = num2str(handles.richt.B_1z_max);

    % Nutfuellfaktor (Stator)
    handles.edit_phi_1n.String = num2str(handles.richt.phi_1n);

    % Wicklungsfaktor Grundwelle (Stator)
    handles.edit_xi_1p.String = num2str(handles.richt.xi_1p);

    % min. Nutteilung (Stator)
    handles.edit_tau_1n_min.String = num2str(handles.richt.tau_1n_min);

    % Eisenfuellfaktor (Stator)
    handles.edit_phi_1Fe.String = num2str(handles.richt.phi_1Fe);
    
    % max. Rueckeninduktion (Rotor)
    handles.edit_B_2r_max.String = num2str(handles.richt.B_2r_max);
    
    % Eisenfuellfaktor (Rotor)
    handles.edit_phi_2Fe.String = num2str(handles.richt.phi_2Fe);
    
    % Maschinenausfuehrung
    idx = find(strcmp(handles.popupmenu_Maschinenausfuehrung.String,handles.opt.Maschinenausfuehrung));
    handles.popupmenu_Maschinenausfuehrung.Value = idx;

    % Schaltung
    idx = find(strcmp(handles.popupmenu_Schaltung.String,handles.opt.Schaltung));
    handles.popupmenu_Schaltung.Value = idx;

    % Spulenform Stator
    idx = find(strcmp(handles.popupmenu_Spulenform_Stator.String,handles.opt.Spulenform_Stator));
    handles.popupmenu_Spulenform_Stator.Value = idx;

    % Nutform Stator
    idx = find(strcmp(handles.popupmenu_Nutform_Stator.String,handles.opt.Nutform_Stator));
    handles.popupmenu_Nutform_Stator.Value = idx;

    % Kuehlungsart
    idx = find(strcmp(handles.popupmenu_Kuehlungsart.String,handles.opt.Kuehlungsart));
    handles.popupmenu_Kuehlungsart.Value = idx;

    % Eisenmaterial Stator
    idx = find(strcmp(handles.popupmenu_Stator_Eisenmaterial.String,handles.opt.Stator_Eisenmaterial.String));
    handles.popupmenu_Stator_Eisenmaterial.Value = idx;

    % Leitermaterial Stator
    idx = find(strcmp(handles.popupmenu_Stator_Leitermaterial.String,handles.opt.Stator_Leitermaterial.String));
    handles.popupmenu_Stator_Leitermaterial.Value = idx;

    % Eisenmaterial Rotor
    idx = find(strcmp(handles.popupmenu_Rotor_Eisenmaterial.String,handles.opt.Rotor_Eisenmaterial.String));
    handles.popupmenu_Rotor_Eisenmaterial.Value = idx;

    % Modus Wicklungsauslegung
    idx = find(strcmp(handles.popupmenu_Mode_Wicklung.String,handles.opt.Mode_Wicklung));
    handles.popupmenu_Mode_Wicklung.Value = idx;

    % Wicklungstyp
    if(strcmp(handles.opt.Mode_Wicklung,'Klassisch'))
        handles.popupmenu_Wicklungstyp.String = handles.default_opt.Wicklungstyp;
        idx = find(strcmp(handles.popupmenu_Wicklungstyp.String,handles.opt.Wicklungstyp));
        handles.popupmenu_Wicklungstyp.Value = idx;
    else
        idx = find(strcmp(handles.popupmenu_Wicklungstyp.String,handles.opt.Wicklungstyp));
        handles.popupmenu_Wicklungstyp.Value = idx;
        handles.popupmenu_Wicklungstyp.String = {'-'};
    end
    
    if(strcmp(handles.opt.Maschinentyp,'ASM'))
        % Stabstromdichte (Rotor)
        handles.edit_S_2s.String = num2str(handles.richt.S_2s);
        
        % Ringstromdichte (Rotor)
        handles.edit_S_2r.String = num2str(handles.richt.S_2r);

        % Stromdichte (Rotor)
        handles.edit_S_2.String = num2str(handles.richt.S_2);

        % max. Zahninduktion (Rotor)
        handles.edit_B_2z_max.String = num2str(handles.richt.B_2z_max);

        % Nutfuellfaktor (Rotor)
        handles.edit_phi_2n.String = num2str(handles.richt.phi_2n);

        % Wicklungsfaktor Grundwelle (Rotor)
        handles.edit_xi_2p.String = num2str(handles.richt.xi_2p);

        % min. Nutteilung (Rotor)
        handles.edit_tau_2n_min.String = num2str(handles.richt.tau_2n_min);

        % Nutform Rotor
        idx = find(strcmp(handles.popupmenu_Nutform_Rotor.String,handles.opt.Nutform_Rotor));
        handles.popupmenu_Nutform_Rotor.Value = idx;
        
        % Spulenform Rotor
        idx = find(strcmp(handles.popupmenu_Spulenform_Rotor.String,handles.opt.Spulenform_Rotor));
        handles.popupmenu_Spulenform_Rotor.Value = idx;
        
        % Leitermaterial Rotor
        idx = find(strcmp(handles.popupmenu_Rotor_Leitermaterial.String,handles.opt.Rotor_Leitermaterial.String));
        handles.popupmenu_Rotor_Leitermaterial.Value = idx;
    elseif(strcmp(handles.opt.Maschinentyp,'PMSM'))
        % Leistungsfaktor
        handles.edit_cos_phi_N.String = num2str(handles.rated.cos_phi_N);
        
        % Magnetmaterial Rotor
        idx = find(strcmp(handles.popupmenu_Rotor_Magnetmaterial.String,handles.opt.Rotor_Magnetmaterial.String));
        handles.popupmenu_Rotor_Magnetmaterial.Value = idx;
    else
        error('Ungueltige Eingabe bei Variable "handles.opt.Maschinentyp"');
    end

% Set Ergebnisse
function handles = setErgebnisse(handles)
    % Innendurchmesser Stator
    handles.edit_D_1i.String = num2str(handles.Entwurf.Geometrie.D_1i);

    % Aussendurchmesser Stator
    handles.edit_D_1a.String = num2str(handles.Entwurf.Geometrie.D_1a);

    % Ideelle Laenge
    handles.edit_l.String = num2str(handles.Entwurf.Geometrie.l);

    % Bohrungsvolumen
    handles.edit_V_Bohrung.String = num2str(handles.Entwurf.Geometrie.V_Bohrung);

    % Laenge geometrischer Luftspalt
    handles.edit_delta.String = num2str(handles.Entwurf.Geometrie.delta);
    
    if(strcmp(handles.opt.Maschinentyp,'ASM'))
        % Ungesaettigte Hauptinduktivitaet
        handles.edit_L_1h_ASM.String = num2str(handles.Entwurf.EMAG.L_1h);

        % Streuinduktivitaet Stator
        handles.edit_L_1sigma_ASM.String = num2str(handles.Entwurf.EMAG.L_1sigma);
        
        % Streuinduktivitaet Rotor
        handles.edit_L_2sigma_ASM.String = num2str(handles.Entwurf.EMAG.L_2sigma);
        
        % Statorwiderstand
        handles.edit_R_1_ASM.String = num2str(handles.Entwurf.EMAG.R_1);
        
        % Rotorwiderstand
        handles.edit_R_2_ASM.String = num2str(handles.Entwurf.EMAG.R_2_trans);
    elseif(strcmp(handles.opt.Maschinentyp,'PMSM'))
        % Ungesaettigte Hauptinduktivitaet
        handles.edit_L_1h_PMSM.String = num2str(handles.Entwurf.EMAG.L_1h);

        % Streuinduktivitaet Stator
        handles.edit_L_1sigma_PMSM.String = num2str(handles.Entwurf.EMAG.L_1sigma);
        
        % Verketteter Fluss der PM
        handles.edit_Psi_PM_PMSM.String = num2str(handles.Entwurf.EMAG.Psi_PM);
        
        % Induktivitaet d-Achse
        handles.edit_L_d_PMSM.String = num2str(handles.Entwurf.EMAG.L_d*1e3);

        % Induktivitaet q-Achse
        handles.edit_L_q_PMSM.String = num2str(handles.Entwurf.EMAG.L_q*1e3);
        
        if(strcmp(handles.opt.Maschinenausfuehrung,'IPMSM (eingelassen)'))
            handles.edit_b_PM.String = num2str(handles.Entwurf.Geometrie.b_PM);
            handles.edit_h_PM.String = num2str(handles.Entwurf.Geometrie.h_PM);
            handles.edit_Abstand_PM_Rotoroberflaeche.String = '-';
            handles.edit_Abstand_PM_Rotoroberflaeche.Enable = 'off';
            handles.edit_alpha_PM.String = '-';
            handles.edit_alpha_PM.Enable = 'off';
            handles.edit_Abstand_PM_unten.String = '-';
            handles.edit_Abstand_PM_unten.Enable = 'off';
        elseif(strcmp(handles.opt.Maschinenausfuehrung,'IPMSM (tangential)'))
            handles.edit_b_PM.String = num2str(handles.Entwurf.Geometrie.b_PM);
            handles.edit_h_PM.String = num2str(handles.Entwurf.Geometrie.h_PM);
            handles.edit_Abstand_PM_Rotoroberflaeche.String = num2str(handles.Entwurf.Geometrie.Abstand_PM_Rotoroberflaeche);
            handles.edit_Abstand_PM_Rotoroberflaeche.Enable = 'on';
            handles.edit_alpha_PM.String = '-';
            handles.edit_alpha_PM.Enable = 'off';
            handles.edit_Abstand_PM_unten.String = '-';
            handles.edit_Abstand_PM_unten.Enable = 'off';
        elseif(strcmp(handles.opt.Maschinenausfuehrung,'IPMSM (V-Form)'))
            handles.edit_b_PM.String = num2str(handles.Entwurf.Geometrie.b_PM);
            handles.edit_h_PM.String = num2str(handles.Entwurf.Geometrie.h_PM);
            handles.edit_Abstand_PM_Rotoroberflaeche.String = num2str(handles.Entwurf.Geometrie.Abstand_PM_Rotoroberflaeche);
            handles.edit_Abstand_PM_Rotoroberflaeche.Enable = 'on';
            handles.edit_alpha_PM.String = num2str(handles.Entwurf.Geometrie.alpha_PM);
            handles.edit_alpha_PM.Enable = 'on';
            handles.edit_Abstand_PM_unten.String = num2str(handles.Entwurf.Geometrie.Abstand_PM_unten);
            handles.edit_Abstand_PM_unten.Enable = 'on';
        end
        
    else
        error('Ungueltige Eingabe bei Variable "handles.opt.Maschinentyp"');
    end
    
% Plot B-H Kurve
function plotBH(var)
    figure('Name','B-H Kurve')
    title(var.String,'interpreter','latex','FontSize', 20)
    hold on
    grid on
    ax = gca;
    ax.YMinorGrid = 'on';
    ax.XScale = 'log';
    xlim([10 max(var.H)*1.1])
    xlabel('Feldstaerke $H$ in A/m','interpreter','latex','FontSize', 15)

    if(isfield(var, 'mu_r'))
        yyaxis left
        h1 = plot(var.H, var.B, 'o');
        h2 = plot(var.H, var.B);
        ylim([0 2.5])
        ylabel('Magnetische Induktion $B$ in T','interpreter','latex','FontSize', 15)

        yyaxis right
        h3 = plot(var.H, var.mu_r, 'o');
        h4 = plot(var.H, var.mu_r);
        ylim([0 4e4])
        ylabel('Relative Permeabilitaet $\mu_r$','interpreter','latex','FontSize', 15)

        legend([h2 h4], {'$B(H)$','$\mu(H)$'}, 'Location','northwest','interpreter','latex')
    else
        h1 = plot(var.H, var.B, 'o');
        h2 = plot(var.H, var.B);
        ylim([0 2.5])
        ylabel('Magnetische Induktion $B$ in T','interpreter','latex','FontSize', 15)

        legend(h2, {'$B(H)$'}, 'Location','northwest','interpreter','latex')
    end

% Plot Eisenverlust Kurve
function plotVerluste(var)
    figure('Name','Eisenverlust Kurve')
    title(var.String,'interpreter','latex','FontSize', 20)
    hold on
    grid on
    ax = gca;
    ax.YMinorGrid = 'on';
    xlim([0 max(var.p_vFe_f_vec)*1.1])
    xlabel('Frequenz $f$ in Hz','interpreter','latex','FontSize', 15)
    
    legendString = '';
    for i = 1:length(var.p_vFe_B_vec)
        plot(var.p_vFe_f_vec, var.p_vFe(i,:));
        legendString = [legendString; {['$B = ' num2str(var.p_vFe_B_vec(i)) ' T$']}];
    end
    
    ylim([0 max(max(var.p_vFe))*1.2])
    ylabel('Spezifische Eisenverluste $p_{Fe}$ in W/kg','interpreter','latex','FontSize', 15)
    
    
    legend(legendString, 'Location','northwest','interpreter','latex')

% Plot Maschine
function plotMaschine(handles,var)
    hold(handles.axes_Ergebnis_Plot,'on')
    box(handles.axes_Ergebnis_Plot,'on')
    
    % Plot
    plot(handles.axes_Ergebnis_Plot,var.Stator_aussen_x,var.Stator_aussen_y,'k')
    plot(handles.axes_Ergebnis_Plot,var.Stator_innen_x,var.Stator_innen_y,'k')
    plot(handles.axes_Ergebnis_Plot,var.Rotor_aussen_x,var.Rotor_aussen_y,'k')
    plot(handles.axes_Ergebnis_Plot,var.Rotor_innen_x,var.Rotor_innen_y,'k')
    for i = 1:handles.Entwurf.Wicklung.N_1
        plot(handles.axes_Ergebnis_Plot,var.Fuellung_Stator_x((length(var.Fuellung_Stator_x)/handles.Entwurf.Wicklung.N_1)*(i-1)+1:(length(var.Fuellung_Stator_x)/handles.Entwurf.Wicklung.N_1)*i,1),var.Fuellung_Stator_y((length(var.Fuellung_Stator_y)/handles.Entwurf.Wicklung.N_1)*(i-1)+1:(length(var.Fuellung_Stator_y)/handles.Entwurf.Wicklung.N_1)*i,1),'Color',[0 115 189]./255)
    end
    if(strcmp(handles.opt.Maschinentyp,'ASM'))
        for i = 1:handles.Entwurf.Wicklung.N_2
            fill(handles.axes_Ergebnis_Plot,var.Fuellung_Rotor_x((length(var.Fuellung_Rotor_x)/handles.Entwurf.Wicklung.N_2)*(i-1)+1:(length(var.Fuellung_Rotor_x)/handles.Entwurf.Wicklung.N_2)*i,1),var.Fuellung_Rotor_y((length(var.Fuellung_Rotor_y)/handles.Entwurf.Wicklung.N_2)*(i-1)+1:(length(var.Fuellung_Rotor_y)/handles.Entwurf.Wicklung.N_2)*i,1),[100 160 200]./255)
        end
    elseif(strcmp(handles.opt.Maschinentyp,'PMSM'))
        if(strcmp(handles.opt.Maschinenausfuehrung,'SPMSM'))
            for i = 1:(2*handles.Entwurf.Bemessungswerte.p)
                fill(handles.axes_Ergebnis_Plot,var.Magnet_Rotor_x(var.points_Magnet(i):var.points_Magnet(i+1)-1,1),var.Magnet_Rotor_y(var.points_Magnet(i):var.points_Magnet(i+1)-1,1),[100 160 200]./255)
            end
        elseif(strcmp(handles.opt.Maschinenausfuehrung,'IPMSM (eingelassen)'))
            for i = 1:(2*handles.Entwurf.Bemessungswerte.p)
                fill(handles.axes_Ergebnis_Plot,var.Magnet_Rotor_x(var.points_Magnet(i):var.points_Magnet(i+1)-1,1),var.Magnet_Rotor_y(var.points_Magnet(i):var.points_Magnet(i+1)-1,1),[100 160 200]./255)
            end
        elseif(strcmp(handles.opt.Maschinenausfuehrung,'IPMSM (tangential)'))
            for i = 1:(2*handles.Entwurf.Bemessungswerte.p)
                fill(handles.axes_Ergebnis_Plot,var.Magnet_Rotor_x((length(var.Magnet_Rotor_x)/(2*handles.Entwurf.Bemessungswerte.p))*(i-1)+1:(length(var.Magnet_Rotor_x)/(2*handles.Entwurf.Bemessungswerte.p))*i,1),var.Magnet_Rotor_y((length(var.Magnet_Rotor_y)/(2*handles.Entwurf.Bemessungswerte.p))*(i-1)+1:(length(var.Magnet_Rotor_y)/(2*handles.Entwurf.Bemessungswerte.p))*i,1),[100 160 200]./255)
            end
        elseif(strcmp(handles.opt.Maschinenausfuehrung,'IPMSM (V-Form)'))
            for i = 1:(4*handles.Entwurf.Bemessungswerte.p)
                fill(handles.axes_Ergebnis_Plot,var.Magnet_Rotor_x((length(var.Magnet_Rotor_x)/(4*handles.Entwurf.Bemessungswerte.p))*(i-1)+1:(length(var.Magnet_Rotor_x)/(4*handles.Entwurf.Bemessungswerte.p))*i,1),var.Magnet_Rotor_y((length(var.Magnet_Rotor_y)/(4*handles.Entwurf.Bemessungswerte.p))*(i-1)+1:(length(var.Magnet_Rotor_y)/(4*handles.Entwurf.Bemessungswerte.p))*i,1),[100 160 200]./255)
            end
        else
            error('Ungueltige Eingabe bei Variable "opt.Maschinenausfuehrung"')
        end
    else
        error('Ungueltige Eingabe bei Variable "handles.opt.Maschinentyp"');
    end
    
    % Wicklung Stator
    %var1 = handles.Entwurf.Wicklung.N_1/handles.Entwurf.Bemessungswerte.m;
    [~,idx] = max(var.Fuellung_Stator_y);
    coord_2 = [0 var.Fuellung_Stator_y(idx)*0.97; 0 var.Fuellung_Stator_y(idx-1)*1.03];
    theta = linspace(0,2*pi,handles.Entwurf.Wicklung.N_1+1);
    theta(end) = [];
    color_3 = [{'r'},{'b'},{'g'}];
    for k = 1:handles.Entwurf.Wicklung.n_lay
        %M6 = handles.Entwurf.Wicklung.Matrix(:,var1*(k-1)+1:var1*k);
        %M7 = sort(M6,2,'ComparisonMethod','abs');
        coord = coord_2(k,:);
        for j = 1:handles.Entwurf.Bemessungswerte.m
            if(k==1)
                M6 = handles.Entwurf.Wicklung.Matrix(:,1:length(handles.Entwurf.Wicklung.Matrix_lay{j,k}));
            elseif(k==2)
                M6 = handles.Entwurf.Wicklung.Matrix(:,length(handles.Entwurf.Wicklung.Matrix_lay{j,k-1})+1:end);
            else
                error('')
            end
            M7 = sort(M6,2,'ComparisonMethod','abs');
            vec = M7(j,:);
            color = color_3{j};
            for i = 1:length(vec)
                theta_1 = theta(abs(vec(i)));
                M = [cos(theta_1) -sin(theta_1); sin(theta_1) cos(theta_1)];
                coord_trans = coord * M;
                if(sign(vec(i)<1))
                    plot(handles.axes_Ergebnis_Plot,coord_trans(1,1),coord_trans(1,2),'Color',color,'Marker','x')
                else
                    plot(handles.axes_Ergebnis_Plot,coord_trans(1,1),coord_trans(1,2),'Color',color,'Marker','o')
                end
            end
        end
    end

    % Anpassung figure
    axis(handles.axes_Ergebnis_Plot,'equal') 
    ylim(handles.axes_Ergebnis_Plot,[-((var.D_1a*1e3)/2)*1.1 ((var.D_1a*1e3)/2)*1.1]);
    xlim(handles.axes_Ergebnis_Plot,[-((var.D_1a*1e3)/2)*1.1 ((var.D_1a*1e3)/2)*1.1]);
    grid(handles.axes_Ergebnis_Plot,'off')
    yticks(handles.axes_Ergebnis_Plot,'auto')
    xticks(handles.axes_Ergebnis_Plot,'auto')
    ylabel(handles.axes_Ergebnis_Plot,'y in mm','interpreter','latex','FontSize', 15)
    xlabel(handles.axes_Ergebnis_Plot,'x in mm','interpreter','latex','FontSize', 15)
    title(handles.axes_Ergebnis_Plot,'Abmessungen Maschine','interpreter','latex','FontSize', 20)
    hold(handles.axes_Ergebnis_Plot,'off')
    
    axtoolbar(handles.axes_Ergebnis_Plot,{'pan','zoomin','zoomout','restoreview'});
    handles.axes_Ergebnis_Plot.Toolbar.Visible = 'on';

% Plot Stator-Nut
function plotNut_1(handles,var)
    hold(handles.axes_Ergebnis_Plot,'on')
    box(handles.axes_Ergebnis_Plot,'on')
    
    % Nutform 'Trapezform (eckig)'
    % Trapeznut-Koordinaten
    % x-Nullpunkt liegt genau in der Mitte der Nut
    % y-Nullpunkt liegt in Nutoeffnung
    % Skizze in SA
    x1_Nut = var.b_1ns / 2;
    x2_Nut = var.b_1n_u / 2;
    x3_Nut = var.b_1n_o / 2;
    y1_Nut = 0;
    y2_Nut = var.h_1ns;
    y3_Nut = var.h_1ns + tan(var.alpha_1nk) * ((var.b_1n_u-var.b_1ns) / 2);
    y4_Nut = var.h_1n;
    Nut = [x1_Nut,x1_Nut,x2_Nut,x3_Nut,-x3_Nut,-x2_Nut,-x1_Nut,-x1_Nut; y1_Nut,y2_Nut,y3_Nut,y4_Nut,y4_Nut,y3_Nut,y2_Nut,y1_Nut]';
    
    % Nutfuellung
    % Fuell-Koordinaten
    x1_Fuellung = x2_Nut - var.d_1iso;
    x2_Fuellung = x3_Nut - var.d_1iso;
    y1_Fuellung = var.h_1nk;
    y2_Fuellung = var.h_1n - var.d_1iso;
    Fuellung = [x1_Fuellung,x2_Fuellung,-x2_Fuellung,-x1_Fuellung; y1_Fuellung,y2_Fuellung,y2_Fuellung,y1_Fuellung]';
    
    % Plot
    plot(handles.axes_Ergebnis_Plot,Nut(:,1),Nut(:,2),'k')
    fill(handles.axes_Ergebnis_Plot,Fuellung(:,1),Fuellung(:,2),[0 115 189]./255)

    % Anpassung figure
    if(y4_Nut>=var.b_1n_o)
        ylim(handles.axes_Ergebnis_Plot,[-0.5,y4_Nut+0.5]);
        xlim(handles.axes_Ergebnis_Plot,[-y4_Nut/2-0.5,y4_Nut/2+0.5]);
    else
        ylim(handles.axes_Ergebnis_Plot,[-0.5,x3_Nut*2+0.5]);
        xlim(handles.axes_Ergebnis_Plot,[-x3_Nut-0.5,x3_Nut+0.5]);
    end
    axis(handles.axes_Ergebnis_Plot,'equal') 
    grid(handles.axes_Ergebnis_Plot,'on')
    yticks(handles.axes_Ergebnis_Plot,0:1:y4_Nut+0.5)
    xticks(handles.axes_Ergebnis_Plot,fix(min(handles.axes_Ergebnis_Plot.XLim)):1:fix(max(handles.axes_Ergebnis_Plot.XLim)))
    ylabel(handles.axes_Ergebnis_Plot,'y in mm','interpreter','latex','FontSize', 15)
    xlabel(handles.axes_Ergebnis_Plot,'x in mm','interpreter','latex','FontSize', 15)
    title(handles.axes_Ergebnis_Plot,'Abmessungen Statornut','interpreter','latex','FontSize', 20)
    hold(handles.axes_Ergebnis_Plot, 'off')
    
    axtoolbar(handles.axes_Ergebnis_Plot,{'pan','zoomin','zoomout','restoreview'});
    handles.axes_Ergebnis_Plot.Toolbar.Visible = 'on';
    
% Plot Rotor-Nut
function plotNut_2(handles,var)
    hold(handles.axes_Ergebnis_Plot,'on')
    box(handles.axes_Ergebnis_Plot,'on')
    
    % Nutform 'Trapezform (eckig)'
    % Trapeznut-Koordinaten
    % x-Nullpunkt liegt genau in der Mitte der Nut
    % y-Nullpunkt liegt in Nutoeffnung
    % Skizze in SA
    x1_Nut = var.b_2ns / 2;
    x2_Nut = var.b_2n_u / 2;
    x3_Nut = var.b_2n_o / 2;
    y1_Nut = 0;
    y2_Nut = var.h_2ns;
    y3_Nut = var.h_2ns + tan(var.alpha_2nk) * ((var.b_2n_u-var.b_2ns) / 2);
    y4_Nut = var.h_2n;
    Nut = [x1_Nut,x1_Nut,x2_Nut,x3_Nut,-x3_Nut,-x2_Nut,-x1_Nut,-x1_Nut; y1_Nut,y2_Nut,y3_Nut,y4_Nut,y4_Nut,y3_Nut,y2_Nut,y1_Nut]';
    
    % Nutfuellung
    % Fuell-Koordinaten
    x1_Fuellung = x2_Nut - var.d_2iso;
    x2_Fuellung = x3_Nut - var.d_2iso;
    y1_Fuellung = var.h_2nk;
    y2_Fuellung = var.h_2n - var.d_2iso;
    Fuellung = [x1_Fuellung,x2_Fuellung,-x2_Fuellung,-x1_Fuellung; y1_Fuellung,y2_Fuellung,y2_Fuellung,y1_Fuellung]';
    
    % Plot
    plot(handles.axes_Ergebnis_Plot,Nut(:,1),Nut(:,2),'k')
    fill(handles.axes_Ergebnis_Plot,Fuellung(:,1),Fuellung(:,2),[100 160 200]./255)

    % Anpassung figure
    if(y4_Nut>=var.b_2n_o)
        ylim(handles.axes_Ergebnis_Plot,[-0.5,y4_Nut+0.5]);
        xlim(handles.axes_Ergebnis_Plot,[-y4_Nut/2-0.5,y4_Nut/2+0.5]);
    else
        ylim(handles.axes_Ergebnis_Plot,[-0.5,x3_Nut*2+0.5]);
        xlim(handles.axes_Ergebnis_Plot,[-x3_Nut-0.5,x3_Nut+0.5]);
    end
    axis(handles.axes_Ergebnis_Plot,'equal') 
    grid(handles.axes_Ergebnis_Plot,'on')
    yticks(handles.axes_Ergebnis_Plot,0:1:y4_Nut+0.5)
    xticks(handles.axes_Ergebnis_Plot,fix(min(handles.axes_Ergebnis_Plot.XLim)):1:fix(max(handles.axes_Ergebnis_Plot.XLim)))
    ylabel(handles.axes_Ergebnis_Plot,'y in mm','interpreter','latex','FontSize', 15)
    xlabel(handles.axes_Ergebnis_Plot,'x in mm','interpreter','latex','FontSize', 15)
    title(handles.axes_Ergebnis_Plot,'Abmessungen Rotornut','interpreter','latex','FontSize', 20)
    hold(handles.axes_Ergebnis_Plot, 'off')
    
    axtoolbar(handles.axes_Ergebnis_Plot,{'pan','zoomin','zoomout','restoreview'});
    handles.axes_Ergebnis_Plot.Toolbar.Visible = 'on';

% Plot Generic Wicklungslayout
function plotWicklungslayout(handles, wick, rated)
    hold(handles.axes_Auswahl_Wicklung,'on')
    box(handles.axes_Auswahl_Wicklung,'on')

    D_1a = 1.3 * 2;
    D_1i = 0.9 * 2;
    b_1ns = 0.03;
    h_1ns = 0.02;
    if(wick.N_1>42)
        b_1n_u = 0.07 * (42/wick.N_1);
        b_1n_o = 0.1 * (42/wick.N_1);
    else
        b_1n_u = 0.07;
        b_1n_o = 0.1;
    end
    alpha_1nk = pi/6;
    h_1n = 0.3;

     % Stator aussen
    theta = linspace(0,2*pi,4e2);
    Stator_aussen_x = (D_1a)/2 * cos(theta) + 0;
    Stator_aussen_y = (D_1a)/2 * sin(theta) + 0;

    Stator_aussen_x = Stator_aussen_x';
    Stator_aussen_y = Stator_aussen_y';
    plot(handles.axes_Auswahl_Wicklung,Stator_aussen_x,Stator_aussen_y,'k','HandleVisibility','off')

    % Stator innen
    % Nut
    Segmenthoehe = ((D_1i)/2) * (1 - sqrt(1 - (b_1ns/(D_1i))^2));
    x1_Nut = b_1ns / 2;
    x2_Nut = b_1n_u / 2;
    x3_Nut = b_1n_o / 2;
    y1_Nut = (D_1i)/2 - Segmenthoehe;
    y2_Nut = h_1ns + (D_1i)/2 - Segmenthoehe;
    y3_Nut = h_1ns + tan(alpha_1nk) * ((b_1n_u-b_1ns) / 2) + (D_1i)/2 - Segmenthoehe;
    y4_Nut = h_1n + (D_1i)/2 - Segmenthoehe;
    Nut = [x1_Nut,x1_Nut,x2_Nut,x3_Nut,-x3_Nut,-x2_Nut,-x1_Nut,-x1_Nut; y1_Nut,y2_Nut,y3_Nut,y4_Nut,y4_Nut,y3_Nut,y2_Nut,y1_Nut]';

    % Kreissegment
    theta = linspace((asin(b_1ns/(D_1i))+(pi/2)),(((2*pi)/wick.N_1)-(asin(b_1ns/(D_1i)))+(pi/2)),1e1);
    Kreissegment = [(D_1i)/2 * cos(theta); (D_1i)/2 * sin(theta)]';

    % Koordinatentransformation
    Stator_innen_x = [];
    Stator_innen_y = [];
    Fuellung_Stator_x = [];
    Fuellung_Stator_y = [];
    for i = 1:wick.N_1
        theta = -((2*(i-1))/wick.N_1)*pi;
        M = [cos(theta) -sin(theta); sin(theta) cos(theta)];

        % Nut
        Nut_trans = Nut * M;

        % Kreissegment
        Kreissegment_trans = Kreissegment * M;

        % gesamt
        Stator_innen_x = [Stator_innen_x; Nut_trans(2:end-1,1); Kreissegment_trans(:,1)];
        Stator_innen_y = [Stator_innen_y; Nut_trans(2:end-1,2); Kreissegment_trans(:,2)];
    end
    Stator_innen_x(end+1,1) = Stator_innen_x(1,1);
    Stator_innen_y(end+1,1) = Stator_innen_y(1,1);

    plot(handles.axes_Auswahl_Wicklung,Stator_innen_x,Stator_innen_y,'k','HandleVisibility','off')

    % Wicklung
    coord_2 = [0 1.1; 0 1];
    coord_text = [0 0.85];
    theta = linspace(0,2*pi,wick.N_1+1);
    theta(end) = [];
    color_3 = [{'r'},{'b'},{'g'}];
    name_3 = [{'A'},{'C'},{'B'}; {'a'},{'c'},{'b'}];
    reset_minus = ones(rated.m,1);
    reset_plus = ones(rated.m,1);
    for k = 1:wick.n_lay
        coord = coord_2(k,:);
        for j = 1:rated.m
            if(k==1)
                M6 = wick.Matrix(:,1:length(wick.Matrix_lay{j,k}));
            elseif(k==2)
                M6 = wick.Matrix(:,length(wick.Matrix_lay{j,k-1})+1:end);
            else
                error('')
            end
            M7 = sort(M6,2,'ComparisonMethod','abs');
            vec = M7(j,:);
            color = color_3{j};
            for i = 1:length(vec)
                theta_1 = theta(abs(vec(i)));
                M = [cos(theta_1) -sin(theta_1); sin(theta_1) cos(theta_1);];
                coord_trans = coord * M;
                coord_text_trans = coord_text * M;
                text(handles.axes_Auswahl_Wicklung,coord_text_trans(1),coord_text_trans(2),num2str(abs(vec(i))))
                if(sign(vec(i)<1))
                    if(reset_minus(j))
                        plot(handles.axes_Auswahl_Wicklung,coord_trans(1,1),coord_trans(1,2),'Color',color,'Marker','x','DisplayName',name_3{2,j})
                    else
                        plot(handles.axes_Auswahl_Wicklung,coord_trans(1,1),coord_trans(1,2),'Color',color,'Marker','x','DisplayName',name_3{2,j},'HandleVisibility','off')
                    end
                    reset_minus(j) = 0;
                else
                    if(reset_plus(j))
                        plot(handles.axes_Auswahl_Wicklung,coord_trans(1,1),coord_trans(1,2),'Color',color,'Marker','o','DisplayName',name_3{1,j})
                    else
                        plot(handles.axes_Auswahl_Wicklung,coord_trans(1,1),coord_trans(1,2),'Color',color,'Marker','o','DisplayName',name_3{1,j},'HandleVisibility','off')
                    end
                    reset_plus(j) = 0;
                end
            end
        end
    end

    legend(handles.axes_Auswahl_Wicklung,'Location', 'eastoutside');
    
    axis(handles.axes_Auswahl_Wicklung,'equal')
    axis(handles.axes_Auswahl_Wicklung,'off')
    grid(handles.axes_Auswahl_Wicklung,'off')
    title(handles.axes_Auswahl_Wicklung,'Wicklungslayout Maschine','interpreter','latex','FontSize', 20)
    hold(handles.axes_Auswahl_Wicklung, 'off')
    
% Save Geometrie in DXF file
function saveDXF(handles,var)
    % Open DXF file
    FID = dxf_open(['3_Ergebnisse/',handles.opt.folder_id,'/1_Entwurf/Geometrie_',handles.opt.file_id,'.dxf']);

    % Produce Stator aussen
    dxf_polyline(FID,var.Stator_aussen_x,var.Stator_aussen_y,zeros(length(var.Stator_aussen_x),1));

    % Produce Stator innen
    dxf_polyline(FID,var.Stator_innen_x,var.Stator_innen_y,zeros(length(var.Stator_innen_x),1));

    % Produce Rotor aussen
    dxf_polyline(FID,var.Rotor_aussen_x,var.Rotor_aussen_y,zeros(length(var.Rotor_aussen_x),1));

    % Produce Rotor innen
    dxf_polyline(FID,var.Rotor_innen_x,var.Rotor_innen_y,zeros(length(var.Rotor_innen_x),1));
    
    % Produce Stator Nutfuellung
    for i = 1:handles.Entwurf.Wicklung.N_1
        dxf_polyline(FID,var.Fuellung_Stator_x((length(var.Fuellung_Stator_x)/handles.Entwurf.Wicklung.N_1)*(i-1)+1:(length(var.Fuellung_Stator_x)/handles.Entwurf.Wicklung.N_1)*i,1),var.Fuellung_Stator_y((length(var.Fuellung_Stator_y)/handles.Entwurf.Wicklung.N_1)*(i-1)+1:(length(var.Fuellung_Stator_y)/handles.Entwurf.Wicklung.N_1)*i,1),zeros((length(var.Fuellung_Stator_y)/handles.Entwurf.Wicklung.N_1),1));
    end
    
    if(strcmp(handles.opt.Maschinentyp,'ASM'))
        % Produce Rotor Nutfuellung
        for i = 1:handles.Entwurf.Wicklung.N_2
            dxf_polyline(FID,var.Fuellung_Rotor_x((length(var.Fuellung_Rotor_x)/handles.Entwurf.Wicklung.N_2)*(i-1)+1:(length(var.Fuellung_Rotor_x)/handles.Entwurf.Wicklung.N_2)*i,1),var.Fuellung_Rotor_y((length(var.Fuellung_Rotor_y)/handles.Entwurf.Wicklung.N_2)*(i-1)+1:(length(var.Fuellung_Rotor_y)/handles.Entwurf.Wicklung.N_2)*i,1),zeros((length(var.Fuellung_Rotor_y)/handles.Entwurf.Wicklung.N_2),1));
        end
    elseif(strcmp(handles.opt.Maschinentyp,'PMSM'))
        % Produce Magnets
        if(strcmp(handles.opt.Maschinenausfuehrung,'SPMSM'))
            for i = 1:(2*handles.Entwurf.Bemessungswerte.p)
                dxf_polyline(FID,var.Magnet_Rotor_x(var.points_Magnet(i):var.points_Magnet(i+1)-1,1),var.Magnet_Rotor_y(var.points_Magnet(i):var.points_Magnet(i+1)-1,1),zeros(length(var.points_Magnet(i):var.points_Magnet(i+1)-1),1));
            end
        elseif(strcmp(handles.opt.Maschinenausfuehrung,'IPMSM (eingelassen)'))
            for i = 1:(2*handles.Entwurf.Bemessungswerte.p)
                dxf_polyline(FID,var.Magnet_Rotor_x(var.points_Magnet(i):var.points_Magnet(i+1)-1,1),var.Magnet_Rotor_y(var.points_Magnet(i):var.points_Magnet(i+1)-1,1),zeros(length(var.points_Magnet(i):var.points_Magnet(i+1)-1),1));
            end
        elseif(strcmp(handles.opt.Maschinenausfuehrung,'IPMSM (tangential)'))
            for i = 1:(2*handles.Entwurf.Bemessungswerte.p)
                dxf_polyline(FID,var.Magnet_Rotor_x((length(var.Magnet_Rotor_x)/(2*handles.Entwurf.Bemessungswerte.p))*(i-1)+1:(length(var.Magnet_Rotor_x)/(2*handles.Entwurf.Bemessungswerte.p))*i,1),var.Magnet_Rotor_y((length(var.Magnet_Rotor_y)/(2*handles.Entwurf.Bemessungswerte.p))*(i-1)+1:(length(var.Magnet_Rotor_y)/(2*handles.Entwurf.Bemessungswerte.p))*i,1),zeros((length(var.Magnet_Rotor_y)/(2*handles.Entwurf.Bemessungswerte.p)),1));
            end    
        elseif(strcmp(handles.opt.Maschinenausfuehrung,'IPMSM (V-Form)'))
            for i = 1:(4*handles.Entwurf.Bemessungswerte.p)
                dxf_polyline(FID,var.Magnet_Rotor_x((length(var.Magnet_Rotor_x)/(4*handles.Entwurf.Bemessungswerte.p))*(i-1)+1:(length(var.Magnet_Rotor_x)/(4*handles.Entwurf.Bemessungswerte.p))*i,1),var.Magnet_Rotor_y((length(var.Magnet_Rotor_y)/(4*handles.Entwurf.Bemessungswerte.p))*(i-1)+1:(length(var.Magnet_Rotor_y)/(4*handles.Entwurf.Bemessungswerte.p))*i,1),zeros((length(var.Magnet_Rotor_y)/(4*handles.Entwurf.Bemessungswerte.p)),1));
            end
        else
            error('Ungueltige Eingabe bei Variable "opt.Maschinenausfuehrung"')
        end
    else
        error('Ungueltige Eingabe bei Variable "handles.opt.Maschinentyp"');
    end
    
    % Close DXF file
    dxf_close(FID);

% Export Excel
function export_table = export_excel(input)
    var = fieldnames(input);
    
    for i = 1:length(var)
        export_struct{i,1} = input.(var{i});
    end
    
    export_table = table(export_struct,'VariableNames',{'var'},'RowNames',var);
    
% Save Excel
function saveExcel(handles, Entwurf)
    copyfile('3_Ergebnisse/1_Misc/Vorlage_Ergebnisse.xlsx',['3_Ergebnisse/',handles.opt.folder_id,'/1_Entwurf/Entwurf_',handles.opt.file_id,'.xlsx']);
    
    rated_export = export_excel(Entwurf.Bemessungswerte);
    richt_export = export_excel(Entwurf.Richtwerte);
    
    Entwurf.Optionen = rmfield(Entwurf.Optionen,{'Content','Locked','Saved','file_id','folder_id','unique_id'});
    var1 = fieldnames(Entwurf.Optionen);
    for i = 1:length(var1)
        if(isstruct(Entwurf.Optionen.(var1{i})))
            Entwurf.Optionen.(var1{i}) = Entwurf.Optionen.(var1{i}).String;
        end
    end
    opt_export = export_excel(Entwurf.Optionen);
    
    range_oben = 2;
    range_unten = range_oben + height(rated_export) - 1;
    range = ['B' num2str(range_oben) ':C' num2str(range_unten)];
    writetable(rated_export,['3_Ergebnisse/',handles.opt.folder_id,'/1_Entwurf/Entwurf_',handles.opt.file_id,'.xlsx'],'Sheet','Bemessungswerte','Range',range,'WriteVariableNames',0,'WriteRowNames',1)
    
    range_oben = 2;
    range_unten = range_oben + height(richt_export) - 1;
    range = ['B' num2str(range_oben) ':C' num2str(range_unten)];
    writetable(richt_export,['3_Ergebnisse/',handles.opt.folder_id,'/1_Entwurf/Entwurf_',handles.opt.file_id,'.xlsx'],'Sheet','Richtwerte','Range',range,'WriteVariableNames',0,'WriteRowNames',1)
    
    range_oben = 2;
    range_unten = range_oben + height(opt_export) - 1;
    range = ['B' num2str(range_oben) ':C' num2str(range_unten)];
    writetable(opt_export,['3_Ergebnisse/',handles.opt.folder_id,'/1_Entwurf/Entwurf_',handles.opt.file_id,'.xlsx'],'Sheet','Optionen','Range',range,'WriteVariableNames',0,'WriteRowNames',1)
    
    
    if(handles.opt.Locked)
        if(~strcmp(Entwurf.Optionen.Mode_Wicklung,'Manuell'))
            Entwurf.Wicklung = rmfield(Entwurf.Wicklung,{'Matrix','Matrix_lay','table_1all','table_1pos','table_1sel','xi_1sp','xi_1gr','xi_1n','xi_1spgr','xi_1','nu_1'});
        else
            Entwurf.Wicklung = rmfield(Entwurf.Wicklung,{'Matrix','Matrix_lay','table_1all','table_1pos','table_1sel','xi_1n','xi_1spgr','xi_1','nu_1'});
        end
        wick_export = export_excel(Entwurf.Wicklung);
        
        Entwurf.Geometrie.Nut_1 = rmfield(Entwurf.Geometrie.Nut_1,{'iter','maxIter'});
        geo_Nut_1_export = export_excel(Entwurf.Geometrie.Nut_1);
        
        if(strcmp(handles.opt.Maschinentyp,'ASM'))
            Entwurf.Geometrie.Nut_2 = rmfield(Entwurf.Geometrie.Nut_2,{'iter','maxIter'});
            geo_Nut_2_export = export_excel(Entwurf.Geometrie.Nut_2);
            Entwurf.Geometrie = rmfield(Entwurf.Geometrie,{'Nut_1','Nut_2'});
        elseif(strcmp(handles.opt.Maschinentyp,'PMSM'))
            Entwurf.Geometrie = rmfield(Entwurf.Geometrie,{'Nut_1'});
        else
            error('Ungueltige Eingabe bei Variable "Entwurf.Optionen.Maschinentyp"');
        end
        
        var2 = [fieldnames(Entwurf.Geometrie)' fieldnames(Entwurf.Geometrie.misc)'; struct2cell(Entwurf.Geometrie)' struct2cell(Entwurf.Geometrie.misc)'];
        Entwurf.Geometrie = struct(var2{:});
        Entwurf.Geometrie = orderfields(Entwurf.Geometrie);
        Entwurf.Geometrie = rmfield(Entwurf.Geometrie,{'misc'});
        var1 = fieldnames(Entwurf.Geometrie);
        for i = 1:length(var1)
            if(length(Entwurf.Geometrie.(var1{i}))>1)
                Entwurf.Geometrie = rmfield(Entwurf.Geometrie,var1(i));
            end
        end
        geo_export = export_excel(Entwurf.Geometrie);
        
        var2 = [fieldnames(Entwurf.EMAG)' fieldnames(Entwurf.EMAG.misc)'; struct2cell(Entwurf.EMAG)' struct2cell(Entwurf.EMAG.misc)'];
        Entwurf.EMAG = struct(var2{:});
        Entwurf.EMAG = orderfields(Entwurf.EMAG);
        if(strcmp(handles.opt.Maschinentyp,'ASM'))
        elseif(strcmp(handles.opt.Maschinentyp,'PMSM'))
            Entwurf.EMAG = rmfield(Entwurf.EMAG,{'misc'});
        else
            error('Ungueltige Eingabe bei Variable "Entwurf.Optionen.Maschinentyp"');
        end
        emag_export = export_excel(Entwurf.EMAG);
        
        range_oben = 2;
        range_unten = range_oben + height(wick_export) - 1;
        range = ['B' num2str(range_oben) ':C' num2str(range_unten)];
        writetable(wick_export,['3_Ergebnisse/',handles.opt.folder_id,'/1_Entwurf/Entwurf_',handles.opt.file_id,'.xlsx'],'Sheet','Wicklung','Range',range,'WriteVariableNames',0,'WriteRowNames',1)

        range_oben = 2;
        range_unten = range_oben + height(geo_export) - 1;
        range = ['B' num2str(range_oben) ':C' num2str(range_unten)];
        writetable(geo_export,['3_Ergebnisse/',handles.opt.folder_id,'/1_Entwurf/Entwurf_',handles.opt.file_id,'.xlsx'],'Sheet','Geometrie','Range',range,'WriteVariableNames',0,'WriteRowNames',1)

        range_oben = range_unten + 2;
        range_unten = range_oben + height(geo_Nut_1_export) - 1;
        range = ['B' num2str(range_oben) ':C' num2str(range_unten)];
        writetable(geo_Nut_1_export,['3_Ergebnisse/',handles.opt.folder_id,'/1_Entwurf/Entwurf_',handles.opt.file_id,'.xlsx'],'Sheet','Geometrie','Range',range,'WriteVariableNames',0,'WriteRowNames',1)

        if(strcmp(handles.opt.Maschinentyp,'ASM'))
            range_oben = range_unten + 2;
            range_unten = range_oben + height(geo_Nut_2_export) - 1;
            range = ['B' num2str(range_oben) ':C' num2str(range_unten)];
            writetable(geo_Nut_2_export,['3_Ergebnisse/',handles.opt.folder_id,'/1_Entwurf/Entwurf_',handles.opt.file_id,'.xlsx'],'Sheet','Geometrie','Range',range,'WriteVariableNames',0,'WriteRowNames',1) 
        elseif(strcmp(handles.opt.Maschinentyp,'PMSM'))
        else
            error('Ungueltige Eingabe bei Variable "Entwurf.Optionen.Maschinentyp"');
        end
        
        range_oben = 2;
        range_unten = range_oben + height(emag_export) - 1;
        range = ['B' num2str(range_oben) ':C' num2str(range_unten)];
        writetable(emag_export,['3_Ergebnisse/',handles.opt.folder_id,'/1_Entwurf/Entwurf_',handles.opt.file_id,'.xlsx'],'Sheet','EMAG','Range',range,'WriteVariableNames',0,'WriteRowNames',1)
    end

% Tingley Algorithmus
function wick = TingleyAlg(rated, wick)

    % Nuten auf die Straenge verteilen
    M1 = zeros(1,wick.q_1n*wick.N_1);
    for i = 1:wick.N_1
        M1(1,(wick.q_1n*(i-1))+1) = i;
    end

    % Anzahl Zeilen
    n_z = 2*rated.p;

    % Anzahl Spalten
    n_s = wick.q_1n * wick.N_1/(2*rated.p);

    % Tingley-Matrix formen
    M2 = reshape(M1,n_s,n_z)';

    % Vorzeichen anpassen
    vz_odd = ones(n_z,1);
    vz_odd(2:2:end,1) = -1;
    vz_even = ones(n_z,1);
    vz_even(1:2:end,1) = -1;
    for i = 1:rated.m
        if(mod(i,2)) %odd
            M3(:,(i-1)*(n_s/rated.m)+1:i*(n_s/rated.m)) = M2(:,(i-1)*(n_s/rated.m)+1:i*(n_s/rated.m)) .* vz_odd;
        else %even
            M3(:,(i-1)*(n_s/rated.m)+1:i*(n_s/rated.m)) = M2(:,(i-1)*(n_s/rated.m)+1:i*(n_s/rated.m)) .* vz_even;
        end
    end

    % Strang zusammenfuegen und Null entfernen 
    M4 = M3(M3~=0);
    M4 = reshape(M4,[],rated.m);
    M5 = num2cell(M4,1)';
    
    % Oberschicht hinzufuegen
    switch wick.Wicklungstyp_1
        case {'1SGL', '1SBL'}
            M7 = M4';
            wick.Matrix_lay = M5;
        case {'2SGL', '2SBL'}
            for i = 1:rated.m
                M6(:,i) = sign(M4(:,i)).*(abs(M4(:,i))+wick.y_1).*(-1);
                M7(i,:) = [M4(:,i); M6(:,i)];
            end
            M8 = num2cell(M6,1)';
            M9 = [M5 M8];
            wick.Matrix_lay = M9;
    end
    M7(abs(M7)>wick.N_1) = M7(abs(M7)>wick.N_1)-sign(M7(abs(M7)>wick.N_1))*wick.N_1;
    wick.Matrix = M7;
    
    clear M1 M2 M3 M4 M5 n_z n_s vz_odd vz_even i
    
% Berechnung x-y-Koordinaten fuer Plot und DXF-Export
function geo = Update_x_y_Koordinaten(rated, geo, wick, opt)
    % Stator aussen
    theta = linspace(0,2*pi,4e2);
    geo.Stator_aussen_x = (geo.D_1a*1e3)/2 * cos(theta) + 0;
    geo.Stator_aussen_y = (geo.D_1a*1e3)/2 * sin(theta) + 0;

    geo.Stator_aussen_x = geo.Stator_aussen_x';
    geo.Stator_aussen_y = geo.Stator_aussen_y';
    
    % Stator innen
    % Nut
    if(strcmp(opt.Nutform_Stator,'Trapezform (eckig)'))
        Segmenthoehe = ((geo.D_1i*1e3)/2) * (1 - sqrt(1 - (geo.Nut_1.b_1ns/(geo.D_1i*1e3))^2));
        x1_Nut = geo.Nut_1.b_1ns / 2;
        x2_Nut = geo.Nut_1.b_1n_u / 2;
        x3_Nut = geo.Nut_1.b_1n_o / 2;
        y1_Nut = (geo.D_1i*1e3)/2 - Segmenthoehe;
        y2_Nut = geo.Nut_1.h_1ns + (geo.D_1i*1e3)/2 - Segmenthoehe;
        y3_Nut = geo.Nut_1.h_1ns + tan(geo.Nut_1.alpha_1nk) * ((geo.Nut_1.b_1n_u-geo.Nut_1.b_1ns) / 2) + (geo.D_1i*1e3)/2 - Segmenthoehe;
        y4_Nut = geo.Nut_1.h_1n + (geo.D_1i*1e3)/2 - Segmenthoehe;
        Nut = [x1_Nut,x1_Nut,x2_Nut,x3_Nut,-x3_Nut,-x2_Nut,-x1_Nut,-x1_Nut; y1_Nut,y2_Nut,y3_Nut,y4_Nut,y4_Nut,y3_Nut,y2_Nut,y1_Nut]';

        % Nutfuellung
        x1_Fuellung = x2_Nut - geo.Nut_1.d_1iso;
        x2_Fuellung = x3_Nut - geo.Nut_1.d_1iso;
        y1_Fuellung = geo.Nut_1.h_1nk + (geo.D_1i*1e3)/2 - Segmenthoehe;
        y2_Fuellung = geo.Nut_1.h_1n - geo.Nut_1.d_1iso + (geo.D_1i*1e3)/2 - Segmenthoehe;
        Fuellung = [x1_Fuellung,x2_Fuellung,-x2_Fuellung,-x1_Fuellung; y1_Fuellung,y2_Fuellung,y2_Fuellung,y1_Fuellung]';
    else
        error('Ungueltige Eingabe bei Variable "opt.Nutform_Stator"');
    end
    
    % Kreissegment
    theta = linspace((asin(geo.Nut_1.b_1ns/(geo.D_1i*1e3))+(pi/2)),(((2*pi)/wick.N_1)-(asin(geo.Nut_1.b_1ns/(geo.D_1i*1e3)))+(pi/2)),1e1);
    Kreissegment = [(geo.D_1i*1e3)/2 * cos(theta); (geo.D_1i*1e3)/2 * sin(theta)]';
    
    % Koordinatentransformation
    geo.Stator_innen_x = [];
    geo.Stator_innen_y = [];
    geo.Fuellung_Stator_x = [];
    geo.Fuellung_Stator_y = [];
    for i = 1:wick.N_1
        theta = -((2*(i-1))/wick.N_1)*pi;
        M = [cos(theta) -sin(theta); sin(theta) cos(theta)];
        
        % Nut
        Nut_trans = Nut * M;
        
        % Nutfuellung
        Fuellung_trans = Fuellung * M;
        
        % Kreissegment
        Kreissegment_trans = Kreissegment * M;

        % gesamt
        geo.Stator_innen_x = [geo.Stator_innen_x; Nut_trans(2:end-1,1); Kreissegment_trans(:,1)];
        geo.Stator_innen_y = [geo.Stator_innen_y; Nut_trans(2:end-1,2); Kreissegment_trans(:,2)];
        geo.Fuellung_Stator_x = [geo.Fuellung_Stator_x; Fuellung_trans(:,1)];
        geo.Fuellung_Stator_y = [geo.Fuellung_Stator_y; Fuellung_trans(:,2)];
        geo.Fuellung_Stator_x(end+1,1) = Fuellung_trans(1,1);
        geo.Fuellung_Stator_y(end+1,1) = Fuellung_trans(1,2);
    end
    geo.Stator_innen_x(end+1,1) = geo.Stator_innen_x(1,1);
    geo.Stator_innen_y(end+1,1) = geo.Stator_innen_y(1,1);
    
    % Magnet
    if(strcmp(opt.Maschinenausfuehrung,'SPMSM'))
        % Rotor aussen
        theta = [0:((2*geo.b_PM*1e-3)/(geo.D_2a + 2*geo.h_PM*1e-3))/1e2:2*pi, 0];%linspace(0,2*pi,1e3);
        geo.Rotor_aussen_x = (geo.D_2a*1e3)/2 * cos(theta) + 0;
        geo.Rotor_aussen_y = (geo.D_2a*1e3)/2 * sin(theta) + 0;

        geo.Rotor_aussen_x = geo.Rotor_aussen_x';
        geo.Rotor_aussen_y = geo.Rotor_aussen_y';
        
        theta = ((2*geo.b_PM*1e-3)/(geo.D_2a + 2*geo.h_PM*1e-3))/1e2:((2*geo.b_PM*1e-3)/(geo.D_2a + 2*geo.h_PM*1e-3))/1e2:2*pi;
        geo.Magnet_Rotor_x = [];
        geo.Magnet_Rotor_y = [];
        geo.points_Magnet = 1;
        for i = 1:(2*rated.p)
            % Kreissegment
            theta_1 = theta(1,(theta-(((2*geo.b_PM*1e-3)/(geo.D_2a + 2*geo.h_PM*1e-3))*i + (((2*geo.tau_2p)/geo.D_2a)-((2*geo.b_PM*1e-3)/(geo.D_2a + 2*geo.h_PM*1e-3)))*(i-1)))<=1e-2) + (pi/2);
            theta((theta-((2*geo.tau_2p)/geo.D_2a)*i)<=0) = [];
            Kreissegment_aussen = [((geo.D_2a*1e3)/2 + geo.h_PM) * cos(theta_1); ((geo.D_2a*1e3)/2 + geo.h_PM) * sin(theta_1)]';
            Kreissegment_innen = [((geo.D_2a*1e3)/2) * cos(theta_1); ((geo.D_2a*1e3)/2) * sin(theta_1)]';
            Kreissegment = [Kreissegment_innen(1,:); Kreissegment_aussen; flipud(Kreissegment_innen)];
            
            geo.points_Magnet(i+1) = geo.points_Magnet(i) + length(Kreissegment(:,1));
            % gesamt
            geo.Magnet_Rotor_x = [geo.Magnet_Rotor_x; Kreissegment(:,1)];
            geo.Magnet_Rotor_y = [geo.Magnet_Rotor_y; Kreissegment(:,2)];
        end
    elseif(strcmp(opt.Maschinenausfuehrung,'IPMSM (eingelassen)'))
        % Korrektur der Rueckenhoehe (Rotor) h_2r [mm]
        geo.h_2r = (((geo.D_2a)/2) * 1e3) - geo.h_PM;
        
        % Rotor aussen
        theta = [0:((2*geo.b_PM*1e-3)/geo.D_2a)/1e2:2*pi, 0];%linspace(0,2*pi,1e3);
        geo.Rotor_aussen_x = (geo.D_2a*1e3)/2 * cos(theta) + 0;
        geo.Rotor_aussen_y = (geo.D_2a*1e3)/2 * sin(theta) + 0;

        geo.Rotor_aussen_x = geo.Rotor_aussen_x';
        geo.Rotor_aussen_y = geo.Rotor_aussen_y';
        
        theta = ((2*geo.b_PM*1e-3)/geo.D_2a)/1e2:((2*geo.b_PM*1e-3)/geo.D_2a)/1e2:2*pi;
        geo.Magnet_Rotor_x = [];
        geo.Magnet_Rotor_y = [];
        geo.points_Magnet = 1;
        for i = 1:(2*rated.p)
            % Kreissegment
            theta_1 = theta(1,(theta-(((2*geo.b_PM*1e-3)/geo.D_2a)*i + (((2*geo.tau_2p)/geo.D_2a)-((2*geo.b_PM*1e-3)/geo.D_2a))*(i-1)))<=1e-2);
            theta((theta-((2*geo.tau_2p)/geo.D_2a)*i)<=0) = [];
            Kreissegment_innen = [((geo.D_2a*1e3)/2 - geo.h_PM) * cos(theta_1); ((geo.D_2a*1e3)/2 - geo.h_PM) * sin(theta_1)]';
            Kreissegment_aussen = [((geo.D_2a*1e3)/2) * cos(theta_1); ((geo.D_2a*1e3)/2) * sin(theta_1)]';
            Kreissegment = [Kreissegment_innen(1,:); Kreissegment_aussen; flipud(Kreissegment_innen)];
            
            geo.points_Magnet(i+1) = geo.points_Magnet(i) + length(Kreissegment(:,1));
            % gesamt
            geo.Magnet_Rotor_x = [geo.Magnet_Rotor_x; Kreissegment(:,1)];
            geo.Magnet_Rotor_y = [geo.Magnet_Rotor_y; Kreissegment(:,2)];
        end
    elseif(strcmp(opt.Maschinenausfuehrung,'IPMSM (tangential)'))
        % Hoehe der PM-Konstruktion [mm]
        geo.h_1_PM = ((geo.D_2a)/2 *1e3) - 0.5*sqrt((geo.D_2a*1e3)^2 - geo.b_PM^2);
        
        % Korrektur der Rueckenhoehe (Rotor) h_2r [mm]
        geo.h_2r = ((geo.D_2a)/2 *1e3) - (geo.h_1_PM + geo.Abstand_PM_Rotoroberflaeche + geo.h_PM);
        
        % Rotor aussen
        theta = linspace(0,2*pi,4e2);
        geo.Rotor_aussen_x = (geo.D_2a*1e3)/2 * cos(theta) + 0;
        geo.Rotor_aussen_y = (geo.D_2a*1e3)/2 * sin(theta) + 0;

        geo.Rotor_aussen_x = geo.Rotor_aussen_x';
        geo.Rotor_aussen_y = geo.Rotor_aussen_y';
        
        x1_Magnet = geo.b_PM/2;
        x2_Magnet = -geo.b_PM/2;
        y1_Magnet = geo.h_2r;
        y2_Magnet = geo.h_2r + geo.h_PM;
        Magnet = [x1_Magnet,x1_Magnet,x2_Magnet,x2_Magnet; y1_Magnet,y2_Magnet,y2_Magnet,y1_Magnet]';
        
        geo.Magnet_Rotor_x = [];
        geo.Magnet_Rotor_y = [];
        for i = 1:(2*rated.p)
            theta = -((2*(i-1))/(2*rated.p))*pi;
            M = [cos(theta) -sin(theta); sin(theta) cos(theta)];

            % Magnet
            Magnet_trans = Magnet * M;

            % gesamt
            geo.Magnet_Rotor_x = [geo.Magnet_Rotor_x; Magnet_trans(:,1)];
            geo.Magnet_Rotor_y = [geo.Magnet_Rotor_y; Magnet_trans(:,2)];
            geo.Magnet_Rotor_x(end+1,1) = Magnet_trans(1,1);
            geo.Magnet_Rotor_y(end+1,1) = Magnet_trans(1,2);
        end
        
    elseif(strcmp(opt.Maschinenausfuehrung,'IPMSM (V-Form)'))
        % Hoehe der PM-Konstruktion [mm]
        geo.h_1_PM = ((geo.D_2a)/2 *1e3) - 0.5*sqrt((geo.D_2a*1e3)^2 - ((geo.b_PM * sin(geo.alpha_PM/2)) + geo.Abstand_PM_unten)^2);
        geo.h_2_PM = (geo.b_PM/2) * cos(geo.alpha_PM/2);
        geo.h_3_PM = 0.5 * geo.h_PM * sqrt(4-cos(geo.alpha_PM/2)^2);
        
        % Korrektur der Rueckenhoehe (Rotor) h_2r [mm]
        geo.h_2r = ((geo.D_2a)/2 *1e3) - (geo.h_1_PM + geo.Abstand_PM_Rotoroberflaeche + geo.h_2_PM + geo.h_3_PM);
        
        % Rotor aussen
        theta = linspace(0,2*pi,4e2);
        geo.Rotor_aussen_x = (geo.D_2a*1e3)/2 * cos(theta) + 0;
        geo.Rotor_aussen_y = (geo.D_2a*1e3)/2 * sin(theta) + 0;

        geo.Rotor_aussen_x = geo.Rotor_aussen_x';
        geo.Rotor_aussen_y = geo.Rotor_aussen_y';
        
        x1_Magnet = 0;
        x2_Magnet = -geo.b_PM/2;
        y1_Magnet = 0;
        y2_Magnet = -geo.h_PM;
        Magnet = [x1_Magnet,x2_Magnet,x2_Magnet,x1_Magnet,x1_Magnet; y1_Magnet,y1_Magnet,y2_Magnet,y2_Magnet,y1_Magnet]';
        
        theta = pi/2 - geo.alpha_PM/2;
        M = [cos(theta) -sin(theta); sin(theta) cos(theta)];
        
        Magnet = Magnet * M;
        
        MovePM_y = ((geo.D_2a)/2 *1e3) - (geo.h_1_PM + geo.Abstand_PM_Rotoroberflaeche + geo.h_2_PM);
        MovePM_x = -(geo.Abstand_PM_unten/2);
        MovePM = [MovePM_x MovePM_x MovePM_x MovePM_x MovePM_x; MovePM_y MovePM_y MovePM_y MovePM_y MovePM_y]';
        Magnet = Magnet + MovePM;
        Magnet = [Magnet; [Magnet(:,1).*-1  Magnet(:,2)]];
        
        geo.Magnet_Rotor_x = [];
        geo.Magnet_Rotor_y = [];
        for i = 1:(2*rated.p)
            theta = -((2*(i-1))/(2*rated.p))*pi;
            M = [cos(theta) -sin(theta); sin(theta) cos(theta)];

            % Magnet
            Magnet_trans = Magnet * M;

            % gesamt
            geo.Magnet_Rotor_x = [geo.Magnet_Rotor_x; Magnet_trans(:,1)];
            geo.Magnet_Rotor_y = [geo.Magnet_Rotor_y; Magnet_trans(:,2)];
        end
    else
        error('Ungueltige Eingabe bei Variable "opt.Maschinenausfuehrung"')
    end

    % Rotor innen
    theta = linspace(0,2*pi,4e2);
    geo.Rotor_innen_x = (geo.D_2i*1e3)/2 * cos(theta) + 0;
    geo.Rotor_innen_y = (geo.D_2i*1e3)/2 * sin(theta) + 0;

    geo.Rotor_innen_x = geo.Rotor_innen_x';
    geo.Rotor_innen_y = geo.Rotor_innen_y';
   
% Font size optimization for Windows
function decreaseFontSizesIfReq(handles)
% make all fonts smaller on a non-mac-osx computer
persistent fontSizeDecreased
fontSizeDecreased = [];
if ~ismac()
    % No MAC OSX detected; decrease font sizes
    if isempty(fontSizeDecreased)
        for afield = fieldnames(handles)'
            afield = afield{1}; %#ok<FXSET>
            try %#ok<TRYNC>
                set(handles.(afield),'FontSize',get(handles.(afield),'FontSize')*0.75); % decrease font size
            end
        end
        fontSizeDecreased=1; % do not perform this step again.
    end
end
