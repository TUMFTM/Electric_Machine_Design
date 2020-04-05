% -------------------------------------------------------------------------
% TU Muenchen - Lehrstuhl fuer Fahrzeugtechnik (FTM)
% -------------------------------------------------------------------------
% Modell fuer den Entwurf und die Analyse einer PMSM oder ASM (MEAPA)
% -------------------------------------------------------------------------
% Autor: Svenja Kalt (kalt@ftm.mw.tum.de)
%        Jonathan Erhard
% -------------------------------------------------------------------------

function varargout = GUI_Analyse(varargin)
% GUI_ANALYSE MATLAB code for GUI_Analyse.fig
%      GUI_ANALYSE, by itself, creates a new GUI_ANALYSE or raises the existing
%      singleton*.
%
%      H = GUI_ANALYSE returns the handle to a new GUI_ANALYSE or the handle to
%      the existing singleton*.
%
%      GUI_ANALYSE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_ANALYSE.M with the given input arguments.
%
%      GUI_ANALYSE('Property','Value',...) creates a new GUI_ANALYSE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_Analyse_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_Analyse_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_Analyse

% Last Modified by GUIDE v2.5 13-Mar-2019 12:13:27

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_Analyse_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_Analyse_OutputFcn, ...
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

% --- Executes just before GUI_Analyse is made visible.
function GUI_Analyse_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_Analyse (see VARARGIN)

% Font size optimization for Windows
decreaseFontSizesIfReq(handles);

% Choose default command line output for GUI_Analyse
handles.output = hObject;

% Center GUI_Analyse
set(handles.figure1,'Position',[842 100 740 885]);

% Load Logos
FTM_Logo = imread('FTM_Logo.png');
axes(handles.axes_FTM_Logo);
imshow(FTM_Logo);

TUM_Logo = imread('TUM_Logo.jpeg');
axes(handles.axes_TUM_Logo);
imshow(TUM_Logo);

% Save input
handles.Entwurf = varargin{1};

% Check if there is any Analyse data
listing = dir(['3_Ergebnisse/',handles.Entwurf.Optionen.folder_id,'/2_Analyse']);
listing = listing(~[listing.isdir] & contains({listing.name},['Analyse_',handles.Entwurf.Optionen.file_id]));
[~,idx] = min(now*ones(size(listing)) - [listing.datenum]');
listing = listing(idx);

if(~isempty(listing))
    load([listing.folder '/' listing.name]);
    handles.Analyse = Analyse;
    handles.opt = handles.Analyse.Optionen;
    handles.default_opt = handles.Analyse.Optionen.default;
    
    handles = setOptionen(handles);
else
    handles.opt.Locked = 0;
    handles.opt.Count = 1;
    
    % Set default opt
    handles.default_opt.P_vw = 1;
    handles.default_opt.P_vu = 1;
    handles.default_opt.P_vu_Modell = 'Modellansatz Jordan';
    handles.default_opt.P_vme = 1;
    handles.default_opt.P_vzus = 1;
    handles.default_opt.Generator = 0;
    handles.default_opt.n_max = 10000;
    handles.default_opt.n_tics = 60;
    handles.default_opt.M_tics = 60;
    handles.default_opt.u_1max = handles.Entwurf.EMAG.U_1Str * sqrt(2);
    handles.default_opt.i_1max = handles.Entwurf.EMAG.I_1Str * sqrt(2);
    handles = setDefaultOptionen(handles);
    
    if(strcmp(handles.Entwurf.Optionen.Maschinentyp,'PMSM'))
        omega_k_max = (str2double(handles.edit_n_max.String)/60)*(2*pi)*handles.Entwurf.Bemessungswerte.p;
        i_f = handles.Entwurf.EMAG.Psi_PM / handles.Entwurf.EMAG.L_d;
        if(i_f >= str2double(handles.edit_i_1max.String))
            omega_k_max2 = sqrt(((str2double(handles.edit_u_1max.String))^2 - handles.Entwurf.EMAG.R_1^2*(str2double(handles.edit_i_1max.String))^2) / (handles.Entwurf.EMAG.Psi_PM - handles.Entwurf.EMAG.L_d*(str2double(handles.edit_i_1max.String)))^2) - 5;
            if(omega_k_max > omega_k_max2)
                handles.default_opt.n_max = (omega_k_max2 / ((2*pi)*handles.Entwurf.Bemessungswerte.p)) * 60;
            end
        end
    end

    handles = setDefaultOptionen(handles);

    % Read Optionen
    handles = readOptionen(handles);
end

% Set GUI
handles.text_ID.String = ['ID: Analyse_' handles.Entwurf.Optionen.file_id '_v' num2str(handles.opt.Count)];
idx = find(strcmp(handles.popupmenu_Maschinentyp.String, handles.Entwurf.Optionen.Maschinentyp));
handles.popupmenu_Maschinentyp.Value = idx;
handles = toggleLocked(handles);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI_Analyse wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = GUI_Analyse_OutputFcn(hObject, eventdata, handles) 
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

button = questdlg('Wollen Sie die Analyse beenden? Alle nicht gespeicherten Aenderungen gehen verloren.', 'Beenden','Ja','Nein','Nein'); 
switch button 
    case 'Ja'
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
        var = handles.opt.Count;
        if(handles.opt.Locked)
            if(isfield(handles,'Analyse'))
                handles = rmfield(handles,'Analyse');
            end
        end
        if(isfield(handles,'opt'))
            handles = rmfield(handles,'opt');
        end
        if(isfield(handles,'default_opt'))
            handles = rmfield(handles,'default_opt');
        end
        
        handles.opt.Locked = 0;
        handles.opt.Count = var + 1;

        % Set default opt
        handles.default_opt.P_vw = 1;
        handles.default_opt.P_vu = 1;
        handles.default_opt.P_vu_Modell = 'Modellansatz Jordan';
        handles.default_opt.P_vme = 1;
        handles.default_opt.P_vzus = 1;
        handles.default_opt.Generator = 0;
        handles.default_opt.n_max = 10000;
        handles.default_opt.n_tics = 60;
        handles.default_opt.M_tics = 60;
        handles.default_opt.u_1max = handles.Entwurf.EMAG.U_1Str * sqrt(2);
        handles.default_opt.i_1max = handles.Entwurf.EMAG.I_1Str * sqrt(2);
        handles = setDefaultOptionen(handles);
        
        if(strcmp(handles.Entwurf.Optionen.Maschinentyp,'PMSM'))
            omega_k_max = (str2double(handles.edit_n_max.String)/60)*(2*pi)*handles.Entwurf.Bemessungswerte.p;
            i_f = handles.Entwurf.EMAG.Psi_PM / handles.Entwurf.EMAG.L_d;
            if(i_f >= str2double(handles.edit_i_1max.String))
                omega_k_max2 = sqrt(((str2double(handles.edit_u_1max.String))^2 - handles.Entwurf.EMAG.R_1^2*(str2double(handles.edit_i_1max.String))^2) / (handles.Entwurf.EMAG.Psi_PM - handles.Entwurf.EMAG.L_d*(str2double(handles.edit_i_1max.String)))^2) - 5;
                if(omega_k_max > omega_k_max2)
                    handles.default_opt.n_max = (omega_k_max2 / ((2*pi)*handles.Entwurf.Bemessungswerte.p)) * 60;
                end
            end 
        end
        % Set default data
        handles = setDefaultOptionen(handles);
        
        handles = toggleLocked(handles);

        % Read data
        handles = readOptionen(handles);
        
        handles.text_ID.String = ['ID: Analyse_' handles.Entwurf.Optionen.file_id '_v' num2str(handles.opt.Count)];
    case 'Nein'
end

% Update handles structure
guidata(hObject, handles);

% --------------------------------------------------------------------
function toolbar_OpenFile_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to toolbar_OpenFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[file,path] = uigetfile('*.mat',['3_Ergebnisse/',handles.Entwurf.Optionen.folder_id,'/2_Analyse/']);

if isequal(file,0)
    % disp('Vorgang abgebrochen');
else
    %disp(['User selected ', fullfile(path,file)]);
    if(~strcmp(path,[pwd,'/3_Ergebnisse/',handles.Entwurf.Optionen.folder_id,'/2_Analyse/']))
        uiwait(msgbox({'Laden fehlgeschlagen. Datei muss aus dem Ordner des geoeffneten Entwurfs (ID) geladen werden.'},'Datei laden','help','error'))
        return;
    end
    
    if(handles.opt.Locked)
        if(isfield(handles,'Analyse'))
            handles = rmfield(handles,'Analyse');
        end
    end
    if(isfield(handles,'opt'))
        handles = rmfield(handles,'opt');
    end
    if(isfield(handles,'default_opt'))
        handles = rmfield(handles,'default_opt');
    end
    load([path file]);
    if(Analyse.Optionen.Locked)
        handles.Analyse = Analyse;
    end
    handles.opt = Analyse.Optionen;
    handles.default_opt = Analyse.Optionen.default;
    
    handles = toggleLocked(handles);
    
    % Set from file
    handles = setOptionen(handles);
    
    handles.text_ID.String = ['ID: Analyse_' handles.Entwurf.Optionen.file_id '_v' num2str(handles.opt.Count)];
end

% Update handles structure
guidata(hObject, handles);

% --------------------------------------------------------------------
function toolbar_SaveFile_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to toolbar_SaveFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if(isfile(['3_Ergebnisse/',handles.Entwurf.Optionen.folder_id,'/2_Analyse/Analyse_',handles.Entwurf.Optionen.file_id,'_v',num2str(handles.opt.Count),'.mat']))
    button = questdlg('Speichern als ... ?', 'Datei speichern','Neue Datei','Ueberschreiben','Neue Datei'); 
    switch button 
        case 'Neue Datei'
            handles.opt.Count = handles.opt.Count + 1;
            handles.text_ID.String = ['ID: Analyse_' handles.Entwurf.Optionen.file_id '_v' num2str(handles.opt.Count)];
        case 'Ueberschreiben'
    end
else
end

if(handles.opt.Locked)
    Analyse = handles.Analyse;
else
    Analyse.Optionen = handles.opt;
end
Analyse.Optionen.default = handles.default_opt;

% Save struct Analyse
try
    save(['3_Ergebnisse/',handles.Entwurf.Optionen.folder_id,'/2_Analyse/Analyse_',handles.Entwurf.Optionen.file_id,'_v',num2str(handles.opt.Count),'.mat'],'Analyse');
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
    button = questdlg('Wollen Sie die Datei entsperren? Es wird eine neue Analyse-Datei mit den bisherigen Eingabedaten erzeugt. Alle Ergebnisse der Analyse werden dabei zurueckgesetzt.', 'Datei entsperren','Ja','Nein','Nein'); 
    switch button 
        case 'Ja'
            if(handles.opt.Locked)
                if(isfield(handles,'Analyse'))
                    handles = rmfield(handles,'Analyse');
                end
            end

            handles.opt.Locked = 0;
            handles = toggleLocked(handles);
            handles.opt.Count = handles.opt.Count + 1;
            
            handles.text_ID.String = ['ID: Analyse_' handles.Entwurf.Optionen.file_id '_v' num2str(handles.opt.Count)];
        case 'Nein'
    end
end

handles = toggleLocked(handles);

% Update handles structure
guidata(hObject, handles);

% #########################################################################
% #   B) EINGABE                                                          #
% #########################################################################

% --- Executes on button press in pushbutton_Analyse.
function pushbutton_Analyse_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Analyse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.pushbutton_Analyse.Enable = 'off';
set(findall(handles.uipanel_Optionen_Verluste,'-property','Enable'),'Enable','off')
set(findall(handles.uipanel_Optionen_Berechnung,'-property','Enable'),'Enable','off')
handles.pushbutton_Reset.Enable = 'off';
handles.text_Wait.String = [{'Bitte warten ...'}; {'Die Berechnungen werden ausgefuehrt.'}];

pause(0.1)

if(strcmp(handles.Entwurf.Optionen.Maschinentyp,'ASM'))
    [handles.Analyse] = Analyse_ASM(handles);
elseif(strcmp(handles.Entwurf.Optionen.Maschinentyp,'PMSM'))
    [handles.Analyse] = Analyse_PMSM(handles);
else
    error('Ungueltige Eingabe bei Variable "handles.Entwurf.Optionen.Maschinentyp"');
end

% Berechnung erfolgreich
handles.opt.Locked = 1;
handles.text_Wait.String = [{'Success ...'}; {'Die Berechnungen waren erfolgreich.'}];
pause(1.0)

handles = toggleLocked(handles);

% Update handles structure
guidata(hObject, handles);

% #########################################################################
% #   B.1) OPTIONEN VERLUSTE                                              #
% #########################################################################

% --- Executes on button press in checkbox_P_vw.
function checkbox_P_vw_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_P_vw (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = readOptionen(handles);

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in checkbox_P_vu.
function checkbox_P_vu_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_P_vu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = readOptionen(handles);

% Update handles structure
guidata(hObject, handles);

% --- Executes on selection change in popupmenu_P_vu_Modell.
function popupmenu_P_vu_Modell_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_P_vu_Modell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = readOptionen(handles);

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in pushbutton_P_vu_Modell.
function pushbutton_P_vu_Modell_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_P_vu_Modell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

uiwait(msgbox({'Modellansatz Jordan: Berechnung der Ummagnetisierungsverluste ueber Induktionen und Elektroblechdatenblaetter mit der Erweiterung der Steinmetz-Gleichungen nach Jordan' '' 'Abschaetzung ueber verketteten Fluss: Berechnung der Ummagnetisierungsverluste ueber verketteten Fluss (nicht empfohlen, sehr ungenau)'},'Optionen','help','modal'))

% --- Executes on button press in checkbox_P_vme.
function checkbox_P_vme_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_P_vme (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = readOptionen(handles);

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in checkbox_P_vzus.
function checkbox_P_vzus_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_P_vzus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = readOptionen(handles);

% Update handles structure
guidata(hObject, handles);

% #########################################################################
% #   B.2) OPTIONEN BERECHNUNG                                            #
% #########################################################################

% --- Executes on button press in checkbox_Generator.
function checkbox_Generator_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_Generator (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = readOptionen(handles);

% Update handles structure
guidata(hObject, handles);

function edit_n_max_Callback(hObject, eventdata, handles)
% hObject    handle to edit_n_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = readOptionen(handles);

% Update handles structure
guidata(hObject, handles);

function edit_n_tics_Callback(hObject, eventdata, handles)
% hObject    handle to edit_n_tics (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = readOptionen(handles);

% Update handles structure
guidata(hObject, handles);

function edit_M_tics_Callback(hObject, eventdata, handles)
% hObject    handle to edit_M_tics (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = readOptionen(handles);

% Update handles structure
guidata(hObject, handles);

function edit_u_1max_Callback(hObject, eventdata, handles)
% hObject    handle to edit_u_1max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = readOptionen(handles);

% Update handles structure
guidata(hObject, handles);

function edit_i_1max_Callback(hObject, eventdata, handles)
% hObject    handle to edit_i_1max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = readOptionen(handles);

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in pushbutton_Reset.
function pushbutton_Reset_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = setDefaultOptionen(handles);
handles = readOptionen(handles);

% Update handles structure
guidata(hObject, handles);

% #########################################################################
% #   C) ERGEBNISSE                                                       #
% #########################################################################

% --- Executes on selection change in popupmenu_Ergebnis_Plot.
function popupmenu_Ergebnis_Plot_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_Ergebnis_Plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

contents = cellstr(get(hObject,'String'));
handles.Plot.Auswahl_Plot = contents{get(hObject,'Value')};

handles = plotErgebnisse(handles);

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in pushbutton_Save_Plot.
function pushbutton_Save_Plot_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Save_Plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

f_new = figure;
ax_new = copyobj(handles.axes_Ergebnis_Plot, f_new);
set(ax_new,'Position','default')
saveas(f_new, ['3_Ergebnisse/',handles.Entwurf.Optionen.folder_id,'/2_Analyse/',handles.Plot.Auswahl_Plot,'_',handles.Entwurf.Optionen.file_id,'_v',num2str(handles.opt.Count)], 'fig');
close(f_new);

% #########################################################################
% #   D) HILFSFUNKTIONEN                                                  #
% #########################################################################

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

% Toggle Locked
function handles = toggleLocked(handles)
    if(~handles.opt.Locked)
        handles.popupmenu_Ergebnis_Plot.Value = 1;
        handles.popupmenu_Ergebnis_Plot.String = {'-'};
        
        cla(handles.axes_Ergebnis_Plot,'reset');
        handles.axes_Ergebnis_Plot.Visible = 'off';
        handles.text_Wait.Visible = 'on';
        handles.text_Wait.String = 'Noch keine Berechnungen ausgefuehrt.';
        handles.popupmenu_Ergebnis_Plot.Enable = 'off';
        handles.pushbutton_Save_Plot.Enable = 'off';
        handles.pushbutton_Analyse.Enable = 'on';
        set(findall(handles.uipanel_Optionen_Verluste,'-property','Enable'),'Enable','on')
        set(findall(handles.uipanel_Optionen_Berechnung,'-property','Enable'),'Enable','on')
        handles.pushbutton_Reset.Enable = 'on';
        if(strcmp(handles.Entwurf.Optionen.Maschinentyp,'PMSM'))
        else
            handles.popupmenu_P_vu_Modell.Enable = 'off';
        end
    else
        handles.popupmenu_Ergebnis_Plot.Value = 1;
        var = {''};
        if(any(handles.opt.P_vw==1 | handles.opt.P_vu==1 |...
               handles.opt.P_vme==1 |  handles.opt.P_vzus==1))
           var = [var, 'Wirkungsgrad gesamt', 'Verluste gesamt'];
           if(handles.opt.P_vw==1)
                var = [var, 'Wirkungsgrad Wicklungsverluste','Wicklungsverluste'];
            end
            if(handles.opt.P_vu==1)
                var = [var, 'Wirkungsgrad Ummagnetisierungsverluste', 'Ummagnetisierungsverluste'];
            end
            if(handles.opt.P_vme==1)
                var = [var, 'Wirkungsgrad mechanische Verluste', 'mechanische Verluste'];
            end
            if(handles.opt.P_vzus==1)
                var = [var, 'Wirkungsgrad Zusatzverluste', 'Zusatzverluste'];
            end
            var = [var, 'i_1d', 'i_1q', 'u_1d', 'u_1q'];
        else
            var = [var,'i_1d', 'i_1q', 'u_1d', 'u_1q'];
        end
        var(1) = [];
        handles.popupmenu_Ergebnis_Plot.String = var;
        handles.Plot.Auswahl_Plot = handles.popupmenu_Ergebnis_Plot.String{handles.popupmenu_Ergebnis_Plot.Value};
        handles = plotErgebnisse(handles);
        
        handles.text_Wait.Visible = 'off';
        handles.text_Wait.String = [{'Success ...'}; {'Die Berechnungen waren erfolgreich.'}];
        handles.popupmenu_Ergebnis_Plot.Enable = 'on';
        handles.pushbutton_Save_Plot.Enable = 'on';
        handles.pushbutton_Analyse.Enable = 'off';
        set(findall(handles.uipanel_Optionen_Verluste,'-property','Enable'),'Enable','off')
        set(findall(handles.uipanel_Optionen_Berechnung,'-property','Enable'),'Enable','off')
        handles.pushbutton_Reset.Enable = 'off';
    end

% Set default Optionen
function handles = setDefaultOptionen(handles)
    % Wicklungsverluste
    handles.checkbox_P_vw.Value = handles.default_opt.P_vw;
    
    % Ummagnetisierungsverluste
    handles.checkbox_P_vu.Value = handles.default_opt.P_vu;
    
    % Eisenverlustmodell
    idx = find(strcmp(handles.popupmenu_P_vu_Modell.String,handles.default_opt.P_vu_Modell));
    handles.popupmenu_P_vu_Modell.Value = idx;
    
    % mechanische Verluste
    handles.checkbox_P_vme.Value = handles.default_opt.P_vme;
    
    % Zusatzverluste
    handles.checkbox_P_vzus.Value = handles.default_opt.P_vzus;
    
    % Generatorbereich
    handles.checkbox_Generator.Value = handles.default_opt.Generator;
    
    % max. Drehzahl
    handles.edit_n_max.String = num2str(handles.default_opt.n_max);
    
    % Aufloesung Drehzahl
    handles.edit_n_tics.String = num2str(handles.default_opt.n_tics);
    
    % Aufloesung Drehmoment
    handles.edit_M_tics.String = num2str(handles.default_opt.M_tics);
    
    % Momentensteuerung max. Spannung
    handles.edit_u_1max.String = num2str(handles.default_opt.u_1max);
    
    % Momentensteuerung max. Strom
    handles.edit_i_1max.String = num2str(handles.default_opt.i_1max);

% Set Optionen
function handles = setOptionen(handles)
    % Wicklungsverluste
    handles.checkbox_P_vw.Value = handles.opt.P_vw;
    
    % Ummagnetisierungsverluste
    handles.checkbox_P_vu.Value = handles.opt.P_vu;
    
    % Eisenverlustmodell
    idx = find(strcmp(handles.popupmenu_P_vu_Modell.String,handles.opt.P_vu_Modell));
    handles.popupmenu_P_vu_Modell.Value = idx;
    
    % mechanische Verluste
    handles.checkbox_P_vme.Value = handles.opt.P_vme;
    
    % Zusatzverluste
    handles.checkbox_P_vzus.Value = handles.opt.P_vzus;
    
    % Generatorbereich
    handles.checkbox_Generator.Value = handles.opt.Generator;
    
    % max. Drehzahl
    handles.edit_n_max.String = num2str(handles.opt.n_max);
    
    % Aufloesung Drehzahl
    handles.edit_n_tics.String = num2str(handles.opt.n_tics);
    
    % Aufloesung Drehmoment
    handles.edit_M_tics.String = num2str(handles.opt.M_tics);
    
    % Momentensteuerung max. Spannung
    handles.edit_u_1max.String = num2str(handles.opt.u_1max);
    
    % Momentensteuerung max. Strom
    handles.edit_i_1max.String = num2str(handles.opt.i_1max);

% Read Optionen
% function handles = readOptionen(handles)
%     % Wicklungsverluste
%     handles.opt.P_vw = handles.checkbox_P_vw.Value;
%     
%     % Ummagnetisierungsverluste
%     handles.opt.P_vu = handles.checkbox_P_vu.Value;
%     
%     % Eisenverlustmodell
%     handles.opt.P_vu_Modell = handles.popupmenu_P_vu_Modell.String{handles.popupmenu_P_vu_Modell.Value};
%     
%     % mechanische Verluste
%     handles.opt.P_vme = handles.checkbox_P_vme.Value;
%     
%     % Zusatzverluste
%     handles.opt.P_vzus = handles.checkbox_P_vzus.Value;
%     
%     % Generatorbereich
%     handles.opt.Generator = handles.checkbox_Generator.Value;
%     
%     % max. Drehzahl
%     handles.opt.n_max = str2double(handles.edit_n_max.String);
%     
%     % Aufloesung Drehzahl
%     handles.opt.n_tics = str2double(handles.edit_n_tics.String);
%     
%     % Aufloesung Drehmoment
%     handles.opt.M_tics = str2double(handles.edit_M_tics.String);
%     
%     % Momentensteuerung max. Spannung
%     handles.opt.u_1max = str2double(handles.edit_u_1max.String);
%     
%     % Momentensteuerung max. Strom
%     handles.opt.i_1max = str2double(handles.edit_i_1max.String);

% Read Optionen
function handles = readOptionen(handles)
    % Wicklungsverluste
    handles.opt.P_vw = handles.checkbox_P_vw.Value;
    
    % Ummagnetisierungsverluste
    handles.opt.P_vu = handles.checkbox_P_vu.Value;
    
    % Eisenverlustmodell
    handles.opt.P_vu_Modell = handles.popupmenu_P_vu_Modell.String{handles.popupmenu_P_vu_Modell.Value};
    
    % mechanische Verluste
    handles.opt.P_vme = handles.checkbox_P_vme.Value;
    
    % Zusatzverluste
    handles.opt.P_vzus = handles.checkbox_P_vzus.Value;
    
    % Generatorbereich
    handles.opt.Generator = handles.checkbox_Generator.Value;
    
    % max. Drehzahl
    handles.opt.n_max = str2double(handles.edit_n_max.String);
    
    % Aufloesung Drehzahl
    handles.opt.n_tics = str2double(handles.edit_n_tics.String);
    
    % Aufloesung Drehmoment
    handles.opt.M_tics = str2double(handles.edit_M_tics.String);
    
    % Momentensteuerung max. Spannung
    handles.opt.u_1max = str2double(handles.edit_u_1max.String);
    
    % Momentensteuerung max. Strom
    handles.opt.i_1max = str2double(handles.edit_i_1max.String);

% Plot Kennfeld
function handles = plotKennfeld(handles)
    hold(handles.axes_Ergebnis_Plot,'off');
    handles.axes_Ergebnis_Plot.Visible = 'on';
    plot(handles.axes_Ergebnis_Plot,handles.Analyse.Betriebsdaten.mot.n_m_vec,handles.Analyse.Betriebsdaten.mot.M_max_vec);
    hold(handles.axes_Ergebnis_Plot,'on');
    if(handles.opt.Generator)
        plot(handles.axes_Ergebnis_Plot,handles.Analyse.Betriebsdaten.gen.n_m_vec,handles.Analyse.Betriebsdaten.gen.M_max_vec);
    end
    colorbar(handles.axes_Ergebnis_Plot,'EastOutside');
    grid(handles.axes_Ergebnis_Plot,'on');
    handles.axes_Ergebnis_Plot.Layer = 'top';
    xlabel(handles.axes_Ergebnis_Plot,'Drehzahl n in $\frac{U}{min}$','interpreter','latex','FontSize', 15);
    ylabel(handles.axes_Ergebnis_Plot,'Drehmoment M in Nm','interpreter','latex','FontSize', 15);
    
% Plot Ergebnisse
function handles = plotErgebnisse(handles)
    if(strcmp(handles.Plot.Auswahl_Plot,'Wirkungsgrad gesamt'))
        handles = plotKennfeld(handles);
        title(handles.axes_Ergebnis_Plot,'Wirkungsgrad gesamt $\eta_{V,ges}$ in \%','interpreter','latex','FontSize', 20);
        Plot_Label_Vektor_eta = [0.7,0.75,0.8,0.84,0.86,0.88,0.9,0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.98,0.99,0.995];
        [C,h]=contourf(handles.axes_Ergebnis_Plot,handles.Analyse.Betriebsdaten.n_m_mesh,handles.Analyse.Betriebsdaten.M_max_mesh,handles.Analyse.Wirkungsgrad.eta_ges_mesh,Plot_Label_Vektor_eta);
        
        % Save Kennfeld neu
        % Save Kennfeld
        KennfeldMaschine = struct;
        KennfeldMaschineSpeichern = handles;
        KennfeldMaschine.n = KennfeldMaschineSpeichern.Analyse.Betriebsdaten.n_m_mesh;
        KennfeldMaschine.M = KennfeldMaschineSpeichern.Analyse.Betriebsdaten.M_max_mesh;
        KennfeldMaschine.etages = KennfeldMaschineSpeichern.Analyse.Wirkungsgrad.eta_ges_mesh;
        save (['3_Ergebnisse/',handles.Entwurf.Optionen.folder_id,'/','2_Analyse/', 'KennfeldMaschine_Gesamtwirkungsgrad','.mat'],'KennfeldMaschine');
        folder_id = [KennfeldMaschineSpeichern.Entwurf.Optionen.folder_id]; %file_id, '_data'];
        file_id = [KennfeldMaschineSpeichern.Entwurf.Optionen.file_id  ];
        
        % Konvertieren der Daten in richtige Form
        % Bemessungsgroessen
        Kennfeld_Gesamtwirkungsgrad=struct;
        Kennfeld_Gesamtwirkungsgrad.type = KennfeldMaschineSpeichern.Entwurf.Optionen.Maschinentyp;
        Kennfeld_Gesamtwirkungsgrad.power = KennfeldMaschineSpeichern.Entwurf.Bemessungswerte.P_N/1000;
        Kennfeld_Gesamtwirkungsgrad.n_n = KennfeldMaschineSpeichern.Entwurf.Bemessungswerte.n_N;
        Kennfeld_Gesamtwirkungsgrad.n_max = KennfeldMaschineSpeichern.Analyse.Optionen.n_max;
        Kennfeld_Gesamtwirkungsgrad.U_n = KennfeldMaschineSpeichern.Entwurf.Bemessungswerte.U_N;
        
        %Drehzahlachse - X-Achse Kennfeld
        Kennfeld_Gesamtwirkungsgrad.eff_n_axis = transpose(KennfeldMaschineSpeichern.Analyse.Betriebsdaten.mot.n_m_vec);

        %Drehmoment - Y-Achse
        KennfeldMaschineSpeichern.Analyse.Betriebsdaten.mot.M_vec=transpose(fliplr(KennfeldMaschineSpeichern.Analyse.Betriebsdaten.mot.M_vec));
        Kennfeld_Gesamtwirkungsgrad.eff_T_axis = [(fliplr(transpose(KennfeldMaschineSpeichern.Analyse.Betriebsdaten.mot.M_vec(2:end)))) * (-1), transpose(KennfeldMaschineSpeichern.Analyse.Betriebsdaten.mot.M_vec)];
        %Kennfeld_Gesamtwirkungsgrad.eff_T_axis = [fliplr(transpose(KennfeldMaschineSpeichern.Analyse.Betriebsdaten.mot.M_vec)),(transpose(KennfeldMaschineSpeichern.Analyse.Betriebsdaten.mot.M_vec(2:end))) * (-1)];

       %Diagrammbereich                                                             
        % Manipulate efficiency map
        KennfeldMaschineSpeichern.Analyse.Wirkungsgrad.eta_ges_mesh=flipud(KennfeldMaschineSpeichern.Analyse.Wirkungsgrad.eta_ges_mesh); %hier Unterschied Verlustarten
        KennfeldMaschineSpeichern.Analyse.Wirkungsgrad.eta_ges_mesh(KennfeldMaschineSpeichern.Analyse.Wirkungsgrad.eta_ges_mesh<0.1) = 0.1; %Hier Unterschied Verlustarten %set small efficiency values to 0.1
        KennfeldMaschineSpeichern.Analyse.Wirkungsgrad.eta_ges_mesh(1,:) = 0; %first colomn is zero %hier Unterschied Verlustarten
        KennfeldMaschineSpeichern.Analyse.Wirkungsgrad.eta_ges_mesh(:,1) = 0.1; %first row is zero %hier Unterschied Verlustarten
        % Create one "larger" eff. map that includes acceleration and recuperation
        Kennfeld_Gesamtwirkungsgrad.eff_recu = 1./KennfeldMaschineSpeichern.Analyse.Wirkungsgrad.eta_ges_mesh; %recuperation efficiency %hier Unterschied Verlustarten
        Kennfeld_Gesamtwirkungsgrad.eff_recu(Kennfeld_Gesamtwirkungsgrad.eff_recu==inf)=100;
        Kennfeld_Gesamtwirkungsgrad.eff_recu = Kennfeld_Gesamtwirkungsgrad.eff_recu(2:end,:);
        Kennfeld_Gesamtwirkungsgrad.eff = [flipud(Kennfeld_Gesamtwirkungsgrad.eff_recu); KennfeldMaschineSpeichern.Analyse.Wirkungsgrad.eta_ges_mesh]; %hier Unterschied Verlustarten
    
        %Volllastkennlinie
        %motor.T_max = transpose (Maschinendaten.Analyse.Betriebsdaten.mot.M_max_vec);                                                            %noch auf 2000 interpolieren?
        KennfeldMaschineSpeichern.Analyse.Betriebsdaten.mot.M_max_vec=transpose(KennfeldMaschineSpeichern.Analyse.Betriebsdaten.mot.M_max_vec);
        Kennfeld_Gesamtwirkungsgrad.T_max = transpose(KennfeldMaschineSpeichern.Analyse.Betriebsdaten.mot.M_max_vec);  
        %Drehzahlachse skaliert
        Kennfeld_Gesamtwirkungsgrad.T_max_n_axis = transpose(KennfeldMaschineSpeichern.Analyse.Betriebsdaten.mot.M_vec);                                                                          %noch auf 2000 interpolieren?

         save (['3_Ergebnisse/',handles.Entwurf.Optionen.folder_id,'/','2_Analyse/', KennfeldMaschineSpeichern.Entwurf.Optionen.Maschinentyp,'_',num2str(KennfeldMaschineSpeichern.Entwurf.Bemessungswerte.P_N/1000),'_',num2str(KennfeldMaschineSpeichern.Entwurf.Bemessungswerte.n_N),'_',num2str(KennfeldMaschineSpeichern.Analyse.Optionen.n_max),'_',num2str(KennfeldMaschineSpeichern.Entwurf.Bemessungswerte.U_N),'_','Kennfeld_Gesamtwirkungsgrad','.mat'],'Kennfeld_Gesamtwirkungsgrad');
         save (['LDSimulation/Vehicle/Para_Powertrain/Motor_efficiency/MOT_memory//', KennfeldMaschineSpeichern.Entwurf.Optionen.Maschinentyp,'_',num2str(KennfeldMaschineSpeichern.Entwurf.Bemessungswerte.P_N/1000),'_',num2str(KennfeldMaschineSpeichern.Entwurf.Bemessungswerte.n_N),'_',num2str(KennfeldMaschineSpeichern.Analyse.Optionen.n_max),'_',num2str(KennfeldMaschineSpeichern.Entwurf.Bemessungswerte.U_N),'_', 'Kennfeld_Gesamtwirkungsgrad' ,'.mat'],'Kennfeld_Gesamtwirkungsgrad');
        
        clabel(C,h,Plot_Label_Vektor_eta);
        datacursormode on
        
    elseif(strcmp(handles.Plot.Auswahl_Plot,'Verluste gesamt'))
        handles = plotKennfeld(handles);
        title(handles.axes_Ergebnis_Plot,'Verluste gesamt $P_{V,ges}$ in W','interpreter','latex','FontSize', 20);
        contourf(handles.axes_Ergebnis_Plot,handles.Analyse.Betriebsdaten.n_m_mesh,handles.Analyse.Betriebsdaten.M_max_mesh,handles.Analyse.Verluste.P_vges_mesh);
        datacursormode on
        
        % Save Kennfeld neu
        KennfeldMaschine_Verluste_gesamt = struct;
        KennfeldMaschine_Verluste_gesamt_Speichern = handles;
        KennfeldMaschine_Verluste_gesamt.n = KennfeldMaschine_Verluste_gesamt_Speichern.Analyse.Betriebsdaten.n_m_mesh;
        KennfeldMaschine_Verluste_gesamt.M = KennfeldMaschine_Verluste_gesamt_Speichern.Analyse.Betriebsdaten.M_max_mesh;
        KennfeldMaschine_Verluste_gesamt.verlges = KennfeldMaschine_Verluste_gesamt_Speichern.Analyse.Verluste.P_vges_mesh;
         save (['3_Ergebnisse/',handles.Entwurf.Optionen.folder_id,'/','2_Analyse/', 'KennfeldMaschine_Verluste_gesamt','.mat'],'KennfeldMaschine_Verluste_gesamt');
        folder_id = [KennfeldMaschine_Verluste_gesamt_Speichern.Entwurf.Optionen.folder_id]; %file_id, '_data'];
        file_id = [KennfeldMaschine_Verluste_gesamt_Speichern.Entwurf.Optionen.file_id  ];
        
        % Konvertieren der Daten in richtige Form
        % Bemessungsgroessen
        Kennfeld_Gesamtverluste=struct;
        Kennfeld_Gesamtverluste.type = KennfeldMaschine_Verluste_gesamt_Speichern.Entwurf.Optionen.Maschinentyp;
        Kennfeld_Gesamtverluste.power = KennfeldMaschine_Verluste_gesamt_Speichern.Entwurf.Bemessungswerte.P_N/1000;
        Kennfeld_Gesamtverluste.n_n = KennfeldMaschine_Verluste_gesamt_Speichern.Entwurf.Bemessungswerte.n_N;
        Kennfeld_Gesamtverluste.n_max = KennfeldMaschine_Verluste_gesamt_Speichern.Analyse.Optionen.n_max;
        Kennfeld_Gesamtverluste.U_n = KennfeldMaschine_Verluste_gesamt_Speichern.Entwurf.Bemessungswerte.U_N;
        
        %Drehzahlachse - X-Achse Kennfeld
        Kennfeld_Gesamtverluste.eff_n_axis = transpose(KennfeldMaschine_Verluste_gesamt_Speichern.Analyse.Betriebsdaten.mot.n_m_vec);

        %Drehmoment - Y-Achse
        KennfeldMaschine_Verluste_gesamt_Speichern.Analyse.Betriebsdaten.mot.M_vec=transpose(fliplr(KennfeldMaschine_Verluste_gesamt_Speichern.Analyse.Betriebsdaten.mot.M_vec));
        Kennfeld_Gesamtverluste.eff_T_axis = [(fliplr(transpose(KennfeldMaschine_Verluste_gesamt_Speichern.Analyse.Betriebsdaten.mot.M_vec(2:end)))) * (-1), transpose(KennfeldMaschine_Verluste_gesamt_Speichern.Analyse.Betriebsdaten.mot.M_vec)];
        %Kennfeld_Gesamtwirkungsgrad.eff_T_axis = [fliplr(transpose(KennfeldMaschineSpeichern.Analyse.Betriebsdaten.mot.M_vec)),(transpose(KennfeldMaschineSpeichern.Analyse.Betriebsdaten.mot.M_vec(2:end))) * (-1)];

        %Diagrammbereich                                                             
        % Manipulate efficiency map
        KennfeldMaschine_Verluste_gesamt_Speichern.Analyse.Verluste.P_vges_mesh=flipud(KennfeldMaschine_Verluste_gesamt_Speichern.Analyse.Verluste.P_vges_mesh); %hier Unterschied Verlustarten
        KennfeldMaschine_Verluste_gesamt_Speichern.Analyse.Verluste.P_vges_mesh(KennfeldMaschine_Verluste_gesamt_Speichern.Analyse.Verluste.P_vges_mesh<0.1) = 0.1; %Hier Unterschied Verlustarten %set small efficiency values to 0.1
        KennfeldMaschine_Verluste_gesamt_Speichern.Analyse.Verluste.P_vges_mesh(1,:) = 0; %first colomn is zero %hier Unterschied Verlustarten
        KennfeldMaschine_Verluste_gesamt_Speichern.Analyse.Verluste.P_vges_mesh(:,1) = 0.1; %first row is zero %hier Unterschied Verlustarten
        % Create one "larger" eff. map that includes acceleration and recuperation
        Kennfeld_Gesamtverluste.eff_recu = 1./KennfeldMaschine_Verluste_gesamt_Speichern.Analyse.Verluste.P_vges_mesh; %recuperation efficiency %hier Unterschied Verlustarten
        Kennfeld_Gesamtverluste.eff_recu(Kennfeld_Gesamtverluste.eff_recu==inf)=100;
        Kennfeld_Gesamtverluste.eff_recu = Kennfeld_Gesamtverluste.eff_recu(2:end,:);
        Kennfeld_Gesamtverluste.eff = [flipud(Kennfeld_Gesamtverluste.eff_recu); KennfeldMaschine_Verluste_gesamt_Speichern.Analyse.Verluste.P_vges_mesh]; %hier Unterschied Verlustarten
    
        %Volllastkennlinie
        %motor.T_max = transpose (Maschinendaten.Analyse.Betriebsdaten.mot.M_max_vec);                                                            %noch auf 2000 interpolieren?
        KennfeldMaschine_Verluste_gesamt_Speichern.Analyse.Betriebsdaten.mot.M_max_vec=transpose(KennfeldMaschine_Verluste_gesamt_Speichern.Analyse.Betriebsdaten.mot.M_max_vec);
        Kennfeld_Gesamtverluste.T_max = transpose(KennfeldMaschine_Verluste_gesamt_Speichern.Analyse.Betriebsdaten.mot.M_max_vec);  
        %Drehzahlachse skaliert
        Kennfeld_Gesamtverluste.T_max_n_axis = transpose(KennfeldMaschine_Verluste_gesamt_Speichern.Analyse.Betriebsdaten.mot.M_vec);                                                                          %noch auf 2000 interpolieren?
        
        save (['3_Ergebnisse/',handles.Entwurf.Optionen.folder_id,'/','2_Analyse/', KennfeldMaschine_Verluste_gesamt_Speichern.Entwurf.Optionen.Maschinentyp,'_',num2str(KennfeldMaschine_Verluste_gesamt_Speichern.Entwurf.Bemessungswerte.P_N/1000),'_',num2str(KennfeldMaschine_Verluste_gesamt_Speichern.Entwurf.Bemessungswerte.n_N),'_',num2str(KennfeldMaschine_Verluste_gesamt_Speichern.Analyse.Optionen.n_max),'_',num2str(KennfeldMaschine_Verluste_gesamt_Speichern.Entwurf.Bemessungswerte.U_N),'_','Kennfeld_Gesamtverluste','.mat'],'Kennfeld_Gesamtverluste');
        save (['LDSimulation/Vehicle/Para_Powertrain/Motor_efficiency/MOT_memory//', KennfeldMaschine_Verluste_gesamt_Speichern.Entwurf.Optionen.Maschinentyp,'_',num2str(KennfeldMaschine_Verluste_gesamt_Speichern.Entwurf.Bemessungswerte.P_N/1000),'_',num2str(KennfeldMaschine_Verluste_gesamt_Speichern.Entwurf.Bemessungswerte.n_N),'_',num2str(KennfeldMaschine_Verluste_gesamt_Speichern.Analyse.Optionen.n_max),'_',num2str(KennfeldMaschine_Verluste_gesamt_Speichern.Entwurf.Bemessungswerte.U_N),'_', 'Kennfeld_Gesamtverluste' ,'.mat'],'Kennfeld_Gesamtverluste');
        
    elseif(strcmp(handles.Plot.Auswahl_Plot,'Wirkungsgrad Wicklungsverluste'))
        handles = plotKennfeld(handles);
        title(handles.axes_Ergebnis_Plot,'Wirkungsgrad Wicklungsverluste $\eta_{V,w}$ in \%','interpreter','latex','FontSize', 20);
        Plot_Label_Vektor_eta = [0.7,0.75,0.8,0.84,0.86,0.88,0.9,0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.98,0.99,0.995];
        [C,h]=contourf(handles.axes_Ergebnis_Plot,handles.Analyse.Betriebsdaten.n_m_mesh,handles.Analyse.Betriebsdaten.M_max_mesh,handles.Analyse.Wirkungsgrad.eta_vw_mesh,Plot_Label_Vektor_eta);
        clabel(C,h,Plot_Label_Vektor_eta);
        datacursormode on
        
        % Save Kennfeld neu
        KennfeldMaschine_Wirkungsgrad_Wicklungsverluste = struct;
        KennfeldMaschine_Wirkungsgrad_Wicklungsverluste_Speichern = handles;
        KennfeldMaschine_Wirkungsgrad_Wicklungsverluste.n = KennfeldMaschine_Wirkungsgrad_Wicklungsverluste_Speichern.Analyse.Betriebsdaten.n_m_mesh;
        KennfeldMaschine_Wirkungsgrad_Wicklungsverluste.M = KennfeldMaschine_Wirkungsgrad_Wicklungsverluste_Speichern.Analyse.Betriebsdaten.M_max_mesh;
        KennfeldMaschine_Wirkungsgrad_Wicklungsverluste.wirkwick = KennfeldMaschine_Wirkungsgrad_Wicklungsverluste_Speichern.Analyse.Wirkungsgrad.eta_vw_mesh;
         save (['3_Ergebnisse/',handles.Entwurf.Optionen.folder_id,'/','2_Analyse/', 'KennfeldMaschine_Wirkungsgrad_Wicklungsverluste','.mat'],'KennfeldMaschine_Wirkungsgrad_Wicklungsverluste');
        folder_id = [KennfeldMaschine_Wirkungsgrad_Wicklungsverluste_Speichern.Entwurf.Optionen.folder_id]; %file_id, '_data'];
        file_id = [KennfeldMaschine_Wirkungsgrad_Wicklungsverluste_Speichern.Entwurf.Optionen.file_id];
        
        % Konvertieren der Daten in richtige Form
        % Bemessungsgroessen
        Kennfeld_Wirkungsgrad_Wicklungsverluste = struct;
        Kennfeld_Wirkungsgrad_Wicklungsverluste.type = KennfeldMaschine_Wirkungsgrad_Wicklungsverluste_Speichern.Entwurf.Optionen.Maschinentyp;
        Kennfeld_Wirkungsgrad_Wicklungsverluste.power = KennfeldMaschine_Wirkungsgrad_Wicklungsverluste_Speichern.Entwurf.Bemessungswerte.P_N/1000;
        Kennfeld_Wirkungsgrad_Wicklungsverluste.n_n = KennfeldMaschine_Wirkungsgrad_Wicklungsverluste_Speichern.Entwurf.Bemessungswerte.n_N;
        Kennfeld_Wirkungsgrad_Wicklungsverluste.n_max = KennfeldMaschine_Wirkungsgrad_Wicklungsverluste_Speichern.Analyse.Optionen.n_max;
        Kennfeld_Wirkungsgrad_Wicklungsverluste.U_n = KennfeldMaschine_Wirkungsgrad_Wicklungsverluste_Speichern.Entwurf.Bemessungswerte.U_N;
        
        % Drehzahlachse - X-Achse Kennfeld
        Kennfeld_Wirkungsgrad_Wicklungsverluste.eff_n_axis = transpose(KennfeldMaschine_Wirkungsgrad_Wicklungsverluste_Speichern.Analyse.Betriebsdaten.mot.n_m_vec);

        % Drehmoment - Y-Achse
        KennfeldMaschine_Wirkungsgrad_Wicklungsverluste_Speichern.Analyse.Betriebsdaten.mot.M_vec=transpose(fliplr(KennfeldMaschine_Wirkungsgrad_Wicklungsverluste_Speichern.Analyse.Betriebsdaten.mot.M_vec));
        Kennfeld_Wirkungsgrad_Wicklungsverluste.eff_T_axis = [(fliplr(transpose(KennfeldMaschine_Wirkungsgrad_Wicklungsverluste_Speichern.Analyse.Betriebsdaten.mot.M_vec(2:end)))) * (-1), transpose(KennfeldMaschine_Wirkungsgrad_Wicklungsverluste_Speichern.Analyse.Betriebsdaten.mot.M_vec)];
        
        % Diagrammbereich                                                             
        % Manipulate efficiency map
        KennfeldMaschine_Wirkungsgrad_Wicklungsverluste_Speichern.Analyse.Wirkungsgrad.eta_vw_mesh=flipud(KennfeldMaschine_Wirkungsgrad_Wicklungsverluste_Speichern.Analyse.Wirkungsgrad.eta_vw_mesh); %hier Unterschied Verlustarten
        KennfeldMaschine_Wirkungsgrad_Wicklungsverluste_Speichern.Analyse.Wirkungsgrad.eta_vw_mesh(KennfeldMaschine_Wirkungsgrad_Wicklungsverluste_Speichern.Analyse.Wirkungsgrad.eta_vw_mesh<0.1) = 0.1; %Hier Unterschied Verlustarten %set small efficiency values to 0.1
        KennfeldMaschine_Wirkungsgrad_Wicklungsverluste_Speichern.Analyse.Wirkungsgrad.eta_vw_mesh(1,:) = 0; %first colomn is zero %hier Unterschied Verlustarten
        KennfeldMaschine_Wirkungsgrad_Wicklungsverluste_Speichern.Analyse.Wirkungsgrad.eta_vw_mesh(:,1) = 0.1; %first row is zero %hier Unterschied Verlustarten
        % Create one "larger" eff. map that includes acceleration and recuperation
        Kennfeld_Wirkungsgrad_Wicklungsverluste.eff_recu = 1./KennfeldMaschine_Wirkungsgrad_Wicklungsverluste_Speichern.Analyse.Wirkungsgrad.eta_vw_mesh; %recuperation efficiency %hier Unterschied Verlustarten
        Kennfeld_Wirkungsgrad_Wicklungsverluste.eff_recu(Kennfeld_Wirkungsgrad_Wicklungsverluste.eff_recu==inf)=100;
        Kennfeld_Wirkungsgrad_Wicklungsverluste.eff_recu = Kennfeld_Wirkungsgrad_Wicklungsverluste.eff_recu(2:end,:);
        Kennfeld_Wirkungsgrad_Wicklungsverluste.eff = [flipud(Kennfeld_Wirkungsgrad_Wicklungsverluste.eff_recu); KennfeldMaschine_Wirkungsgrad_Wicklungsverluste_Speichern.Analyse.Wirkungsgrad.eta_vw_mesh]; %hier Unterschied Verlustarten
    
        % Volllastkennlinie
        KennfeldMaschine_Wirkungsgrad_Wicklungsverluste_Speichern.Analyse.Betriebsdaten.mot.M_max_vec=transpose(KennfeldMaschine_Wirkungsgrad_Wicklungsverluste_Speichern.Analyse.Betriebsdaten.mot.M_max_vec);
        Kennfeld_Wirkungsgrad_Wicklungsverluste.T_max = transpose(KennfeldMaschine_Wirkungsgrad_Wicklungsverluste_Speichern.Analyse.Betriebsdaten.mot.M_max_vec);  
        
        % Drehzahlachse skaliert
        Kennfeld_Wirkungsgrad_Wicklungsverluste.T_max_n_axis = transpose(KennfeldMaschine_Wirkungsgrad_Wicklungsverluste_Speichern.Analyse.Betriebsdaten.mot.M_vec);                                                                          %noch auf 2000 interpolieren?
        
         save (['3_Ergebnisse/',handles.Entwurf.Optionen.folder_id,'/','2_Analyse/', KennfeldMaschine_Wirkungsgrad_Wicklungsverluste_Speichern.Entwurf.Optionen.Maschinentyp,'_',num2str(KennfeldMaschine_Wirkungsgrad_Wicklungsverluste_Speichern.Entwurf.Bemessungswerte.P_N/1000),'_',num2str(KennfeldMaschine_Wirkungsgrad_Wicklungsverluste_Speichern.Entwurf.Bemessungswerte.n_N),'_',num2str(KennfeldMaschine_Wirkungsgrad_Wicklungsverluste_Speichern.Analyse.Optionen.n_max),'_',num2str(KennfeldMaschine_Wirkungsgrad_Wicklungsverluste_Speichern.Entwurf.Bemessungswerte.U_N),'_','Kennfeld_Wirkungsgrad_Wicklungsverluste','.mat'],'Kennfeld_Wirkungsgrad_Wicklungsverluste');
         save (['LDSimulation/Vehicle/Para_Powertrain/Motor_efficiency/MOT_memory//', KennfeldMaschine_Wirkungsgrad_Wicklungsverluste_Speichern.Entwurf.Optionen.Maschinentyp,'_',num2str(KennfeldMaschine_Wirkungsgrad_Wicklungsverluste_Speichern.Entwurf.Bemessungswerte.P_N/1000),'_',num2str(KennfeldMaschine_Wirkungsgrad_Wicklungsverluste_Speichern.Entwurf.Bemessungswerte.n_N),'_',num2str(KennfeldMaschine_Wirkungsgrad_Wicklungsverluste_Speichern.Analyse.Optionen.n_max),'_',num2str(KennfeldMaschine_Wirkungsgrad_Wicklungsverluste_Speichern.Entwurf.Bemessungswerte.U_N),'_', 'Kennfeld_Wirkungsgrad_Wicklungsverluste' ,'.mat'],'Kennfeld_Wirkungsgrad_Wicklungsverluste');
        
    elseif(strcmp(handles.Plot.Auswahl_Plot,'Wicklungsverluste'))
        handles = plotKennfeld(handles);
        title(handles.axes_Ergebnis_Plot,'Wicklungsverluste $P_{V,w}$ in W','interpreter','latex','FontSize', 20);
        contourf(handles.axes_Ergebnis_Plot,handles.Analyse.Betriebsdaten.n_m_mesh,handles.Analyse.Betriebsdaten.M_max_mesh,handles.Analyse.Verluste.P_vw_mesh);
        datacursormode on
        
        % Save Kennfeld neu
        KennfeldMaschine_Wicklungsverluste = struct;
        KennfeldMaschine_Wicklungsverluste_Speichern = handles;
        KennfeldMaschine_Wicklungsverluste.n = KennfeldMaschine_Wicklungsverluste_Speichern.Analyse.Betriebsdaten.n_m_mesh;
        KennfeldMaschine_Wicklungsverluste.M = KennfeldMaschine_Wicklungsverluste_Speichern.Analyse.Betriebsdaten.M_max_mesh;
        KennfeldMaschine_Wicklungsverluste.verlwick = KennfeldMaschine_Wicklungsverluste_Speichern.Analyse.Verluste.P_vw_mesh;
         save (['3_Ergebnisse/',handles.Entwurf.Optionen.folder_id,'/','2_Analyse/', 'KennfeldMaschine_Wicklungsverluste','.mat'],'KennfeldMaschine_Wicklungsverluste');
        folder_id = [KennfeldMaschine_Wicklungsverluste_Speichern.Entwurf.Optionen.folder_id]; %file_id, '_data'];
        file_id = [KennfeldMaschine_Wicklungsverluste_Speichern.Entwurf.Optionen.file_id];
        
        % Konvertieren der Daten in richtige Form
        % Bemessungsgroessen
        Kennfeld_Wicklungsverluste = struct;
        Kennfeld_Wicklungsverluste.type = KennfeldMaschine_Wicklungsverluste_Speichern.Entwurf.Optionen.Maschinentyp;
        Kennfeld_Wicklungsverluste.power = KennfeldMaschine_Wicklungsverluste_Speichern.Entwurf.Bemessungswerte.P_N/1000;
        Kennfeld_Wicklungsverluste.n_n = KennfeldMaschine_Wicklungsverluste_Speichern.Entwurf.Bemessungswerte.n_N;
        Kennfeld_Wicklungsverluste.n_max = KennfeldMaschine_Wicklungsverluste_Speichern.Analyse.Optionen.n_max;
        Kennfeld_Wicklungsverluste.U_n = KennfeldMaschine_Wicklungsverluste_Speichern.Entwurf.Bemessungswerte.U_N;
        
        % Drehzahlachse - X-Achse Kennfeld
        Kennfeld_Wicklungsverluste.eff_n_axis = transpose(KennfeldMaschine_Wicklungsverluste_Speichern.Analyse.Betriebsdaten.mot.n_m_vec);

        % Drehmoment - Y-Achse
        KennfeldMaschine_Wicklungsverluste_Speichern.Analyse.Betriebsdaten.mot.M_vec=transpose(fliplr(KennfeldMaschine_Wicklungsverluste_Speichern.Analyse.Betriebsdaten.mot.M_vec));
        Kennfeld_Wicklungsverluste.eff_T_axis = [(fliplr(transpose(KennfeldMaschine_Wicklungsverluste_Speichern.Analyse.Betriebsdaten.mot.M_vec(2:end)))) * (-1), transpose(KennfeldMaschine_Wicklungsverluste_Speichern.Analyse.Betriebsdaten.mot.M_vec)];
        
        % Diagrammbereich                                                             
        % Manipulate efficiency map
        KennfeldMaschine_Wicklungsverluste_Speichern.Analyse.Verluste.P_vw_mesh=flipud(KennfeldMaschine_Wicklungsverluste_Speichern.Analyse.Verluste.P_vw_mesh); %hier Unterschied Verlustarten
        KennfeldMaschine_Wicklungsverluste_Speichern.Analyse.Verluste.P_vw_mesh(KennfeldMaschine_Wicklungsverluste_Speichern.Analyse.Verluste.P_vw_mesh<0.1) = 0.1; %Hier Unterschied Verlustarten %set small efficiency values to 0.1
        KennfeldMaschine_Wicklungsverluste_Speichern.Analyse.Verluste.P_vw_mesh(1,:) = 0; %first colomn is zero %hier Unterschied Verlustarten
        KennfeldMaschine_Wicklungsverluste_Speichern.Analyse.Verluste.P_vw_mesh(:,1) = 0.1; %first row is zero %hier Unterschied Verlustarten
        
        % Create one "larger" eff. map that includes acceleration and recuperation
        Kennfeld_Wicklungsverluste.eff_recu = 1./KennfeldMaschine_Wicklungsverluste_Speichern.Analyse.Verluste.P_vw_mesh; %recuperation efficiency %hier Unterschied Verlustarten
        Kennfeld_Wicklungsverluste.eff_recu(Kennfeld_Wicklungsverluste.eff_recu==inf)=100;
        Kennfeld_Wicklungsverluste.eff_recu = Kennfeld_Wicklungsverluste.eff_recu(2:end,:);
        Kennfeld_Wicklungsverluste.eff = [flipud(Kennfeld_Wicklungsverluste.eff_recu); KennfeldMaschine_Wicklungsverluste_Speichern.Analyse.Verluste.P_vw_mesh]; %hier Unterschied Verlustarten
    
        % Volllastkennlinie
        KennfeldMaschine_Wicklungsverluste_Speichern.Analyse.Betriebsdaten.mot.M_max_vec=transpose(KennfeldMaschine_Wicklungsverluste_Speichern.Analyse.Betriebsdaten.mot.M_max_vec);
        Kennfeld_Wicklungsverluste.T_max = transpose(KennfeldMaschine_Wicklungsverluste_Speichern.Analyse.Betriebsdaten.mot.M_max_vec);  
        
        % Drehzahlachse skaliert
        Kennfeld_Wicklungsverluste.T_max_n_axis = transpose(KennfeldMaschine_Wicklungsverluste_Speichern.Analyse.Betriebsdaten.mot.M_vec);                                                          %noch auf 2000 interpolieren?
        
         save (['3_Ergebnisse/',handles.Entwurf.Optionen.folder_id,'/', '2_Analyse/', KennfeldMaschine_Wicklungsverluste_Speichern.Entwurf.Optionen.Maschinentyp,'_',num2str(KennfeldMaschine_Wicklungsverluste_Speichern.Entwurf.Bemessungswerte.P_N/1000),'_',num2str(KennfeldMaschine_Wicklungsverluste_Speichern.Entwurf.Bemessungswerte.n_N),'_',num2str(KennfeldMaschine_Wicklungsverluste_Speichern.Analyse.Optionen.n_max),'_',num2str(KennfeldMaschine_Wicklungsverluste_Speichern.Entwurf.Bemessungswerte.U_N),'_','Kennfeld_Wicklungsverluste','.mat'],'Kennfeld_Wicklungsverluste');
         save (['LDSimulation/Vehicle/Para_Powertrain/Motor_efficiency/MOT_memory//', KennfeldMaschine_Wicklungsverluste_Speichern.Entwurf.Optionen.Maschinentyp,'_',num2str(KennfeldMaschine_Wicklungsverluste_Speichern.Entwurf.Bemessungswerte.P_N/1000),'_',num2str(KennfeldMaschine_Wicklungsverluste_Speichern.Entwurf.Bemessungswerte.n_N),'_',num2str(KennfeldMaschine_Wicklungsverluste_Speichern.Analyse.Optionen.n_max),'_',num2str(KennfeldMaschine_Wicklungsverluste_Speichern.Entwurf.Bemessungswerte.U_N),'_', 'Kennfeld_Wicklungsverluste' ,'.mat'],'Kennfeld_Wicklungsverluste');
        
    elseif(strcmp(handles.Plot.Auswahl_Plot,'Wirkungsgrad Ummagnetisierungsverluste'))
        handles = plotKennfeld(handles);
        title(handles.axes_Ergebnis_Plot,'Wirkungsgrad Ummagnetisierungsverluste $\eta_{V,u}$ in \%','interpreter','latex','FontSize', 20);
        Plot_Label_Vektor_eta = [0.7,0.75,0.8,0.84,0.86,0.88,0.9,0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.98,0.99,0.995];
        [C,h]=contourf(handles.axes_Ergebnis_Plot,handles.Analyse.Betriebsdaten.n_m_mesh,handles.Analyse.Betriebsdaten.M_max_mesh,handles.Analyse.Wirkungsgrad.eta_fe_mesh,Plot_Label_Vektor_eta);
        clabel(C,h,Plot_Label_Vektor_eta);
        datacursormode on
        
        % Save Kennfeld neu
        KennfeldMaschine_Wirkungsgrad_Ummagnetisierungsverluste = struct;
        KennfeldMaschine_WirkUmmagnetisierungsverluste_Speichern = handles;
        KennfeldMaschine_Wirkungsgrad_Ummagnetisierungsverluste.n = KennfeldMaschine_WirkUmmagnetisierungsverluste_Speichern.Analyse.Betriebsdaten.n_m_mesh;
        KennfeldMaschine_Wirkungsgrad_Ummagnetisierungsverluste.M = KennfeldMaschine_WirkUmmagnetisierungsverluste_Speichern.Analyse.Betriebsdaten.M_max_mesh;
        KennfeldMaschine_Wirkungsgrad_Ummagnetisierungsverluste.wirkummagn = KennfeldMaschine_WirkUmmagnetisierungsverluste_Speichern.Analyse.Wirkungsgrad.eta_fe_mesh;
         save (['3_Ergebnisse/',handles.Entwurf.Optionen.folder_id,'/','2_Analyse/','KennfeldMaschine_Wirkungsgrad_Ummagnetisierungsverluste','.mat'],'KennfeldMaschine_Wirkungsgrad_Ummagnetisierungsverluste');
        folder_id = [KennfeldMaschine_WirkUmmagnetisierungsverluste_Speichern.Entwurf.Optionen.folder_id]; %file_id, '_data'];
        file_id = [KennfeldMaschine_WirkUmmagnetisierungsverluste_Speichern.Entwurf.Optionen.file_id];
        
        % Konvertieren der Daten in richtige Form
        % Bemessungsgroessen
        Kennfeld_Wirkungsgrad_Ummagnetisierungsverluste = struct;
        Kennfeld_Wirkungsgrad_Ummagnetisierungsverluste.type = KennfeldMaschine_WirkUmmagnetisierungsverluste_Speichern.Entwurf.Optionen.Maschinentyp;
        Kennfeld_Wirkungsgrad_Ummagnetisierungsverluste.power = KennfeldMaschine_WirkUmmagnetisierungsverluste_Speichern.Entwurf.Bemessungswerte.P_N/1000;
        Kennfeld_Wirkungsgrad_Ummagnetisierungsverluste.n_n = KennfeldMaschine_WirkUmmagnetisierungsverluste_Speichern.Entwurf.Bemessungswerte.n_N;
        Kennfeld_Wirkungsgrad_Ummagnetisierungsverluste.n_max = KennfeldMaschine_WirkUmmagnetisierungsverluste_Speichern.Analyse.Optionen.n_max;
        Kennfeld_Wirkungsgrad_Ummagnetisierungsverluste.U_n = KennfeldMaschine_WirkUmmagnetisierungsverluste_Speichern.Entwurf.Bemessungswerte.U_N;
        
        % Drehzahlachse - X-Achse Kennfeld
        Kennfeld_Wirkungsgrad_Ummagnetisierungsverluste.eff_n_axis = transpose(KennfeldMaschine_WirkUmmagnetisierungsverluste_Speichern.Analyse.Betriebsdaten.mot.n_m_vec);

        % Drehmoment - Y-Achse
        KennfeldMaschine_WirkUmmagnetisierungsverluste_Speichern.Analyse.Betriebsdaten.mot.M_vec=transpose(fliplr(KennfeldMaschine_WirkUmmagnetisierungsverluste_Speichern.Analyse.Betriebsdaten.mot.M_vec));
        Kennfeld_Wirkungsgrad_Ummagnetisierungsverluste.eff_T_axis = [(fliplr(transpose(KennfeldMaschine_WirkUmmagnetisierungsverluste_Speichern.Analyse.Betriebsdaten.mot.M_vec(2:end)))) * (-1), transpose(KennfeldMaschine_WirkUmmagnetisierungsverluste_Speichern.Analyse.Betriebsdaten.mot.M_vec)];
        
        % Diagrammbereich                                                             
        % Manipulate efficiency map
        KennfeldMaschine_WirkUmmagnetisierungsverluste_Speichern.Analyse.Wirkungsgrad.eta_fe_mesh=flipud(KennfeldMaschine_WirkUmmagnetisierungsverluste_Speichern.Analyse.Wirkungsgrad.eta_fe_mesh); %hier Unterschied Verlustarten
        KennfeldMaschine_WirkUmmagnetisierungsverluste_Speichern.Analyse.Wirkungsgrad.eta_fe_mesh(KennfeldMaschine_WirkUmmagnetisierungsverluste_Speichern.Analyse.Wirkungsgrad.eta_fe_mesh<0.1) = 0.1; %Hier Unterschied Verlustarten %set small efficiency values to 0.1
        KennfeldMaschine_WirkUmmagnetisierungsverluste_Speichern.Analyse.Wirkungsgrad.eta_fe_mesh(1,:) = 0; %first colomn is zero %hier Unterschied Verlustarten
        KennfeldMaschine_WirkUmmagnetisierungsverluste_Speichern.Analyse.Wirkungsgrad.eta_fe_mesh(:,1) = 0.1; %first row is zero %hier Unterschied Verlustarten
        
        % Create one "larger" eff. map that includes acceleration and recuperation
        Kennfeld_Wirkungsgrad_Ummagnetisierungsverluste.eff_recu = 1./KennfeldMaschine_WirkUmmagnetisierungsverluste_Speichern.Analyse.Wirkungsgrad.eta_fe_mesh; %recuperation efficiency %hier Unterschied Verlustarten
        Kennfeld_Wirkungsgrad_Ummagnetisierungsverluste.eff_recu(Kennfeld_Wirkungsgrad_Ummagnetisierungsverluste.eff_recu==inf)=100;
        Kennfeld_Wirkungsgrad_Ummagnetisierungsverluste.eff_recu = Kennfeld_Wirkungsgrad_Ummagnetisierungsverluste.eff_recu(2:end,:);
        Kennfeld_Wirkungsgrad_Ummagnetisierungsverluste.eff = [flipud(Kennfeld_Wirkungsgrad_Ummagnetisierungsverluste.eff_recu); KennfeldMaschine_WirkUmmagnetisierungsverluste_Speichern.Analyse.Wirkungsgrad.eta_fe_mesh]; %hier Unterschied Verlustarten
    
        % Volllastkennlinie
        KennfeldMaschine_WirkUmmagnetisierungsverluste_Speichern.Analyse.Betriebsdaten.mot.M_max_vec=transpose(KennfeldMaschine_WirkUmmagnetisierungsverluste_Speichern.Analyse.Betriebsdaten.mot.M_max_vec);
        Kennfeld_Wirkungsgrad_Ummagnetisierungsverluste.T_max = transpose(KennfeldMaschine_WirkUmmagnetisierungsverluste_Speichern.Analyse.Betriebsdaten.mot.M_max_vec);  
        
        % Drehzahlachse skaliert
        Kennfeld_Wirkungsgrad_Ummagnetisierungsverluste.T_max_n_axis = transpose(KennfeldMaschine_WirkUmmagnetisierungsverluste_Speichern.Analyse.Betriebsdaten.mot.M_vec);                                                          %noch auf 2000 interpolieren?
        
         save (['3_Ergebnisse/',handles.Entwurf.Optionen.folder_id,'/','2_Analyse/',KennfeldMaschine_WirkUmmagnetisierungsverluste_Speichern.Entwurf.Optionen.Maschinentyp,'_',num2str(KennfeldMaschine_WirkUmmagnetisierungsverluste_Speichern.Entwurf.Bemessungswerte.P_N/1000),'_',num2str(KennfeldMaschine_WirkUmmagnetisierungsverluste_Speichern.Entwurf.Bemessungswerte.n_N),'_',num2str(KennfeldMaschine_WirkUmmagnetisierungsverluste_Speichern.Analyse.Optionen.n_max),'_',num2str(KennfeldMaschine_WirkUmmagnetisierungsverluste_Speichern.Entwurf.Bemessungswerte.U_N),'_','Kennfeld_Wirkungsgrad_Ummagnetisierungsverluste','.mat'],'Kennfeld_Wirkungsgrad_Ummagnetisierungsverluste');
         save (['LDSimulation/Vehicle/Para_Powertrain/Motor_efficiency/MOT_memory//', KennfeldMaschine_WirkUmmagnetisierungsverluste_Speichern.Entwurf.Optionen.Maschinentyp,'_',num2str(KennfeldMaschine_WirkUmmagnetisierungsverluste_Speichern.Entwurf.Bemessungswerte.P_N/1000),'_',num2str(KennfeldMaschine_WirkUmmagnetisierungsverluste_Speichern.Entwurf.Bemessungswerte.n_N),'_',num2str(KennfeldMaschine_WirkUmmagnetisierungsverluste_Speichern.Analyse.Optionen.n_max),'_',num2str(KennfeldMaschine_WirkUmmagnetisierungsverluste_Speichern.Entwurf.Bemessungswerte.U_N),'_', 'Kennfeld_Wirkungsgrad_Ummagnetisierungsverluste' ,'.mat'],'Kennfeld_Wirkungsgrad_Ummagnetisierungsverluste');
        
    elseif(strcmp(handles.Plot.Auswahl_Plot,'Ummagnetisierungsverluste'))
        handles = plotKennfeld(handles);
        title(handles.axes_Ergebnis_Plot,'Ummagnetisierungsverluste $P_{V,u}$ in W','interpreter','latex','FontSize', 20);
        contourf(handles.axes_Ergebnis_Plot,handles.Analyse.Betriebsdaten.n_m_mesh,handles.Analyse.Betriebsdaten.M_max_mesh,handles.Analyse.Verluste.P_vu_mesh);
        datacursormode on
        
        % Save Kennfeld neu
        KennfeldMaschine_Ummagnetisierungsverluste = struct;
        KennfeldMaschine_Ummagnetisierungsverluste_Speichern = handles;
        KennfeldMaschine_Ummagnetisierungsverluste.n = KennfeldMaschine_Ummagnetisierungsverluste_Speichern.Analyse.Betriebsdaten.n_m_mesh;
        KennfeldMaschine_Ummagnetisierungsverluste.M = KennfeldMaschine_Ummagnetisierungsverluste_Speichern.Analyse.Betriebsdaten.M_max_mesh;
        KennfeldMaschine_Ummagnetisierungsverluste.verlummagn = KennfeldMaschine_Ummagnetisierungsverluste_Speichern.Analyse.Verluste.P_vu_mesh;
         save (['3_Ergebnisse/', handles.Entwurf.Optionen.folder_id,'/','2_Analyse/', 'KennfeldMaschine_Ummagnetisierungsverluste','.mat'],'KennfeldMaschine_Ummagnetisierungsverluste');
        folder_id = [KennfeldMaschine_Ummagnetisierungsverluste_Speichern.Entwurf.Optionen.folder_id]; %file_id, '_data'];
        file_id = [KennfeldMaschine_Ummagnetisierungsverluste_Speichern.Entwurf.Optionen.file_id];
        
        % Konvertieren der Daten in richtige Form
        % Bemessungsgroessen
        Kennfeld_Ummagnetisierungsverluste = struct;
        Kennfeld_Ummagnetisierungsverluste.type = KennfeldMaschine_Ummagnetisierungsverluste_Speichern.Entwurf.Optionen.Maschinentyp;
        Kennfeld_Ummagnetisierungsverluste.power = KennfeldMaschine_Ummagnetisierungsverluste_Speichern.Entwurf.Bemessungswerte.P_N/1000;
        Kennfeld_Ummagnetisierungsverluste.n_n = KennfeldMaschine_Ummagnetisierungsverluste_Speichern.Entwurf.Bemessungswerte.n_N;
        Kennfeld_Ummagnetisierungsverluste.n_max = KennfeldMaschine_Ummagnetisierungsverluste_Speichern.Analyse.Optionen.n_max;
        Kennfeld_Ummagnetisierungsverluste.U_n = KennfeldMaschine_Ummagnetisierungsverluste_Speichern.Entwurf.Bemessungswerte.U_N;
        
        % Drehzahlachse - X-Achse Kennfeld
        Kennfeld_Ummagnetisierungsverluste.eff_n_axis = transpose(KennfeldMaschine_Ummagnetisierungsverluste_Speichern.Analyse.Betriebsdaten.mot.n_m_vec);

        % Drehmoment - Y-Achse
        KennfeldMaschine_Ummagnetisierungsverluste_Speichern.Analyse.Betriebsdaten.mot.M_vec=transpose(fliplr(KennfeldMaschine_Ummagnetisierungsverluste_Speichern.Analyse.Betriebsdaten.mot.M_vec));
        Kennfeld_Ummagnetisierungsverluste.eff_T_axis = [(fliplr(transpose(KennfeldMaschine_Ummagnetisierungsverluste_Speichern.Analyse.Betriebsdaten.mot.M_vec(2:end)))) * (-1), transpose(KennfeldMaschine_Ummagnetisierungsverluste_Speichern.Analyse.Betriebsdaten.mot.M_vec)];
        
        % Diagrammbereich                                                             
        % Manipulate efficiency map
        KennfeldMaschine_Ummagnetisierungsverluste_Speichern.Analyse.Verluste.P_vu_mesh=flipud(KennfeldMaschine_Ummagnetisierungsverluste_Speichern.Analyse.Verluste.P_vu_mesh); %hier Unterschied Verlustarten
        KennfeldMaschine_Ummagnetisierungsverluste_Speichern.Analyse.Verluste.P_vu_mesh(KennfeldMaschine_Ummagnetisierungsverluste_Speichern.Analyse.Verluste.P_vu_mesh<0.1) = 0.1; %Hier Unterschied Verlustarten %set small efficiency values to 0.1
        KennfeldMaschine_Ummagnetisierungsverluste_Speichern.Analyse.Verluste.P_vu_mesh(1,:) = 0; %first colomn is zero %hier Unterschied Verlustarten
        KennfeldMaschine_Ummagnetisierungsverluste_Speichern.Analyse.Verluste.P_vu_mesh(:,1) = 0.1; %first row is zero %hier Unterschied Verlustarten
        
        % Create one "larger" eff. map that includes acceleration and recuperation
        Kennfeld_Ummagnetisierungsverluste.eff_recu = 1./KennfeldMaschine_Ummagnetisierungsverluste_Speichern.Analyse.Verluste.P_vu_mesh; %recuperation efficiency %hier Unterschied Verlustarten
        Kennfeld_Ummagnetisierungsverluste.eff_recu(Kennfeld_Ummagnetisierungsverluste.eff_recu==inf)=100;
        Kennfeld_Ummagnetisierungsverluste.eff_recu = Kennfeld_Ummagnetisierungsverluste.eff_recu(2:end,:);
        Kennfeld_Ummagnetisierungsverluste.eff = [flipud(Kennfeld_Ummagnetisierungsverluste.eff_recu); KennfeldMaschine_Ummagnetisierungsverluste_Speichern.Analyse.Verluste.P_vu_mesh]; %hier Unterschied Verlustarten
    
        % Volllastkennlinie
        KennfeldMaschine_Ummagnetisierungsverluste_Speichern.Analyse.Betriebsdaten.mot.M_max_vec=transpose(KennfeldMaschine_Ummagnetisierungsverluste_Speichern.Analyse.Betriebsdaten.mot.M_max_vec);
        Kennfeld_Ummagnetisierungsverluste.T_max = transpose(KennfeldMaschine_Ummagnetisierungsverluste_Speichern.Analyse.Betriebsdaten.mot.M_max_vec);  
        
        % Drehzahlachse skaliert
        Kennfeld_Ummagnetisierungsverluste.T_max_n_axis = transpose(KennfeldMaschine_Ummagnetisierungsverluste_Speichern.Analyse.Betriebsdaten.mot.M_vec);                                                          %noch auf 2000 interpolieren?

         save (['3_Ergebnisse/',handles.Entwurf.Optionen.folder_id,'/', '2_Analyse/', KennfeldMaschine_Ummagnetisierungsverluste_Speichern.Entwurf.Optionen.Maschinentyp,'_',num2str(KennfeldMaschine_Ummagnetisierungsverluste_Speichern.Entwurf.Bemessungswerte.P_N/1000),'_',num2str(KennfeldMaschine_Ummagnetisierungsverluste_Speichern.Entwurf.Bemessungswerte.n_N),'_',num2str(KennfeldMaschine_Ummagnetisierungsverluste_Speichern.Analyse.Optionen.n_max),'_',num2str(KennfeldMaschine_Ummagnetisierungsverluste_Speichern.Entwurf.Bemessungswerte.U_N),'_','Kennfeld_Ummagnetisierungsverluste','.mat'],'Kennfeld_Ummagnetisierungsverluste');
         save (['LDSimulation/Vehicle/Para_Powertrain/Motor_efficiency/MOT_memory//', KennfeldMaschine_Ummagnetisierungsverluste_Speichern.Entwurf.Optionen.Maschinentyp,'_',num2str(KennfeldMaschine_Ummagnetisierungsverluste_Speichern.Entwurf.Bemessungswerte.P_N/1000),'_',num2str(KennfeldMaschine_Ummagnetisierungsverluste_Speichern.Entwurf.Bemessungswerte.n_N),'_',num2str(KennfeldMaschine_Ummagnetisierungsverluste_Speichern.Analyse.Optionen.n_max),'_',num2str(KennfeldMaschine_Ummagnetisierungsverluste_Speichern.Entwurf.Bemessungswerte.U_N),'_', 'Kennfeld_Ummagnetisierungsverluste' ,'.mat'],'Kennfeld_Ummagnetisierungsverluste');
                
    elseif(strcmp(handles.Plot.Auswahl_Plot,'Wirkungsgrad mechanische Verluste'))
        handles = plotKennfeld(handles);
        title(handles.axes_Ergebnis_Plot,'Wirkungsgrad mechanische Verluste $\eta_{V,mech}$ in \%','interpreter','latex','FontSize', 20);
        Plot_Label_Vektor_eta = [0.7,0.75,0.8,0.84,0.86,0.88,0.9,0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.98,0.99,0.995];
        [C,h]=contourf(handles.axes_Ergebnis_Plot,handles.Analyse.Betriebsdaten.n_m_mesh,handles.Analyse.Betriebsdaten.M_max_mesh,handles.Analyse.Wirkungsgrad.eta_vme_mesh,Plot_Label_Vektor_eta);
        clabel(C,h,Plot_Label_Vektor_eta);
        datacursormode on
        
        % Save Kennfeld neu
        KennfeldMaschine_WirkmechVerluste = struct;
        KennfeldMaschine_WirkmechVerluste_Speichern = handles;
        KennfeldMaschine_WirkmechVerluste.n = KennfeldMaschine_WirkmechVerluste_Speichern.Analyse.Betriebsdaten.n_m_mesh;
        KennfeldMaschine_WirkmechVerluste.M = KennfeldMaschine_WirkmechVerluste_Speichern.Analyse.Betriebsdaten.M_max_mesh;
        KennfeldMaschine_WirkmechVerluste.wirkmechverl = KennfeldMaschine_WirkmechVerluste_Speichern.Analyse.Wirkungsgrad.eta_vme_mesh;
         save (['3_Ergebnisse/',handles.Entwurf.Optionen.folder_id,'/','2_Analyse/', 'KennfeldMaschine_WirkmechVerluste','.mat'],'KennfeldMaschine_WirkmechVerluste');
        folder_id = [KennfeldMaschine_WirkmechVerluste_Speichern.Entwurf.Optionen.folder_id]; %file_id, '_data'];
        file_id = [KennfeldMaschine_WirkmechVerluste_Speichern.Entwurf.Optionen.file_id];
        
        % Konvertieren der Daten in richtige Form
        % Bemessungsgroessen
        Kennfeld_WirkmechVerluste = struct;
        Kennfeld_WirkmechVerluste.type = KennfeldMaschine_WirkmechVerluste_Speichern.Entwurf.Optionen.Maschinentyp;
        Kennfeld_WirkmechVerluste.power = KennfeldMaschine_WirkmechVerluste_Speichern.Entwurf.Bemessungswerte.P_N/1000;
        Kennfeld_WirkmechVerluste.n_n = KennfeldMaschine_WirkmechVerluste_Speichern.Entwurf.Bemessungswerte.n_N;
        Kennfeld_WirkmechVerluste.n_max = KennfeldMaschine_WirkmechVerluste_Speichern.Analyse.Optionen.n_max;
        Kennfeld_WirkmechVerluste.U_n = KennfeldMaschine_WirkmechVerluste_Speichern.Entwurf.Bemessungswerte.U_N;
        
        % Drehzahlachse - X-Achse Kennfeld
        Kennfeld_WirkmechVerluste.eff_n_axis = transpose(KennfeldMaschine_WirkmechVerluste_Speichern.Analyse.Betriebsdaten.mot.n_m_vec);

        % Drehmoment - Y-Achse
        KennfeldMaschine_WirkmechVerluste_Speichern.Analyse.Betriebsdaten.mot.M_vec=transpose(fliplr(KennfeldMaschine_WirkmechVerluste_Speichern.Analyse.Betriebsdaten.mot.M_vec));
        Kennfeld_WirkmechVerluste.eff_T_axis = [(fliplr(transpose(KennfeldMaschine_WirkmechVerluste_Speichern.Analyse.Betriebsdaten.mot.M_vec(2:end)))) * (-1), transpose(KennfeldMaschine_WirkmechVerluste_Speichern.Analyse.Betriebsdaten.mot.M_vec)];
        
        % Diagrammbereich                                                             
        % Manipulate efficiency map
        KennfeldMaschine_WirkmechVerluste_Speichern.Analyse.Wirkungsgrad.eta_vme_mesh=flipud(KennfeldMaschine_WirkmechVerluste_Speichern.Analyse.Wirkungsgrad.eta_vme_mesh); %hier Unterschied Verlustarten
        KennfeldMaschine_WirkmechVerluste_Speichern.Analyse.Wirkungsgrad.eta_vme_mesh(KennfeldMaschine_WirkmechVerluste_Speichern.Analyse.Wirkungsgrad.eta_vme_mesh<0.1) = 0.1; %Hier Unterschied Verlustarten %set small efficiency values to 0.1
        KennfeldMaschine_WirkmechVerluste_Speichern.Analyse.Wirkungsgrad.eta_vme_mesh(1,:) = 0; %first colomn is zero %hier Unterschied Verlustarten
        KennfeldMaschine_WirkmechVerluste_Speichern.Analyse.Wirkungsgrad.eta_vme_mesh(:,1) = 0.1; %first row is zero %hier Unterschied Verlustarten
        
        % Create one "larger" eff. map that includes acceleration and recuperation
        Kennfeld_WirkmechVerluste.eff_recu = 1./KennfeldMaschine_WirkmechVerluste_Speichern.Analyse.Wirkungsgrad.eta_vme_mesh; %recuperation efficiency %hier Unterschied Verlustarten
        Kennfeld_WirkmechVerluste.eff_recu(Kennfeld_WirkmechVerluste.eff_recu==inf)=100;
        Kennfeld_WirkmechVerluste.eff_recu = Kennfeld_WirkmechVerluste.eff_recu(2:end,:);
        Kennfeld_WirkmechVerluste.eff = [flipud(Kennfeld_WirkmechVerluste.eff_recu); KennfeldMaschine_WirkmechVerluste_Speichern.Analyse.Wirkungsgrad.eta_vme_mesh]; %hier Unterschied Verlustarten
    
        % Volllastkennlinie
        KennfeldMaschine_WirkmechVerluste_Speichern.Analyse.Betriebsdaten.mot.M_max_vec=transpose(KennfeldMaschine_WirkmechVerluste_Speichern.Analyse.Betriebsdaten.mot.M_max_vec);
        Kennfeld_WirkmechVerluste.T_max = transpose(KennfeldMaschine_WirkmechVerluste_Speichern.Analyse.Betriebsdaten.mot.M_max_vec);  
        
        % Drehzahlachse skaliert
        Kennfeld_WirkmechVerluste.T_max_n_axis = transpose(KennfeldMaschine_WirkmechVerluste_Speichern.Analyse.Betriebsdaten.mot.M_vec);                                                          %noch auf 2000 interpolieren?
        
         save (['3_Ergebnisse/',handles.Entwurf.Optionen.folder_id,'/', '2_Analyse/', KennfeldMaschine_WirkmechVerluste_Speichern.Entwurf.Optionen.Maschinentyp,'_',num2str(KennfeldMaschine_WirkmechVerluste_Speichern.Entwurf.Bemessungswerte.P_N/1000),'_',num2str(KennfeldMaschine_WirkmechVerluste_Speichern.Entwurf.Bemessungswerte.n_N),'_',num2str(KennfeldMaschine_WirkmechVerluste_Speichern.Analyse.Optionen.n_max),'_',num2str(KennfeldMaschine_WirkmechVerluste_Speichern.Entwurf.Bemessungswerte.U_N),'_','Kennfeld_WirkmechVerluste','.mat'],'Kennfeld_WirkmechVerluste');
         save (['LDSimulation/Vehicle/Para_Powertrain/Motor_efficiency/MOT_memory//', KennfeldMaschine_WirkmechVerluste_Speichern.Entwurf.Optionen.Maschinentyp,'_',num2str(KennfeldMaschine_WirkmechVerluste_Speichern.Entwurf.Bemessungswerte.P_N/1000),'_',num2str(KennfeldMaschine_WirkmechVerluste_Speichern.Entwurf.Bemessungswerte.n_N),'_',num2str(KennfeldMaschine_WirkmechVerluste_Speichern.Analyse.Optionen.n_max),'_',num2str(KennfeldMaschine_WirkmechVerluste_Speichern.Entwurf.Bemessungswerte.U_N),'_', 'Kennfeld_WirkmechVerluste' ,'.mat'],'Kennfeld_WirkmechVerluste');
               
    elseif(strcmp(handles.Plot.Auswahl_Plot,'mechanische Verluste'))
        handles = plotKennfeld(handles);
        title(handles.axes_Ergebnis_Plot,'mechanische Verluste $P_{V,mech}$ in W','interpreter','latex','FontSize', 20);
        contourf(handles.axes_Ergebnis_Plot,handles.Analyse.Betriebsdaten.n_m_mesh,handles.Analyse.Betriebsdaten.M_max_mesh,handles.Analyse.Verluste.P_vme_mesh);
        datacursormode on
        
        % Save Kennfeld neu
        KennfeldMaschine_mechVerluste = struct;
        KennfeldMaschine_mechVerluste_Speichern = handles;
        KennfeldMaschine_mechVerluste.n = KennfeldMaschine_mechVerluste_Speichern.Analyse.Betriebsdaten.n_m_mesh;
        KennfeldMaschine_mechVerluste.M = KennfeldMaschine_mechVerluste_Speichern.Analyse.Betriebsdaten.M_max_mesh;
        KennfeldMaschine_mechVerluste.mechverl = KennfeldMaschine_mechVerluste_Speichern.Analyse.Verluste.P_vme_mesh;
         save (['3_Ergebnisse/', handles.Entwurf.Optionen.folder_id,'/','2_Analyse/','KennfeldMaschine_mechVerluste','.mat'],'KennfeldMaschine_mechVerluste');
        folder_id = [KennfeldMaschine_mechVerluste_Speichern.Entwurf.Optionen.folder_id]; %file_id, '_data'];
        file_id = [KennfeldMaschine_mechVerluste_Speichern.Entwurf.Optionen.file_id];
        
        % Konvertieren der Daten in richtige Form
        % Bemessungsgroessen
        Kennfeld_mechVerluste = struct;
        Kennfeld_mechVerluste.type = KennfeldMaschine_mechVerluste_Speichern.Entwurf.Optionen.Maschinentyp;
        Kennfeld_mechVerluste.power = KennfeldMaschine_mechVerluste_Speichern.Entwurf.Bemessungswerte.P_N/1000;
        Kennfeld_mechVerluste.n_n = KennfeldMaschine_mechVerluste_Speichern.Entwurf.Bemessungswerte.n_N;
        Kennfeld_mechVerluste.n_max = KennfeldMaschine_mechVerluste_Speichern.Analyse.Optionen.n_max;
        Kennfeld_mechVerluste.U_n = KennfeldMaschine_mechVerluste_Speichern.Entwurf.Bemessungswerte.U_N;
        
        % Drehzahlachse - X-Achse Kennfeld
        Kennfeld_mechVerluste.eff_n_axis = transpose(KennfeldMaschine_mechVerluste_Speichern.Analyse.Betriebsdaten.mot.n_m_vec);

        % Drehmoment - Y-Achse
        KennfeldMaschine_mechVerluste_Speichern.Analyse.Betriebsdaten.mot.M_vec=transpose(fliplr(KennfeldMaschine_mechVerluste_Speichern.Analyse.Betriebsdaten.mot.M_vec));
        Kennfeld_mechVerluste.eff_T_axis = [(fliplr(transpose(KennfeldMaschine_mechVerluste_Speichern.Analyse.Betriebsdaten.mot.M_vec(2:end)))) * (-1), transpose(KennfeldMaschine_mechVerluste_Speichern.Analyse.Betriebsdaten.mot.M_vec)];
        
        % Diagrammbereich                                                             
        % Manipulate efficiency map
        KennfeldMaschine_mechVerluste_Speichern.Analyse.Verluste.P_vme_mesh=flipud(KennfeldMaschine_mechVerluste_Speichern.Analyse.Verluste.P_vme_mesh); %hier Unterschied Verlustarten
        KennfeldMaschine_mechVerluste_Speichern.Analyse.Verluste.P_vme_mesh(KennfeldMaschine_mechVerluste_Speichern.Analyse.Verluste.P_vme_mesh<0.1) = 0.1; %Hier Unterschied Verlustarten %set small efficiency values to 0.1
        KennfeldMaschine_mechVerluste_Speichern.Analyse.Verluste.P_vme_mesh(1,:) = 0; %first colomn is zero %hier Unterschied Verlustarten
        KennfeldMaschine_mechVerluste_Speichern.Analyse.Verluste.P_vme_mesh(:,1) = 0.1; %first row is zero %hier Unterschied Verlustarten
        
        % Create one "larger" eff. map that includes acceleration and recuperation
        Kennfeld_mechVerluste.eff_recu = 1./KennfeldMaschine_mechVerluste_Speichern.Analyse.Verluste.P_vme_mesh; %recuperation efficiency %hier Unterschied Verlustarten
        Kennfeld_mechVerluste.eff_recu(Kennfeld_mechVerluste.eff_recu==inf)=100;
        Kennfeld_mechVerluste.eff_recu = Kennfeld_mechVerluste.eff_recu(2:end,:);
        Kennfeld_mechVerluste.eff = [flipud(Kennfeld_mechVerluste.eff_recu); KennfeldMaschine_mechVerluste_Speichern.Analyse.Verluste.P_vme_mesh]; %hier Unterschied Verlustarten
    
        % Volllastkennlinie
        KennfeldMaschine_mechVerluste_Speichern.Analyse.Betriebsdaten.mot.M_max_vec=transpose(KennfeldMaschine_mechVerluste_Speichern.Analyse.Betriebsdaten.mot.M_max_vec);
        Kennfeld_mechVerluste.T_max = transpose(KennfeldMaschine_mechVerluste_Speichern.Analyse.Betriebsdaten.mot.M_max_vec);  
        
        % Drehzahlachse skaliert
        Kennfeld_mechVerluste.T_max_n_axis = transpose(KennfeldMaschine_mechVerluste_Speichern.Analyse.Betriebsdaten.mot.M_vec);                                                          %noch auf 2000 interpolieren?

         save (['3_Ergebnisse/',handles.Entwurf.Optionen.folder_id,'/', '2_Analyse/', KennfeldMaschine_mechVerluste_Speichern.Entwurf.Optionen.Maschinentyp,'_',num2str(KennfeldMaschine_mechVerluste_Speichern.Entwurf.Bemessungswerte.P_N/1000),'_',num2str(KennfeldMaschine_mechVerluste_Speichern.Entwurf.Bemessungswerte.n_N),'_',num2str(KennfeldMaschine_mechVerluste_Speichern.Analyse.Optionen.n_max),'_',num2str(KennfeldMaschine_mechVerluste_Speichern.Entwurf.Bemessungswerte.U_N),'_','Kennfeld_mechVerluste','.mat'],'Kennfeld_mechVerluste');
         save (['LDSimulation/Vehicle/Para_Powertrain/Motor_efficiency/MOT_memory//', KennfeldMaschine_mechVerluste_Speichern.Entwurf.Optionen.Maschinentyp,'_',num2str(KennfeldMaschine_mechVerluste_Speichern.Entwurf.Bemessungswerte.P_N/1000),'_',num2str(KennfeldMaschine_mechVerluste_Speichern.Entwurf.Bemessungswerte.n_N),'_',num2str(KennfeldMaschine_mechVerluste_Speichern.Analyse.Optionen.n_max),'_',num2str(KennfeldMaschine_mechVerluste_Speichern.Entwurf.Bemessungswerte.U_N),'_', 'Kennfeld_mechVerluste' ,'.mat'],'Kennfeld_mechVerluste');
               
    elseif(strcmp(handles.Plot.Auswahl_Plot,'Wirkungsgrad Zusatzverluste'))
        handles = plotKennfeld(handles);
        title(handles.axes_Ergebnis_Plot,'Wirkungsgrad Zusatzverluste $\eta_{V,zus}$ in \%','interpreter','latex','FontSize', 20);
        Plot_Label_Vektor_eta = [0.7,0.75,0.8,0.84,0.86,0.88,0.9,0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.98,0.99,0.995];
        [C,h]=contourf(handles.axes_Ergebnis_Plot,handles.Analyse.Betriebsdaten.n_m_mesh,handles.Analyse.Betriebsdaten.M_max_mesh,handles.Analyse.Wirkungsgrad.eta_zus_mesh,Plot_Label_Vektor_eta);
        clabel(C,h,Plot_Label_Vektor_eta);
        datacursormode on
        
        % Save Kennfeld neu
        KennfeldMaschine_WirkZusVerl = struct;
        KennfeldMaschine_WirkZusVerl_Speichern = handles;
        KennfeldMaschine_WirkZusVerl.n = KennfeldMaschine_WirkZusVerl_Speichern.Analyse.Betriebsdaten.n_m_mesh;
        KennfeldMaschine_WirkZusVerl.M = KennfeldMaschine_WirkZusVerl_Speichern.Analyse.Betriebsdaten.M_max_mesh;
        KennfeldMaschine_WirkZusVerl.wirkzusverl = KennfeldMaschine_WirkZusVerl_Speichern.Analyse.Wirkungsgrad.eta_zus_mesh;
         save (['3_Ergebnisse/',handles.Entwurf.Optionen.folder_id,'/','2_Analyse/', 'KennfeldMaschine_WirkZusVerl','.mat'],'KennfeldMaschine_WirkZusVerl');
        folder_id = [KennfeldMaschine_WirkZusVerl_Speichern.Entwurf.Optionen.folder_id]; %file_id, '_data'];
        file_id = [KennfeldMaschine_WirkZusVerl_Speichern.Entwurf.Optionen.file_id];
        
        % Konvertieren der Daten in richtige Form
        % Bemessungsgroessen
        Kennfeld_WirkZusVerl = struct;
        Kennfeld_WirkZusVerl.type = KennfeldMaschine_WirkZusVerl_Speichern.Entwurf.Optionen.Maschinentyp;
        Kennfeld_WirkZusVerl.power = KennfeldMaschine_WirkZusVerl_Speichern.Entwurf.Bemessungswerte.P_N/1000;
        Kennfeld_WirkZusVerl.n_n = KennfeldMaschine_WirkZusVerl_Speichern.Entwurf.Bemessungswerte.n_N;
        Kennfeld_WirkZusVerl.n_max = KennfeldMaschine_WirkZusVerl_Speichern.Analyse.Optionen.n_max;
        Kennfeld_WirkZusVerl.U_n = KennfeldMaschine_WirkZusVerl_Speichern.Entwurf.Bemessungswerte.U_N;
        
        % Drehzahlachse - X-Achse Kennfeld
        Kennfeld_WirkZusVerl.eff_n_axis = transpose(KennfeldMaschine_WirkZusVerl_Speichern.Analyse.Betriebsdaten.mot.n_m_vec);

        % Drehmoment - Y-Achse
        KennfeldMaschine_WirkZusVerl_Speichern.Analyse.Betriebsdaten.mot.M_vec=transpose(fliplr(KennfeldMaschine_WirkZusVerl_Speichern.Analyse.Betriebsdaten.mot.M_vec));
        Kennfeld_WirkZusVerl.eff_T_axis = [(fliplr(transpose(KennfeldMaschine_WirkZusVerl_Speichern.Analyse.Betriebsdaten.mot.M_vec(2:end)))) * (-1), transpose(KennfeldMaschine_WirkZusVerl_Speichern.Analyse.Betriebsdaten.mot.M_vec)];
        
        % Diagrammbereich                                                             
        % Manipulate efficiency map
        KennfeldMaschine_WirkZusVerl_Speichern.Analyse.Wirkungsgrad.eta_zus_mesh=flipud(KennfeldMaschine_WirkZusVerl_Speichern.Analyse.Wirkungsgrad.eta_zus_mesh); %hier Unterschied Verlustarten
        KennfeldMaschine_WirkZusVerl_Speichern.Analyse.Wirkungsgrad.eta_zus_mesh(KennfeldMaschine_WirkZusVerl_Speichern.Analyse.Wirkungsgrad.eta_zus_mesh<0.1) = 0.1; %Hier Unterschied Verlustarten %set small efficiency values to 0.1
        KennfeldMaschine_WirkZusVerl_Speichern.Analyse.Wirkungsgrad.eta_zus_mesh(1,:) = 0; %first colomn is zero %hier Unterschied Verlustarten
        KennfeldMaschine_WirkZusVerl_Speichern.Analyse.Wirkungsgrad.eta_zus_mesh(:,1) = 0.1; %first row is zero %hier Unterschied Verlustarten
        
        % Create one "larger" eff. map that includes acceleration and recuperation
        Kennfeld_WirkZusVerl.eff_recu = 1./KennfeldMaschine_WirkZusVerl_Speichern.Analyse.Wirkungsgrad.eta_zus_mesh; %recuperation efficiency %hier Unterschied Verlustarten
        Kennfeld_WirkZusVerl.eff_recu(Kennfeld_WirkZusVerl.eff_recu==inf)=100;
        Kennfeld_WirkZusVerl.eff_recu = Kennfeld_WirkZusVerl.eff_recu(2:end,:);
        Kennfeld_WirkZusVerl.eff = [flipud(Kennfeld_WirkZusVerl.eff_recu); KennfeldMaschine_WirkZusVerl_Speichern.Analyse.Wirkungsgrad.eta_zus_mesh]; %hier Unterschied Verlustarten
    
        % Volllastkennlinie
        KennfeldMaschine_WirkZusVerl_Speichern.Analyse.Betriebsdaten.mot.M_max_vec=transpose(KennfeldMaschine_WirkZusVerl_Speichern.Analyse.Betriebsdaten.mot.M_max_vec);
        Kennfeld_WirkZusVerl.T_max = transpose(KennfeldMaschine_WirkZusVerl_Speichern.Analyse.Betriebsdaten.mot.M_max_vec);  
        
        % Drehzahlachse skaliert
        Kennfeld_WirkZusVerl.T_max_n_axis = transpose(KennfeldMaschine_WirkZusVerl_Speichern.Analyse.Betriebsdaten.mot.M_vec);                                                          %noch auf 2000 interpolieren?
        
         save (['3_Ergebnisse/',handles.Entwurf.Optionen.folder_id,'/','2_Analyse/', KennfeldMaschine_WirkZusVerl_Speichern.Entwurf.Optionen.Maschinentyp,'_',num2str(KennfeldMaschine_WirkZusVerl_Speichern.Entwurf.Bemessungswerte.P_N/1000),'_',num2str(KennfeldMaschine_WirkZusVerl_Speichern.Entwurf.Bemessungswerte.n_N),'_',num2str(KennfeldMaschine_WirkZusVerl_Speichern.Analyse.Optionen.n_max),'_',num2str(KennfeldMaschine_WirkZusVerl_Speichern.Entwurf.Bemessungswerte.U_N),'_','Kennfeld_WirkZusVerl','.mat'],'Kennfeld_WirkZusVerl');
         save (['LDSimulation/Vehicle/Para_Powertrain/Motor_efficiency/MOT_memory//', KennfeldMaschine_WirkZusVerl_Speichern.Entwurf.Optionen.Maschinentyp,'_',num2str(KennfeldMaschine_WirkZusVerl_Speichern.Entwurf.Bemessungswerte.P_N/1000),'_',num2str(KennfeldMaschine_WirkZusVerl_Speichern.Entwurf.Bemessungswerte.n_N),'_',num2str(KennfeldMaschine_WirkZusVerl_Speichern.Analyse.Optionen.n_max),'_',num2str(KennfeldMaschine_WirkZusVerl_Speichern.Entwurf.Bemessungswerte.U_N),'_', 'Kennfeld_WirkZusVerl' ,'.mat'],'Kennfeld_WirkZusVerl');
               
    elseif(strcmp(handles.Plot.Auswahl_Plot,'Zusatzverluste'))
        handles = plotKennfeld(handles);
        title(handles.axes_Ergebnis_Plot,'Zusatzverluste $P_{V,zus}$ in W','interpreter','latex','FontSize', 20);
        contourf(handles.axes_Ergebnis_Plot,handles.Analyse.Betriebsdaten.n_m_mesh,handles.Analyse.Betriebsdaten.M_max_mesh,handles.Analyse.Verluste.P_vzus_mesh);
        datacursormode on
        
        % Save Kennfeld neu
        KennfeldMaschine_ZusVerl = struct;
        KennfeldMaschine_ZusVerl_Speichern = handles;
        KennfeldMaschine_ZusVerl.n = KennfeldMaschine_ZusVerl_Speichern.Analyse.Betriebsdaten.n_m_mesh;
        KennfeldMaschine_ZusVerl.M = KennfeldMaschine_ZusVerl_Speichern.Analyse.Betriebsdaten.M_max_mesh;
        KennfeldMaschine_ZusVerl.zusverl = KennfeldMaschine_ZusVerl_Speichern.Analyse.Verluste.P_vzus_mesh;
         save (['3_Ergebnisse/', handles.Entwurf.Optionen.folder_id,'/','2_Analyse/', 'KennfeldMaschine_ZusVerl','.mat'],'KennfeldMaschine_ZusVerl');
        folder_id = [KennfeldMaschine_ZusVerl_Speichern.Entwurf.Optionen.folder_id]; %file_id, '_data'];
        file_id = [KennfeldMaschine_ZusVerl_Speichern.Entwurf.Optionen.file_id];
        
        % Konvertieren der Daten in richtige Form
        % Bemessungsgroessen
        Kennfeld_ZusVerl = struct;
        Kennfeld_ZusVerl.type = KennfeldMaschine_ZusVerl_Speichern.Entwurf.Optionen.Maschinentyp;
        Kennfeld_ZusVerl.power = KennfeldMaschine_ZusVerl_Speichern.Entwurf.Bemessungswerte.P_N/1000;
        Kennfeld_ZusVerl.n_n = KennfeldMaschine_ZusVerl_Speichern.Entwurf.Bemessungswerte.n_N;
        Kennfeld_ZusVerl.n_max = KennfeldMaschine_ZusVerl_Speichern.Analyse.Optionen.n_max;
        Kennfeld_ZusVerl.U_n = KennfeldMaschine_ZusVerl_Speichern.Entwurf.Bemessungswerte.U_N;
        
        % Drehzahlachse - X-Achse Kennfeld
        Kennfeld_ZusVerl.eff_n_axis = transpose(KennfeldMaschine_ZusVerl_Speichern.Analyse.Betriebsdaten.mot.n_m_vec);

        % Drehmoment - Y-Achse
        KennfeldMaschine_ZusVerl_Speichern.Analyse.Betriebsdaten.mot.M_vec=transpose(fliplr(KennfeldMaschine_ZusVerl_Speichern.Analyse.Betriebsdaten.mot.M_vec));
        Kennfeld_ZusVerl.eff_T_axis = [(fliplr(transpose(KennfeldMaschine_ZusVerl_Speichern.Analyse.Betriebsdaten.mot.M_vec(2:end)))) * (-1), transpose(KennfeldMaschine_ZusVerl_Speichern.Analyse.Betriebsdaten.mot.M_vec)];
        
        % Diagrammbereich                                                             
        % Manipulate efficiency map
        KennfeldMaschine_ZusVerl_Speichern.Analyse.Verluste.P_vzus_mesh=flipud(KennfeldMaschine_ZusVerl_Speichern.Analyse.Verluste.P_vzus_mesh); %hier Unterschied Verlustarten
        KennfeldMaschine_ZusVerl_Speichern.Analyse.Verluste.P_vzus_mesh(KennfeldMaschine_ZusVerl_Speichern.Analyse.Verluste.P_vzus_mesh<0.1) = 0.1; %Hier Unterschied Verlustarten %set small efficiency values to 0.1
        KennfeldMaschine_ZusVerl_Speichern.Analyse.Verluste.P_vzus_mesh(1,:) = 0; %first colomn is zero %hier Unterschied Verlustarten
        KennfeldMaschine_ZusVerl_Speichern.Analyse.Verluste.P_vzus_mesh(:,1) = 0.1; %first row is zero %hier Unterschied Verlustarten
        
        % Create one "larger" eff. map that includes acceleration and recuperation
        Kennfeld_ZusVerl.eff_recu = 1./KennfeldMaschine_ZusVerl_Speichern.Analyse.Verluste.P_vzus_mesh; %recuperation efficiency %hier Unterschied Verlustarten
        Kennfeld_ZusVerl.eff_recu(Kennfeld_ZusVerl.eff_recu==inf)=100;
        Kennfeld_ZusVerl.eff_recu = Kennfeld_ZusVerl.eff_recu(2:end,:);
        Kennfeld_ZusVerl.eff = [flipud(Kennfeld_ZusVerl.eff_recu); KennfeldMaschine_ZusVerl_Speichern.Analyse.Verluste.P_vzus_mesh]; %hier Unterschied Verlustarten
    
        % Volllastkennlinie
        KennfeldMaschine_ZusVerl_Speichern.Analyse.Betriebsdaten.mot.M_max_vec=transpose(KennfeldMaschine_ZusVerl_Speichern.Analyse.Betriebsdaten.mot.M_max_vec);
        Kennfeld_ZusVerl.T_max = transpose(KennfeldMaschine_ZusVerl_Speichern.Analyse.Betriebsdaten.mot.M_max_vec);  
        
        % Drehzahlachse skaliert
        Kennfeld_ZusVerl.T_max_n_axis = transpose(KennfeldMaschine_ZusVerl_Speichern.Analyse.Betriebsdaten.mot.M_vec);                                                          %noch auf 2000 interpolieren?

         save (['3_Ergebnisse/',handles.Entwurf.Optionen.folder_id,'/','2_Analyse/', KennfeldMaschine_ZusVerl_Speichern.Entwurf.Optionen.Maschinentyp,'_',num2str(KennfeldMaschine_ZusVerl_Speichern.Entwurf.Bemessungswerte.P_N/1000),'_',num2str(KennfeldMaschine_ZusVerl_Speichern.Entwurf.Bemessungswerte.n_N),'_',num2str(KennfeldMaschine_ZusVerl_Speichern.Analyse.Optionen.n_max),'_',num2str(KennfeldMaschine_ZusVerl_Speichern.Entwurf.Bemessungswerte.U_N),'_','Kennfeld_ZusVerl','.mat'],'Kennfeld_ZusVerl');
         save (['LDSimulation/Vehicle/Para_Powertrain/Motor_efficiency/MOT_memory//', KennfeldMaschine_ZusVerl_Speichern.Entwurf.Optionen.Maschinentyp,'_',num2str(KennfeldMaschine_ZusVerl_Speichern.Entwurf.Bemessungswerte.P_N/1000),'_',num2str(KennfeldMaschine_ZusVerl_Speichern.Entwurf.Bemessungswerte.n_N),'_',num2str(KennfeldMaschine_ZusVerl_Speichern.Analyse.Optionen.n_max),'_',num2str(KennfeldMaschine_ZusVerl_Speichern.Entwurf.Bemessungswerte.U_N),'_', 'Kennfeld_ZusVerl' ,'.mat'],'Kennfeld_ZusVerl');
                
    elseif(strcmp(handles.Plot.Auswahl_Plot,'i_1d'))
        handles = plotKennfeld(handles);
        title(handles.axes_Ergebnis_Plot,'$i_1d$ in A','interpreter','latex','FontSize', 20);
        contourf(handles.axes_Ergebnis_Plot,handles.Analyse.Betriebsdaten.n_m_mesh,handles.Analyse.Betriebsdaten.M_max_mesh,handles.Analyse.Momentensteuerung.i_1d_mesh);
        datacursormode on
        
        % Save Kennfeld
        KennfeldMaschine_i_ld = handles.Analyse.Momentensteuerung.i_1d_mesh;
         save (['3_Ergebnisse/',handles.Entwurf.Optionen.folder_id,'/','2_Analyse/', 'KennfeldMaschine_i_ld','.mat'],'KennfeldMaschine_i_ld');
        
    elseif(strcmp(handles.Plot.Auswahl_Plot,'i_1q'))
        handles = plotKennfeld(handles);
        title(handles.axes_Ergebnis_Plot,'$i_1q$ in A','interpreter','latex','FontSize', 20);
        contourf(handles.axes_Ergebnis_Plot,handles.Analyse.Betriebsdaten.n_m_mesh,handles.Analyse.Betriebsdaten.M_max_mesh,handles.Analyse.Momentensteuerung.i_1q_mesh);
        datacursormode on
        
        % Save Kennfeld
        KennfeldMaschine_i_lq = handles.Analyse.Momentensteuerung.i_1q_mesh;
         save (['3_Ergebnisse/',handles.Entwurf.Optionen.folder_id,'/','2_Analyse/', 'KennfeldMaschine_i_lq','.mat'],'KennfeldMaschine_i_lq');
        
    elseif(strcmp(handles.Plot.Auswahl_Plot,'u_1d'))
        handles = plotKennfeld(handles);
        title(handles.axes_Ergebnis_Plot,'$u_1d$ in V','interpreter','latex','FontSize', 20);
        contourf(handles.axes_Ergebnis_Plot,handles.Analyse.Betriebsdaten.n_m_mesh,handles.Analyse.Betriebsdaten.M_max_mesh,handles.Analyse.Momentensteuerung.u_1d_mesh);
        datacursormode on
        
        % Save Kennfeld
        KennfeldMaschine_u_ld = handles.Analyse.Momentensteuerung.u_1d_mesh;
         save (['3_Ergebnisse/',handles.Entwurf.Optionen.folder_id,'/','2_Analyse/', 'KennfeldMaschine_u_ld','.mat'],'KennfeldMaschine_u_ld');
        
    elseif(strcmp(handles.Plot.Auswahl_Plot,'u_1q'))
        handles = plotKennfeld(handles);
        title(handles.axes_Ergebnis_Plot,'$u_1q$ in V','interpreter','latex','FontSize', 20);
        contourf(handles.axes_Ergebnis_Plot,handles.Analyse.Betriebsdaten.n_m_mesh,handles.Analyse.Betriebsdaten.M_max_mesh,handles.Analyse.Momentensteuerung.u_1q_mesh);
        datacursormode on
        
        % Save Kennfeld
        KennfeldMaschine_u_lq = handles.Analyse.Momentensteuerung.u_1q_mesh;
         save (['3_Ergebnisse/',handles.Entwurf.Optionen.folder_id,'/','2_Analyse/', 'KennfeldMaschine_u_lq','.mat'],'KennfeldMaschine_u_lq');
        
    end

%% Save Analyse.mat
ZWiSP=handles.Analyse;
save (['3_Ergebnisse/',handles.Entwurf.Optionen.folder_id,'/','2_Analyse/', 'Analyse','.mat'],'ZWiSP');


% Plot Kennfeld
% function handles = plotKennfeld(handles)
%     hold(handles.axes_Ergebnis_Plot,'off');
%     handles.axes_Ergebnis_Plot.Visible = 'on';
%     plot(handles.axes_Ergebnis_Plot,handles.Analyse.Betriebsdaten.mot.n_m_vec,handles.Analyse.Betriebsdaten.mot.M_max_vec);
%     hold(handles.axes_Ergebnis_Plot,'on');
%     if(handles.opt.Generator)
%         plot(handles.axes_Ergebnis_Plot,handles.Analyse.Betriebsdaten.gen.n_m_vec,handles.Analyse.Betriebsdaten.gen.M_max_vec);
%     end
%     colorbar(handles.axes_Ergebnis_Plot,'EastOutside');
%     grid(handles.axes_Ergebnis_Plot,'on');
%     handles.axes_Ergebnis_Plot.Layer = 'top';
%     xlabel(handles.axes_Ergebnis_Plot,'Drehzahl n in $\frac{U}{min}$','interpreter','latex','FontSize', 15);
%     ylabel(handles.axes_Ergebnis_Plot,'Drehmoment M in Nm','interpreter','latex','FontSize', 15);
    
% Plot Ergebnisse
% function handles = plotErgebnisse(handles)
%     if(strcmp(handles.Plot.Auswahl_Plot,'Wirkungsgrad gesamt'))
%         handles = plotKennfeld(handles);
%         title(handles.axes_Ergebnis_Plot,'Wirkungsgrad gesamt $\eta_{V,ges}$ in \%','interpreter','latex','FontSize', 20);
%         Plot_Label_Vektor_eta = [0.7,0.75,0.8,0.84,0.86,0.88,0.9,0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.98,0.99,0.995];
%         [C,h]=contourf(handles.axes_Ergebnis_Plot,handles.Analyse.Betriebsdaten.n_m_mesh,handles.Analyse.Betriebsdaten.M_max_mesh,handles.Analyse.Wirkungsgrad.eta_ges_mesh,Plot_Label_Vektor_eta);
%         clabel(C,h,Plot_Label_Vektor_eta);
%         datacursormode on
%     elseif(strcmp(handles.Plot.Auswahl_Plot,'Verluste gesamt'))
%         handles = plotKennfeld(handles);
%         title(handles.axes_Ergebnis_Plot,'Verluste gesamt $P_{V,ges}$ in W','interpreter','latex','FontSize', 20);
%         contourf(handles.axes_Ergebnis_Plot,handles.Analyse.Betriebsdaten.n_m_mesh,handles.Analyse.Betriebsdaten.M_max_mesh,handles.Analyse.Verluste.P_vges_mesh);
%         datacursormode on
%     elseif(strcmp(handles.Plot.Auswahl_Plot,'Wirkungsgrad Wicklungsverluste'))
%         handles = plotKennfeld(handles);
%         title(handles.axes_Ergebnis_Plot,'Wirkungsgrad Wicklungsverluste $\eta_{V,w}$ in \%','interpreter','latex','FontSize', 20);
%         Plot_Label_Vektor_eta = [0.7,0.75,0.8,0.84,0.86,0.88,0.9,0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.98,0.99,0.995];
%         [C,h]=contourf(handles.axes_Ergebnis_Plot,handles.Analyse.Betriebsdaten.n_m_mesh,handles.Analyse.Betriebsdaten.M_max_mesh,handles.Analyse.Wirkungsgrad.eta_vw_mesh,Plot_Label_Vektor_eta);
%         clabel(C,h,Plot_Label_Vektor_eta);
%         datacursormode on
%     elseif(strcmp(handles.Plot.Auswahl_Plot,'Wicklungsverluste'))
%         handles = plotKennfeld(handles);
%         title(handles.axes_Ergebnis_Plot,'Wicklungsverluste $P_{V,w}$ in W','interpreter','latex','FontSize', 20);
%         contourf(handles.axes_Ergebnis_Plot,handles.Analyse.Betriebsdaten.n_m_mesh,handles.Analyse.Betriebsdaten.M_max_mesh,handles.Analyse.Verluste.P_vw_mesh);
%         datacursormode on
%     elseif(strcmp(handles.Plot.Auswahl_Plot,'Wirkungsgrad Ummagnetisierungsverluste'))
%         handles = plotKennfeld(handles);
%         title(handles.axes_Ergebnis_Plot,'Wirkungsgrad Ummagnetisierungsverluste $\eta_{V,u}$ in \%','interpreter','latex','FontSize', 20);
%         Plot_Label_Vektor_eta = [0.7,0.75,0.8,0.84,0.86,0.88,0.9,0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.98,0.99,0.995];
%         [C,h]=contourf(handles.axes_Ergebnis_Plot,handles.Analyse.Betriebsdaten.n_m_mesh,handles.Analyse.Betriebsdaten.M_max_mesh,handles.Analyse.Wirkungsgrad.eta_fe_mesh,Plot_Label_Vektor_eta);
%         clabel(C,h,Plot_Label_Vektor_eta);
%         datacursormode on
%     elseif(strcmp(handles.Plot.Auswahl_Plot,'Ummagnetisierungsverluste'))
%         handles = plotKennfeld(handles);
%         title(handles.axes_Ergebnis_Plot,'Ummagnetisierungsverluste $P_{V,u}$ in W','interpreter','latex','FontSize', 20);
%         contourf(handles.axes_Ergebnis_Plot,handles.Analyse.Betriebsdaten.n_m_mesh,handles.Analyse.Betriebsdaten.M_max_mesh,handles.Analyse.Verluste.P_vu_mesh);
%         datacursormode on
%     elseif(strcmp(handles.Plot.Auswahl_Plot,'Wirkungsgrad mechanische Verluste'))
%         handles = plotKennfeld(handles);
%         title(handles.axes_Ergebnis_Plot,'Wirkungsgrad mechanische Verluste $\eta_{V,mech}$ in \%','interpreter','latex','FontSize', 20);
%         Plot_Label_Vektor_eta = [0.7,0.75,0.8,0.84,0.86,0.88,0.9,0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.98,0.99,0.995];
%         [C,h]=contourf(handles.axes_Ergebnis_Plot,handles.Analyse.Betriebsdaten.n_m_mesh,handles.Analyse.Betriebsdaten.M_max_mesh,handles.Analyse.Wirkungsgrad.eta_vme_mesh,Plot_Label_Vektor_eta);
%         clabel(C,h,Plot_Label_Vektor_eta);
%         datacursormode on
%     elseif(strcmp(handles.Plot.Auswahl_Plot,'mechanische Verluste'))
%         handles = plotKennfeld(handles);
%         title(handles.axes_Ergebnis_Plot,'mechanische Verluste $P_{V,mech}$ in W','interpreter','latex','FontSize', 20);
%         contourf(handles.axes_Ergebnis_Plot,handles.Analyse.Betriebsdaten.n_m_mesh,handles.Analyse.Betriebsdaten.M_max_mesh,handles.Analyse.Verluste.P_vme_mesh);
%         datacursormode on
%     elseif(strcmp(handles.Plot.Auswahl_Plot,'Wirkungsgrad Zusatzverluste'))
%         handles = plotKennfeld(handles);
%         title(handles.axes_Ergebnis_Plot,'Wirkungsgrad Zusatzverluste $\eta_{V,zus}$ in \%','interpreter','latex','FontSize', 20);
%         Plot_Label_Vektor_eta = [0.7,0.75,0.8,0.84,0.86,0.88,0.9,0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.98,0.99,0.995];
%         [C,h]=contourf(handles.axes_Ergebnis_Plot,handles.Analyse.Betriebsdaten.n_m_mesh,handles.Analyse.Betriebsdaten.M_max_mesh,handles.Analyse.Wirkungsgrad.eta_zus_mesh,Plot_Label_Vektor_eta);
%         clabel(C,h,Plot_Label_Vektor_eta);
%         datacursormode on
%     elseif(strcmp(handles.Plot.Auswahl_Plot,'Zusatzverluste'))
%         handles = plotKennfeld(handles);
%         title(handles.axes_Ergebnis_Plot,'Zusatzverluste $P_{V,zus}$ in W','interpreter','latex','FontSize', 20);
%         contourf(handles.axes_Ergebnis_Plot,handles.Analyse.Betriebsdaten.n_m_mesh,handles.Analyse.Betriebsdaten.M_max_mesh,handles.Analyse.Verluste.P_vzus_mesh);
%         datacursormode on
%     elseif(strcmp(handles.Plot.Auswahl_Plot,'i_1d'))
%         handles = plotKennfeld(handles);
%         title(handles.axes_Ergebnis_Plot,'$i_1d$ in A','interpreter','latex','FontSize', 20);
%         contourf(handles.axes_Ergebnis_Plot,handles.Analyse.Betriebsdaten.n_m_mesh,handles.Analyse.Betriebsdaten.M_max_mesh,handles.Analyse.Momentensteuerung.i_1d_mesh);
%         datacursormode on
%     elseif(strcmp(handles.Plot.Auswahl_Plot,'i_1q'))
%         handles = plotKennfeld(handles);
%         title(handles.axes_Ergebnis_Plot,'$i_1q$ in A','interpreter','latex','FontSize', 20);
%         contourf(handles.axes_Ergebnis_Plot,handles.Analyse.Betriebsdaten.n_m_mesh,handles.Analyse.Betriebsdaten.M_max_mesh,handles.Analyse.Momentensteuerung.i_1q_mesh);
%         datacursormode on
%     elseif(strcmp(handles.Plot.Auswahl_Plot,'u_1d'))
%         handles = plotKennfeld(handles);
%         title(handles.axes_Ergebnis_Plot,'$u_1d$ in V','interpreter','latex','FontSize', 20);
%         contourf(handles.axes_Ergebnis_Plot,handles.Analyse.Betriebsdaten.n_m_mesh,handles.Analyse.Betriebsdaten.M_max_mesh,handles.Analyse.Momentensteuerung.u_1d_mesh);
%         datacursormode on
%     elseif(strcmp(handles.Plot.Auswahl_Plot,'u_1q'))
%         handles = plotKennfeld(handles);
%         title(handles.axes_Ergebnis_Plot,'$u_1q$ in V','interpreter','latex','FontSize', 20);
%         contourf(handles.axes_Ergebnis_Plot,handles.Analyse.Betriebsdaten.n_m_mesh,handles.Analyse.Betriebsdaten.M_max_mesh,handles.Analyse.Momentensteuerung.u_1q_mesh);
%         datacursormode on
%     end
