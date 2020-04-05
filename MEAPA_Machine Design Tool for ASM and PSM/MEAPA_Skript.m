% -------------------------------------------------------------------------
% TU Muenchen - Lehrstuhl fuer Fahrzeugtechnik (FTM)
% -------------------------------------------------------------------------
% Modell fuer den Entwurf und die Analyse einer PMSM oder ASM (MEAPA)
% -------------------------------------------------------------------------
% Autor: Svenja Kalt (kalt@ftm.mw.tum.de)
%        Jonathan Erhard
% -------------------------------------------------------------------------

% Hinweise zur Verwendung dieses Skripts:
% (1) Die Eingangsparameter koennen beliebig festgelegt werden, z.B. kann
%     rated dem Skript uebergeben werden, entsprechend muss der rated-Teil
%     unten auskommentiert werden. Es muss allerdings nicht zwingend etwas
%     uebergeben werden.
% (2) Wird das Skript ohne die Main-Datei benutzt, muss darauf geachtet 
%     werden, dass die Pfade korrekt gesetzt werden.
% (3) Da keine Benutzeroberflaeche zur Auswahl der Wicklung aufgerufen
%     wird, kann nur der klassische Entwurf durchgefuehrt werden (keine
%     Benutzerinteraktion notwendig)

function [Entwurf, Analyse] = MEAPA_Skript(rated)

%% Maschinentyp auswaehlen
opt.Maschinentyp = 'ASM';                                                  % 'ASM', 'PMSM'

%% Eingabeparameter Entwurf
if(strcmp(opt.Maschinentyp,'ASM'))
%     %% BEMESSUNGSWERTE ASM
%     % Nennleistung P_N [W]
%     rated.P_N   = 3000;
%     
%     % Nenndrehzahl n_N [U/min]
%     rated.n_N   = 1500;
%     
%     % Nennspannung U_N [V]
%     rated.U_N   = 400;
%     
%     % Polpaarzahl p [-]
%     rated.p     = 2;
%     
%     % Nennfrequenz f_N [-]
%     rated.f_N   = (rated.p * rated.n_N) / 60;
%     
%     % Strangzahl m [-]
%     rated.m     = 3;

    %% OPTIONEN ASM
    % Maschinenausfuehrung
    opt.Maschinenausfuehrung            = 'Kaefiglaeufer';                 % 'Kaefiglaeufer'
    
    % Schaltung
    opt.Schaltung                       = 'Dreieck';                       % 'Stern', 'Dreieck'
    
    % Spulenform Stator
    opt.Spulenform_Stator               = 'Runddraht';                     % 'Runddraht'
    
    % Spulenform Rotor
    opt.Spulenform_Rotor                = 'Runddraht';                     % 'Runddraht'
    
    % Nutform Stator
    opt.Nutform_Stator                  = 'Trapezform (eckig)';            % 'Trapezform (eckig)'
    
    % Nutform Rotor
    opt.Nutform_Rotor                   = 'Trapezform (eckig)';            % 'Trapezform (eckig)'
    
    % Kuehlung Stator
    opt.Kuehlungsart                    = 'Oberflaechenkuehlung';          % 'Oberflaechenkuehlung', 'Innen- oder Kreislaufkuehlung'
    
    % Eisenmaterial Stator
    opt.Stator_Eisenmaterial.String     = 'VACOFLUX 48';                   % 'M250-35A', 'M800-50A', 'VACOFLUX 48', 'VACOFLUX 50'
    opt.Stator_Eisenmaterial = loadMaterial(opt.Stator_Eisenmaterial,'Elektroblech');
    
    % Leitermaterial Stator
    opt.Stator_Leitermaterial.String    = 'Kupfer';                        % 'Aluminiumdraht', 'Aluminiumguss', 'Kupfer'
    opt.Stator_Leitermaterial = loadMaterial(opt.Stator_Leitermaterial,'Leiter');
    
    % Temperatur Leitermaterial Stator [°C]
    opt.theta_1                         = 90;
    
    % Eisenmaterial Rotor
    opt.Rotor_Eisenmaterial.String      = 'VACOFLUX 48';                   % 'M250-35A', 'M800-50A', 'VACOFLUX 48', 'VACOFLUX 50'
    opt.Rotor_Eisenmaterial = loadMaterial(opt.Stator_Eisenmaterial,'Elektroblech');
    
    % Leitermaterial Rotor
    opt.Rotor_Leitermaterial.String     = 'Aluminiumguss';                 % 'Aluminiumdraht', 'Aluminiumguss', 'Kupfer'
    opt.Rotor_Leitermaterial = loadMaterial(opt.Rotor_Leitermaterial,'Leiter');
    
    % Temperatur Leitermaterial Rotor [°C]
    opt.theta_2                         = 115;
    
    % Modus Wicklungsauslegung
    opt.Mode_Wicklung                   = 'Klassisch';                     % 'Klassisch'
    
    % Optimierungsziel Wicklung
    opt.Wicklungstyp                    = 'A';                             % 'A','B','C'

    %% RICHTWERTE ASM
    % Richtwert fuer die relative Ankerlaenge lambda [-]
    % Quelle: [Mueller08, S.577 - Tabelle 9.1.3], [Meyer18, S.117]
    richt.lambda = 2.5;                                                   % zwischen 0.6 und 1.0 fuer p=1, zwischen 1.0 und 4.0 fuer p>1

    % Richtwert fuer die Kanalbreite der Ventilationskanaele l_v [m]
    % Quelle: [Meyer09, S.41], [Mueller08, S.585]
    richt.l_v = 0.01;                                                      % zwischen 0.006 und 0.01

    % Richtwert fuer den Mittelwert der Luftspaltinduktion B_m [T]
    % Quelle: [Mueller08, S.582 - Tabelle 9.1.5]
    richt.B_m = 0.58;                                                      % zwischen 0.4 und 0.65

    % Richtwert fuer den Strombelag A [A/mm]
    % Quelle: [Mueller08, S.580 - Tabelle 9.1.4]
    % ACHTUNG: entweder B_m oder A angeben
    %richt.A = 50.0;                                                       % zwischen 20.0 und 120.0

    % Richtwert fuer die Stromdichte (Stator) S_1 [A/mm^2]
    % Quelle: [Mueller08, S.580 - Tabelle 9.1.4]
    richt.S_1 = 7.0;                                                       % zwischen 3.0 und 8.0

    % Richtwert fuer die max. zulaessige Induktion im Ruecken (Stator) B_1r_max [T]
    % Quelle: [Mueller08, S.582 - Tabelle 9.1.5]
    richt.B_1r_max = 1.4;                                                  % zwischen 1.3 und 1.65

    % Richtwert fuer die max. zulaessige Induktion in den Zaehnen (Stator) B_1z_max [T]
    % Quelle: [Mueller08, S.582 - Tabelle 9.1.5]
    richt.B_1z_max = 1.8;                                                 % zwischen 1.4 und 2.1

    % Richtwert fuer den Nutfuellfaktor (Stator) phi_1n [-]
    % Quelle: [Mueller08, S.586 - Tabelle 9.1.6]
    % V/A: Niederspannung
    richt.phi_1n = 0.5;                                                    % zwischen 0.3 und 0.5 fuer Runddraht, zwischen 0.35 und 0.6 fuer Formspule- oder Stab

    % Richtwert fuer den Wicklungsfaktor (Stator) xi_1p [-]
    % Quelle: [Mueller08, S.596]
    % V/A: keine Unterscheidung zwischen Einschicht- und Zweischichtwicklung,
    % da lediglich fuer initiale Abschaetzung benoetigt (wird im weiteren
    % Entwurfsverlauf genau berechnet)
    richt.xi_1p = 0.96;                                                    % zwischen 0.92 und 0.96

    % Richtwert fuer die minimale Nutteilung (Stator) tau_1n_min [m]
    % Quelle: [Meyer09, S.46], [Pyr14]
    richt.tau_1n_min = 0.007;                                              % zwischen 0.007 und 0.07

    % Richtwert fuer den Eisenfuellfaktor (Stator) phi_1Fe [-]
    % Quelle: [Mueller08, S.599]
    richt.phi_1Fe = 0.95;                                                  % zwischen 0.9 und 1.0

    % Richtwert fuer die Stromdichte (Rotor) S_2 [A/mm^2]
    % Quelle: [Mueller08, S.580 - Tabelle 9.1.4]
    richt.S_2 = 5.0;                                                       % zwischen 3.0 und 8.0 fuer Kupfer, zwischen 3.0 und 6.5 fuer Aluminiumguss, zwischen 3.0 und 6.5 fuer Aluminiumdraht

    % Richtwert fuer die Stabstromdichte (Rotor) S_2s [A/mm^2]
    % Quelle: [Mueller08, S.580 - Tabelle 9.1.4]
    richt.S_2s = 4.0;                                                      % zwischen 3.0 und 8.0 fuer Kupfer, zwischen 3.0 und 6.5 fuer Aluminiumguss, zwischen 3.0 und 6.5 fuer Aluminiumdraht

    % Richtwert fuer die Ringstromdichte (Rotor) S_2r [A/mm^2]
    % Quelle: [Mueller08, S.580 - Tabelle 9.1.4]
    richt.S_2r = 5.0;                                                      % zwischen 3.0 und 8.0 fuer Kupfer, zwischen 3.0 und 6.5 fuer Aluminiumguss, zwischen 3.0 und 6.5 fuer Aluminiumdraht

    % Richtwert fuer die max. zulaessige Induktion im Ruecken (Rotor) B_2r_max [T]
    % Quelle: [Mueller08, S.582 - Tabelle 9.1.5]
    richt.B_2r_max = 1.5;                                                  % zwischen 0.4 und 1.6

    % Richtwert fuer die max. zulaessige Induktion in den Zaehnen (Rotor) B_2z_max [T]
    % Quelle: [Mueller08, S.582 - Tabelle 9.1.5]
    richt.B_2z_max = 1.9;                                                  % zwischen 1.5 und 2.2

    % Richtwert fuer den Nutfuellfaktor (Rotor) phi_2n [-]
    % Quelle: [Mueller08, S.586 - Tabelle 9.1.6]
    % V/A: Niederspannung
    richt.phi_2n = 0.5;                                                    % zwischen 0.3 und 0.5 fuer Runddraht, zwischen 0.35 und 0.6 fuer Formspule- oder Stab

    % Richtwert fuer den Wicklungsfaktor (Rotor) xi_2p [-]
    % Quelle: [Mueller08, S.596]
    % V/A: keine Unterscheidung zwischen Einschicht- und Zweischichtwicklung,
    % da lediglich fuer initiale Abschaetzung benoetigt (wird im weiteren
    % Entwurfsverlauf genau berechnet)
    richt.xi_2p = 0.96;                                                    % zwischen 0.92 und 0.96

    % Richtwert fuer die minimale Nutteilung (Rotor) tau_2n_min [m]
    % Quelle: [Meyer09, S.46], [Pyr14]
    richt.tau_2n_min = 0.007;                                              % zwischen 0.007 und 0.07

    % Richtwert fuer den Eisenfuellfaktor (Rotor) phi_2Fe [-]
    % Quelle: [Mueller08, S.599]
    richt.phi_2Fe = 0.95;                                                  % zwischen 0.9 und 1.0

elseif(strcmp(opt.Maschinentyp,'PMSM'))
%     %% BEMESSUNGSWERTE PMSM
%     % Nennleistung P_N [W]
%     rated.P_N           = 10000;
%     
%     % Nenndrehzahl n_N [U/min]
%     rated.n_N           = 3000;
%     
%     % Nennspannung U_N [V]
%     rated.U_N           = 400;
%     
%     % Polpaarzahl p [-]
%     rated.p             = 2;
%     
%     % Nennfrequenz f_N [-]
%     rated.f_N           = (rated.p * rated.n_N) / 60;
%     
%     % Leistungsfaktor cos_phi_N [-]
%     rated.cos_phi_N     = 0.9;
%     
%     % Strangzahl m [-]
%     rated.m             = 3;

    %% OPTIONEN PMSM
    % Maschinenausfuehrung
    opt.Maschinenausfuehrung            = 'SPMSM';                         % 'SPMSM', 'IPMSM (eingelassen)', 'IPMSM (tangential)', 'IPMSM (V-Form)';
    
    % Schaltung
    opt.Schaltung                       = 'Stern';                         % 'Stern', 'Dreieck'
    
    % Spulenform Stator
    opt.Spulenform_Stator               = 'Runddraht';                     % 'Runddraht'
    
    % Nutform Stator
    opt.Nutform_Stator                  = 'Trapezform (eckig)';            % 'Trapezform (eckig)'
    
    % Kuehlung Stator
    opt.Kuehlungsart                    = 'Luft (indirekt)';               % 'Luft (indirekt)', 'Wasser (direkt)'
    
    % Eisenmaterial Stator
    opt.Stator_Eisenmaterial.String     = 'VACOFLUX 50';                   % 'M250-35A', 'M800-50A', 'VACOFLUX 48', 'VACOFLUX 50'
    opt.Stator_Eisenmaterial = loadMaterial(opt.Stator_Eisenmaterial,'Elektroblech');
    
    % Leitermaterial Stator
    opt.Stator_Leitermaterial.String    = 'Kupfer';                        % 'Aluminiumdraht', 'Aluminiumguss', 'Kupfer'
    opt.Stator_Leitermaterial = loadMaterial(opt.Stator_Leitermaterial,'Leiter');
    
    % Temperatur Leitermaterial Stator [°C]
    opt.theta_1                         = 90;
    
    % Eisenmaterial Rotor
    opt.Rotor_Eisenmaterial.String      = 'VACOFLUX 50';                   % 'M250-35A', 'M800-50A', 'VACOFLUX 48', 'VACOFLUX 50'
    opt.Rotor_Eisenmaterial = loadMaterial(opt.Stator_Eisenmaterial,'Elektroblech');
    
    % Magnetmaterial Rotor
    opt.Rotor_Magnetmaterial.String     = 'VACODYM 238 TP';                % 'VACODYM 238 TP', 'VACODYM 225 TP'
    opt.Rotor_Magnetmaterial = loadMaterial(opt.Rotor_Magnetmaterial,'Magnet');
    
    % Modus Wicklungsauslegung
    opt.Mode_Wicklung                   = 'Klassisch';                     % 'Klassisch'
    
    % Optimierungsziel Wicklung
    opt.Wicklungstyp                    = 'A';                             % 'A','B','C'

    %% RICHTWERTE PMSM
    % Richtwert fuer die relative Ankerlaenge lambda [-]
    % Quelle: [Mueller08, S.577 - Tabelle 9.1.3], [Meyer18, S.117]
    richt.lambda = 2.5;                                                    % zwischen 0.6 und 1.0 fuer p=1, zwischen 1.0 und 4.0 fuer p>1

    % Richtwert fuer die Kanalbreite der Ventilationskanaele l_v [m]
    % Quelle: [Meyer09, S.41], [Mueller08, S.585]
    richt.l_v = 0.01;                                                      % zwischen 0.006 und 0.01

    % Richtwert fuer die Amplitude der Luftspaltinduktion B_p [T]
    % Quelle: [Mueller08, S.582 - Tabelle 9.1.5]
    richt.B_p = 0.85;                                                      % zwischen 0.75 und 1.05

    % Richtwert fuer den Strombelag A [A/mm]
    % Quelle: [Mueller08, S.580 - Tabelle 9.1.4]
    % ACHTUNG: entweder B_p oder A angeben
    %richt.A = 50.0;                                                       % zwischen 30.0 und 120.0 fuer Luft (indirekt), zwischen 160.0 und 300.0 fuer Wasser (direkt)

    % Richtwert fuer die Stromdichte (Stator) S_1 [A/mm^2]
    % Quelle: [Mueller08, S.580 - Tabelle 9.1.4]
    richt.S_1 = 7.0;                                                       % zwischen 3.0 und 7.0 fuer Luft (indirekt), zwischen 13.0 und 18.0 fuer Wasser (direkt)

    % Richtwert fuer die max. zulaessige Induktion im Ruecken (Stator) B_1r_max [T]
    % Quelle: [Mueller08, S.582 - Tabelle 9.1.5]
    richt.B_1r_max = 1.4;                                                  % zwischen 1.0 und 1.5

    % Richtwert fuer die max. zulaessige Induktion in den Zaehnen (Stator) B_1z_max [T]
    % Quelle: [Mueller08, S.582 - Tabelle 9.1.5]
    richt.B_1z_max = 1.8;                                                  % zwischen 1.6 und 2.0

    % Richtwert fuer den Nutfuellfaktor (Stator) phi_1n [-]
    % Quelle: [Mueller08, S.586 - Tabelle 9.1.6]
    % V/A: Niederspannung
    richt.phi_1n = 0.5;                                                    % zwischen 0.3 und 0.5 fuer Runddraht, zwischen 0.35 und 0.6 fuer Formspule- oder Stab

    % Richtwert fuer den Wicklungsfaktor (Stator) xi_1p [-]
    % Quelle: [Mueller08, S.596]
    % V/A: keine Unterscheidung zwischen Einschicht- und Zweischichtwicklung,
    % da lediglich fuer initiale Abschaetzung benoetigt (wird im weiteren
    % Entwurfsverlauf genau berechnet)
    richt.xi_1p = 0.96;                                                    % zwischen 0.92 und 0.96

    % Richtwert fuer die minimale Nutteilung (Stator) tau_1n_min [m]
    % Quelle: [Meyer09, S.46], [Pyr14]
    richt.tau_1n_min = 0.007;                                              % zwischen 0.007 und 0.07

    % Richtwert fuer den Eisenfuellfaktor (Stator) phi_1Fe [-]
    % Quelle: [Mueller08, S.599]
    richt.phi_1Fe = 0.95;                                                  % zwischen 0.9 und 1.0

    % Richtwert fuer die max. zulaessige Induktion im Ruecken (Rotor) B_2r_max [T]
    % Quelle: [Mueller08, S.582 - Tabelle 9.1.5]
    richt.B_2r_max = 1.4;                                                  % zwischen 1.0 und 1.5

    % Richtwert fuer den Eisenfuellfaktor (Rotor) phi_2Fe [-]
    % Quelle: [Mueller08, S.599]
    richt.phi_2Fe = 0.95;                                                  % zwischen 0.9 und 1.0
else
    error('Ungueltige Eingabe bei Variable "opt.Maschinentyp"');
end

%% Umspeichern
handles.rated = rated;
handles.richt = richt;
handles.opt = opt;
clear rated richt opt

%% Entwurf starten
if(strcmp(handles.opt.Maschinentyp,'ASM'))
    [handles.Entwurf] = Entwurf_ASM(handles);
elseif(strcmp(handles.opt.Maschinentyp,'PMSM'))
    [handles.Entwurf] = Entwurf_PMSM(handles);
else
    error('Ungueltige Eingabe bei Variable "Entwurf.Optionen.Maschinentyp"');
end
disp('Entwurf abgeschlossen');

%% Eingabeparameter Analyse
    %% OPTIONEN VERLUSTE
    % Wicklungsverluste
    handles.opt.P_vw = 1;                                                  % 0, 1
    
    % Ummagnetisierungsverluste
    handles.opt.P_vu = 1;                                                  % 0, 1
    
    % Eisenverlustmodell
    handles.opt.P_vu_Modell = '1';                                         % '1'
    
    % mechanische Verluste
    handles.opt.P_vme = 1;                                                 % 0, 1
    
    % Zusatzverluste
    handles.opt.P_vzus = 1;                                                % 0, 1

    %% OPTIONEN BERECHNUNG
    % Generatorbereich
    handles.opt.Generator = 0;                                             % 0, 1
    
    % max. Drehzahl [U/min]
    handles.opt.n_max = 10000;
    
    % Aufloesung Drehzahl
    handles.opt.n_tics = 60;
    
    % Aufloesung Drehmoment
    handles.opt.M_tics = 60;
    
    % Ansteuerung max. Spannung [V]
    handles.opt.u_1max = handles.Entwurf.EMAG.U_1Str * sqrt(2);
    
    % Ansteuerung max. Strom [A]
    handles.opt.i_1max = handles.Entwurf.EMAG.I_1Str * sqrt(2);

%% Analyse starten
if(strcmp(handles.Entwurf.Optionen.Maschinentyp,'ASM'))
    [handles.Analyse] = Analyse_ASM(handles);
elseif(strcmp(handles.Entwurf.Optionen.Maschinentyp,'PMSM'))
    [handles.Analyse] = Analyse_PMSM(handles);
else
    error('Ungueltige Eingabe bei Variable "handles.Entwurf.Optionen.Maschinentyp"');
end
disp('Analyse abgeschlossen');

Entwurf = handles.Entwurf;
Analyse = handles.Analyse;

end

%% Zusatzfunktion
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
end
