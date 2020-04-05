% -------------------------------------------------------------------------
% TU Muenchen - Lehrstuhl fuer Fahrzeugtechnik (FTM)
% -------------------------------------------------------------------------
% Modell fuer den Entwurf und die Analyse einer PMSM oder ASM (MEAPA)
% -------------------------------------------------------------------------
% Autor: Svenja Kalt (kalt@ftm.mw.tum.de)
%        Jonathan Erhard
% -------------------------------------------------------------------------

%% Inhaltsverzeichnis
% A) Preprocessing
% B) Ermittlung der Hauptabmessungen
%   B.1) Statorinnendurchmesser
%   B.2) Geometrischer Luftspalt und Rotoraussendurchmesser
%   B.3) Axiale Laenge der aktiven Teile
% C) Statorentwurf
%   C.1) Wicklungsschema
%   C.2) Verschaltung
%   C.3) Nutform und magnetischer Kreis
% D) Rotorentwurf
%   D.1) Wicklungsschema
%   D.2) Nutform und magnetischer Kreis
% E) Nachrechnung
%   E.1) Nachrechnung des magnetischen Kreises
%   E.2) Induktivitaeten
%   E.3) Widerstaende
% F) Postprocessing
% G) Hilfsfunktionen

function [Entwurf] = Entwurf_ASM(handles)
%Entwurf Diese Funktion legt eine ASM aus.
%   In der Funktion Entwurf wird eine ASM entsprechend der
%   Bemessungsgroessen ausgelegt. Vereinfachungen oder Annahmen, die
%   getroffen werden, sind entsprechend gekennzeichnet.

%% A) Preprocessing
% #########################################################################
% #   A) PREPROCESSING                                                    #
% #########################################################################

warning('off','backtrace')
% Umspeichern der Variablen zur einfacheren Nutzung
rated = handles.rated;
richt = handles.richt;
opt = handles.opt;

% Allgemeine Konstanten
% Magnetische Feldkonstante mu_0 [H/m]
const.mu_0 = pi * 4e-7;

%% B) Ermittlung der Hauptabmessungen
% #########################################################################
% #   B) ERMITTLUNG DER HAUPTABMESSUNGEN                                  #
% #########################################################################

% Berechnung des Nenndrehmoments M_n [Nm]
rated.M_N = rated.P_N * 60 / (2 * pi * rated.n_N);

% Berechnung der Strangspannung (Stator) U_1Str [V]
% Quelle: [Mueller14, S.33 - Formel 0.6.5], [Mueller14, S.35 - Formel 0.6.6]
if(strcmp(opt.Schaltung,'Dreieck'))
    emag.U_1Str = rated.U_N;
elseif(strcmp(opt.Schaltung,'Stern'))
    emag.U_1Str = rated.U_N / sqrt(3);
else
    error('Ungueltige Eingabe bei Variable "opt.Schaltung"');
end

% Abschaetzung der induzierten Spannung E_h [V] aus Erfahrungswerten
% Quelle: [Mueller08, S.595]
if(rated.P_N>(30*1e3))
    emag.E_h = emag.U_1Str;
else
    emag.E_h = emag.U_1Str * 0.95; % zwischen 0.92 und 0.96
end

% Abschaetzung des Nennwirkungsgrads eta_N [-] aus Erfahrungswerten
% Quelle: [Mueller08, S.569 - Bild 9.1.4]
% V/A: eta_N_max zur Abschaetzung heranziehen [Mueller08, S.565]
% eta_N_data = [P_mech {kW}; eta_N_max {%}; eta_N_min {%}]
eta_N_data = [ 1.0  2.0  4.0  6.0  8.0 10.0 20.0 40.0 60.0 80.0 100.0 200.0 400.0 600.0 800.0 1000.0 2000.0 4000.0 5000.0; ...
              80.0 83.0 85.7 87.2 88.2 89.0 91.2 92.7 93.5 93.9  94.2  95.0  95.6  95.8  96.0   96.1   96.3   96.4   96.4; ...
              67.3 71.5 75.8 78.0 79.6 80.6 84.0 86.8 88.0 88.8  89.2  90.7  91.7  92.4  92.8   93.0   93.8   94.5   94.7];
emag.eta_N = interp1(eta_N_data(1,:), eta_N_data(2,:),(rated.P_N/1e3),'linear')/100;
if(any(isnan(emag.eta_N)))
    error('Datenbereich der hinterlegten Kurve zur Abschaetzung des Nennwirkungsgrads zu klein. \n(keine Werte fuer den Nennwirkungsgrad eta_N bei einer Nenneistung von P_N = %.2f W hinterlegt)',rated.P_N(1,end));
end
clear eta_N_data

% Abschaetzung des Nennleistungsfaktors cos_phi_N [-] aus Erfahrungswerten
% Quelle: [Mueller08, S.570 - Bild 9.1.5]
% V/A: cos_phi_max zur Abschaetzung heranziehen [Mueller08, S.565]
% cos_phi_N_data = [P_mech {kW}; cos_phi_max(p=1...3) {%}; cos_phi_min(p=1...3) {%}; cos_phi_max(p=4...8) {%}; cos_phi_min(p=4...8) {%}]
cos_phi_N_data = [ 1.0  2.0  4.0  6.0  8.0 10.0 20.0 40.0 60.0 80.0 100.0 200.0 400.0 600.0 800.0 1000.0 2000.0 4000.0 5000.0; ...
                  78.0 80.9 83.5 85.0 85.9 86.6 88.4 89.4 90.0 90.2  90.3  90.8  91.1  91.2  91.2   91.2   91.2   91.2   91.2; ...
                  76.2 78.5 80.5 81.7 82.4 82.9 83.9 84.5 84.9 85.0  85.0  85.0  85.0  85.0  85.0   84.9   84.5   84.0   83.9; ...
                   NaN  NaN  NaN  NaN  NaN 78.5 80.2 81.8 82.5 83.0  83.3  84.1  84.9  85.1  85.6   85.8   86.2   86.4   86.5; ...
                   NaN  NaN  NaN  NaN  NaN 68.9 71.0 72.4 73.1 73.4  73.6  74.3  74.9  75.0  75.0   75.0   74.9   74.9   74.9];

if(rated.p>=1 &&rated.p<=3)
   emag.cos_phi_N = interp1(cos_phi_N_data(1,:), cos_phi_N_data(2,:),(rated.P_N/1e3),'linear')/100;
elseif(rated.p>=4 && rated.p<=8)
   emag.cos_phi_N = interp1(cos_phi_N_data(1,:), cos_phi_N_data(4,:),(rated.P_N/1e3),'linear')/100;
else
    error('Ungueltige Eingabe bei Variable "rated.p"')
end
if(any(isnan(emag.cos_phi_N)))
    error('Datenbereich der hinterlegten Kurve zur Abschaetzung des Nennleistungsfaktors zu klein. \n(keine Werte fuer den Nennleistungsfaktor cos_phi_N bei einer Nenneistung von P_N = %.2f W hinterlegt)',rated.P_N(1,end));
end
clear cos_phi_N_data

% Berechnung der elektrischen Nennleistung P_el_N [W]
% Quelle: [Mueller08, S.569 - Formel 9.1.12]
emag.P_el_N = rated.P_N / emag.eta_N;

% Berechnung der Scheinleistung P_s [VA] fuer Motoren
% Quelle: [Mueller08, S.572 - Formel 9.1.17a]
emag.P_s = rated.P_N / (emag.eta_N * emag.cos_phi_N);

% Berechnung des Strangstroms (Stator) I_1Str [A]
% Quelle: [Mueller08, S.569 - Formel 9.1.12], [Mueller14, S.34 - Formel 0.6.8]
emag.I_1Str = emag.P_el_N / (rated.m * emag.U_1Str * emag.cos_phi_N);

% Berechnung des Bemessungsstroms I_N [A]
% Quelle: [Mueller14, S.33 - Formel 0.6.4], [Mueller14, S.34 - Formel 0.6.7]
if(strcmp(opt.Schaltung,'Dreieck'))
    emag.I_N = emag.I_1Str * sqrt(3);
elseif(strcmp(opt.Schaltung,'Stern'))
    emag.I_N = emag.I_1Str;
else
    error('Ungueltige Eingabe bei Variable "opt.Schaltung"');
end

% Berechnung der inneren Scheinleistung P_si [VA]
% Quelle: [Mueller08, S.569 - Formel 9.1.12], [Mueller08, S.572 - Formel 9.1.17a]
emag.P_si = (emag.P_el_N / emag.cos_phi_N) * (emag.E_h / emag.U_1Str);
% emag.P_si = emag.P_s * (emag.E_h / emag.U_1Str);

%% B.1) Statorinnendurchmesser
% #########################################################################
% #   B.1) STATORINNENDURCHMESSER                                         #
% #########################################################################

% Abschaetzung des Ausnutzungsfaktors C_mech [kWmin/m^3] aus Erfahrungswerten
% Quelle: [Mueller08, S.571 - Bild 9.1.6]
% V/A: obere Grenze der Bereiche zur Abschaetzung heranziehen [Mueller08, S.565]
% C_mech_data = [P_mech/(2*p) {kW}; C_mech(p=1) {kWmin/m^3}; C_mech(p=2...6) {kWmin/m^3}; C_mech(p=8...12) {kWmin/m^3}, ; C_mech(p=24...28) {kWmin/m^3}]
C_mech_data = [0.10 0.20 0.40 0.60 0.80 1.00 2.00 4.00 6.00 8.00 10.00 20.00 40.00 60.00 80.00 100.00 200.00 400.00 600.00 800.00 1000.00; ...
                NaN  NaN  NaN 1.25 1.28 1.32 1.44 1.56 1.62 1.66  1.70  1.85  2.08  2.25  2.39   2.54   3.01   3.39   3.55   3.65    3.69; ...
               1.46 1.63 1.83 1.97 2.06 2.13 2.37 2.66 2.85 3.02  3.14  3.54  3.94  4.20  4.40   4.55   5.02   5.51   6.00   6.49    6.98; ...
               NaN  NaN  NaN  NaN  NaN  NaN  NaN  1.47 1.74 1.97  2.17  2.93  3.63  3.92  4.07   4.17    NaN    NaN    NaN    NaN     NaN; ...
               NaN  NaN  NaN  NaN  NaN  NaN  NaN  1.21 1.39 1.57  1.74  2.26  2.77  3.01  3.17   3.27    NaN    NaN    NaN    NaN     NaN];
 
if(rated.p==1)
    geo.misc.C_mech = interp1(C_mech_data(1,:), C_mech_data(2,:),((rated.P_N/1e3)/(2*rated.p)),'linear');
elseif(rated.p>=2 && rated.p<=6)
    geo.misc.C_mech = interp1(C_mech_data(1,:), C_mech_data(3,:),((rated.P_N/1e3)/(2*rated.p)),'linear');
elseif(rated.p>=8 && rated.p<=12)
    geo.misc.C_mech = interp1(C_mech_data(1,:), C_mech_data(4,:),((rated.P_N/1e3)/(2*rated.p)),'linear');
elseif(rated.p>=24 && rated.p<=28)
    geo.misc.C_mech = interp1(C_mech_data(1,:), C_mech_data(4,:),((rated.P_N/1e3)/(2*rated.p)),'linear');
else
    error('Ungueltige Eingabe bei Variable "Polpaarzahl"')
end
if(any(isnan(geo.misc.C_mech)))
    error('Datenbereich der hinterlegten Kurve zur Abschaetzung des Ausnutzungsfaktors zu klein. \n(keine Werte fuer den Ausnutzungsfaktors C_mech bei einer Nenneistung von P_N = %.2f W hinterlegt)',rated.P_N(1,end));
end
clear C_mech_data

% Reduktion des Ausnutzungsfaktors von Induktionsmaschinen anhand der Ausfuehrung und Betriebsart k_red [-]
% Quelle: [Mueller08, S.572 - Tabelle 9.1.1]
% V/A: Eintrittstemperatur des Primaerkuehlmittels = 40°C, keine Ex-geschuetzten Maschinen
if(strcmp(opt.Kuehlungsart,'Innen- oder Kreislaufkuehlung'))
    if(strcmp(opt.Maschinenausfuehrung,'Kaefiglaeufer'))
        geo.misc.k_red = 1.0;
    elseif(strcmp(opt.Maschinenausfuehrung,'Schleifringlaeufer'))
        geo.misc.k_red = 0.95; % zwischen 0.9 und 0.95
    else
        error('Ungueltige Eingabe bei Variable "opt.Maschinenausfuehrung"')
    end
elseif(strcmp(opt.Kuehlungsart,'Oberflaechenkuehlung'))
    if(strcmp(opt.Maschinenausfuehrung,'Kaefiglaeufer'))
        geo.misc.k_red = 0.88; % zwischen 0.7 und 0.9
    elseif(strcmp(opt.Maschinenausfuehrung,'Schleifringlaeufer'))
        geo.misc.k_red = 0.8; % zwischen 0.65 und 0.8
    else
        error('Ungueltige Eingabe bei Variable "opt.Maschinenausfuehrung"')
    end
else
    error('Ungueltige Eingabe bei Variable "opt.Kuehlungsart"')
end

% Berechnung des reduzierten Ausnutzungsfaktors C_mech_red [kWmin/m^3]
% Quelle: [Mueller08, S.572]
geo.misc.C_mech_red = geo.misc.k_red * geo.misc.C_mech;

% Berechnung des Ausnutzungsfaktors C [kVAmin/m^3]
% Quelle: [Mueller08, S.570 - Formel 9.1.14]
geo.misc.C = (geo.misc.C_mech_red / (emag.cos_phi_N * emag.eta_N)) * (emag.E_h / emag.U_1Str);

% Richtwert fuer die relative Ankerlaenge lambda [-]
% Quelle: [Mueller08, S.577 - Tabelle 9.1.3], [Meyer18, S.117]
% if(rated.p==1)
%     richt.lambda = 0.8; % zwischen 0.6 und 1.0
% elseif(rated.p>1)
%     richt.lambda = nthroot(rated.p,3);
%     if(richt.lambda>4.0 || richt.lambda<1.0)
%         richt.lambda = 2.5; % zwischen 1.0 und 4.0
%     end
% end

% Berechnung des Bohrungsdurchmessers (Stator) D_1i [m]
% Quelle: [Mueller08, S.576 - Formel 9.1.19a]
geo.D_1i = nthroot(((rated.P_N/1e3) * 2*rated.p) / (rated.n_N * geo.misc.C_mech_red * richt.lambda * pi),3);
D_1i_alt = nthroot(((emag.P_si/1e3) * 2*rated.p) / (rated.n_N * geo.misc.C * richt.lambda * pi),3);

% Kontrolle Bohrungsdurchmesser (kann spaeter auskommentiert werden)
if((geo.D_1i-D_1i_alt)>1e-3)
    error('Anomalie in der Berechnung des Bohrungsdurchmessers "D_1i"')
end
clear D_1i_alt

%% B.2) Geometrischer Luftspalt und Rotoraussendurchmesser
% #########################################################################
% #   B.2) GEOMETRISCHER LUFTSPALT UND ROTORAUSSENDURCHMESSER             #
% #########################################################################

% Berechnung des geometrischen Luftspalts delta [mm]
% Quelle: [Mueller08, S.585 - Formel 9.1.36a-b]
if(rated.p==1)
    geo.delta = 0.4 * nthroot(rated.P_N/1e3,4);
else
    geo.delta = 0.25 * nthroot(rated.P_N/1e3,4);
end

% Kontrolle ob groesser als min Luftspalt -> Fertigungsgruende, Sicherheit
% Quelle: [Mueller08, S.585 - Formel 9.1.36c]
if(geo.delta<0.2)
    geo.delta = 0.2;
end

% Berechnung des Aussendurchmessers (Rotor) D_2a [m]
% Quelle: [Meyer09, S.42 - Formel 5.25]
geo.D_2a = geo.D_1i - 2*(geo.delta/1e3);

%% B.3) Axiale Laenge der aktiven Teile
% #########################################################################
% #   B.3) AXIALE LAENGE DER AKTIVEN TEILE                                #
% #########################################################################

% Berechnung der Polteilung (Stator) tau_1p [m]
% Quelle: [Mueller08, S.577], [Mueller08, S.19 - Formel 1.1.16]
geo.tau_1p = (geo.D_1i * pi) / (2 * rated.p);

% Berechnung der ideellen Laenge l_i [m]
% Quelle: [Mueller08, S.570 - Formel 9.1.13], [Mueller08, S.577]
geo.l_i = (rated.P_N/1e3) / (geo.misc.C_mech_red * geo.D_1i^2 * rated.n_N);
l_i_alt = geo.tau_1p * richt.lambda;

% Kontrolle ideelle Laenge (kann spaeter auskommentiert werden)
if((geo.l_i-l_i_alt)>0.001)
    error('Anomalie in der Berechnung der ideellen Laenge "l_i"')
end
clear l_i_alt

% Richtwert fuer die Kanalbreite der Ventilationskanaele l_v [m]
% Quelle: [Meyer09, S.41], [Mueller08, S.585]
% richt.l_v = 0.01; % zwischen 0.006 und 0.01

% Berechnung der Anzahl der benoetigten Luftschlitze n_v [-]
% Quelle: eigenentwickelt mit Angaben von [Mueller08, S.585]
% V/A: Teilpaketlaenge festgelegt zu 80mm
if(strcmp(opt.Kuehlungsart,'Oberflaechenkuehlung'))
    if(geo.l_i>0.2)
        geo.misc.n_v = floor(geo.l_i / 0.08);
    else
        geo.misc.n_v = 0;
    end
elseif(strcmp(opt.Kuehlungsart,'Innen- oder Kreislaufkuehlung'))
    geo.misc.n_v = 0;
else
    error('Ungueltige Eingabe bei Variable "opt.Kuehlungsart"')
end

% Berechnung Hilfsfaktor gamma_v [-] (Carter-Faktor)
% Quelle: [Mueller08, S.204 - Formel 2.3.20], [Mueller08, S.206]
% Bemerkung: l_v an die Stelle von b_s
geo.misc.gamma_v = (1 / (1 + 5 * (geo.delta/1e3 / richt.l_v)));

% Berechnung der Blechpaketlaenge (Stator und Rotor) l [m]
% Quelle: [Mueller08, S.208 - Formel 2.3.25], [Mueller08, S.596]
% V/A: Rotor und Stator haben die gleiche Laenge, Luftschlitze im Rotor
% mit denen des Stators deckungsgleich
geo.l = geo.l_i + (geo.misc.n_v * geo.misc.gamma_v * richt.l_v) - 2*(geo.delta/1e3);

% Berechnung der reinen Eisenlaenge (Stator und Rotor) l_Fe [m]
% Quelle: [Mueller08, S.210 - Formel 2.3.30]
geo.l_Fe = geo.l - (geo.misc.n_v * richt.l_v);

% Berechnung des Bohrungsvolumens V_Bohrung [m^3]
% Quelle: [Mueller08, S.570 - Formel 9.1.13]
% geo.V_Bohrung = (rated.P_N/1e3) / (rated.n_N * geo.misc.C_mech_red) * (pi/4);
geo.V_Bohrung = geo.D_1i^2 * geo.l * (pi/4);

%% C) Statorentwurf
% #########################################################################
% #   C) STATORENTWURF                                                    #
% #########################################################################

% Richtwert fuer den Mittelwert der Luftspaltinduktion B_m [T]
% Quelle: [Mueller08, S.582 - Tabelle 9.1.5]
% richt.B_m = 0.58; % zwischen 0.4 und 0.65

if(isfield(richt, 'B_m'))
    % Berechnung des Strombelags A [A/mm]
    % Quelle: [Mueller08, S.570 - Formel 9.1.14]
    emag.A = (geo.misc.C * 2*sqrt(2)) / (pi^3 * richt.B_m * richt.xi_1p) * 60;

    % Ueberpruefung ob A im zulaessigen Bereich
    % Quelle: [Mueller08, S.580 - Tabelle 9.1.4]
    if(emag.A<20 || emag.A>120) % zwischen 20.0 und 120.0
        warning('Variable "emag.A" ausserhalb der Grenzen')
    end
elseif(isfield(richt, 'A'))
    % Berechnung des Mittelwerts der Luftspaltinduktion B_m [T]
    % Quelle: [Mueller08, S.572 - Formel 9.1.14]
    emag.B_m = (geo.misc.C * 2*sqrt(2)) / (pi^3 * richt.A * richt.xi_1p) * 60;
    
    % Ueberpruefung ob B_m im zulaessigen Bereich
    % Quelle: [Mueller08, S.582 - Tabelle 9.1.5]
    if(emag.B_m<0.4 || emag.B_m>0.65) % zwischen 0.4 und 0.65
        warning('Variable "emag.B_m" ausserhalb der Grenzen')
    end
    
    richt.B_m = emag.B_m;
    emag.A = richt.A;
else
    error('Kein Wert fuer den Mittelwert der Luftspaltinduktion "B_m" oder den Strombelag "A" gefunden.')
end

% Berechnung der Amplitude der Luftspaltinduktion B_p [T]
% Quelle: [Mueller08, S.570]
emag.B_p = richt.B_m / (2/pi);

% Berechnung des Luftspaltflusses Phi_delta [Wb]
% Quelle: [Mueller08, S.198 - Formel 2.3.5]
emag.Phi_delta = richt.B_m * geo.tau_1p * geo.l_i;

% Berechnung des Hauptwellenflusses Phi_h [Wb]
% Quelle: [Mueller08, S.596], [Mueller08, S.199 - Formel 2.3.11]
% V/A: Hauptwellenfluss kann fuer Wicklungsentwurf mit genuegender
% Genauigkeit gleich dem Luftspaltfluss gesetzt werden
emag.Phi_h = emag.Phi_delta;
Phi_h_alt = 2/pi * emag.B_p * geo.tau_1p * geo.l_i;

% Kontrolle Hauptwellenfluss (kann spaeter auskommentiert werden)
if((emag.Phi_h-Phi_h_alt)>0.1)
    error('Anomalie in der Berechnung des Hauptwellenflusses "Phi_h"')
end
clear Phi_h_alt

% Berechnung des maximalen Flusses im Ruecken (Stator) Phi_1r_max [Wb]
% Quelle: [Mueller08, S.599]
emag.Phi_1r_max = (1.05*emag.Phi_delta) / 2;

% Richtwert fuer die max. zulaessige Induktion im Ruecken (Stator) B_1r_max [T]
% Quelle: [Mueller08, S.582 - Tabelle 9.1.5]
% V/A: Begrenzung auf 1.8T wegen Diagramm fuer Rueckenreduktionsfaktor
% richt.B_1r_max = 1.65; % zwischen 1.3 und 1.65 (...2)

% Richtwert fuer die max. zulaessige Induktion in den Zaehnen (Stator) B_1z_max [T]
% Quelle: [Mueller08, S.582 - Tabelle 9.1.5]
% richt.B_1z_max = 1.65; % zwischen 1.4 und 2.1

% Richtwert fuer den Nutfuellfaktor (Stator) phi_1n [-]
% Quelle: [Mueller08, S.586 - Tabelle 9.1.6]
% V/A: Niederspannung
% if(strcmp(opt.Spulenform_Stator,'Runddrahtwicklungen'))
%     richt.phi_1n = 0.5; % zwischen 0.3 und 0.5
% elseif(strcmp(opt.Spulenform_Stator,'Formspulen- oder Stabwicklungen'))
%     richt.phi_1n = 0.6; % zwischen 0.35 und 0.6
% else
%     error('Ungueltige Eingabe bei Variable "opt.Spulenform_Stator"')
% end

% Richtwert fuer den Wicklungsfaktor (Stator) xi_1p [-]
% Quelle: [Mueller08, S.596]
% if(strcmp(opt.Schichten,'Einschichtwicklung'))
%     richt.xi_1p = 0.92;
% elseif(strcmp(opt.Schichten,'Zweischichtwicklung'))
%     richt.xi_1p = 0.96;
% else
%     error('Ungueltige Eingabe bei Variable "opt.Schichten"')
% end

% Richtwert fuer den Eisenfuellfaktor (Stator) phi_1Fe [-]
% Quelle: [Mueller08, S.599]
% richt.phi_1Fe = 0.95;

% Richtwert fuer die Stromdichte (Stator) S_1 [A/mm^2]
% Quelle: [Mueller08, S.580 - Tabelle 9.1.4]
% richt.S_1 = 7.0; % zwischen 3.0 und 8.0

% Ueberpruefung ob Produkt (Strombelag * Stromdichte) [(A/mm)*(A/mm^2)] im zulaessigen Bereich
% Quelle: [Mueller08, S.581 - Formel 9.1.28]
if(emag.A*richt.S_1<100 || emag.A*richt.S_1>350) % zwischen 100 und 350 je nach Maschinengroesse und Intensitaet der Kuehlung
    warning('Produkt (Strombelag * Stromdichte) ausserhalb der Grenzen')
end

% Abschaetzung der mittleren Windungslaenge der Wicklungen (Stator) l_1m [m]
% Quelle: [Mueller08, S.586 - Formel 9.1.39]
wick.l_1m = 2.0 * (geo.l + 1.3 * geo.tau_1p + (0.03 + 0.02*(emag.U_1Str/1000)));

% Berechnung der minimal notwendigen Rueckenhoehe (Stator) h_1r_min [m]
% Quelle: [Mueller08, S.581 - Formel 9.1.30]
geo.Nut.h_1r_min = (emag.Phi_1r_max) / (geo.l_Fe * richt.phi_1Fe * richt.B_1r_max);

% Berechnung der maximalen Rueckenhoehe (Stator) h_1r_max [m]
% Quelle: [Mueller08, S.581 - Formel 9.1.30]
% V/A: Saettigungseffekte und nicht sinusfoermige Verteilung --> Faktor (2/pi) [Meyer09, S.40]
geo.Nut.h_1r_max = (emag.Phi_1r_max / (2/pi)) / (geo.l_Fe * richt.phi_1Fe * richt.B_1r_max);

% Berechnung der minimal notwendigen Nuthoehe (Stator) h_1n_min [m]
% Quelle: [Mueller08, S.597 - Formel 9.2.7b, Mueller08, S.581, Mueller08, S.215 - Formel 2.4.3]
% V/A: parallelflankige offene Nuten
% 4 Gleichungen umformen: h_n=A_n/b_n, S=(A*tau_n)/(A_n*phi_n), tau_n=b_n+b_z, b_z=(B_m*1.4*tau_n*l_i)/(B_zmax*phi_1Fe*l_Fe)
geo.Nut.h_1n_min = emag.A / (richt.S_1 * richt.phi_1n * (1-(((richt.B_m*1.4)*geo.l_i)/(richt.B_1z_max*richt.phi_1Fe*geo.l_Fe))))/1e3;

% Berechnung der maximal notwendigen Nuthoehe (Stator) h_1n_max [m]
% Quelle: [Meyer09, S.41 - Formel 5.20]
geo.Nut.h_1n_max = 1.25 * geo.Nut.h_1n_min;

% Berechnung des maximal erforderlichen Aussendurchmessers (Stator) D_1a_max [m]
% Quelle: [Mueller, S.599 - Formel 9.2.15a]
geo.D_1a_max = geo.D_1i + 2*geo.Nut.h_1n_max + 2*geo.Nut.h_1r_max;

%% C.1) Wicklungsschema
% #########################################################################
% #   C.1) WICKLUNGSSCHEMA                                                #
% #########################################################################

% Richtwert fuer die minimale Nutteilung (Stator) tau_1n_min [m]
% Quelle: [Meyer09, S.46], [Pyr14] -> evtl auch kleiner (bis zu 0.007m)
% richt.tau_1n_min = 0.02; % zwischen 0.01 und 0.07

% Berechnung der maximalen Nutzahl (Stator) N_1max [-]
% Quelle: [Mueller08, S.19 - Formel 1.1.14], [Meyer09, S.46 - Formel 5.37]
wick.N_1max = floor((geo.D_1i * pi) / richt.tau_1n_min);

% Vektor fuer Nutenzahl (Stator) N_1 --> 1. Randbedingung
N_1 = (1:wick.N_1max)';

% Berechnung des Lochzahl-Vektors (Stator) q_1 [-]
% Quelle: [Mueller08, S.21 - Formel 1.2.2], [Meyer09, S.45 - Formel 5.35]
q_1 = N_1 ./ (2.*rated.p.*rated.m);

% Zaehler und Nenner von Lochzahl (Stator) ermitteln [-]
[q_1z, q_1n] = numden(sym(q_1));
q_1z = double(q_1z);
q_1n = double(q_1n);

% Groessten gemeinsamen Teiler von Lochzahl-Nenner und Strangzahl ermitteln [-]
for i = 1:length(q_1n)
    ggT_q_1n_m(i,1) = ggT_fun(q_1n(i,1),rated.m);
end

% Berechnung der Nutteilung (Stator) tau_1n [m]
% Quelle: [Mueller08, S.19 - Formel 1.1.14]
tau_1n = (geo.D_1i.*pi) ./ N_1;

% Durchmesserschritt Wicklung (Stator) y_Durchmesser [-]
% Quelle: [Mueller08, S.19 - Formel 1.1.15]
y_1Durchmesser = N_1 ./ (2.*rated.p);

% Zusammenfuegen in eine Tabelle
wick.table_1all = table(N_1, q_1, q_1z, q_1n, ggT_q_1n_m, tau_1n,  y_1Durchmesser);
clear N_1 q_1 q_1z q_1n ggT_q_1n_m tau_1n y_1Durchmesser

% Richtwert fuer die minimale Lochzahl q_min fuer Ganzlochwicklungen [-]
% Quelle: [Mueller08, S.121], [Meyer09, S.45]
wick.q_min_GL = 2;

% Berechnung der minimalen Statornutzahl N_min fuer Ganzlochwicklung [-]
% Quelle: [Mueller08, S.21 - Formel 1.2.2]
wick.N_min_GL = (2 * rated.p * rated.m * wick.q_min_GL);

% Wicklungstypen
% Wicklung 1SGL: Einschicht, Ganzloch
% Wicklung 1SBL: Einschicht, Bruchloch
% Wicklung 2SGL: Zweischicht, Ganzloch
% Wicklung 2SBL: Zweischicht, Bruchloch

% Modus Wicklungsauslegung
if(strcmp(opt.Mode_Wicklung,'Klassisch'))
    % Wicklung A: Ganzlochwicklung, Einschicht, ungesehnt, gezont
    % Wicklung B: Ganzlochwicklung, Zweischicht, gesehnt, gezont
    % Wicklung C: Bruchlochwicklung, Zweischicht, gesehnt, gezont, Zahnspulenwicklung (q<1) wird dann umgesetzt, wenn der Durchmesser q>1 nicht zulaesst
    
    % Ueberpruefung ob Ganzlochwicklung moeglich ist, ansonsten wird Benutzereingabe ueberschrieben, weil fuer Realisierung eine Bruchlochwicklung verwendet werden muss!
    if(wick.N_min_GL>wick.N_1max && ~strcmp(opt.Wicklungstyp,'C'))
        warning('Wicklungstyp A und B nicht moeglich! Wicklungstyp C wird verwendet!');
        opt.Wicklungstyp = 'C';
    end
    
    switch opt.Wicklungstyp
        case 'A'
            Wicklungstypen = {'1SGL'};
        case 'B'
            Wicklungstypen = {'2SGL'};
        case 'C'
            Wicklungstypen = {'2SBL'};
    end
elseif(strcmp(opt.Mode_Wicklung,'Optimierung'))
    Wicklungstypen = {'1SGL', '1SBL', '2SGL', '2SBL'};
    
elseif(strcmp(opt.Mode_Wicklung,'Manuell'))
    Wicklungstypen = {'1SGL', '2SGL', '2SBL'}; % 1SBL nur als Zahnspulenwicklung keine Strangverschachtelung o.ae. moeglich
    
else
    error('Ungueltige Eingabe bei Variable "opt.Mode_Wicklung"');
end

% Schleife ueber alle Wicklungen
for i_Wicklung = 1:length(Wicklungstypen)
    % Ganzloch- bzw. Bruchlochwicklungen heraussuchen
    switch Wicklungstypen{i_Wicklung}
        case {'1SGL', '2SGL'}
            % moegliche Wicklungen fuer ganzzahlige Lochzahl suchen
            % -> q_1n == 1
            W = wick.table_1all(wick.table_1all.q_1n==1,:);
        case {'1SBL', '2SBL'}
            % moegliche Wicklungen fuer nicht ganzzahlige Lochzahl suchen
            % -> q_1n ~= 1
            W = wick.table_1all(wick.table_1all.q_1n~=1,:);
    end

    % Speichern des Wicklungstyps
    W.Wicklungstyp_1(:) = Wicklungstypen(i_Wicklung);
    W = [W(:,end) W(:,1:end-1)];

    % Anzahl Schichten hinzufuegen
    switch Wicklungstypen{i_Wicklung}
        case {'1SGL', '1SBL'}
            W.n_lay = 1 * ones(length(W.N_1),1);
        case {'2SGL', '2SBL'}
            W.n_lay = 2 * ones(length(W.N_1),1);
    end

    % Symmetrie --> 2. Randbedingung
    % nur fuer Bruchlochwicklungen relevant, da Ganzlochwicklungen stets symmetrisch
    % Quelle: [Mueller08, S.30]
    switch Wicklungstypen{i_Wicklung}
        case '1SBL'
            % Symmetriebedingung: p/n ist eine natuerliche Zahl & ggT{n,m} == 1
            % Quelle: [Mueller08, S.33 - Tabelle 1.2.3]
            Symmetrie = (~mod((rated.p),W.q_1n) & W.ggT_q_1n_m==1);

            % Tabelle aktualisieren
            W = W(Symmetrie,:);

        case '2SBL'
            % Symmetriebedingung: 2p/n ist eine natuerliche Zahl & ggT{n,m} == 1
            % Quelle: [Mueller08, S.33 - Tabelle 1.2.3]
            Symmetrie = (~mod((2*rated.p),W.q_1n) & W.ggT_q_1n_m==1);

            % Tabelle aktualisieren
            W = W(Symmetrie,:);
    end
    clear Symmetrie

    % distributed 1SBL aussortieren (q>1), bei 1SBL werden nur Zahnspulenwicklungen betrachtet
    % Quelle: [Mueller08, S.51]
    switch Wicklungstypen{i_Wicklung}
        case '1SBL'
            W = W(W.q_1<1,:);
    end
    if(isempty(W))
        continue;
    end

    % ungezonte GL aussortieren (q<2)
    % Quelle: [Mueller08, S.121], [Meyer09, S.45]
    switch Wicklungstypen{i_Wicklung}
        case {'1SGL', '2SGL'}
            W = W(W.q_1>=wick.q_min_GL,:);
    end
    if(isempty(W))
        continue;
    end

    % Wicklungsoptimierung durch Schrittverkuerzung y_v [-]
    % Wicklungsschritt y = y_Durchmesser - y_v [Mueller08, S.43 - Formel 1.2.11]
    % Berechnung der moeglichen sinnvollen Sehnungen Sehnung = y/y_Durchmesser
    for i = 1:length(W.N_1)
        count_Sehnung = 0;
        if(W.q_1(i)<1) % Zahnspulenwicklung
            % Zahnspulenwicklungen haben immer den Wicklungsschritt y=1
            % Quelle: [Mueller08, S.75]
            y_pos = 1;
        else % keine Zahnspulenwicklung
            switch Wicklungstypen{i_Wicklung}
                case {'1SGL'}
                    % Einschichtwicklung immer Durchmesserwicklungen
                    % Quelle: [Mueller08, S.40]
                    y_pos = W.y_1Durchmesser(i);
                case {'1SBL'}
                    %muss nicht betrachtet werden, weil aussortiert (s.o.)
                case {'2SGL'}
                    y_pos = 1:1:W.y_1Durchmesser(i)-1; % ungesehnte 2SGL werden nicht betrachtet, da gleiches elektromagnetisches Verhalten wie ungesehnte 1SGL
                case {'2SBL'}
                    y_pos = 1:1:(W.y_1Durchmesser(i)*2);
            end
        end
        y_v_pos = W.y_1Durchmesser(i) - y_pos;
        for j = 1:length(y_v_pos)
            % Sinnvolle Sehnung zwischen 2/3<=y/y_Durchmesser<=4/3
            % Quelle: [Mueller08, S.76]
            if((y_pos(j)/W.y_1Durchmesser(i))>=(2/3) && (y_pos(j)/W.y_1Durchmesser(i))<=(4/3))
                if(count_Sehnung>0)
                    W(end+1,:) = W(i,:);
                    W.y_1(end) = y_pos(j);
                    W.y_1v(end) = y_v_pos(j);
                else
                    W.y_1(i) = y_pos(j);
                    W.y_1v(i) = y_v_pos(j);
                end
                count_Sehnung = count_Sehnung + 1;
            end

        end
        clear count_Sehnung y_pos y_v_pos
    end
    W.Sehnung_1 = W.y_1 ./ W.y_1Durchmesser;

    % Nicht sinnvoll sehnbare Wicklungen loeschen (betrifft Zahnspulenwicklungen)
    W(W.Sehnung_1==0,:) = [];
    if(isempty(W))
        continue;
    end

    % Sortieren der Tabelle
    W = sortrows(W,{'N_1','Sehnung_1'},{'ascend','ascend'});

    % Richtwert fuer die Nutschlitzbreite (Stator) b_1ns [mm]
    % Quelle: [Meyer18, S.94 - Tabelle 5.3]
    b_ns_data = [0.05 0.2 0.3 0.45; 1.5 3 3 5];
    b_ns_fun = polyfit(b_ns_data(1,:),b_ns_data(2,:),1);
    NG.b_1ns = polyval(b_ns_fun,geo.D_1i);
    clear b_ns_data b_ns_fun
    
    if(~strcmp(opt.Mode_Wicklung,'Manuell'))
        for i = 1:length(W.N_1)
            % Berechnung der auftretenden Oberwellen nu [-]
            % Quelle: [Mueller08, S.54 - Formel 1.2.36]
            g = [0,1:100000,-1:-1:-100000];
            nu_1s = rated.p * (1 + ((2 .* rated.m .* g) ./ W.q_1n(i)));
            nu_1s = sort(nu_1s,2,'ComparisonMethod','abs');
            nu_1 = nu_1s./rated.p;

            % Vorabberechnungen Zonenfaktor
            % Quelle: [Binder12, S.128-135]
            if(mod(W.q_1n(i),2)) %ungerade
                p_u = W.q_1n(i);
                Q_u = 2 * rated.m * W.q_1z(i);
            else %gerade
                p_u = W.q_1n(i) / 2;
                Q_u = rated.m * W.q_1z(i);
            end

            if(mod(Q_u,2)) %ungerade
                q_1 = (Q_u + rated.m) / (2 * rated.m);
                q_2 = q_1 - 1;
            else %gerade
                q_1 = Q_u / (2 * rated.m);
                q_2 = q_1;
            end

            alpha_Q = (2 * pi * p_u) / Q_u;

            g_min = 0;
            Y = ((g_min * Q_u) + 1) / p_u;
            while(mod(Y,1))
                g_min = g_min + 1;
                Y = ((g_min * Q_u) + 1) / p_u;
            end

            % Berechnung des Gruppenfaktors / Zonenfaktors (Stator) xi_1gr [-]
            % Quelle: [Binder12, S.136 - Formel 3.3-25]
            xi_1gr = (sin(nu_1 .* alpha_Q .* Y .* (q_1./2)) - cos(nu_1 .* p_u .* pi .* Y) .* sin(nu_1 .* alpha_Q .* Y .* (q_2./2))) ./ ((q_1 + q_2) .* sin(nu_1 .* alpha_Q .* (Y./2)));

            % Berechnung des Spulenfaktors / Sehnungsfaktors (Stator) xi_1sp [-]
            % Quelle: [Binder12, S.136 - Formel 3.3-25]
            xi_1sp = sin(nu_1 .* (pi./2) .* W.Sehnung_1(i));

            % Berechnung des Spulenfaktors / Sehnungsfaktors (Stator) xi_1sp [-]
            % Quelle: [Mueller08, S.88 - Formel 1.2.67]
            % xi_1sp = sin((nu_1s.*W.y_1.*pi)./(rated.p.*W.y_1Durchmesser.*2));

            % Berechnung des Gruppenfaktors / Zonenfaktors (Stator) xi_1gr [-]
            % Quelle: [Mueller08, S.88 - Formel 1.2.67]
            % xi_1gr = (sin((nu_1s.*W.q_1.*pi)./(rated.p.*W.y_1Durchmesser.*2)))./(W.q_1.*sin((nu_1s.*1.*pi)./(rated.p.*W.y_1Durchmesser.*2)));

            % Berechnung des Nutschlitzbreitenfaktors (Stator) xi_1n [-]
            % Quelle: [Mueller08, S.87 - Formel 1.2.62]
            xi_1n = (sin((nu_1s.*(NG.b_1ns./1e3))./geo.D_1i))./((nu_1s.*(NG.b_1ns./1e3))./geo.D_1i);

            % Berechnung des Produkts aus Gruppen- und Spulenfaktor (Stator) xi_1sz [-]
            % Quelle: [Mueller08, S.88 - Formel 1.2.67]
            xi_1spgr = xi_1sp .* xi_1gr;

            % Berechnung des gesamten Wicklungsfaktors (Stator) xi [-]
            % Quelle: [Mueller08, S.88 - Formel 1.2.67]
            xi_1 = xi_1sp .* xi_1gr .* xi_1n;

            % Hauptwellen Wicklungsfaktor (Stator) xi_1p [-]
            xi_1p = xi_1(nu_1==1);

            % Berechnung des Streukoeffizienten der Oberwellenstreuung sigma_1o (Stator) [-]
            % Quelle: [Mueller08, S.337 - Formel 3.7.32]
            % V/A: nur die ersten 200000 Oberwellen zur Berechnung herangezogen
            var1 = nu_1s(nu_1~=1);
            var2 = xi_1spgr(nu_1~=1);
            W.sigma_1o(i) = sum(((rated.p./var1) .* (var2./xi_1spgr(nu_1==1))).^2);

            % Speichern des Faktoren
            W.xi_1sp(i) = {xi_1sp(1:15)'};
            W.xi_1gr(i) = {xi_1gr(1:15)'};
            W.xi_1n(i) = {xi_1n(1:15)'};
            W.xi_1spgr(i) = {xi_1spgr(1:15)'};
            W.xi_1(i) = {xi_1(1:15)'};
            W.xi_1p(i) = abs(xi_1p);
            W.nu_1(i) = {nu_1(1:15)'};

            clear g nu_1s nu_1 p_u Q_u q_1 q_2 alpha_Q g_min Y xi_1sp xi_1gr xi_1n xi_1spgr xi_1 xi_1p var1 var2
        end
    end

%% C.2) Verschaltung
% #########################################################################
% #   C.2) VERSCHALTUNG                                                   #
% #########################################################################

    % Berechnen der maximalen Anzahl an Parallelschaltungen a_max [-]
    % und der moeglichen Kombinationen von Parallelschaltungen a_pos [-]
    switch Wicklungstypen{i_Wicklung}
        case '1SGL'
            % Maximalzahl paralleler Wicklungszweige a_max [-]
            % Quelle: [Mueller08, S.115]
            a_max = rated.p;

            % moegliche Kombinationen von Parallelschaltungen
            % Bedingung: p/a ist eine natuerliche Zahl
            % Quelle: [Mueller08, S.115]
            a_pos = 1;
            for a = 2:a_max
                if(~mod((rated.p/a),1))
                    a_pos = [a_pos, a];
                end
            end

        case '2SGL'
            % Maximalzahl paralleler Wicklungszweige a_max [-]
            % Quelle: [Mueller08, S.115]
            a_max = 2 * rated.p;

            % moegliche Kombinationen von Parallelschaltungen
            % Bedingung: p/a ist eine natuerliche Zahl
            % Quelle: [Mueller08, S.115]
            a_pos = 1;
            for a = 2:a_max
                if(~mod((rated.p/a),1))
                    a_pos = [a_pos, a];
                end
            end
        case {'1SBL', '2SBL'}
            W_copy = W;
            for i = 1:length(W.N_1)
                if(mod(W.q_1n(i),2)) % ungerade
                    % Zahl der Urverteilungen [-]
                    % Quelle: [Mueller08, S.36 - Tabelle 1.2.4]
                    t = rated.p / W.q_1n(i);

                    % Maximalzahl paralleler Wicklungszweige a_max [-]
                    % Quelle: [Mueller08, S.36 - Tabelle 1.2.4]
                    a_max = 2*t;

                else % gerade
                    % Zahl der Urverteilungen [-]
                    % Quelle: [Mueller08, S.36 - Tabelle 1.2.4]
                    t = 2 * (rated.p / W.q_1n(i));

                    % Maximalzahl paralleler Wicklungszweige a_max [-]
                    % Quelle: [Mueller08, S.36 - Tabelle 1.2.4]
                    a_max = t;

                end

                % moegliche Kombinationen von Parallelschaltungen
                % Bedingung: t/a ist eine natuerliche Zahl
                % Quelle: [Mueller08, S.36 - Tabelle 1.2.4]
                a_pos = 1;
                for a = 2:a_max
                    if(~mod((t/a),1))
                        a_pos = [a_pos, a];
                    end
                end

                % Wicklungen zu Tabelle hinzufuegen
                count_Parallel = 0;
                for j = 1:length(a_pos)
                    if(count_Parallel>0)
                        W_copy(end+1,:) = W_copy(i,:);
                        W_copy.a_1(end) = a_pos(j);
                    else
                        W_copy.a_1(i) = a_pos(j);
                    end
                    count_Parallel = count_Parallel + 1;
                end

            end
            W = W_copy; 
    end

    switch Wicklungstypen{i_Wicklung}
        case {'1SGL','2SGL'}
            % Wicklungen zu Tabelle hinzufuegen
            for i = 1:length(W.N_1)
                count_Parallel = 0;
                for j = 1:length(a_pos)
                    if(count_Parallel>0)
                        W(end+1,:) = W(i,:);
                        W.a_1(end) = a_pos(j);
                    else
                        W.a_1(i) = a_pos(j);
                    end
                    count_Parallel = count_Parallel + 1;
                end
            end
    end
    clear a_max a_pos a count_Parallel i j t W_copy

    % Sortieren der Tabelle
    W = sortrows(W,{'N_1','Sehnung_1','a_1'},{'ascend','ascend','ascend'});
    
    if(strcmp(opt.Mode_Wicklung,'Manuell') && length(Wicklungstypen)~=i_Wicklung)
    elseif(strcmp(opt.Mode_Wicklung,'Manuell') && length(Wicklungstypen)==i_Wicklung)
        if(i_Wicklung==1)
            wick.table_1pos = W;
        else
            wick.table_1pos = [wick.table_1pos; W];
        end
        clear W saveNG
        
        popup_String = '';
        for i = 1:length(wick.table_1pos.N_1)
            popup_String = [popup_String; {[wick.table_1pos.Wicklungstyp_1{i} ', N_1=' num2str(wick.table_1pos.N_1(i)) ', q_1=' num2str(wick.table_1pos.q_1(i)) ', Sehnung_1=' num2str(wick.table_1pos.Sehnung_1(i)) ', a_1=' num2str(wick.table_1pos.a_1(i))]}];
        end

        handles.popupmenu_Auswahl_Wicklung.String = ['-'; popup_String];
        handles.popupmenu_Auswahl_Wicklung.Value = 1;
        set(handles.uipanel_Auswahl_Wicklung,'Position',[0,0,742,826])
        set(handles.uipanel_Auswahl_Wicklung,'Visible','on')
        set(handles.uipanel_Ergebnisse,'Visible','off')
        handles.text_Wicklungslayout.Visible = 'on';
        handles.edit_Wicklungslayout.Visible = 'on';
        handles.text_y_1v.Visible = 'on';
        handles.edit_y_1v.Visible = 'on';
        handles.pushbutton_Wicklungslayout.Visible = 'on';
        handles.axes_Auswahl_Wicklung.Visible = 'off';
        handles.toolbar.Children(1).Enable = 'off';
        handles.toolbar.Children(2).Enable = 'off';
        handles.toolbar.Children(3).Enable = 'off';
        handles.toolbar.Children(4).Enable = 'off';
        save(['3_Ergebnisse/1_Misc/Entwurf_',opt.file_id,'_temp.mat'],'wick','rated');

        uiwait(gcf);

        load(['3_Ergebnisse/1_Misc/Entwurf_',opt.file_id,'_temp.mat']);
        delete(['3_Ergebnisse/1_Misc/Entwurf_',opt.file_id,'_temp.mat']);

        handles.toolbar.Children(1).Enable = 'on';
        handles.toolbar.Children(2).Enable = 'on';
        handles.toolbar.Children(3).Enable = 'on';
        handles.toolbar.Children(4).Enable = 'on';

        if(handles.popupmenu_Auswahl_Wicklung.Value~=1)
            wick.table_1sel = wick.table_1pos(handles.popupmenu_Auswahl_Wicklung.Value-1,:);
        else
            % Parallelverschaltung
            % V/A: bei Modifikation der Wicklung kann keine Aussage ueber die
            % Moeglichkeiten der Parallelverschaltung getroffen werden
            wick.table_1sel.a_1 = 1;
        end
        W = wick.table_1sel;

        set(handles.uipanel_Auswahl_Wicklung,'Position',[1484,0,742,826])
        set(handles.uipanel_Auswahl_Wicklung,'Visible','off')
        set(handles.uipanel_Ergebnisse,'Visible','on')
        clear popup_String
        pause(0.1)

        cla(handles.axes_Auswahl_Wicklung,'reset')
        handles.popupmenu_Auswahl_Wicklung.Value = 1;
        handles.popupmenu_Auswahl_Wicklung.String = {'-'};
        handles.edit_y_1v.String = {'-'};

        Wicklungstypen = W.Wicklungstyp_1;
        i_Wicklung = 1 ;

        % Berechnung der auftretenden Oberwellen nu [-]
        % Quelle: [Mueller08, S.54 - Formel 1.2.36]
        g = [0,1:10000,-1:-1:-10000];
        nu_1s = rated.p * (1 + ((2 .* rated.m .* g) ./ W.q_1n));
        nu_1s = sort(nu_1s,2,'ComparisonMethod','abs');
        nu_1 = nu_1s./rated.p;

        % Windungfaktor mittels Nutenstern berechnen (geometrische Basis)
        h = 1; % da Wicklung symmetrisch ist, reicht es aus den Windungsfaktor fuer einen Strang zu berechnen (hier fuer Strang h=1)
        x=zeros(length(nu_1),1,length(wick.Matrix(h,:)));
        y=zeros(length(nu_1),1,length(wick.Matrix(h,:)));
        z=zeros(length(nu_1),1);
        for i = 1:length(nu_1)
            % x-Koordinate von Nutenpfeil
            x(i,h,:) = sin(((2.*pi)./W.N_1) .* rated.p .* nu_1(i) .* (abs(wick.Matrix(h,:))-1)) .* sign(wick.Matrix(h,:));
            % y-Koordinate von Nutenpfeil
            y(i,h,:) = cos(((2.*pi)./W.N_1) .* rated.p .* nu_1(i) .* (abs(wick.Matrix(h,:))-1)) .* sign(wick.Matrix(h,:));
            % Laenge von Nutenpfeil
            z(i,h) = sqrt((sum(x(i,h,:))).^2 + (sum(y(i,h,:))).^2);
        end

        % Windungsfaktor
        xi_1spgr = (z ./ ((W.N_1*W.n_lay)/rated.m))';

        clear x y z

        % Richtwert fuer die Nutschlitzbreite (Stator) b_1ns [mm]
        % Quelle: [Meyer18, S.94 - Tabelle 5.3]
        b_ns_data = [0.05 0.2 0.3 0.45; 1.5 3 3 5];
        b_ns_fun = polyfit(b_ns_data(1,:),b_ns_data(2,:),1);
        NG.b_1ns = polyval(b_ns_fun,geo.D_1i);
        clear b_ns_data b_ns_fun

        % Berechnung des Nutschlitzbreitenfaktors (Stator) xi_1n [-]
        % Quelle: [Mueller08, S.87 - Formel 1.2.62]
        xi_1n = (sin((nu_1s.*(NG.b_1ns./1e3))./geo.D_1i))./((nu_1s.*(NG.b_1ns./1e3))./geo.D_1i);

        % Berechnung des gesamten Wicklungsfaktors (Stator) xi [-]
        % Quelle: [Mueller08, S.88 - Formel 1.2.67]
        xi_1 = xi_1spgr .* xi_1n;

        % Hauptwellen Wicklungsfaktor (Stator) xi_1p [-]
        xi_1p = xi_1(nu_1==1);

        % Berechnung des Streukoeffizienten der Oberwellenstreuung sigma_1o (Stator) [-]
        % Quelle: [Mueller08, S.337 - Formel 3.7.32]
        % V/A: nur die ersten 200000 Oberwellen zur Berechnung herangezogen
        var1 = nu_1s(nu_1~=1);
        var2 = xi_1spgr(nu_1~=1);
        W.sigma_1o = sum(((rated.p./var1) .* (var2./xi_1spgr(nu_1==1))).^2);

        % Speichern des Faktoren
        W.xi_1n = {xi_1n(1:15)'};
        W.xi_1spgr = {xi_1spgr(1:15)'};
        W.xi_1 = {xi_1(1:15)'};
        W.xi_1p = abs(xi_1p);
        W.nu_1 = {nu_1(1:15)'};

        clear i h g nu_1s nu_1 xi_1n xi_1spgr xi_1 xi_1p var1 var2
    
    else
    end
    
    if(strcmp(opt.Mode_Wicklung,'Manuell') && length(Wicklungstypen)~=i_Wicklung)
    else
        % Berechnung der optimalen Strangwindungszahl der Wicklung (Stator) w_1Str_opt [-]
        % Quelle: [Mueller08, S.114 - Formel 1.2.89]
        W.w_1Str_opt = (sqrt(2) .* emag.E_h) ./ (2 .* pi .* rated.f_N .* W.xi_1p .* emag.Phi_h);

        % moegliche ganzzahlige Leiterzahlen je Nut (Stator) z_1n [-]
        % Quelle: [Mueller08, S.115], [Mueller08, S.114]
        % schaerfere Bedingung fuer Zweischichtwicklungen
        % Leiterzahlen je Nut z_n muessen zudem geradzahlig sein
        % Quelle: [Mueller08, S.115]
        switch Wicklungstypen{i_Wicklung}
            case {'1SGL', '1SBL'}
                W.z_1n = round((W.w_1Str_opt .* 2 .* W.a_1 .* rated.m) ./ W.N_1);
            case {'2SGL', '2SBL'}
                var = floor((W.w_1Str_opt .* 2 .* W.a_1 .* rated.m) ./ W.N_1);
                W.z_1n = ceil((W.w_1Str_opt .* 2 .* W.a_1 .* rated.m) ./ W.N_1);
                W.z_1n(logical(~mod(var,2))) = var(logical(~mod(var,2)));
                clear var
        end

        % Nicht umsetzbare Wicklungen (z_n < 1) loeschen
        W(W.z_1n<1,:) = [];
        if(isempty(W))
            continue;
        end

        % moegliche resultierende Strangwindungszahlen der Wicklung (Stator)
        % Quelle: [Mueller08, S.114]
        W.w_1Str = (W.N_1 .* W.z_1n) ./ (2 .* W.a_1 .* rated.m);

        % Abweichung von optimaler Strangwindungszahl berechnen
        W.err_1 = W.w_1Str_opt - W.w_1Str;

        % Korrektur des Hauptwellenflusses Phi_h [Wb]
        % Quelle: [Mueller08, S.114 - Formel 1.2.89]
        W.Phi_h = (sqrt(2) .* emag.E_h) ./ (2 .* pi .* rated.f_N .* W.xi_1p .* W.w_1Str);

        % Korrektur des Luftspaltflusses Phi_delta [Wb]
        % Quelle: [Mueller08, S.596], [Mueller08, S.199 - Formel 2.3.11]
        % V/A: Luftspaltfluss kann fuer Wicklungsentwurf mit genuegender
        % Genauigkeit gleich dem Hauptwellenfluss gesetzt werden
        W.Phi_delta = W.Phi_h;

        % Korrektur des Mittelwerts der Luftspaltinduktion B_m [T]
        % Quelle: [Mueller08, S.198 - Formel 2.3.5]
        W.B_m = W.Phi_delta ./ (geo.tau_1p .* geo.l_i);

        % Korrektur der Amplitude der Luftspaltinduktion B_p [T]
        % Quelle: [Mueller08, S.198-199], [Mueller08, S.244-245], [Mueller08, S.199 - Formel 2.3.11]
        % Hinweis: B_p hat reine Sinus-Form (Hauptwelle), und damit ein alpha_p von pi/2
        % W.B_p = W.B_m ./ (2./pi);
        W.B_p = (W.Phi_h .* pi)  ./ (2 .* geo.tau_1p .* geo.l_i);

        % Korrektur des Strombelags A [A/mm]
        % Quelle: [Mueller08, S.579 - Formel 9.1.23c]
        W.A = (2 .* W.w_1Str .* rated.m .* emag.I_1Str) ./ (pi .* geo.D_1i.*1e3);

        % Korrektur der Stromdichte S_1 (Stator) [A/mm]
        % Produkt S*A konstant halten
        W.S_1 = (emag.A .* richt.S_1) ./ W.A;

%% C.3) Nutform und magnetischer Kreis
% #########################################################################
% #   C.3) NUTFORM UND MAGNETISCHER KREIS                                 #
% #########################################################################

        % Berechnung des Zweigstroms (Stator) I_1zw [A]
        % Quelle: [Mueller08, S.597]
        W.I_1zw = emag.I_1Str ./ W.a_1;

        % Berechnung des Leiterquerschnitts (Stator) A_1L [mm^2]
        % Quelle: [Mueller08, S.581 - Formel 9.1.28]
        W.A_1L = W.I_1zw ./ W.S_1;

        % Berechnung des Durchmessers des Leiters (Stator) d_1L [mm]
        W.d_1L = sqrt((4 .* W.A_1L) ./ pi);

        % Berechnung der Leitermasse (Stator) m_1Le [kg]
        W.V_1Le = ((W.A_1L.*1e-6) .* wick.l_1m .* W.w_1Str .* W.a_1 .* rated.m);
        W.m_1Le = W.V_1Le .* opt.Stator_Leitermaterial.rho_Le;

        % Berechnung des erforderlichen Nutquerschnitts (Stator) A_1n [mm^2]
        % Quelle: [Mueller08, S.597], [Mueller08, S.586 - Formel 9.1.37]
        W.A_1n = (W.z_1n .* W.A_1L) ./ richt.phi_1n;

        % Berechnung des maximalen Flusses im Ruecken (Stator) Phi_1r_max [Wb]
        % Quelle: [Mueller08, S.599]
        W.Phi_1r_max = (1.05.*W.Phi_delta) ./ 2;

        % Annahme Abplattungsfaktor alpha_p [-]
        % Quelle: [Mueller08, S.258]
        W.alpha_p = 1.4 * ones(length(W.N_1),1);

        % Berechnung der maximalen Induktion der Feldkurve B_max [T]
        % Quelle: [Mueller08, S.244 - Formel 2.5.21]
        W.B_max = W.alpha_p .* W.B_m;

        % Richtwert fuer die Nutschlitzhoehe (Stator) h_1ns [mm]
        % Quelle: [Meyer18, S.94 - Tabelle 5.3]
        h_ns_data = [0.05 0.2 0.3 0.45; 0.5 1 1 2];
        h_ns_fun = polyfit(h_ns_data(1,:),h_ns_data(2,:),1);
        NG.h_1ns = polyval(h_ns_fun,geo.D_1i);
        clear h_ns_data h_ns_fun

        % Obergrenze Iterationen
        NG.maxIter = 1e6;
        NG.iter = 0;

        % Zwischenspeichern
        saveNG = NG;

        % Iterative Generierung der Nut- und Zahnform (Stator)
        if(strcmp(opt.Nutform_Stator,'Trapezform (eckig)'))
            for i = 1:length(W.N_1)
                % Zwischenspeichern
                NG = saveNG;

                % Abschaetzung der Nutbreite (Stator) b_1n [mm]
                % Quelle: [Mueller08, S.597 - Formel 9.2.7a]
                % NG.b_1n_alt = W.tau_1n(i)/2 * 1e3;

                % Abschaetzung der Nuthoehe (Stator) h_1n [mm]
                % Quelle: [Mueller08, S.597 - Formel 9.2.7b]
                % NG.h_1n_alt = W.A_1n(i) / NG.b_1n_alt;

                % Definition der Parameter zur eindeutigen Beschreibung der Nutform
                % Quelle: [Meyer09, S.53 - Tabelle 5.3]

                % Vorgabe einer minimalen Nutbreite (Stator) b_1n_o, b_1n_u, b_1n_m (oben, unten, mitte) [mm]
                % Quelle: [Meyer09, S.55 - Formeln 5.49,5.50,5.51]
                if((0.2 * W.tau_1n(i)*1e3)<NG.b_1ns)
                    NG.b_1n_u = NG.b_1ns;                                    
                else
                    NG.b_1n_u = 0.2 * W.tau_1n(i)*1e3;
                end
                NG.b_1n_o = NG.b_1n_u;
                NG.b_1n_m = (NG.b_1n_o + NG.b_1n_u) / 2;

                % Festlegung der Dicke der Nutisolierung (Stator) d_1iso [mm]
                % Quelle: [Mueller Berechnungsbeispiel], [Mueller06, S.224 - Abbildung 2.3.10]
                NG.d_1iso = 0.3;

                % Festlegung des Nutkeilwinkels (Stator) alpha_1nk [rad]
                % Quelle: Annahme
                NG.alpha_1nk = pi/6;

                % Berechnung der Hoehe des Nutkeilwinkels (Stator) h_1k [mm]
                NG.h_1k = tan(NG.alpha_1nk) * 0.5 * (NG.b_1n_u - NG.b_1ns);

                % Berechnung der Hoehe bis Wicklungsanfang (Stator) h_1nk [mm]
                % V/A: Nutkeil zu h_k+0.5mm angenommen
                NG.h_1nk = NG.h_1ns + NG.h_1k + 0.5;

                % Vorgabe einer minimalen Nuthoehe (Stator) h_1n [mm]
                % Quelle: Annahme
                NG.h_1n = NG.h_1nk + 1;

                % Berechnung der Zahnbreite oben (Stator) b_1z_o [mm]
                NG.b_1z_o = (geo.D_1i*1e3 + 2*NG.h_1n) * pi / W.N_1(i) - NG.b_1n_o;

                % Berechnung der Zahnbreite unten (Stator) b_1z_u [mm]
                NG.b_1z_u = (geo.D_1i*1e3 + 2*NG.h_1ns + 2*NG.h_1k) * pi / W.N_1(i) - NG.b_1n_u;

                % Berechnung der Zahnbreite mitte (Stator) b_1z_m [mm]
                NG.b_1z_m = (NG.b_1z_o + NG.b_1z_u) / 2;

                % Berechnung der Rueckenhoehe (Stator) h_1r [mm]
                NG.h_1r = (geo.D_1a_max - geo.D_1i)/2*1e3 - NG.h_1n;

                % Berechnung der tatsaechlichen Nutflaeche (Stator) A_1n_tat [mm^2]
                % Quelle: Flaechenberechnung Trapez
                W.A_1n_tat(i) = 0.5 * (NG.b_1n_u + NG.b_1n_o) * (NG.h_1n - NG.h_1nk);

                while(W.A_1n_tat(i)<W.A_1n(i) && NG.iter<NG.maxIter)
                    % Berechnung der Zahninduktion oben (Stator) B_1z_o [T]
                    % Quelle: [Mueller08, S.215 - Formel 2.4.3]
                    W.B_1z_o(i) = (W.B_max(i) * W.tau_1n(i) * geo.l_i) / (NG.b_1z_o/1e3 * richt.phi_1Fe * geo.l_Fe);

                    % Berechnung der Zahninduktion unten (Stator) B_1z_u [T]
                    % Quelle: [Mueller08, S.215 - Formel 2.4.3]
                    W.B_1z_u(i) = (W.B_max(i) * W.tau_1n(i) * geo.l_i) / (NG.b_1z_u/1e3 * richt.phi_1Fe * geo.l_Fe);

                    % Berechnung der Zahninduktion mitte (Stator) B_1z_m [T]
                    % Quelle: [Mueller08, S.215 - Formel 2.4.3]
                    W.B_1z_m(i) = (W.B_max(i) * W.tau_1n(i) * geo.l_i) / (NG.b_1z_m/1e3 * richt.phi_1Fe * geo.l_Fe);

                    % Berechnung der Rueckeninduktion (Stator) B_1r [T]
                    % Quelle: [Mueller08, S.581 - Formel 9.1.30]
                    W.B_1r(i) = W.Phi_1r_max(i) / (NG.h_1r*1e-3 * richt.phi_1Fe * geo.l_Fe);

                    % Nutform-Algorithmus
                    % Quelle: [Meyer09, S.56 - Abbildung 5.8]
                    if(W.B_1z_o(i) < richt.B_1z_max || W.B_1z_u(i) < richt.B_1z_max)
                        if(W.B_1z_o(i) < W.B_1z_u(i))
                            NG.b_1n_o = NG.b_1n_o + 0.1;
                        else
                            NG.b_1n_u = NG.b_1n_u + 0.1;
                        end
                    elseif(W.B_1r(i) < W.B_1z_o(i) && W.B_1r(i) < W.B_1z_u(i))
                        NG.h_1n = NG.h_1n + 0.1;
                    elseif(W.B_1z_o(i) < W.B_1z_u(i))
                        NG.b_1n_o = NG.b_1n_o + 0.1;
                    else
                        NG.b_1n_u = NG.b_1n_u + 0.1;
                    end

                    % Update Statorgeometrie
                    NG.b_1n_m = (NG.b_1n_o + NG.b_1n_u) / 2;
                    NG.b_1z_o = (geo.D_1i*1e3 + 2*NG.h_1n) * pi / W.N_1(i) - NG.b_1n_o;
                    NG.b_1z_u = (geo.D_1i*1e3 + 2*NG.h_1ns + 2*NG.h_1k) * pi / W.N_1(i) - NG.b_1n_u;
                    NG.b_1z_m = (NG.b_1z_o + NG.b_1z_u) / 2;
                    NG.h_1r = (geo.D_1a_max - geo.D_1i)/2*1e3 - NG.h_1n/1e3;
                    NG.h_1k = tan(NG.alpha_1nk) * 0.5 * (NG.b_1n_u - NG.b_1ns);
                    NG.h_1nk = NG.h_1ns + NG.h_1k + 0.5;

                    % Ueberpruefen ob entartet
                    if NG.b_1z_o<=0 || NG.b_1z_u<=0 || NG.b_1z_m<=0 || NG.h_1r<=0
                        error('Geometrie ist entartet');
                    end

                    % Berechnung der tatsaechlichen Nutflaeche (Stator) A_1n_tat [mm^2]
                    % Quelle: Flaechenberechnung Trapez
                    W.A_1n_tat(i) = 0.5 * (NG.b_1n_u + NG.b_1n_o) * (NG.h_1n - NG.h_1nk);

                    NG.iter = NG.iter + 1;
                end

                % Kontrolle Iterationen
                if(NG.iter>NG.maxIter)
                    error('Kontrolle Nutraumbilanz notwendig!')
                end

                % Anpassen der Rueckenhoehe h_1r an die maximale Rueckeninduktion (Stator) [mm]
                % Quelle: [Mueller08, S.581 - Formel 9.1.30]
                if(W.B_1r(i)<richt.B_1r_max)
                    NG.h_1r = W.Phi_1r_max(i) / (geo.l_Fe * richt.phi_1Fe * richt.B_1r_max) * 1e3;
                end

                % Berechnung der Rueckeninduktion (Stator) B_1r [T]
                % Quelle: [Mueller08, S.581 - Formel 9.1.30]
                W.B_1r(i) = W.Phi_1r_max(i) / (NG.h_1r*1e-3 * richt.phi_1Fe * geo.l_Fe);

                % Berechnung der relativen Rueckenhoehe (Stator) h_1r_rel [-]
                % Quelle: [Mueller08, S.586]
                NG.h_1r_rel = NG.h_1r*1e-3 / geo.tau_1p;

                % Berechnung des Aussendurchmessers (Stator) D_1a [m]
                % Quelle: [Mueller, S.599 - Formel 9.2.15a]
                W.D_1a(i) = geo.D_1i + 2 * (NG.h_1r*1e-3 + NG.h_1n*1e-3);

                % Speichern der Nutgeometrie
                W.Nut_1(i) = NG;

                clear NG
            end
        else
            error('Ungueltige Eingabe bei Variable "opt.Nutform_Stator"');
        end
    end

    if(i_Wicklung==1)
        wick.table_1pos = W;
    else
        wick.table_1pos = [wick.table_1pos; W];
    end
    clear W saveNG

end

clear Wicklungstypen i_Wicklung i j

if(~isfield(wick,'table_1pos'))
    error('Keine gueltige Wicklung gefunden"');
end

% Modus Wicklungsauslegung
if(strcmp(opt.Mode_Wicklung,'Klassisch'))
    % Wicklung A: Ganzlochwicklung, Einschicht, ungesehnt, gezont
    % Wicklung B: Ganzlochwicklung, Zweischicht, gesehnt, gezont
    % Wicklung C: Bruchlochwicklung, Zweischicht, gesehnt, gezont, Zahnspulenwicklung (q<1) wird dann umgesetzt, wenn der Durchmesser q>1 nicht zulaesst
    
    switch opt.Wicklungstyp
        case 'A'
            wick.table_1sel = wick.table_1pos(wick.table_1pos.Sehnung_1==1,:);
        case {'B', 'C'}
            % Oberwellenunterdrueckung -> 5/6-Sehnung zur optimalen Oberwellenunterdrueckung
            % Quelle: [Meyer09, S.47]
            Sehnung_opt = 5/6;
            wick.table_1sel = wick.table_1pos(wick.table_1pos.Sehnung_1==Sehnung_opt,:);
            if(isempty(wick.table_1sel))
                Sehnung_opt = 2/3;
                wick.table_1sel = wick.table_1pos(wick.table_1pos.Sehnung_1==Sehnung_opt,:);
                if(isempty(wick.table_1sel))
                    error('Keine gueltige Wicklung gefunden"');
                end
            end
    end
    
    % Auswahl der Nutzahl N [-]
    if(max(wick.table_1pos.q_1)<1) % Zahnspulenwicklung
        % weitere Bedingung fuer die Nutzahl (Ruecksicht auf den Wicklungsfaktor der Hauptwelle)
        % Quelle: [Mueller08, S.76]
        wick.table_1sel = wick.table_1sel(:,wick.table_1sel.N_1<=(3*rated.p) & wick.table_1sel.N_1(1,:)>=(1.5*rated.p));

        % wenn eine Zahnspulenwicklung (q<1) erforderlich ist, dann soll q
        % moeglichst nahe an 1/m sein
        % Quelle: [Binder12, S.94]
        [~,idx] = min(abs((1/rated.m) - wick.table_1sel.q_1));
        wick.table_1sel = wick.table_1sel(idx,:);
    else
        % grosse Lochzahl gewuenscht zur optimalen Oberwellenunterdrueckung -> N muss maximal werden
        % Quelle: [Meyer09, S.45-46]
        wick.table_1sel = wick.table_1sel(max(wick.table_1sel.N_1)==wick.table_1sel.N_1,:);
        
        % minimale Abweichung auswaehlen und Parameter entsprechend uebernehmen,
        % falls es mehrere Minima gibt, wird die Loesung mit der geringsten
        % Parallelverschaltung gewaehlt
        [~, idx] = min(wick.table_1sel.err_1);
        wick.table_1sel = wick.table_1sel(idx,:);
    end
    var1 = table2struct(wick.table_1sel);
    var2 = [fieldnames(wick)' fieldnames(var1)'; struct2cell(wick)' struct2cell(var1)'];
    wick = struct(var2{:});

    clear Sehnung_opt idx var1 var2
    
    wick = TingleyAlg(rated, wick);
    
elseif(strcmp(opt.Mode_Wicklung,'Optimierung'))
    % Wicklungen aussortieren, die in mindestens einem Kriterium
    % (Oberwellenstreuung, Hauptwellenwicklungsfaktor, Aussendurchmesser und
    % Abweichung von der optimalen Strangwindungszahl) schlechter sind als
    % eine andere Wicklung und in den anderen Kriterien hoechstens
    % aequivalent zur Vergleichswicklung sind
    del=0;
    for i = 1:length(wick.table_1pos.N_1)
        for j = 1:length(wick.table_1pos.N_1)
            if(i==j)
            else
            if(wick.table_1pos.sigma_1o(i)>=wick.table_1pos.sigma_1o(j) && wick.table_1pos.xi_1p(i)<=wick.table_1pos.xi_1p(j) && abs(wick.table_1pos.err_1(i))>=abs(wick.table_1pos.err_1(j)) && wick.table_1pos.D_1a(i)>=wick.table_1pos.D_1a(j) && wick.table_1pos.m_1Le(i)>=wick.table_1pos.m_1Le(j))
                if(wick.table_1pos.sigma_1o(i)==wick.table_1pos.sigma_1o(j) && wick.table_1pos.xi_1p(i)==wick.table_1pos.xi_1p(j) && abs(wick.table_1pos.err_1(i))==abs(wick.table_1pos.err_1(j)) && wick.table_1pos.D_1a(i)==wick.table_1pos.D_1a(j) && wick.table_1pos.m_1Le(i)==wick.table_1pos.m_1Le(j))
                else
                    del(i) = i;
                    break;
                end
            end
            end
        end
    end
    del = unique(del);
    del = del(del~=0);
    del = flipud(del');
    wick.table_1sel = wick.table_1pos;
    for i = 1:length(del)
        wick.table_1sel(del(i),:) = [];
    end    
    
    % Geringe Abweichung zur optimalen Strangwindungszahl
    % V/A: max. 10% Abweichung von optimaler Strangwindungszahl
    wick.table_1sel = wick.table_1sel(abs(wick.table_1sel.err_1)<=((wick.table_1sel.w_1Str_opt*1.1)-wick.table_1sel.w_1Str_opt),:);
    
    if(length(wick.table_1sel.N_1)>7)
        % % Hoher Wicklungsfaktor der Hauptwelle
        % V/A: Wicklungsfaktor Hauptwelle >= 90 %
        if(max(wick.table_1sel.xi_1p)>=0.9)
            wick.table_1sel = wick.table_1sel(wick.table_1sel.xi_1p>=0.9,:);
        end
    end
    
    if(length(wick.table_1sel.N_1)>7)
        % Geringe Oberwellenstreuung
        % V/A: Oberwellenstreuung <= 10 % bei verteilten Wicklungen
        if(min(wick.table_1sel.sigma_1o)<=0.1 && max(wick.table_1sel.q_1)>1)
            wick.table_1sel = wick.table_1sel(wick.table_1sel.sigma_1o<=0.1,:);
        end
    end
    
    % Kleiner Aussendurchmesser
    %[counts,centers] = hist(W_1.D_1a);
    %W_1 = W_1(W_1.D_1a<=centers(1,2),:);
    %clear counts centers
    
    set(handles.figure1, 'currentaxes', handles.axes_Auswahl_Wicklung);
    cla(handles.axes_Auswahl_Wicklung)
    for i = 1:length(wick.table_1sel.N_1)
        P(i,:) = [wick.table_1sel.xi_1p(i) wick.table_1sel.sigma_1o(i)*1e2 wick.table_1sel.err_1(i) wick.table_1sel.B_m(i) wick.table_1sel.D_1a(i)*1e3 wick.table_1sel.m_1Le(i)];
    end
    P_labels = [{'xi_{1p}'};{'sigma_{1o} in %'};{'err_1 (abs)'};{'B_m in T'};{'D_{1a} in mm'};{'m_{1Le} in kg'}];
    axes_interval = 5;
	axes_precision = 3;
    spider_plot(P, P_labels, axes_interval, axes_precision,...
                'Marker', 'o',...
                'LineWidth', 2,...
                'MarkerSize', 5);
            
    % Title properties
    title(handles.axes_Auswahl_Wicklung,'Wicklungen','interpreter','latex','FontSize', 20)
    
    % Legend properties
    lgd_String = '';
    for i = 1:length(wick.table_1sel.N_1)
        lgd_String = [lgd_String; {[wick.table_1sel.Wicklungstyp_1{i} ', N_1=' num2str(wick.table_1sel.N_1(i)) ', q_1=' num2str(wick.table_1sel.q_1(i)) ', Sehnung_1=' num2str(wick.table_1sel.Sehnung_1(i)) ', a_1=' num2str(wick.table_1sel.a_1(i))]}];
    end
    lgd = legend(lgd_String, 'Location', 'southoutside');
    if(length(wick.table_1sel.N_1)>10)
        lgd.NumColumns = 2;
    end
    
    % Wicklung auswaehlen
    handles.popupmenu_Auswahl_Wicklung.String = ['-'; lgd_String];
    handles.popupmenu_Auswahl_Wicklung.Value = 1;
    set(handles.uipanel_Auswahl_Wicklung,'Position',[0,0,742,826])
    set(handles.uipanel_Auswahl_Wicklung,'Visible','on')
    set(handles.uipanel_Ergebnisse,'Visible','off')
    handles.text_Wicklungslayout.Visible = 'off';
    handles.edit_Wicklungslayout.Visible = 'off';
    handles.pushbutton_Wicklungslayout.Visible = 'off';
    handles.text_y_1v.Visible = 'off';
    handles.edit_y_1v.Visible = 'off';
    handles.axes_Auswahl_Wicklung.Visible = 'off';
    handles.toolbar.Children(1).Enable = 'off';
    handles.toolbar.Children(2).Enable = 'off';
    handles.toolbar.Children(3).Enable = 'off';
    handles.toolbar.Children(4).Enable = 'off';

    uiwait(gcf);
    
    wick.table_1sel = wick.table_1sel(handles.popupmenu_Auswahl_Wicklung.Value-1,:);
    var1 = table2struct(wick.table_1sel);
    var2 = [fieldnames(wick)' fieldnames(var1)'; struct2cell(wick)' struct2cell(var1)'];
    wick = struct(var2{:});
    
    set(handles.uipanel_Auswahl_Wicklung,'Position',[1484,0,742,826])
    set(handles.uipanel_Auswahl_Wicklung,'Visible','off')
    set(handles.uipanel_Ergebnisse,'Visible','on')
    handles.toolbar.Children(1).Enable = 'on';
    handles.toolbar.Children(2).Enable = 'on';
    handles.toolbar.Children(3).Enable = 'on';
    handles.toolbar.Children(4).Enable = 'on';
    
    cla(handles.axes_Auswahl_Wicklung,'reset')
    handles.popupmenu_Auswahl_Wicklung.Value = 1;
    handles.popupmenu_Auswahl_Wicklung.String = {'-'};
    
    clear i j del P P_labels axes_interval axes_precision lgd_String var1 var2
    
    wick = TingleyAlg(rated, wick);
    
elseif(strcmp(opt.Mode_Wicklung,'Manuell'))
    var1 = table2struct(wick.table_1pos);
    temp = wick.Matrix_lay;
    wick = rmfield(wick,'Matrix_lay');
    var2 = [fieldnames(wick)' fieldnames(var1)'; struct2cell(wick)' struct2cell(var1)'];
    wick = struct(var2{:});
    wick.Matrix_lay = temp;
    
    clear var1 var2 temp
else
    error('Ungueltige Eingabe bei Variable "opt.Mode_Wicklung"');
end

% Data handling
emag.w_1Str_opt = wick.w_1Str_opt;
wick = rmfield(wick,'w_1Str_opt');
emag.w_1Str = wick.w_1Str;
wick = rmfield(wick,'w_1Str');
emag.Phi_h = wick.Phi_h;
wick = rmfield(wick,'Phi_h');
emag.Phi_delta = wick.Phi_delta;
wick = rmfield(wick,'Phi_delta');
emag.misc.alpha_p = wick.alpha_p;
wick = rmfield(wick,'alpha_p');
emag.B_m = wick.B_m;
wick = rmfield(wick,'B_m');
emag.B_p = wick.B_p;
wick = rmfield(wick,'B_p');
emag.A = wick.A;
wick = rmfield(wick,'A');
emag.I_1zw = wick.I_1zw;
wick = rmfield(wick,'I_1zw');
emag.Phi_1r_max = wick.Phi_1r_max;
wick = rmfield(wick,'Phi_1r_max');
emag.B_max = wick.B_max;
wick = rmfield(wick,'B_max');
emag.B_1z_o = wick.B_1z_o;
wick = rmfield(wick,'B_1z_o');
emag.B_1z_u = wick.B_1z_u;
wick = rmfield(wick,'B_1z_u');
emag.B_1z_m = wick.B_1z_m;
wick = rmfield(wick,'B_1z_m');
emag.B_1r = wick.B_1r;
wick = rmfield(wick,'B_1r');
geo.Nut.A_1n = wick.A_1n;
wick = rmfield(wick,'A_1n');
geo.Nut.A_1n_tat = wick.A_1n_tat;
wick = rmfield(wick,'A_1n_tat');
geo.D_1a = wick.D_1a;
wick = rmfield(wick,'D_1a');
geo.Nut_1 = wick.Nut_1;
wick = rmfield(wick,'Nut_1');

var2 = [fieldnames(geo.Nut_1)' fieldnames(geo.Nut)'; struct2cell(geo.Nut_1)' struct2cell(geo.Nut)'];
geo.Nut_1 = struct(var2{:});
geo = rmfield(geo,'Nut');
clear var2

%% D) Rotorentwurf
% #########################################################################
% #   D) ROTORENTWURF                                                     #
% #########################################################################

%% D.1) Wicklungsschema
% #########################################################################
% #   D.1) WICKLUNGSSCHEMA                                                #
% #########################################################################

if(strcmp(opt.Maschinenausfuehrung,'Kaefiglaeufer'))
    % Bei der Wahl der Nutzahl im Rotor muss folgendes beachtet werden:
    % Quelle: [Mueller08, S.122]
    % Unterdrueckung asynchroner Oberwellenmomente:
    % N_2 <= 1.25*N_1
    % Unterdrueckung synchroner Oberwellen:
    % N_2 ~= N_1 && N_2 ~= N_1+-2p
    % Unterdrueckung von Geraeuschen bzw. Ruettelkraeften:
    % N_2 ~= N_1 && N_2 ~= N_1+-1 && N_2~= N_1+- (2p+-1)
    
    % Berechnung der Nutzahl (Rotor) N_2 [-]
    % Quelle: [Mueller08, S.122 - Formel 1.2.102]
    N_2max = wick.N_1 + 4*rated.p;
    N_2min = wick.N_1 - 4*rated.p;
    
    % N_2 < N_1 fuer kleinere zusaetzliche Verluste
    % Quelle: [Mueller08, S.123]
    if(N_2min>0)
        wick.N_2 = N_2min;
    else
        wick.N_2 = N_2max;
    end
    clear N_2max N_2min
    
    % Berechnung des Stroms (Rotor) I_2 [A]
    % Quelle: [Mueller08, S.598 - Formel 9.2.8b]
    emag.I_2 = (rated.m * emag.w_1Str * wick.xi_1p)/((wick.N_2/(2*rated.p)) * 1) * emag.I_1Str * emag.cos_phi_N;

    % Berechnung des Stabstroms (Rotor) I_2s [A]
    % Quelle: [Mueller08, S.598 - Formel 9.2.9]
    emag.I_2s = emag.I_2 / rated.p;

    % Berechnung des Ringstroms (Rotor) I_2r [A]
    % Quelle: [Mueller08, S.598 - Formel 9.2.10]
    emag.I_2r = emag.I_2 / (2 * rated.p * sin((pi*rated.p)/wick.N_2));
    
    % Berechnung des Stabquerschnitts (Rotor) A_2s [mm^2]
    % Quelle: [Mueller08, S.581 - Formel 9.1.28]
    wick.A_2s = emag.I_2s / richt.S_2s;

    % Berechnung des Ringquerschnitts (Rotor) A_2r [mm^2]
    % Quelle: [Mueller08, S.581 - Formel 9.1.28]
    wick.A_2r = emag.I_2r ./ richt.S_2r;
    
    % Berechnung des erforderlichen Nutquerschnitts (Rotor) A_2n [mm^2]
    % Quelle: [Mueller08, S.603], [Mueller08, S.586 - Formel 9.1.37]
    geo.Nut_2.A_2n = wick.A_2s;
    
elseif(strcmp(opt.Maschinenausfuehrung,'Schleifringlaeufer'))
    % Bei der Wahl der Nutzahl im Rotor muss folgendes beachtet werden:
    % Quelle: [Mueller08, S.122]
    % Unterdrueckung asynchroner Oberwellenmomente:
    % N_2 <= 1.25*N_1
    % Unterdrueckung synchroner Oberwellen:
    % N_2 ~= N_1 && N_2 ~= N_1+-2p
    % Unterdrueckung von Geraeuschen bzw. Ruettelkraeften:
    % N_2 ~= N_1 && N_2 ~= N_1+-1 && N_2~= N_1+- (2p+-1)
    
    % Berechnung der Nutzahl (Rotor) N_2 [-]
    % Quelle: [Mueller08, S.122 - Formel 1.2.101]
    q_2max = q_1 + 1;
    q_2min = q_1 - 1;
    N_2max = 2 .* rated.m .* rated.p .* q_2max;
    N_2min = 2 .* rated.m .* rated.p .* q_2min;
    
    % N_2 < N_1 fuer kleinere zusaetzliche Verluste
    % Quelle: [Mueller08, S.123]
    wick.N_2 = N_2min;
    wick.q_2 = q_2min;
    
    clear N_2_max N_2_min q_2_max q_2_min
    
   	% Berechnung des Stroms (Rotor) I_2 [A]
    % Quelle: [Mueller08, S.598 - Formel 9.2.8b]
    emag.I_2 = (rated.m .* emag.w_1Str .* wick.xi_1p)./((wick.N_2./(2.*rated.p)) .* 1) .* emag.I_1Str .* emag.cos_phi_N;
    
    %to-do
else
    error('Ungueltige Eingabe bei Variable "opt.Maschinenausfuehrung"')
end

%% D.2) Nutform und magnetischer Kreis
% #########################################################################
% #   D.2) NUTFORM UND MAGNETISCHER KREIS                                 #
% #########################################################################

% Berechnung der Nutteilung (Rotor) tau_2n [m]
% Quelle: [Mueller08, S.19 - Formel 1.1.14]
wick.tau_2n = (geo.D_2a.*pi) ./ wick.N_2;

% Ueberpruefung ob tau_2n im zulaessigen Bereich
% Quelle: [Meyer09, S.46], [Pyr14] -> evtl auch kleiner (bis zu 0.007m)
%if(wick.tau_2n>=richt.tau_2n_min)
%    warning('tau_2n zu klein. Korrektur notwendig')
%end

% Berechnung der Polteilung (Rotor) tau_2p [m]
% Quelle: [Mueller08, S.577], [Mueller08, S.19 - Formel 1.1.16]
geo.tau_2p = (geo.D_2a * pi) / (2 * rated.p);

% Richtwert fuer die max. zulaessige Induktion im Ruecken (Rotor) B_2r_max [T]
% Quelle: [Mueller08, S.582 - Tabelle 9.1.5]
% richt.B_2r_max = 1.6; % zwischen 0.4 und 1.6 (...2)

% Berechnung des maximalen Flusses im Ruecken (Rotor) Phi_2r_max [Wb]
% Quelle: [Mueller08, S.599]
emag.Phi_2r_max = emag.Phi_delta / 2;

% Richtwert fuer die max. zulaessige Induktion in den Zaehnen (Rotor) B_2z_max [T]
% Quelle: [Mueller08, S.582 - Tabelle 9.1.5]
% richt.B_2z_max = 1.7; % zwischen 1.5 und 2.2

% Richtwert fuer die Nutschlitzbreite (Rotor) b_2ns [mm]
% Quelle: [Meyer18, S.94 - Tabelle 5.3]
b_ns_data = [0.05 0.2 0.3 0.45; 1.5 3 3 5];
b_ns_fun = polyfit(b_ns_data(1,:),b_ns_data(2,:),1);
geo.Nut_2.b_2ns = polyval(b_ns_fun,geo.D_2a);
clear b_ns_data b_ns_fun

% Richtwert fuer die Nutschlitzhoehe (Rotor) h_2ns [mm]
% Quelle: [Meyer18, S.94 - Tabelle 5.3]
h_ns_data = [0.05 0.2 0.3 0.45; 0.5 1 1 2];
h_ns_fun = polyfit(h_ns_data(1,:),h_ns_data(2,:),1);
geo.Nut_2.h_2ns = polyval(h_ns_fun,geo.D_2a);
clear h_ns_data h_ns_fun

% Definition der Parameter zur eindeutigen Beschreibung der Nutform
% Quelle: [Meyer09, S.53 - Tabelle 5.3]

% Vorgabe einer minimalen Nutbreite (Rotor) b_2n_o, b_2n_u, b_2n_m (oben, unten, mitte) [mm]
% Quelle: [Meyer09, S.55 - Formeln 5.49,5.50,5.51]
if((0.2 * wick.tau_2n*1e3)<geo.Nut_2.b_2ns)
    geo.Nut_2.b_2n_u = geo.Nut_2.b_2ns;                                    
else
    geo.Nut_2.b_2n_u = 0.2 * wick.tau_2n*1e3;
end
geo.Nut_2.b_2n_o = geo.Nut_2.b_2n_u;
geo.Nut_2.b_2n_m = (geo.Nut_2.b_2n_o + geo.Nut_2.b_2n_u) / 2;

% Festlegung der Dicke der Nutisolierung (Rotor) d_2iso [mm]
% Quelle: [Mueller Berechnungsbeispiel], [Mueller06, S.224 - Abbildung 2.3.10]
geo.Nut_2.d_2iso = 0.3;

% Festlegung des Nutkeilwinkels (Rotor) alpha_2nk [rad]
% Quelle: Annahme
geo.Nut_2.alpha_2nk = pi/6;

% Berechnung der Hoehe des Nutkeilwinkels (Rotor) h_2k [mm]
geo.Nut_2.h_2k = tan(geo.Nut_2.alpha_2nk) * 0.5 * (geo.Nut_2.b_2n_u - geo.Nut_2.b_2ns);

% Berechnung der Hoehe bis Wicklungsanfang (Rotor) h_2nk [mm]
% V/A: Nutkeil zu h_k+0.5mm angenommen
geo.Nut_2.h_2nk = geo.Nut_2.h_2ns + geo.Nut_2.h_2k + 0.5;

% Vorgabe einer minimalen Nuthoehe (Rotor) h_2n [mm]
% Quelle: Annahme
geo.Nut_2.h_2n = geo.Nut_2.h_2nk + 1;

% Berechnung der Zahnbreite oben (Rotor) b_2z_o [mm]
geo.Nut_2.b_2z_o = (geo.D_2a*1e3 - 2*geo.Nut_2.h_2n) * pi / wick.N_2 - geo.Nut_2.b_2n_o;

% Berechnung der Zahnbreite unten (Rotor) b_2z_u [mm]
geo.Nut_2.b_2z_u = (geo.D_2a*1e3 - 2*geo.Nut_2.h_2ns - 2*geo.Nut_2.h_2k) * pi / wick.N_2 - geo.Nut_2.b_2n_u;

% Berechnung der Zahnbreite mitte (Rotor) b_2z_m [mm]
geo.Nut_2.b_2z_m = (geo.Nut_2.b_2z_o + geo.Nut_2.b_2z_u) / 2;

% Berechnung der Rueckenhoehe (Rotor) h_2r [mm]
geo.Nut_2.h_2r = geo.D_2a/2*1e3 - geo.Nut_2.h_2n;

% Berechnung der tatsaechlichen Nutflaeche (Rotor) A_2n_tat [mm^2]
% Quelle: Flaechenberechnung Trapez
geo.Nut_2.A_2n_tat = 0.5 * (geo.Nut_2.b_2n_u + geo.Nut_2.b_2n_o) * (geo.Nut_2.h_2n - geo.Nut_2.h_2nk);

% Obergrenze Iterationen
geo.Nut_2.maxIter = 1e6;
geo.Nut_2.iter = 0;

while(geo.Nut_2.A_2n_tat<geo.Nut_2.A_2n && geo.Nut_2.iter<geo.Nut_2.maxIter)
    % Berechnung der Zahninduktion oben (Rotor) B_2z_o [T]
    % Quelle: [Mueller08, S.215 - Formel 2.4.3]
    emag.B_2z_o = (emag.B_max * wick.tau_2n * geo.l_i) / (geo.Nut_2.b_2z_o/1e3 * richt.phi_2Fe * geo.l_Fe);

    % Berechnung der Zahninduktion unten (Rotor) B_2z_u [T]
    % Quelle: [Mueller08, S.215 - Formel 2.4.3]
    emag.B_2z_u = (emag.B_max * wick.tau_2n * geo.l_i) / (geo.Nut_2.b_2z_u/1e3 * richt.phi_2Fe * geo.l_Fe);

    % Berechnung der Zahninduktion mitte (Rotor) B_2z_m [T]
    % Quelle: [Mueller08, S.215 - Formel 2.4.3]
    emag.B_2z_m = (emag.B_max * wick.tau_2n * geo.l_i) / (geo.Nut_2.b_2z_m/1e3 * richt.phi_2Fe * geo.l_Fe);

    % Berechnung der Rueckeninduktion (Rotor) B_2r [T]
    % Quelle: [Mueller08, S.581 - Formel 9.1.30]
    emag.B_2r = emag.Phi_2r_max / (geo.Nut_2.h_2r*1e-3 * richt.phi_2Fe * geo.l_Fe);

    % Nutform-Algorithmus
    % Quelle: [Meyer09, S.56 - Abbildung 5.8]
    if(emag.B_2z_u < richt.B_2z_max || emag.B_2z_o < richt.B_2z_max)
        if(emag.B_2z_u < emag.B_2z_o)
            geo.Nut_2.b_2n_u = geo.Nut_2.b_2n_u + 0.1;
        else
            geo.Nut_2.b_2n_o = geo.Nut_2.b_2n_o + 0.1;
        end
    elseif(emag.B_2r < emag.B_2z_u && emag.B_2r < emag.B_2z_o && geo.Nut_2.b_2n_o > 1.0)
        geo.Nut_2.h_2n = geo.Nut_2.h_2n + 0.1;
        while(emag.B_2z_o > richt.B_2z_max && geo.Nut_2.b_2n_o > 1.0)
            geo.Nut_2.b_2n_o = geo.Nut_2.b_2n_o - 0.1;
            geo.Nut_2.b_2z_o = (geo.D_2a*1e3 - 2*geo.Nut_2.h_2n) * pi / wick.N_2 - geo.Nut_2.b_2n_o;
            emag.B_2z_o = (emag.B_max * wick.tau_2n * geo.l_i) / (geo.Nut_2.b_2z_o/1e3 * richt.phi_2Fe * geo.l_Fe);
        end
    elseif(emag.B_2z_u < emag.B_2z_o)
        geo.Nut_2.b_2n_u = geo.Nut_2.b_2n_u + 0.1;
    else
        geo.Nut_2.b_2n_o = geo.Nut_2.b_2n_o + 0.1;
    end

    % Update Statorgeometrie
    geo.Nut_2.b_2n_m = (geo.Nut_2.b_2n_o + geo.Nut_2.b_2n_u) / 2;
    geo.Nut_2.b_2z_o = (geo.D_2a*1e3 - 2*geo.Nut_2.h_2n) * pi / wick.N_2 - geo.Nut_2.b_2n_o;
    geo.Nut_2.b_2z_u = (geo.D_2a*1e3 - 2*geo.Nut_2.h_2ns - 2*geo.Nut_2.h_2k) * pi / wick.N_2 - geo.Nut_2.b_2n_u;
    geo.Nut_2.b_2z_m = (geo.Nut_2.b_2z_o + geo.Nut_2.b_2z_u) / 2;
    geo.Nut_2.h_2r = geo.D_2a/2*1e3 - geo.Nut_2.h_2n;
    geo.Nut_2.h_2k = tan(geo.Nut_2.alpha_2nk) * 0.5 * (geo.Nut_2.b_2n_u - geo.Nut_2.b_2ns);
    geo.Nut_2.h_2nk = geo.Nut_2.h_2ns + geo.Nut_2.h_2k + 0.5;

    % Ueberpruefen ob entartet
    if geo.Nut_2.b_2z_o<=0 || geo.Nut_2.b_2z_u<=0 || geo.Nut_2.b_2z_m<=0 || geo.Nut_2.h_2r<=0
        error('Geometrie ist entartet');
    end
    
    % Berechnung der tatsaechlichen Nutflaeche (Rotor) A_2n_tat [mm^2]
    % Quelle: Flaechenberechnung Trapez
    geo.Nut_2.A_2n_tat = 0.5 * (geo.Nut_2.b_2n_u + geo.Nut_2.b_2n_o) * (geo.Nut_2.h_2n - geo.Nut_2.h_2nk);

    geo.Nut_2.iter = geo.Nut_2.iter + 1;
end

% Kontrolle Iterationen
if(geo.Nut_2.iter>geo.Nut_2.maxIter)
    error('Kontrolle Nutraumbilanz notwendig!')
end

% Anpassen der Rueckenhoehe h_2r an die maximale Rueckeninduktion (Rotor) [mm]
% Quelle: [Mueller08, S.581 - Formel 9.1.30]
if(emag.B_2r<richt.B_2r_max)
    geo.Nut_2.h_2r = emag.Phi_2r_max / (geo.l_Fe * richt.phi_2Fe * richt.B_2r_max) * 1e3;
end

% Berechnung der Rueckeninduktion (Rotor) B_2r [T]
% Quelle: [Mueller08, S.581 - Formel 9.1.30]
emag.B_2r = emag.Phi_2r_max / (geo.Nut_2.h_2r*1e-3 * richt.phi_2Fe * geo.l_Fe);

% Berechnung der relativen Rueckenhoehe (Rotor) h_2r_rel [-]
% Quelle: [Mueller08, S.586]
geo.Nut_2.h_2r_rel = geo.Nut_2.h_2r*1e-3 / geo.tau_2p;

% Berechnung des Innendurchmessers (Rotor) D_2i [m]
geo.D_2i = geo.D_2a - 2.*geo.Nut_2.h_2n*1e-3 - 2.*geo.Nut_2.h_2r*1e-3;

%% E) Nachrechnung
% #########################################################################
% #   E) NACHRECHNUNG                                                     #
% #########################################################################

%% E.1) Nachrechnung des magnetischen Kreises
% #########################################################################
% #   E.1) NACHRECHNUNG DES MAGNETISCHEN KREISES                          #
% #########################################################################

tol_Nut = 2e-2;
b_1z_o_prev = 1e3;
b_1z_m_prev = 1e3;
b_1z_u_prev = 1e3;
b_2z_o_prev = 1e3;
b_2z_m_prev = 1e3;
b_2z_u_prev = 1e3;
maxIter_Korrektur_Nut = 1e3;
iter_Korrektur_Nut = 1;

while(1)
    if(abs(b_1z_o_prev-geo.Nut_1.b_1z_o)>tol_Nut || abs(b_1z_m_prev-geo.Nut_1.b_1z_m)>tol_Nut || abs(b_1z_u_prev-geo.Nut_1.b_1z_u)>tol_Nut || ...
            abs(b_2z_o_prev-geo.Nut_2.b_2z_o)>tol_Nut || abs(b_2z_m_prev-geo.Nut_2.b_2z_m)>tol_Nut || abs(b_2z_u_prev-geo.Nut_2.b_2z_u)>tol_Nut)
    else
        var1 = [emag.B_1z_o; emag.B_1z_m; emag.B_1z_u; emag.B_2z_o; emag.B_2z_m; emag.B_2z_u];
        var2 = [richt.B_1z_max; richt.B_1z_max; richt.B_1z_max; richt.B_2z_max; richt.B_2z_max; richt.B_2z_max];
        var1_names = {'B_1z_o'; 'B_1z_m'; 'B_1z_u'; 'B_2z_o'; 'B_2z_m'; 'B_2z_u'};
        var2_names = {'B_1z_max'; 'B_1z_max'; 'B_1z_max'; 'B_2z_max'; 'B_2z_max'; 'B_2z_max'};
        for i = 1:length(var1)
            if(var1(i) > var2(i)*1.02)
                warning(['Induktion ' var1_names{i} ' = ' num2str(var1(i)) ' ueberschreitet vorgegebenen Maximalwert ' var2_names{i} ' = ' num2str(var2(i))])
            end
        end
        clear var1 var1_names var2 var2_names
        break;
    end
    if(maxIter_Korrektur_Nut<=iter_Korrektur_Nut)
        warning('Kontrolle Korrektur Nut notwendig!')
        break;
    end
    
    b_1z_o_prev = geo.Nut_1.b_1z_o;
    b_1z_m_prev = geo.Nut_1.b_1z_m;
    b_1z_u_prev = geo.Nut_1.b_1z_u;
    b_2z_o_prev = geo.Nut_2.b_2z_o;
    b_2z_m_prev = geo.Nut_2.b_2z_m;
    b_2z_u_prev = geo.Nut_2.b_2z_u;
    
    tol_alpha_p = 2e-2;
    alpha_p_prev = 1e3;
    maxIter_Korrektur_alpha_p = 1e3;
    iter_Korrektur_alpha_p = 1;
    while(abs(alpha_p_prev-emag.misc.alpha_p)>tol_alpha_p && maxIter_Korrektur_alpha_p>iter_Korrektur_alpha_p)

        alpha_p_prev = emag.misc.alpha_p;

        % Berechnung des Carter-Faktors (Stator) k_1c [-]
        % Quelle: [Mueller08, S.203 - Formel 2.3.19], [Mueller08, S.204 - Formel 2.3.20]
        emag.misc.gamma_1 = 1 / (1 + 5 * geo.delta / (geo.Nut_1.b_1ns));
        emag.misc.k_1c = wick.tau_1n / (wick.tau_1n - emag.misc.gamma_1* geo.Nut_1.b_1ns*1e-3);

        % Berechnung des Carter-Faktors (Rotor) k_2c [-]
        % Quelle: [Mueller08, S.203 - Formel 2.3.19], [Mueller08, S.204 - Formel 2.3.20]
        emag.misc.gamma_2 = 1 / (1 + 5 * geo.delta / (geo.Nut_2.b_2ns));
        emag.misc.k_2c = wick.tau_2n / (wick.tau_2n - emag.misc.gamma_2 * geo.Nut_2.b_2ns*1e-3);

        % Berechnung des Carter-Faktors k [-]
        % Quelle: [Mueller08, S.205 - Formel 2.3.21c]
        emag.misc.k_c = emag.misc.k_1c * emag.misc.k_2c;

        % Berechnung des magnetischen Spannungsabfalls ueber dem Luftspalt V_delta [A]
        % Quelle: [Mueller08, S.209 - Formel 2.3.28]
        emag.V_delta = (1/const.mu_0) * emag.B_max * emag.misc.k_c * geo.delta*1e-3;

        % Berechnung der scheinbaren Zahninduktion (Stator) B_1zs [T]
        % Quelle: [Mueller08, S.217 - Formel 2.4.9]
        emag.B_1zs_o = (emag.B_max * wick.tau_1n * geo.l_i) / (geo.Nut_1.b_1z_o*1e-3 * richt.phi_1Fe * geo.l_Fe);
        emag.B_1zs_m = (emag.B_max * wick.tau_1n * geo.l_i) / (geo.Nut_1.b_1z_m*1e-3 * richt.phi_1Fe * geo.l_Fe);
        emag.B_1zs_u = (emag.B_max * wick.tau_1n * geo.l_i) / (geo.Nut_1.b_1z_u*1e-3 * richt.phi_1Fe * geo.l_Fe);

        % Berechnung der Zahnflaeche (Stator) A_1z [m^2]
        % Quelle: [Mueller08, S.217]
        % A_1z_o = geo.Nut_1.b_1z_o*1e-3 * richt.phi_1Fe * geo.l_Fe;
        % A_1z_m = geo.Nut_1.b_1z_m*1e-3 * richt.phi_1Fe * geo.l_Fe;
        % A_1z_u = geo.Nut_1.b_1z_u*1e-3 * richt.phi_1Fe * geo.l_Fe;

        % Berechnung der Flaeche der parallel liegenden Luftwege (Stator) A_1Luft [m^2]
        % Quelle: [Mueller08, S.217 - Formel 2.4.7]
        % A_1Luft_o = (geo.Nut_1.b_1n_o*1e-3 * geo.l) + (geo.Nut_1.b_1z_o*1e-3 * geo.misc.n_v * richt.l_v) + (geo.Nut_1.b_1z_o*1e-3 * (1-richt.phi_1Fe) * geo.l_Fe);
        % A_1Luft_m = (geo.Nut_1.b_1n_m*1e-3 * geo.l) + (geo.Nut_1.b_1z_m*1e-3 * geo.misc.n_v * richt.l_v) + (geo.Nut_1.b_1z_m*1e-3 * (1-richt.phi_1Fe) * geo.l_Fe);
        % A_1Luft_u = (geo.Nut_1.b_1n_u*1e-3 * geo.l) + (geo.Nut_1.b_1z_u*1e-3 * geo.misc.n_v * richt.l_v) + (geo.Nut_1.b_1z_u*1e-3 * (1-richt.phi_1Fe) * geo.l_Fe);

        % Berechnung der Zahninduktion (Stator) B_1z [T]
        % Quelle: [Mueller08, S.217 - Formel 2.4.8]
        % V/A: scheinbare Zahninduktion wird als Zahninduktion angenommen (Zahninduktion wird dadurch bei Werten oberhalb von 1.7 T evtl zu gross abgeschaetzt)
        % emag.B_1z_o = emag.B_1zs_o - const.mu_0*emag.H_1z_u*(A_1Luft_o/A_1z_o);
        % emag.B_1z_m = emag.B_1zs_m - const.mu_0*emag.H_1z_m*(A_1Luft_m/A_1z_m);
        % emag.B_1z_u = emag.B_1zs_u - const.mu_0*emag.H_1z_o*(A_1Luft_u/A_1z_u);
        emag.B_1z_o = emag.B_1zs_o;
        emag.B_1z_m = emag.B_1zs_m;
        emag.B_1z_u = emag.B_1zs_u;

        % Berechnung der Zahninduktion (Rotor) B_2z [T]
        % Quelle: [Mueller08, S.215 - Formel 2.4.3]
        emag.B_2z_o = (emag.B_max * wick.tau_2n * geo.l_i) / (geo.Nut_2.b_2z_o*1e-3 * richt.phi_2Fe * geo.l_Fe);
        emag.B_2z_m = (emag.B_max * wick.tau_2n * geo.l_i) / (geo.Nut_2.b_2z_m*1e-3 * richt.phi_2Fe * geo.l_Fe);
        emag.B_2z_u = (emag.B_max * wick.tau_2n * geo.l_i) / (geo.Nut_2.b_2z_u*1e-3 * richt.phi_2Fe * geo.l_Fe);

        % Ermittlung der Zahnfeldstaerke (Stator) H_1z [A/m]
        % Quelle: [Datenblatt des gewaehlten Elektrobands]
        emag.H_1z_o = interp1(opt.Stator_Eisenmaterial.B, opt.Stator_Eisenmaterial.H, emag.B_1z_o);
        emag.H_1z_m = interp1(opt.Stator_Eisenmaterial.B, opt.Stator_Eisenmaterial.H, emag.B_1z_m);
        emag.H_1z_u = interp1(opt.Stator_Eisenmaterial.B, opt.Stator_Eisenmaterial.H, emag.B_1z_u);
        if(any(isnan(emag.H_1z_o)))
            error('Datenbereich der hinterlegten BH-Kurve (Elektroband) zu klein. Bitte anderes Eisenmaterial waehlen. \n(keine Werte fuer die Zahnfeldstaerke (Stator) H_1z_o bei einer Zahninduktion (Stator) von B_1z_o = %.2f T hinterlegt)',emag.B_1z_o(1,end));
        end
        if(any(isnan(emag.H_1z_m)))
            error('Datenbereich der hinterlegten BH-Kurve (Elektroband) zu klein. Bitte anderes Eisenmaterial waehlen. \n(keine Werte fuer die Zahnfeldstaerke (Stator) H_1z_m bei einer Zahninduktion (Stator) von B_1z_m = %.2f T hinterlegt)',emag.B_1z_m(1,end));
        end
        if(any(isnan(emag.H_1z_u)))
            error('Datenbereich der hinterlegten BH-Kurve (Elektroband) zu klein. Bitte anderes Eisenmaterial waehlen. \n(keine Werte fuer die Zahnfeldstaerke (Stator) H_1z_u bei einer Zahninduktion (Stator) von B_1z_u = %.2f T hinterlegt)',emag.B_1z_u(1,end));
        end

        % Ermittlung der Zahnfeldstaerke (Rotor) H_2z [A/m]
        % Quelle: [Datenblatt des gewaehlten Elektrobands]
        emag.H_2z_o = interp1(opt.Rotor_Eisenmaterial.B, opt.Rotor_Eisenmaterial.H, emag.B_2z_o);
        emag.H_2z_m = interp1(opt.Rotor_Eisenmaterial.B, opt.Rotor_Eisenmaterial.H, emag.B_2z_m);
        emag.H_2z_u = interp1(opt.Rotor_Eisenmaterial.B, opt.Rotor_Eisenmaterial.H, emag.B_2z_u);
        if(any(isnan(emag.H_2z_o)))
            error('Datenbereich der hinterlegten BH-Kurve (Elektroband) zu klein. Bitte anderes Eisenmaterial waehlen. \n(keine Werte fuer die Zahnfeldstaerke (Rotor) H_2z bei einer Zahninduktion (Rotor) von B_2z_o = %.2f T hinterlegt)',emag.B_2z_o(1,end));
        end
        if(any(isnan(emag.H_2z_m)))
            error('Datenbereich der hinterlegten BH-Kurve (Elektroband) zu klein. Bitte anderes Eisenmaterial waehlen. \n(keine Werte fuer die Zahnfeldstaerke (Rotor) H_2z bei einer Zahninduktion (Rotor) von B_2z_m = %.2f T hinterlegt)',emag.B_2z_m(1,end));
        end
        if(any(isnan(emag.H_2z_u)))
            error('Datenbereich der hinterlegten BH-Kurve (Elektroband) zu klein. Bitte anderes Eisenmaterial waehlen. \n(keine Werte fuer die Zahnfeldstaerke (Rotor) H_2z bei einer Zahninduktion (Rotor) von B_2z_u = %.2f T hinterlegt)',emag.B_2z_u(1,end));
        end

        % Berechnung des magnetischen Spannungsabfalls ueber den Zaehnen (Stator) V_1z [A]
        % Quelle: [Mueller08, S.216 - Formel 2.4.5b]
        % Simpsonsche Regel
        emag.V_1z = (1/6) * (emag.H_1z_o + 4*emag.H_1z_m + emag.H_1z_u) * geo.Nut_1.h_1n*1e-3;

        % Berechnung des magnetischen Spannungsabfalls ueber den Zaehnen (Rotor) V_2z [A]
        % Quelle: [Mueller08, S.216 - Formel 2.4.5b]
        emag.V_2z = (1/6) * (emag.H_2z_o + 4*emag.H_2z_m + emag.H_2z_u) * geo.Nut_2.h_2n*1e-3;

        % Berechnung des magnetischen Spannungsabfalls ueber den Zaehnen V_z [A]
        % Quelle: [Mueller08, S.244 - Formel 2.5.20]
        emag.V_z = emag.V_1z + emag.V_2z;

        % Berechnung des Zahnsaettigungsfaktors k [-]
        % Quelle: [Mueller08, S.244 - Formel 2.5.20]
        emag.misc.k = emag.V_z / emag.V_delta;

        % Ermittlung des Abplattungsfaktors alpha_p [-]
        % Quelle: [Mueller08, S.245 - Bild 2.5.12]
        alpha_p_data = [0:0.1:1.6; ...
            1.57 1.49 1.45 1.41 1.38 1.36 1.34 1.32 1.3 1.28 1.27 1.265 1.26 1.25 1.245 1.24 1.23];
        emag.misc.alpha_p = interp1(alpha_p_data(1,:), alpha_p_data(2,:), emag.misc.k, 'linear');
        if(any(isnan(emag.misc.alpha_p)))
            error('Datenbereich der hinterlegten Kurve zur Abschaetzung des Abplattungsfaktors zu klein. \n(keine Werte fuer den Abplattungsfaktors alpha_p bei einem Zahnsaettigungsfaktors von k = %.2f hinterlegt)',emag.misc.k(1,end));
        end
        clear alpha_p_data

        % Berechnung der maximalen Induktion der Feldkurve B_max [T]
        % Quelle: [Mueller08, S.244 - Formel 2.5.21]
        emag.B_max = emag.B_m * emag.misc.alpha_p;

        iter_Korrektur_alpha_p = iter_Korrektur_alpha_p + 1;
    end

    if(iter_Korrektur_alpha_p>=maxIter_Korrektur_alpha_p)
        warning('Kontrolle Korrektur Abplattungsfaktor notwendig!')
    end
    
    % Korrektur Nutform
    [emag, wick, geo] = Statornut(emag, wick, geo, richt);
    [emag, wick, geo] = Rotornut(emag, wick, geo, richt);
    
    iter_Korrektur_Nut = iter_Korrektur_Nut + 1;
end
clear tol_alpha_p iter_Korrektur_alpha_p maxIter_Korrektur_alpha_p tol_B iter_Korrektur_Nut maxIter_Korrektur_Nut ...
    alpha_p_prev B_1z_o_prev B_1z_m_prev B_1z_u_prev B_2z_o_prev B_2z_m_prev B_2z_u_prev i

% Ueberpruefung ob h_1r_rel im zulaessigen Bereich
% Quelle: [Mueller08, S.586 - Tabelle 9.1.7]
if(geo.Nut_1.h_1r_rel<0.13 || geo.Nut_1.h_1r_rel>0.28)
    warning('h_1r_rel ausserhalb der Grenzen')
end

% Ueberpruefung ob h_2r_rel im zulaessigen Bereich
% Quelle: [Mueller08, S.586 - Tabelle 9.1.7]
if(geo.Nut_2.h_2r_rel<0.13 || geo.Nut_2.h_2r_rel>0.28)
    warning('h_2r_rel ausserhalb der Grenzen')
end

% Berechnung der Flaechen
% Rueckenflaeche (Stator) A_1r [m^2]
geo.A_1r = (pi/4) * (geo.D_1a^2 - (geo.D_1i + (2*geo.Nut_1.h_1n*1e-3))^2);
% Zahnflaeche (Stator) A_1z [m^2]
geo.A_1z = (pi/4) * (geo.D_1a^2 - geo.D_1i^2) - geo.A_1r - wick.N_1*(geo.Nut_1.A_1n_tat + geo.Nut_1.b_1ns*geo.Nut_1.h_1ns + geo.Nut_1.b_1n_u*(geo.Nut_1.h_1nk-geo.Nut_1.h_1k-geo.Nut_1.h_1ns) + 0.5*(geo.Nut_1.b_1ns*geo.Nut_1.b_1n_u)*geo.Nut_1.h_1k)/1000000;
% Rueckenflaeche (Rotor) A_2r [m^2]
geo.A_2r = (pi/4) * (geo.D_2a^2 - geo.D_2i^2);

% Berechnung der Volumina
% Rueckenvolumen (Stator) Vo_1r [m^3]
geo.Vo_1r = geo.A_1r * geo.l_Fe * richt.phi_1Fe; 
% Zahnvolumen (Stator) Vo_1z [m^3]
geo.Vo_1z = geo.A_1z * geo.l_Fe * richt.phi_1Fe;
% Rueckenvolumen (Rotor) Vo_1r [m^3]
geo.Vo_2r = geo.A_2r * geo.l_Fe * richt.phi_2Fe;

% Ermittlung der Rueckenfeldstaerke (Stator) H_1r [A/m]
% Quelle: [Datenblatt des gewaehlten Elektrobands]
emag.H_1r = interp1(opt.Stator_Eisenmaterial.B, opt.Stator_Eisenmaterial.H, emag.B_1r);
if(any(isnan(emag.H_1r)))
    error('Datenbereich der hinterlegten BH-Kurve (Elektroband) zu klein. Bitte anderes Eisenmaterial waehlen. \n(keine Werte fuer die Rueckenfeldstaerke (Stator) H_1r bei einer Rueckeninduktion (Stator) von B_1r = %.2f T hinterlegt)',emag.B_1r(1,end));
end

% Ermittlung der Rueckenfeldstaerke (Rotor) H_2r [A/m]
% Quelle: [Datenblatt des gewaehlten Elektrobands]
emag.H_2r = interp1(opt.Rotor_Eisenmaterial.B, opt.Rotor_Eisenmaterial.H, emag.B_2r);
if(any(isnan(emag.H_2r)))
    error('Datenbereich der hinterlegten BH-Kurve (Elektroband) zu klein. Bitte anderes Eisenmaterial waehlen. \n(keine Werte fuer die Rueckenfeldstaerke (Rotor) H_2r bei einer Rueckeninduktion (Rotor) von B_2r = %.2f T hinterlegt)',emag.B_2r(1,end));
end

% Berechnung der Polteilung (Stator) an der Grenzflaeche zwischen Zahngebiet und Ruecken
% Quelle: [Mueller08, S.227 - Formel 2.4.29]
geo.tau_1r = (geo.D_1i + 2*geo.Nut_1.h_1n*1e-3)* pi / (2*rated.p);

% Berechnung der Polteilung (Rotor) an der Grenzflaeche zwischen Zahngebiet und Ruecken
% Quelle: [Mueller08, S.227 - Formel 2.4.29]
geo.tau_2r = (geo.D_2a - 2*geo.Nut_2.h_2n*1e-3)* pi / (2*rated.p);

% % Ermittlung des Rueckenreduktionsfaktors C_r
% % Quelle: [Mueller08, S.227 - Bild 2.4.11]
% % C_1r_xx_data = [h_1r/tau_1r {-}; C_1r {-}]
% % Rohdaten Statorruecken & p = 1
% C_1r_data_p1.B_08 = [0:0.020:0.260; ...
%                      0.697 0.761 0.813 0.859 0.896 0.928 0.962 0.990 1.015 1.045 1.075 1.109 1.145 1.186];
% C_1r_data_p1.B_10 = [0:0.020:0.320; ...
%                      0.603 0.656 0.699 0.738 0.770 0.801 0.829 0.857 0.880 0.904 0.927 0.957 0.991 1.027 1.074 1.128 1.19];
% C_1r_data_p1.B_12 = [0:0.020:0.320 0.350; ...
%                      0.505 0.550 0.582 0.616 0.648 0.675 0.699 0.720 0.737 0.756 0.779 0.807 0.835 0.869 0.909 0.956 1.02 1.119];
% C_1r_data_p1.B_14 = [0:0.020:0.320 0.350; ...
%                      0.402 0.430 0.453 0.472 0.489 0.504 0.517 0.536 0.553 0.576 0.608 0.638 0.672 0.715 0.756 0.796 0.843 0.945];
% C_1r_data_p1.B_16 = [0:0.020:0.220; ...
%                      0.269 0.287 0.297 0.306 0.314 0.324 0.335 0.350 0.367 0.388 0.414 0.444];
% % Rohdaten Rotorruecken & p = 1
% C_2r_data_p1.B_08 = [0:0.020:0.300 0.315; ...
%                      0.697 0.671 0.652 0.639 0.626 0.617 0.603 0.592 0.581 0.568 0.563 0.561 0.563 0.570 0.591 0.623 0.653];
% C_2r_data_p1.B_10 = [0:0.020:0.300 0.315; ...
%                      0.603 0.588 0.575 0.564 0.553 0.541 0.528 0.515 0.504 0.491 0.486 0.484 0.491 0.512 0.546 0.603 0.653];
% C_2r_data_p1.B_12 = [0:0.020:0.300 0.315; ...
%                      0.507 0.494 0.483 0.470 0.459 0.446 0.429 0.416 0.407 0.399 0.401 0.409 0.426 0.454 0.512 0.582 0.655];
% C_2r_data_p1.B_14 = [0:0.020:0.300 0.315; ...
%                      0.402 0.387 0.376 0.367 0.359 0.352 0.346 0.339 0.335 0.333 0.339 0.351 0.371 0.405 0.464 0.554 0.655];
% C_2r_data_p1.B_16 = [0:0.020:0.300 0.315; ...
%                      0.269 0.256 0.242 0.231 0.220 0.205 0.194 0.183 0.176 0.170 0.170 0.178 0.199 0.238 0.310 0.443 0.653];
% % Rohdaten Stator- und Rotorruecken  & p = Inf
% C_r_data_pInf.B_08 = [0:0.020:0.30; ...
%                       0.695 0.699 0.701 0.708 0.718 0.729 0.739 0.752 0.769 0.786 0.803 0.820 0.837 0.854 0.871 0.888];
% C_r_data_pInf.B_10 = [0:0.020:0.30; ...
%                       0.603 0.607 0.616 0.622 0.631 0.641 0.649 0.658 0.668 0.681 0.696 0.711 0.728 0.751 0.774 0.797];
% C_r_data_pInf.B_12 = [0:0.020:0.30; ...
%                       0.507 0.509 0.511 0.519 0.526 0.532 0.542 0.551 0.564 0.576 0.589 0.608 0.627 0.647 0.667 0.687];
% C_r_data_pInf.B_14 = [0:0.020:0.30; ...
%                       0.404 0.406 0.408 0.412 0.417 0.423 0.431 0.440 0.452 0.461 0.480 0.499 0.518 0.537 0.556 0.575];
% C_r_data_pInf.B_16 = [0:0.020:0.30; ...
%                       0.269 0.272 0.278 0.280 0.284 0.290 0.296 0.303 0.313 0.323 0.333 0.343 0.353 0.363 0.373 0.383];
% 
% % Interpolation
% [~, C_1r_p1] = C_r_data_interp(C_1r_data_p1, geo.Nut_1.h_1r, geo.tau_1r, emag.B_1r, '');%Statorruecken, p=1');
% [~, C_2r_p1] = C_r_data_interp(C_2r_data_p1, geo.Nut_2.h_2r, geo.tau_2r, emag.B_2r, '');%Rotorruecken, p=1');
% [~, C_1r_pInf] = C_r_data_interp(C_r_data_pInf, geo.Nut_1.h_1r, geo.tau_1r, emag.B_1r, '');%Statorruecken, p=Inf');
% [~, C_2r_pInf] = C_r_data_interp(C_r_data_pInf, geo.Nut_2.h_2r, geo.tau_2r, emag.B_2r, '');%Rotorruecken, p=Inf');
% 
% % Berechnung des Rueckenreduktionsfaktors C_r
% % V/A: Werte von C_r fuer p>1 sind Mittelwert der Daten fuer p=1 und p=Inf
% if(rated.p>1)
%     emag.misc.C_1r = (C_1r_p1+C_1r_pInf) / 2;
%     emag.misc.C_2r = (C_2r_p1+C_2r_pInf) / 2;
% else
%     emag.misc.C_1r = C_1r_p1;
%     emag.misc.C_2r = C_2r_p1;
% end
% clear C_1r_data_p1 C_2r_data_p1 C_r_data_pInf C_1r_p1 C_2r_p1 C_1r_pInf C_2r_pInf

% Ermittlung des Rueckenreduktionsfaktors C_r - ALTERNATIV
% Quelle: [Pyr14, S.182 - Bild 3.17]
% C_1r_xx_data = [B_r {T}; C_1r {-}]
% Rohdaten Statorruecken & p = 1
C_r_data = [0:0.1:2.0; ...
            0.72 0.72 0.72 0.72 0.72 0.72 0.71 0.70 0.68 0.63 0.57 0.49 0.40 0.31 0.25 0.20 0.17 0.15 0.14 0.13 0.12];

emag.misc.C_1r = interp1(C_r_data(1,:), C_r_data(2,:),emag.B_1r,'linear');
emag.misc.C_2r = interp1(C_r_data(1,:), C_r_data(2,:),emag.B_2r,'linear');

if(any(isnan(emag.misc.C_1r)))
    error('Datenbereich der hinterlegten Kurve zur Abschaetzung des Rueckenreduktionsfaktors zu klein. \n(keine Werte fuer den Rueckenreduktionsfaktor C_1r bei einem Quotient von h_1r/tau_1r, = %.2f hinterlegt)',geo.Nut_1.h_1r*1e-3/geo.tau_1r);
end
if(any(isnan(emag.misc.C_2r)))
    error('Datenbereich der hinterlegten Kurve zur Abschaetzung des Rueckenreduktionsfaktors zu klein. \n(keine Werte fuer den Rueckenreduktionsfaktor C_2r bei einem Quotient von h_2r/tau_2r, = %.2f hinterlegt)',geo.Nut_2.h_2r*1e-3/geo.tau_2r);
end

% Berechnung des magnetischen Spannungsabfalls ueber den Ruecken (Stator) V_1r [A]
% Quelle: [Mueller08, S.226 - Bild 2.4.28a]
emag.V_1r = emag.misc.C_1r * emag.H_1r * (geo.tau_1r/2);

% Berechnung des magnetischen Spannungsabfalls ueber den Ruecken (Rotor) V_2r [A]
% Quelle: [Mueller08, S.226 - Bild 2.4.28a]
emag.V_2r = emag.misc.C_2r * emag.H_2r * (geo.tau_2r/2);

% Berechnung des magnetischen Spannungsabfalls ueber den Ruecken V_r [A]
emag.V_r = emag.V_1r + emag.V_2r;

% Berechnung der Durchflutung Theta_p [A]
% Quelle: [Mueller08, S.258]
emag.Theta_p = emag.V_delta + emag.V_z + emag.V_r;

% Berechnung des Magnetisierungsstroms [A]
% Quelle: [Mueller08, S.256 - Formel 2.6.4]
emag.I_mu = (emag.Theta_p * pi * rated.p) / (rated.m * emag.w_1Str * wick.xi_1p *sqrt(2));

%% E.2) Induktivitaeten
% #########################################################################
% #   E.2) INDUKTIVITAETEN                                                #
% #########################################################################

% Berechnung der ideelen Luftspaltlaenge unter Beruecksichtigung der Nutung delta_i [mm]
% Quelle: [Mueller08, S. 209]
geo.delta_i = emag.misc.k_c * geo.delta;

% Berechnung der ungesaettigten Hauptinduktivitaet L_1h (Stator) [H]
% Quelle: [Mueller08, S.524 - Formel 8.1.42a]
emag.L_1h = (rated.m/2) * (const.mu_0/(geo.delta_i*1e-3)) * (2/pi) * geo.tau_1p * geo.l_i * (4/pi) * ((emag.w_1Str*wick.xi_1p)^2/(2*rated.p));

% Berechnung der gesaettigten Hauptreaktanz X_1h_ges [Ohm]
% Quelle:  [Mueller11, S.354 - Formel 2.3.49]
emag.X_1h_ges = emag.E_h / emag.I_mu;
emag.L_1h_ges = emag.X_1h_ges / (2 * pi * rated.f_N);

% Berechnung des Ersatzluftspalts delta_i_ss [mm]
% Quelle: [Mueller08, S.531 - Formel 8.1.68]
geo.delta_i_ss = (1+((emag.V_z+emag.V_r)/emag.V_delta)) * geo.delta_i;

% Berechnung der gesaettigten Hautpinduktivitaet L_1h_ges [H]
% Quelle:  [Mueller08, S.532 - Formel 8.1.69]
% emag.L_1h_ges_alt = emag.L_1h * (geo.delta_i/geo.delta_i_ss);

% Richtwert fuer den relativen Streuleitwert der Zahnkopfstreuung (Stator und Rotor) lambda_1z, lambda_2z[-]
% Quelle: [Mueller08, S.324 - Bild 3.7.2]
b_1ns_delta = geo.Nut_1.b_1ns / geo.delta;
b_2ns_delta = geo.Nut_2.b_2ns / geo.delta; 
if(b_1ns_delta<3)
    lambda_z_data = [1.15,1,0.8,0.62,0.5,0.39,0.3,0.26,0.19,0.13,0.1,0.07,0.03,0;...
                    0,0.125,0.3125,0.5,0.75,1,1.3125,1.5,1.75,2,2.25,2.5,2.75,3];
    lambda_z_fun = polyfit(lambda_z_data(2,:),lambda_z_data(1,:),4);
    emag.lambda_1z = polyval(lambda_z_fun,b_1ns_delta);
    emag.lambda_2z = polyval(lambda_z_fun,b_2ns_delta);
else
    lambda_z_data = [0,-0.03,-0.05,-0.075,-0.1,-0.11,-0.12,-0.13,-0.14,-0.15,-0.16,-0.17,-0.18;...
                    3,4,5,6,7,8,9,10,11,12,13,14,16];
    lambda_z_fun = polyfit(lambda_z_data(2,:),lambda_z_data(1,:),4);
    emag.lambda_1z = polyval(lambda_z_fun,b_1ns_delta);
    emag.lambda_2z = polyval(lambda_z_fun,b_2ns_delta);
end
clear b_1ns_delta b_2ns_delta lambda_z_data lambda_z_fun

% Unterteilung der Nut zur Berechnung der Nutstreuung
% Quelle: [Mueller08, S.324 - Bild 3.7.3], [Mueller08, S.325]
geo.Nut_1.h_1ue = geo.Nut_1.h_1nk - geo.Nut_1.h_1ns - geo.Nut_1.h_1k;
geo.Nut_1.h_1l = geo.Nut_1.h_1n - geo.Nut_1.h_1nk - geo.Nut_1.d_1iso;
geo.Nut_1.b_1n = geo.Nut_1.b_1n_u; %(b_n_o + b_n_u) / 2;
geo.Nut_1.b_1k = (geo.Nut_1.b_1n + geo.Nut_1.b_1ns) / 2;

geo.Nut_2.h_2ue = geo.Nut_2.h_2nk - geo.Nut_2.h_2ns - geo.Nut_2.h_2k;
geo.Nut_2.h_2l = geo.Nut_2.h_2n - geo.Nut_2.h_2nk - geo.Nut_2.d_2iso;
geo.Nut_2.b_2n = geo.Nut_2.b_2n_u; %(b_n_o + b_n_u) / 2;
geo.Nut_2.b_2k = (geo.Nut_2.b_2n + geo.Nut_2.b_2ns) / 2;

switch wick.Wicklungstyp_1
    case {'1SGL', '1SBL'}
        % Berechnung des resultierenden relativen Nut-Zahnkopf-Streuleitwert (Stator) lambda_1nz [-]
        % Quelle: [Mueller08, S.326 - Formel 3.7.2a]
        emag.lambda_1nz = (geo.Nut_1.h_1l/(3 * geo.Nut_1.b_1n)) + (geo.Nut_1.h_1ue/geo.Nut_1.b_1n) + (geo.Nut_1.h_1k/geo.Nut_1.b_1k) + (geo.Nut_1.h_1ns/geo.Nut_1.b_1ns) + emag.lambda_1z; 

    case {'2SGL', '2SBL'}
        % Richtwerte fuer Hilfsfaktoren k_1 und k_2 [-]
        % Quelle: [Mueller08, S.331 - Bild 3.7.7]
        % V/A: m=3, einfache Zonenbreite
        k1_data =  [0.82,0.89,0.95,1; 2/3,0.8,0.9,1];
        k1_fun = polyfit(k1_data(2,:),k1_data(1,:),1);
        emag.misc.k_1 = polyval(k1_fun,wick.Sehnung_1);

        k2_data =  [0.75,0.85,0.925,1; 2/3,0.8,0.9,1];
        k2_fun = polyfit(k2_data(2,:),k2_data(1,:),1);
        emag.misc.k_2 = polyval(k2_fun,wick.Sehnung_1);

        % Abstand der beiden Wicklungen bei Zweischichtwicklungen d [mm]
        % Quelle: [Mueller08, S.324 - Bild 3.7.4b]
        % V/A: Abstand entspricht unegfaehr der Isolierdicke
        geo.Nut_1.d = geo.Nut_1.d_1iso;

        % Berechnung des resultierenden relativen Nut-Zahnkopf-Streuleitwert (Stator) lambda_1nz [-]
        % Quelle: [Mueller08, S.330 - Formel 3.7.16]
        emag.lambda_1nz = emag.misc.k_1 * (geo.Nut_1.h_1l/(3 * geo.Nut_1.b_1n)) + emag.misc.k_2 * ((geo.Nut_1.h_1ue/geo.Nut_1.b_1n) + (geo.Nut_1.h_1k/geo.Nut_1.b_1k) + (geo.Nut_1.h_1ns/geo.Nut_1.b_1ns) + emag.lambda_1z) + (geo.Nut_1.d/(4 * geo.Nut_1.b_1n)); 
        
        clear k1_data k1_fun k2_data k2_fun
    otherwise
        error('Ungueltige Eingabe bei Variable Wicklung');
end

% Berechnung der Nut- und Zahnkopfstreuung L_1sigma_nz (Stator) [H]
% Quelle: [Mueller08, S.533 - Formel 8.1.72], [Mueller08, S.533 - Formel 8.1.76]
emag.L_1sigma_nz = 2 * const.mu_0 * geo.l_i * (emag.w_1Str^2/rated.p) * (emag.lambda_1nz/wick.q_1);

% Berechnung des resultierenden relativen Nut-Zahnkopf-Streuleitwert (Rotor) lambda_2nz [-]
% Quelle: [Mueller08, S.326 - Formel 3.7.2a]
if(strcmp(opt.Maschinenausfuehrung,'Kaefiglaeufer'))
    emag.lambda_2nz = (geo.Nut_2.h_2l/(3 * geo.Nut_2.b_2n)) + (geo.Nut_2.h_2ue/geo.Nut_2.b_2n) + (geo.Nut_2.h_2k/geo.Nut_2.b_2k) + (geo.Nut_2.h_2ns/geo.Nut_2.b_2ns) + emag.lambda_2z;
elseif(strcmp(opt.Maschinenausfuehrung,'Schleifringlaeufer'))
    %to-do
    switch wick.Wicklungstyp_2
        case {'1SGL', '1SBL'}
            
        case {'2SGL', '2SBL'}
            
    end
else
    error('Ungueltige Eingabe bei Variable "opt.Maschinenausfuehrung"')
end

% Berechnung der Leiterlaenge im Wicklungskopf (Stator) l_1w [m]
% Quelle: [Geometrie], [Mueller - Berechnungsbeispiel]
wick.l_1w = 0.5 * (wick.l_1m - 2*geo.l);

switch wick.Wicklungstyp_1
    case {'1SGL', '1SBL'}
        % Richtwerte fuer den relativen Streuleitwert der Wickelkopfstreuung (Stator) lambda_1ws [-]
        % Quelle: [Mueller08, S.335 - Tabelle 3.7.2]
        % V/A: m=3
        emag.lambda_1ws = 0.3;
        
        % Berechnung des relativen Streuleitwertes [-]
        % Quelle: [Mueller08, S.534 - Formel 8.1.78]
        emag.lambda_1w = emag.lambda_1ws * (wick.l_1w / (geo.l_i));
        
    case {'2SGL', '2SBL'}
        % Richtwerte fuer den relativen Streuleitwert Wickelkopfstreuung (Stator) lambda_1ws [-]
        % Quelle: [Mueller08, S.335 - Tabelle 3.7.2]
        % V/A: m=3, Zylinderwicklung
        emag.lambda_1ws = 0.25;
        
        % Berechnung des relativen Streuleitwertes
        % Quelle: [Mueller08, S.534 - Formel 8.1.79]
        % V/A: einfache Zonenbreite
        emag.lambda_1w = emag.lambda_1ws * (wick.l_1w / (2 * geo.l_i));
    otherwise
        error('Ungueltige Eingabe bei Variable Wicklung');
end

% Berechnung der Wicklungskopfstreuung L_1sigma_w (Stator) [H]
% Quelle: [Mueller08, S.534 - Formel 8.1.77]
emag.L_1sigma_w = 2 * const.mu_0 * geo.l_i * (emag.w_1Str^2/rated.p) * emag.lambda_1w;

% Berechnung der Oberwellenstreuung L_1sigma_o (Stator) [H]
% Quelle: [Mueller08, S.534 - Formel 8.1.80a]
emag.L_1sigma_o = wick.sigma_1o * emag.L_1h_ges;

% Berechnung des Streukoeffizienten der Schraegungsstreuung sigma_schr [-]
% Quelle: [Mueller08, S.535 - Formel 8.1.82]
% to-do
emag.xi_schr_p = 1;
emag.sigma_schr = 1-emag.xi_schr_p^2;

% Berechnung der Schraegungsstreuung L_sigma_schr [H]
% Quelle: [Mueller08, S.534 - Formel 8.1.81a]
emag.L_sigma_schr = emag.sigma_schr * emag.L_1h_ges;

% Berechnung der Streuinduktivitaet (Stator) L_1sigma [H]
% Quelle: [Mueller08, S.535 - Formel 8.1.86a]
emag.L_1sigma = emag.L_1sigma_nz + emag.L_1sigma_w + emag.L_1sigma_o + emag.L_sigma_schr;

% Berechnung des reellen Uebersetzungsverhaeltnisses ue_h [-]
% Quelle: [Mueller08, S.547 - Formel 8.1.126]
if(strcmp(opt.Maschinenausfuehrung,'Kaefiglaeufer'))
    emag.ue_h = ((emag.w_1Str*wick.xi_1p)/(wick.N_2/6));
elseif(strcmp(opt.Maschinenausfuehrung,'Schleifringlaeufer'))
    %to-do
    emag.ue_h = ((emag.w_1Str*wick.xi_1p)/(emag.w_2Str*wick.xi_2p));
else
    error('Ungueltige Eingabe bei Variable "opt.Maschinenausfuehrung"')
end

% Berechnung der Streuinduktivitaet (Rotor) L_2sigma [H]
if(strcmp(opt.Maschinenausfuehrung,'Kaefiglaeufer'))
    % Berechnung der Streuinduktivitaet eines Stabs (Rotor) L_2sigma_s [H]
    % Quelle: [Mueller08, S.548 - Formel 8.1.130]
    emag.L_2sigma_s = const.mu_0 * geo.l_i * emag.lambda_2nz;
    
    % Berechnung des relativen Streuleitwerts des Rings (Rotor) lambda_2r [-]
    % Quelle: [Mueller08, S.548]
    emag.lambda_2r = 0.35; % zwischen 0.25 und 0.5
    
    % Berechnung des mittleren Ringdurchmessers (Rotor) D_2r [m]
    % Quelle: [Geometrie], [Mueller Berechnungsbeispiel], [Nuernberg52, S.29]
    geo.D_2r = geo.D_2a - geo.Nut_2.h_2n*1e-3;
    
    % Berechnung der Streuinduktivitaet eines Ringsegments (Rotor) L_2sigma_r [H]
    % Quelle: [Mueller08, S.548 - Formel 8.1.131]
    emag.L_2sigma_r = const.mu_0 * ((pi*geo.D_2r)/wick.N_2) * emag.lambda_2r;
    
    % Berechnung des Streukoeffizienten der Oberwellenstreuung sigma_2o (Rotor) [-]
    % Quelle: [Mueller08, S.338 - Formel 3.7.33d]
    emag.sigma_2o = ((pi*rated.p)/wick.N_2)^2/(sin((pi*rated.p)/wick.N_2)^2) - 1;
    
    % Berechnung der Streuinduktivitaet (Rotor) L_2sigma [H]
    % Quelle: [Mueller08, S.548 - Formel 8.1.129b], [Mueller08, S.549 - Formel 8.1.135]
    emag.L_2sigma = (wick.N_2/3) * (emag.L_2sigma_s + (1/(2*sin((pi*rated.p)/wick.N_2)^2))*emag.L_2sigma_r) + ((emag.sigma_2o*emag.L_1h_ges) / emag.ue_h^2);
elseif(strcmp(opt.Maschinenausfuehrung,'Schleifringlaeufer'))
    %to-do
    
    % Berechnung der Streuinduktivitaet (Rotor) L_2sigma [H]
    % Quelle: [Mueller08, S.548 - Formel 8.1.129a]
    emag.L_2sigma = L_2sigma_nz + L_2sigma_w + L_2sigma_o;
else
    error('Ungueltige Eingabe bei Variable "opt.Maschinenausfuehrung"')
end

% Berechnung der gesamten Induktivitaet (Stator) L_11 [H]
% Quelle: [Mueller08, S.547 - Formel 8.1.124] [Mueller14, S.418 - Formel 5.4.8]
emag.L_11 = emag.L_1h_ges + emag.L_1sigma;

% Berechnung der Koppelinduktivitaet [H]
% Quelle: [Mueller08, S.549 - Formel 8.1.133]
% L_12 = M
emag.L_12 = emag.L_1h_ges * (1/emag.ue_h) * emag.xi_schr_p;
emag.L_21 = emag.L_12;

% Berechnung der Hauptinduktivitaet (Rotor) L_2h [H]
% Quelle: [Mueller08, S.549 - Formel 8.1.135]
emag.L_2h = emag.L_1h_ges * (1/emag.ue_h^2);

% Berechnung der gesamten Induktivitaet (Rotor) L_22 [H]
% Quelle: [Mueller08, S.549 - Formel 8.1.134]
emag.L_22 = emag.L_2h + emag.L_2sigma;

% Berechnung des Streuleitfaktors (Stator) sigma_1 [-]
% Quelle: [Mueller08, S.547 - Formel 8.1.124], [Mueller14, S.421 - Formel 5.4.18]
emag.sigma_1 = emag.L_1sigma / emag.L_1h_ges;

% Berechnung der auf die Statorseite transformierte Streuinduktivitaet (Rotor) L_2sigma_trans [H]
% Quelle: [Mueller08, S.547 - Formel 8.1.127], [Mueller14, S.419 - Formel 5.4.12]
emag.L_2sigma_trans = emag.L_2sigma * emag.ue_h^2;

% Berechnung des Streuleitfaktors (Rotor) sigma_2 [-]
% Quelle: [Mueller08, S.549 - Formel 8.1.134], [Mueller14, S.427 - Formel 5.5.7]
emag.sigma_2 = emag.L_2sigma / emag.L_2h;
% emag.sigma_2_alt = emag.L_2sigma_trans / emag.L_1h_ges;

% Berechnung der Blondelschen Streuziffer
% Quelle: [Mueller14, S.420 - Formel 5.4.14], [Schroeder17, S.255 - Formel 5.158]
emag.sigma = 1 - (emag.L_12^2 / (emag.L_11 * emag.L_22));
% emag.sigma_alt = 1 - (1 / ((1+emag.sigma_1) * (1+emag.sigma_2)));

%% E.3) Widerstaende
% #########################################################################
% #   E.3) Widerstaende                                                   #
% #########################################################################

% Spezifischer Widerstand bei theta_1 (Stator) rho_1 [mm^2/S*m]
% Quelle: [Mueller08, S.435 - Formel 6.3.2]
emag.misc.rho_1 = opt.Stator_Leitermaterial.rho_20 * (1 + opt.Stator_Leitermaterial.alpha * (opt.theta_1 - 20));    

% Spezifische Leitfaeheigkeit (Stator) kappa_1 [S*m/mm^2]
% Quelle: [Mueller08, S.434 - Formel 6.3.1]
emag.misc.kappa_1 = 1/emag.misc.rho_1;

% Spezifischer Widerstand bei theta_2 (Rotor) rho_2 [mm^2/S*m]
% Quelle: [Mueller08, S.435 - Formel 6.3.2]
emag.misc.rho_2 = opt.Rotor_Leitermaterial.rho_20 * (1 + opt.Rotor_Leitermaterial.alpha * (opt.theta_2 - 20));

% Spezifische Leitfaeheigkeit (Rotor) kappa_2 [S*m/mm^2]
% Quelle: [Mueller08, S.434 - Formel 6.3.1]
emag.misc.kappa_2 = 1/emag.misc.rho_2;

% Berechnung des Widerstands eines Wicklungsstrangs (Stator) R_1 [Ohm]
% Quelle: [Mueller08, S.437 - Formel 6.3.14]
emag.R_1 = (emag.w_1Str * wick.l_1m) / (wick.a_1 * emag.misc.kappa_1 * wick.A_1L);

if(strcmp(opt.Maschinenausfuehrung,'Kaefiglaeufer'))
    % Berechnung des Stabwiderstands (Rotor) R_2s [Ohm]
    % Quelle: [Mueller08, S.435 - Formel 6.3.6]
    emag.R_2s = geo.l / (emag.misc.kappa_2 * wick.A_2s);
    
    % Berechnung der mittleren Ringsegmentlaenge (Rotor) l_2rm [m]
    % Quelle: [Geometrie], [Mueller Berechnungsbeispiel]
    geo.l_2rm = (pi*geo.D_2r) / wick.N_2;
    
    % Berechnung des Widerstands eines Ringsegments (Rotor) R_2r [Ohm]
    % Quelle: [Mueller08, S.435 - Formel 6.3.6]
    emag.R_2r = geo.l_2rm / (emag.misc.kappa_2 * wick.A_2r);
    
    % Berechnung des Widerstands (Rotor) R_2 [Ohm]
    % Quelle: [Mueller08, S.438 - Formel 6.3.16]
    emag.R_2 = (wick.N_2/3) * (emag.R_2s + (1/(2*sin((pi*rated.p)/wick.N_2)^2))*emag.R_2r);
elseif(strcmp(opt.Maschinenausfuehrung,'Schleifringlaeufer'))
    %to-do
else
    error('Ungueltige Eingabe bei Variable "opt.Maschinenausfuehrung"')
end

% Berechnung des auf die Statorseite transformierten Widerstands (Rotor) R_2_trans [Ohm]
emag.R_2_trans = emag.R_2 * emag.ue_h^2;


%% F) Postprocessing
% #########################################################################
% #   F) POSTPROCESSING                                                   #
% #########################################################################

warning('on','backtrace')

% x-y-Koordinaten fuer Plot und DXF-Export berechnen
geo = x_y_Koordinaten(geo, wick, opt);

% Berechnung erfolgreich
opt.Locked = 1;

% Sortieren der structs
const = orderfields(const);
emag = orderfields(emag);
emag.misc = orderfields(emag.misc);
geo = orderfields(geo);
geo.Nut_1 = orderfields(geo.Nut_1);
geo.Nut_2 = orderfields(geo.Nut_2);
geo.misc = orderfields(geo.misc);
wick = orderfields(wick);

opt = orderfields(opt);
rated = orderfields(rated);
richt = orderfields(richt);

% Speichern in Entwurf
Entwurf.Konstanten = const;
Entwurf.EMAG = emag;
Entwurf.Geometrie = geo;
Entwurf.Wicklung = wick;

Entwurf.Bemessungswerte = rated;
Entwurf.Richtwerte = richt;
Entwurf.Optionen = opt;

end

%% G) Hilfsfunktionen
% #########################################################################
% #   G) HILFSFUNKTIONEN                                                  #
% #########################################################################

% Berechnung groesster gemeinsamer Teiler (ggT)
function x = ggT_fun(a,b)
    while b ~= 0 
        z = mod(a,b); 
        a = b; 
        b = z; 
    end 
    x = a;
end

% Statornutgenerierung
function [emag, wick, geo] = Statornut(emag, wick, geo, richt)
    % Definition der Parameter zur eindeutigen Beschreibung der Nutform
    % Quelle: [Meyer09, S.53 - Tabelle 5.3]

    % Vorgabe einer minimalen Nutbreite (Stator) b_1n_o, b_1n_u, b_1n_m (oben, unten, mitte) [mm]
    % Quelle: [Meyer09, S.55 - Formeln 5.49,5.50,5.51]
    if((0.2 * wick.tau_1n*1e3)<geo.Nut_1.b_1ns)
        geo.Nut_1.b_1n_u = geo.Nut_1.b_1ns;                                    
    else
        geo.Nut_1.b_1n_u = 0.2 * wick.tau_1n*1e3;
    end
    geo.Nut_1.b_1n_o = geo.Nut_1.b_1n_u;
    geo.Nut_1.b_1n_m = (geo.Nut_1.b_1n_o + geo.Nut_1.b_1n_u) / 2;

    % Festlegung der Dicke der Nutisolierung (Stator) d_1iso [mm]
    % Quelle: [Mueller Berechnungsbeispiel], [Mueller06, S.224 - Abbildung 2.3.10]
    geo.Nut_1.d_1iso = 0.3;

    % Festlegung des Nutkeilwinkels (Stator) alpha_1nk [rad]
    % Quelle: Annahme
    geo.Nut_1.alpha_1nk = pi/6;

    % Berechnung der Hoehe des Nutkeilwinkels (Stator) h_1k [mm]
    geo.Nut_1.h_1k = tan(geo.Nut_1.alpha_1nk) * 0.5 * (geo.Nut_1.b_1n_u - geo.Nut_1.b_1ns);

    % Berechnung der Hoehe bis Wicklungsanfang (Stator) h_1nk [mm]
    % V/A: Nutkeil zu h_k+0.5mm angenommen
    geo.Nut_1.h_1nk = geo.Nut_1.h_1ns + geo.Nut_1.h_1k + 0.5;

    % Vorgabe einer minimalen Nuthoehe (Stator) h_1n [mm]
    % Quelle: Annahme
    geo.Nut_1.h_1n = geo.Nut_1.h_1nk + 1;

    % Berechnung der Zahnbreite oben (Stator) b_1z_o [mm]
    geo.Nut_1.b_1z_o = (geo.D_1i*1e3 + 2*geo.Nut_1.h_1n) * pi / wick.N_1 - geo.Nut_1.b_1n_o;

    % Berechnung der Zahnbreite unten (Stator) b_1z_u [mm]
    geo.Nut_1.b_1z_u = (geo.D_1i*1e3 + 2*geo.Nut_1.h_1ns + 2*geo.Nut_1.h_1k) * pi / wick.N_1 - geo.Nut_1.b_1n_u;

    % Berechnung der Zahnbreite mitte (Stator) b_1z_m [mm]
    geo.Nut_1.b_1z_m = (geo.Nut_1.b_1z_o + geo.Nut_1.b_1z_u) / 2;

    % Berechnung der Rueckenhoehe (Stator) h_1r [mm]
    geo.Nut_1.h_1r = (geo.D_1a_max - geo.D_1i)/2*1e3 - geo.Nut_1.h_1n;

    % Berechnung der tatsaechlichen Nutflaeche (Stator) A_1n_tat [mm^2]
    % Quelle: Flaechenberechnung Trapez
    geo.Nut_1.A_1n_tat = 0.5 * (geo.Nut_1.b_1n_u + geo.Nut_1.b_1n_o) * (geo.Nut_1.h_1n - geo.Nut_1.h_1nk);

    % Obergrenze Iterationen
    geo.Nut_1.maxIter = 1e6;
    geo.Nut_1.iter = 0;

    while(geo.Nut_1.A_1n_tat<geo.Nut_1.A_1n && geo.Nut_1.iter<geo.Nut_1.maxIter)
        % Berechnung der Zahninduktion oben (Stator) B_1z_o [T]
        % Quelle: [Mueller08, S.215 - Formel 2.4.3]
        emag.B_1z_o = (emag.B_max * wick.tau_1n * geo.l_i) / (geo.Nut_1.b_1z_o/1e3 * richt.phi_1Fe * geo.l_Fe);

        % Berechnung der Zahninduktion unten (Stator) B_1z_u [T]
        % Quelle: [Mueller08, S.215 - Formel 2.4.3]
        emag.B_1z_u = (emag.B_max * wick.tau_1n * geo.l_i) / (geo.Nut_1.b_1z_u/1e3 * richt.phi_1Fe * geo.l_Fe);

        % Berechnung der Zahninduktion mitte (Stator) B_1z_m [T]
        % Quelle: [Mueller08, S.215 - Formel 2.4.3]
        emag.B_1z_m = (emag.B_max * wick.tau_1n * geo.l_i) / (geo.Nut_1.b_1z_m/1e3 * richt.phi_1Fe * geo.l_Fe);

        % Berechnung der Rueckeninduktion (Stator) B_1r [T]
        % Quelle: [Mueller08, S.581 - Formel 9.1.30]
        emag.B_1r = emag.Phi_1r_max / (geo.Nut_1.h_1r*1e-3 * richt.phi_1Fe * geo.l_Fe);

        % Nutform-Algorithmus
        % Quelle: [Meyer09, S.56 - Abbildung 5.8]
        if(emag.B_1z_o < richt.B_1z_max || emag.B_1z_u < richt.B_1z_max)
            if(emag.B_1z_o < emag.B_1z_u)
                geo.Nut_1.b_1n_o = geo.Nut_1.b_1n_o + 0.1;
            else
                geo.Nut_1.b_1n_u = geo.Nut_1.b_1n_u + 0.1;
            end
        elseif(emag.B_1r < emag.B_1z_o && emag.B_1r < emag.B_1z_u)
            geo.Nut_1.h_1n = geo.Nut_1.h_1n + 0.1;
        elseif(emag.B_1z_o < emag.B_1z_u)
            geo.Nut_1.b_1n_o = geo.Nut_1.b_1n_o + 0.1;
        else
            geo.Nut_1.b_1n_u = geo.Nut_1.b_1n_u + 0.1;
        end

        % Update Statorgeometrie
        geo.Nut_1.b_1n_m = (geo.Nut_1.b_1n_o + geo.Nut_1.b_1n_u) / 2;
        geo.Nut_1.b_1z_o = (geo.D_1i*1e3 + 2*geo.Nut_1.h_1n) * pi / wick.N_1 - geo.Nut_1.b_1n_o;
        geo.Nut_1.b_1z_u = (geo.D_1i*1e3 + 2*geo.Nut_1.h_1ns + 2*geo.Nut_1.h_1k) * pi / wick.N_1 - geo.Nut_1.b_1n_u;
        geo.Nut_1.b_1z_m = (geo.Nut_1.b_1z_o + geo.Nut_1.b_1z_u) / 2;
        geo.Nut_1.h_1r = (geo.D_1a_max - geo.D_1i)/2*1e3 - geo.Nut_1.h_1n/1e3;
        geo.Nut_1.h_1k = tan(geo.Nut_1.alpha_1nk) * 0.5 * (geo.Nut_1.b_1n_u - geo.Nut_1.b_1ns);
        geo.Nut_1.h_1nk = geo.Nut_1.h_1ns + geo.Nut_1.h_1k + 0.5;

        % Ueberpruefen ob entartet
        if geo.Nut_1.b_1z_o<=0 || geo.Nut_1.b_1z_u<=0 || geo.Nut_1.b_1z_m<=0 || geo.Nut_1.h_1r<=0
            error('Geometrie ist entartet');
        end

        % Berechnung der tatsaechlichen Nutflaeche (Stator) A_1n_tat [mm^2]
        % Quelle: Flaechenberechnung Trapez
        geo.Nut_1.A_1n_tat = 0.5 * (geo.Nut_1.b_1n_u + geo.Nut_1.b_1n_o) * (geo.Nut_1.h_1n - geo.Nut_1.h_1nk);

        geo.Nut_1.iter = geo.Nut_1.iter + 1;
    end

    % Berechnung der Zahninduktion oben (Stator) B_1z_o [T]
    % Quelle: [Mueller08, S.215 - Formel 2.4.3]
    emag.B_1z_o = (emag.B_max * wick.tau_1n * geo.l_i) / (geo.Nut_1.b_1z_o/1e3 * richt.phi_1Fe * geo.l_Fe);

    % Berechnung der Zahninduktion unten (Stator) B_1z_u [T]
    % Quelle: [Mueller08, S.215 - Formel 2.4.3]
    emag.B_1z_u = (emag.B_max * wick.tau_1n * geo.l_i) / (geo.Nut_1.b_1z_u/1e3 * richt.phi_1Fe * geo.l_Fe);

    % Berechnung der Zahninduktion mitte (Stator) B_1z_m [T]
    % Quelle: [Mueller08, S.215 - Formel 2.4.3]
    emag.B_1z_m = (emag.B_max * wick.tau_1n * geo.l_i) / (geo.Nut_1.b_1z_m/1e3 * richt.phi_1Fe * geo.l_Fe);

    % Kontrolle Iterationen
    if(geo.Nut_1.iter>geo.Nut_1.maxIter)
        error('Kontrolle Nutraumbilanz notwendig!')
    end

    % Anpassen der Rueckenhoehe h_1r an die maximale Rueckeninduktion (Stator) [mm]
    % Quelle: [Mueller08, S.581 - Formel 9.1.30]
    if(emag.B_1r<richt.B_1r_max)
        geo.Nut_1.h_1r = emag.Phi_1r_max / (geo.l_Fe * richt.phi_1Fe * richt.B_1r_max) * 1e3;
    end

    % Berechnung der Rueckeninduktion (Stator) B_1r [T]
    % Quelle: [Mueller08, S.581 - Formel 9.1.30]
    emag.B_1r = emag.Phi_1r_max / (geo.Nut_1.h_1r*1e-3 * richt.phi_1Fe * geo.l_Fe);

    % Berechnung der relativen Rueckenhoehe (Stator) h_1r_rel [-]
    % Quelle: [Mueller08, S.586]
    geo.Nut_1.h_1r_rel = geo.Nut_1.h_1r*1e-3 / geo.tau_1p;

    % Berechnung des Aussendurchmessers (Stator) D_1a [m]
    % Quelle: [Mueller, S.599 - Formel 9.2.15a]
    geo.D_1a = geo.D_1i + 2 * (geo.Nut_1.h_1r*1e-3 + geo.Nut_1.h_1n*1e-3);
end

% Rotornutgenerierung
function [emag, wick, geo] = Rotornut(emag, wick, geo, richt)
    % Definition der Parameter zur eindeutigen Beschreibung der Nutform
    % Quelle: [Meyer09, S.53 - Tabelle 5.3]

    % Vorgabe einer minimalen Nutbreite (Rotor) b_2n_o, b_2n_u, b_2n_m (oben, unten, mitte) [mm]
    % Quelle: [Meyer09, S.55 - Formeln 5.49,5.50,5.51]
    if((0.2 * wick.tau_2n*1e3)<geo.Nut_2.b_2ns)
        geo.Nut_2.b_2n_u = geo.Nut_2.b_2ns;                                    
    else
        geo.Nut_2.b_2n_u = 0.2 * wick.tau_2n*1e3;
    end
    geo.Nut_2.b_2n_o = geo.Nut_2.b_2n_u;
    geo.Nut_2.b_2n_m = (geo.Nut_2.b_2n_o + geo.Nut_2.b_2n_u) / 2;

    % Festlegung der Dicke der Nutisolierung (Rotor) d_2iso [mm]
    % Quelle: [Mueller Berechnungsbeispiel], [Mueller06, S.224 - Abbildung 2.3.10]
    geo.Nut_2.d_2iso = 0.3;

    % Festlegung des Nutkeilwinkels (Rotor) alpha_2nk [rad]
    % Quelle: Annahme
    geo.Nut_2.alpha_2nk = pi/6;

    % Berechnung der Hoehe des Nutkeilwinkels (Rotor) h_2k [mm]
    geo.Nut_2.h_2k = tan(geo.Nut_2.alpha_2nk) * 0.5 * (geo.Nut_2.b_2n_u - geo.Nut_2.b_2ns);

    % Berechnung der Hoehe bis Wicklungsanfang (Rotor) h_2nk [mm]
    % V/A: Nutkeil zu h_k+0.5mm angenommen
    geo.Nut_2.h_2nk = geo.Nut_2.h_2ns + geo.Nut_2.h_2k + 0.5;

    % Vorgabe einer minimalen Nuthoehe (Rotor) h_2n [mm]
    % Quelle: Annahme
    geo.Nut_2.h_2n = geo.Nut_2.h_2nk + 1;

    % Berechnung der Zahnbreite oben (Rotor) b_2z_o [mm]
    geo.Nut_2.b_2z_o = (geo.D_2a*1e3 - 2*geo.Nut_2.h_2n) * pi / wick.N_2 - geo.Nut_2.b_2n_o;

    % Berechnung der Zahnbreite unten (Rotor) b_2z_u [mm]
    geo.Nut_2.b_2z_u = (geo.D_2a*1e3 - 2*geo.Nut_2.h_2ns - 2*geo.Nut_2.h_2k) * pi / wick.N_2 - geo.Nut_2.b_2n_u;

    % Berechnung der Zahnbreite mitte (Rotor) b_2z_m [mm]
    geo.Nut_2.b_2z_m = (geo.Nut_2.b_2z_o + geo.Nut_2.b_2z_u) / 2;

    % Berechnung der Rueckenhoehe (Rotor) h_2r [mm]
    geo.Nut_2.h_2r = geo.D_2a/2*1e3 - geo.Nut_2.h_2n;

    % Berechnung der tatsaechlichen Nutflaeche (Rotor) A_2n_tat [mm^2]
    % Quelle: Flaechenberechnung Trapez
    geo.Nut_2.A_2n_tat = 0.5 * (geo.Nut_2.b_2n_u + geo.Nut_2.b_2n_o) * (geo.Nut_2.h_2n - geo.Nut_2.h_2nk);

    % Obergrenze Iterationen
    geo.Nut_2.maxIter = 1e6;
    geo.Nut_2.iter = 0;

    while(geo.Nut_2.A_2n_tat<geo.Nut_2.A_2n && geo.Nut_2.iter<geo.Nut_2.maxIter)
        % Berechnung der Zahninduktion oben (Rotor) B_2z_o [T]
        % Quelle: [Mueller08, S.215 - Formel 2.4.3]
        emag.B_2z_o = (emag.B_max * wick.tau_2n * geo.l_i) / (geo.Nut_2.b_2z_o/1e3 * richt.phi_2Fe * geo.l_Fe);

        % Berechnung der Zahninduktion unten (Rotor) B_2z_u [T]
        % Quelle: [Mueller08, S.215 - Formel 2.4.3]
        emag.B_2z_u = (emag.B_max * wick.tau_2n * geo.l_i) / (geo.Nut_2.b_2z_u/1e3 * richt.phi_2Fe * geo.l_Fe);

        % Berechnung der Zahninduktion mitte (Rotor) B_2z_m [T]
        % Quelle: [Mueller08, S.215 - Formel 2.4.3]
        emag.B_2z_m = (emag.B_max * wick.tau_2n * geo.l_i) / (geo.Nut_2.b_2z_m/1e3 * richt.phi_2Fe * geo.l_Fe);

        % Berechnung der Rueckeninduktion (Rotor) B_2r [T]
        % Quelle: [Mueller08, S.581 - Formel 9.1.30]
        emag.B_2r = emag.Phi_2r_max / (geo.Nut_2.h_2r*1e-3 * richt.phi_2Fe * geo.l_Fe);

        % Nutform-Algorithmus
        % Quelle: [Meyer09, S.56 - Abbildung 5.8]
        if(emag.B_2z_u < richt.B_2z_max || emag.B_2z_o < richt.B_2z_max)
            if(emag.B_2z_u < emag.B_2z_o)
                geo.Nut_2.b_2n_u = geo.Nut_2.b_2n_u + 0.1;
            else
                geo.Nut_2.b_2n_o = geo.Nut_2.b_2n_o + 0.1;
            end
        elseif(emag.B_2r < emag.B_2z_u && emag.B_2r < emag.B_2z_o && geo.Nut_2.b_2n_o > 1.0)
            geo.Nut_2.h_2n = geo.Nut_2.h_2n + 0.1;
            while(emag.B_2z_o > richt.B_2z_max && geo.Nut_2.b_2n_o > 1.0)
                geo.Nut_2.b_2n_o = geo.Nut_2.b_2n_o - 0.1;
                geo.Nut_2.b_2z_o = (geo.D_2a*1e3 - 2*geo.Nut_2.h_2n) * pi / wick.N_2 - geo.Nut_2.b_2n_o;
                emag.B_2z_o = (emag.B_max * wick.tau_2n * geo.l_i) / (geo.Nut_2.b_2z_o/1e3 * richt.phi_2Fe * geo.l_Fe);
            end
        elseif(emag.B_2z_u < emag.B_2z_o)
            geo.Nut_2.b_2n_u = geo.Nut_2.b_2n_u + 0.1;
        else
            geo.Nut_2.b_2n_o = geo.Nut_2.b_2n_o + 0.1;
        end

        % Update Statorgeometrie
        geo.Nut_2.b_2n_m = (geo.Nut_2.b_2n_o + geo.Nut_2.b_2n_u) / 2;
        geo.Nut_2.b_2z_o = (geo.D_2a*1e3 - 2*geo.Nut_2.h_2n) * pi / wick.N_2 - geo.Nut_2.b_2n_o;
        geo.Nut_2.b_2z_u = (geo.D_2a*1e3 - 2*geo.Nut_2.h_2ns - 2*geo.Nut_2.h_2k) * pi / wick.N_2 - geo.Nut_2.b_2n_u;
        geo.Nut_2.b_2z_m = (geo.Nut_2.b_2z_o + geo.Nut_2.b_2z_u) / 2;
        geo.Nut_2.h_2r = geo.D_2a/2*1e3 - geo.Nut_2.h_2n;
        geo.Nut_2.h_2k = tan(geo.Nut_2.alpha_2nk) * 0.5 * (geo.Nut_2.b_2n_u - geo.Nut_2.b_2ns);
        geo.Nut_2.h_2nk = geo.Nut_2.h_2ns + geo.Nut_2.h_2k + 0.5;

        % Ueberpruefen ob entartet
        if geo.Nut_2.b_2z_o<=0 || geo.Nut_2.b_2z_u<=0 || geo.Nut_2.b_2z_m<=0 || geo.Nut_2.h_2r<=0
            error('Geometrie ist entartet');
        end

        % Berechnung der tatsaechlichen Nutflaeche (Rotor) A_2n_tat [mm^2]
        % Quelle: Flaechenberechnung Trapez
        geo.Nut_2.A_2n_tat = 0.5 * (geo.Nut_2.b_2n_u + geo.Nut_2.b_2n_o) * (geo.Nut_2.h_2n - geo.Nut_2.h_2nk);

        geo.Nut_2.iter = geo.Nut_2.iter + 1;
    end

    % Berechnung der Zahninduktion oben (Rotor) B_2z_o [T]
    % Quelle: [Mueller08, S.215 - Formel 2.4.3]
    emag.B_2z_o = (emag.B_max * wick.tau_2n * geo.l_i) / (geo.Nut_2.b_2z_o/1e3 * richt.phi_2Fe * geo.l_Fe);

    % Berechnung der Zahninduktion unten (Rotor) B_2z_u [T]
    % Quelle: [Mueller08, S.215 - Formel 2.4.3]
    emag.B_2z_u = (emag.B_max * wick.tau_2n * geo.l_i) / (geo.Nut_2.b_2z_u/1e3 * richt.phi_2Fe * geo.l_Fe);

    % Berechnung der Zahninduktion mitte (Rotor) B_2z_m [T]
    % Quelle: [Mueller08, S.215 - Formel 2.4.3]
    emag.B_2z_m = (emag.B_max * wick.tau_2n * geo.l_i) / (geo.Nut_2.b_2z_m/1e3 * richt.phi_2Fe * geo.l_Fe);

    % Kontrolle Iterationen
    if(geo.Nut_2.iter>geo.Nut_2.maxIter)
        error('Kontrolle Nutraumbilanz notwendig!')
    end

    % Anpassen der Rueckenhoehe h_2r an die maximale Rueckeninduktion (Rotor) [mm]
    % Quelle: [Mueller08, S.581 - Formel 9.1.30]
    if(emag.B_2r<richt.B_2r_max)
        geo.Nut_2.h_2r = emag.Phi_2r_max / (geo.l_Fe * richt.phi_2Fe * richt.B_2r_max) * 1e3;
    end

    % Berechnung der Rueckeninduktion (Rotor) B_2r [T]
    % Quelle: [Mueller08, S.581 - Formel 9.1.30]
    emag.B_2r = emag.Phi_2r_max / (geo.Nut_2.h_2r*1e-3 * richt.phi_2Fe * geo.l_Fe);

    % Berechnung der relativen Rueckenhoehe (Rotor) h_2r_rel [-]
    % Quelle: [Mueller08, S.586]
    geo.Nut_2.h_2r_rel = geo.Nut_2.h_2r*1e-3 / geo.tau_2p;

    % Berechnung des Innendurchmessers (Rotor) D_2i [m]
    geo.D_2i = geo.D_2a - 2.*geo.Nut_2.h_2n*1e-3 - 2.*geo.Nut_2.h_2r*1e-3;
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
    
end

% Berechnung x-y-Koordinaten fuer Plot und DXF-Export
function geo = x_y_Koordinaten(geo, wick ,opt)
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
    
    % Rotor aussen
    % Nut
    if(strcmp(opt.Nutform_Stator,'Trapezform (eckig)'))
        Segmenthoehe = ((geo.D_2a*1e3)/2) * (1-cos((2*asin(geo.Nut_2.b_2ns/(geo.D_2a*1e3)))/2));
        x1_Nut = geo.Nut_2.b_2ns / 2;
        x2_Nut = geo.Nut_2.b_2n_u / 2;
        x3_Nut = geo.Nut_2.b_2n_o / 2;
        y1_Nut = (geo.D_2a*1e3)/2 - Segmenthoehe;
        y2_Nut = -geo.Nut_2.h_2ns + (geo.D_2a*1e3)/2 - Segmenthoehe;
        y3_Nut = -geo.Nut_2.h_2ns - tan(geo.Nut_2.alpha_2nk) * ((geo.Nut_2.b_2n_u-geo.Nut_2.b_2ns) / 2) + (geo.D_2a*1e3)/2 - Segmenthoehe;
        y4_Nut = -geo.Nut_2.h_2n + (geo.D_2a*1e3)/2 - Segmenthoehe;
        Nut = [x1_Nut,x1_Nut,x2_Nut,x3_Nut,-x3_Nut,-x2_Nut,-x1_Nut,-x1_Nut; y1_Nut,y2_Nut,y3_Nut,y4_Nut,y4_Nut,y3_Nut,y2_Nut,y1_Nut]';

        % Nutfuellung
        x1_Fuellung = x2_Nut - geo.Nut_2.d_2iso;
        x2_Fuellung = x3_Nut - geo.Nut_2.d_2iso;
        y1_Fuellung = -geo.Nut_2.h_2nk + (geo.D_2a*1e3)/2 - Segmenthoehe;
        y2_Fuellung = -geo.Nut_2.h_2n + geo.Nut_2.d_2iso + (geo.D_2a*1e3)/2 - Segmenthoehe;
        Fuellung = [x1_Fuellung,x2_Fuellung,-x2_Fuellung,-x1_Fuellung; y1_Fuellung,y2_Fuellung,y2_Fuellung,y1_Fuellung]';
    else
        error('Ungueltige Eingabe bei Variable "opt.Nutform_Stator"');
    end
    % Kreissegment
    theta = linspace((asin(geo.Nut_2.b_2ns/(geo.D_2a*1e3))+(pi/2)),(((2*pi)/wick.N_2)-(asin(geo.Nut_2.b_2ns/(geo.D_2a*1e3)))+(pi/2)),1e1);
    Kreissegment = [(geo.D_2a*1e3)/2 * cos(theta); (geo.D_2a*1e3)/2 * sin(theta)]';
    
    % Koordinatentransformation
    geo.Rotor_aussen_x = [];
    geo.Rotor_aussen_y = [];
    geo.Fuellung_Rotor_x = [];
    geo.Fuellung_Rotor_y = [];
    for i = 1:wick.N_2
        theta = -((2*(i-1))/wick.N_2)*pi;
        M = [cos(theta) -sin(theta); sin(theta) cos(theta);];
        
        % Nut
        Nut_trans = Nut * M;
        
        % Nutfuellung
        Fuellung_trans = Fuellung * M;
        
        % Kreissegment
        Kreissegment_trans = Kreissegment * M;

        % gesamt
        geo.Rotor_aussen_x = [geo.Rotor_aussen_x; Nut_trans(2:end-1,1); Kreissegment_trans(:,1)];
        geo.Rotor_aussen_y = [geo.Rotor_aussen_y; Nut_trans(2:end-1,2); Kreissegment_trans(:,2)];
        geo.Fuellung_Rotor_x = [geo.Fuellung_Rotor_x; Fuellung_trans(:,1)];
        geo.Fuellung_Rotor_y = [geo.Fuellung_Rotor_y; Fuellung_trans(:,2)];
        geo.Fuellung_Rotor_x(end+1,1) = Fuellung_trans(1,1);
        geo.Fuellung_Rotor_y(end+1,1) = Fuellung_trans(1,2);
    end
    geo.Rotor_aussen_x(end+1,1) = geo.Rotor_aussen_x(1,1);
    geo.Rotor_aussen_y(end+1,1) = geo.Rotor_aussen_y(1,1);

    % Rotor innen
    theta = linspace(0,2*pi,4e2);
    geo.Rotor_innen_x = (geo.D_2i*1e3)/2 * cos(theta) + 0;
    geo.Rotor_innen_y = (geo.D_2i*1e3)/2 * sin(theta) + 0;

    geo.Rotor_innen_x = geo.Rotor_innen_x';
    geo.Rotor_innen_y = geo.Rotor_innen_y';
    
end

% Interpolation des Rueckenreduktionsfaktors
function [C_r_data, C_r] = C_r_data_interp(C_r_data, h_r, tau_r, B_r, Plot_C_r_Curve)
    % Vektoren auf gleiche Laenge bringen
    maxLength = max([length(C_r_data.B_08), length(C_r_data.B_10), length(C_r_data.B_12), length(C_r_data.B_14), length(C_r_data.B_16)]);
    C_r_data.B_08(:,length(C_r_data.B_08)+1:maxLength) = NaN;
    C_r_data.B_08(3,:) = ones(1,maxLength).*0.8;
    C_r_data.B_10(:,length(C_r_data.B_10)+1:maxLength) = NaN;
    C_r_data.B_10(3,:) = ones(1,maxLength).*1.0;
    C_r_data.B_12(:,length(C_r_data.B_12)+1:maxLength) = NaN;
    C_r_data.B_12(3,:) = ones(1,maxLength).*1.2;
    C_r_data.B_14(:,length(C_r_data.B_14)+1:maxLength) = NaN;
    C_r_data.B_14(3,:) = ones(1,maxLength).*1.4;
    C_r_data.B_16(:,length(C_r_data.B_16)+1:maxLength) = NaN;
    C_r_data.B_16(3,:) = ones(1,maxLength).*1.6;

    % Differenz zwischen Linien unterschiedlicher B_r
    diff_C_r(1,:) = C_r_data.B_08(2,:) - C_r_data.B_10(2,:);
    diff_C_r(2,:) = C_r_data.B_10(2,:) - C_r_data.B_12(2,:);
    diff_C_r(3,:) = C_r_data.B_12(2,:) - C_r_data.B_14(2,:);
    diff_C_r(4,:) = C_r_data.B_14(2,:) - C_r_data.B_16(2,:);

    % Mittlere Differenz zwischen Linien unterschiedlicher B_r
    diff_m_C_r(1,1) = (sum(diff_C_r(1,~isnan(diff_C_r(1,:))))) / length(diff_C_r(1,~isnan(diff_C_r(1,:))));
    diff_m_C_r(2,1) = (sum(diff_C_r(2,~isnan(diff_C_r(2,:))))) / length(diff_C_r(2,~isnan(diff_C_r(2,:))));
    diff_m_C_r(3,1) = (sum(diff_C_r(3,~isnan(diff_C_r(3,:))))) / length(diff_C_r(3,~isnan(diff_C_r(3,:))));
    diff_m_C_r(4,1) = (sum(diff_C_r(4,~isnan(diff_C_r(4,:))))) / length(diff_C_r(4,~isnan(diff_C_r(4,:))));
    
    % Berechnung der Daten fuer B_r = 0.6T
    % V/A: Daten fuer B_r = 0.6T errechnet durch Berechnen der Differenz
    % zwischen den zwei Linien fuer B_r = 0.8T und B_r = 1.0T
    C_r_data.B_06 = [C_r_data.B_08(1,:);
                      C_r_data.B_08(2,:) + diff_C_r(1,:)];
    C_r_data.B_06(3,:) = ones(1,maxLength).*0.6;

    % Berechnung der Daten fuer B_r = 0.4T
    % V/A: Daten fuer B_r = 0.4T errechnet durch Berechnen der Differenz
    % zwischen den zwei Linien fuer B_r = 0.8T und B_r = 1.0T
    C_r_data.B_04 = [C_r_data.B_06(1,:);
                      C_r_data.B_06(2,:) + diff_C_r(1,:)];
    C_r_data.B_04(3,:) = ones(1,maxLength).*0.4;

    % Berechnung der Daten fuer B_r = 1.8T
    % V/A: Daten fuer B_r = 1.8T errechnet durch Berechnen der Differenz
    % zwischen den zwei Linien fuer B_r = 1.4T und B_r = 1.6T
    C_r_data.B_18 = [C_r_data.B_16(1,:);
                      C_r_data.B_16(2,:) - diff_C_r(4,:)];
    C_r_data.B_18(3,:) = ones(1,maxLength).*1.8;
    
    C_r_data = orderfields(C_r_data);
    
    % Plot
    if(strcmp(Plot_C_r_Curve, ''))
    else
        fig = gcf;
        if(strcmp(fig.Name,'C_r Curve'))
            ax = gca;
            ax(end+1) = copyobj(ax(end),gcf);
            delete(get(ax(end),'Children'))
            set(ax(end), 'Color','none','XTick',[], 'YTick',[], 'Title', [], 'XLabel', [], 'YLabel', []);
        else
            close(gcf)
            figure('Name','C_r Curve')
            title('C_r Curve')
            hold on
            grid on
            ax = gca;
            ax.YMinorGrid = 'on';
            xlim([0 0.5])
            ylim([0 1.4])
            xlabel('h_r / \tau_r')
            ylabel('C_r')
        end
        if(strcmp(Plot_C_r_Curve, 'Statorruecken, p=1'))
            LineStyle = '-';
        elseif(strcmp(Plot_C_r_Curve, 'Rotorruecken, p=1'))
            LineStyle = '--';
        elseif(strcmp(Plot_C_r_Curve, 'Statorruecken, p=Inf'))
            LineStyle = ':';
        elseif(strcmp(Plot_C_r_Curve, 'Rotorruecken, p=Inf'))
            LineStyle = '-.';
        else
            error('Ungueltige Eingabe bei Variable "Plot_C_r_Curve"');
        end
        
        plot(C_r_data.B_04(1,:),C_r_data.B_04(2,:),'Color',[0 0.4470 0.7410],'LineStyle',LineStyle,'Parent',ax(end))
        plot(C_r_data.B_06(1,:),C_r_data.B_06(2,:),'Color',[0.8500 0.3250 0.0980],'LineStyle',LineStyle,'Parent',ax(end))
        plot(C_r_data.B_08(1,:),C_r_data.B_08(2,:),'Color',[0.9290 0.6940 0.1250],'LineStyle',LineStyle,'Parent',ax(end))
        plot(C_r_data.B_10(1,:),C_r_data.B_10(2,:),'Color',[0.4940 0.1840 0.5560],'LineStyle',LineStyle,'Parent',ax(end))
        plot(C_r_data.B_12(1,:),C_r_data.B_12(2,:),'Color',[0.4660 0.6740 0.1880],'LineStyle',LineStyle,'Parent',ax(end))
        plot(C_r_data.B_14(1,:),C_r_data.B_14(2,:),'Color',[0.3010 0.7450 0.9330],'LineStyle',LineStyle,'Parent',ax(end))
        plot(C_r_data.B_16(1,:),C_r_data.B_16(2,:),'Color',[0.6350 0.0780 0.1840],'LineStyle',LineStyle,'Parent',ax(end))
        plot(C_r_data.B_18(1,:),C_r_data.B_18(2,:),'Color',[0.5 0.5 0.5],'LineStyle',LineStyle,'Parent',ax(end))
        
        C_r_data_names = fieldnames(C_r_data);
        if(strcmp(Plot_C_r_Curve, 'Statorruecken, p=1'))
            lgd = legend(flipud(ax(end).Children), C_r_data_names, 'Location','northeast','Interpreter','none','Color','w');
            lgd.Title.String = Plot_C_r_Curve;
        elseif(strcmp(Plot_C_r_Curve, 'Rotorruecken, p=1'))
            lgd = legend(flipud(ax(end).Children), C_r_data_names, 'Location','east','Interpreter','none','Color','w');
            lgd.Title.String = Plot_C_r_Curve;
        elseif(strcmp(Plot_C_r_Curve, 'Statorruecken, p=Inf'))
            lgd = legend(flipud(ax(end).Children), C_r_data_names, 'Location','southeast','Interpreter','none','Color','w');
            lgd.Title.String = Plot_C_r_Curve;
        elseif(strcmp(Plot_C_r_Curve, 'Rotorruecken, p=Inf'))
            lgd = legend(flipud(ax(end).Children), C_r_data_names, 'Location','southeast','Interpreter','none','Color','w');
            lgd.Title.String = Plot_C_r_Curve;
        else
            error('Ungueltige Eingabe bei Variable "Plot_C_r_Curve"');
        end
    end

    % Interpolation der benoetigen Werte fuer B_r und h_r/tau_r
    C_r_data_names = fieldnames(C_r_data);
    for j = 1:length(B_r)
        B_vec = 0.4:0.2:1.8;
        B_vec = B_vec(1,abs(B_r(j) - B_vec)<=0.2);
        count = 1;
        for i = 1:length(C_r_data_names)
            if(any(abs(C_r_data.(C_r_data_names{i})(3,1)-B_vec)<1e-3))
                C_r_data_x(count,:) = C_r_data.(C_r_data_names{i})(1,:);
                C_r_data_y(count,:) = C_r_data.(C_r_data_names{i})(2,:);
                count = count + 1;
            end
        end
        try
            var_y = interp1(B_vec, C_r_data_y(:,(all(~isnan(C_r_data_y)))), B_r(j),'linear');
            var_x = C_r_data_x(1,(all(~isnan(C_r_data_x))));
            C_r(j) = interp1(var_x, var_y, ((h_r*1e-3)/tau_r),'linear');
        catch
            error('Datenbereich der hinterlegten C_r-Kurve (Rueckenreduktionsfaktor) zu klein.');
        end
    end
end
