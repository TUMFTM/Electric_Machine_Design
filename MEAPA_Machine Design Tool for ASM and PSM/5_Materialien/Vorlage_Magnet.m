% -------------------------------------------------------------------------
% TU Muenchen - Lehrstuhl fuer Fahrzeugtechnik (FTM)
% -------------------------------------------------------------------------
% Modell fuer den Entwurf und die Analyse einer PMSM oder ASM (MEAPA)
% -------------------------------------------------------------------------
% Autor: Svenja Kalt (kalt@ftm.mw.tum.de)
%        Jonathan Erhard
% -------------------------------------------------------------------------

% Hinweis: Um neue Dateien fuer Magnete zu erzeugen muss ein struct mit den
%          Eintraegen
%           - Bezeichnung
%           - rho_Mag
%           - B_r
%           - mu_r
%           - rho_el
%           - H_c
%           - TK_B_r
%          angelegt werden (siehe VACOMAX 225 TP). Anschliessend muss das
%          Skript nur noch ausgefuehrt werden.

clear data

%% VACOMAX 225 TP (Samarium-Kobalt-Legierung)
%{%
% Quelle: [https://www.vacuumschmelze.de/de/produkte/dauermagnete-systeme/dauermagnete/sm-co/vacomax/vacomax-225-tp.html]

% Bezeichnung
data.Bezeichnung = 'VACOMAX 225 TP';

% Dichte des Leiterwerkstoffs rho_Mag [kg/m^3]
data.rho_Mag = 7500;

% Magnetische Remanenzinduktion [T]
data.B_r = 1.07;

% Relative Permeabilitaet [-]
data.mu_r = 1.03;

% Spezifischer elektrischer Widerstand [Ohm*mm^2/m]
data.rho_el = 0.65; % zwischen 0.65 und 0.95

% Koerzitivfeldstaerke [A/m]
data.H_c = 879e3;

% Temperaturkoeffizient der magnetischen Remanenzinduktion [%/°C]
% V/A: RT = 100°C
data.TK_B_r = -0.030;

%}

%% VACODYM 238 TP (Nd-Fe-B Basis)
%{
% Quelle: [https://www.vacuumschmelze.de/de/produkte/dauermagnete-systeme/dauermagnete/nd-fe-b/vacodym/vacodym-238-tp.html]

% Bezeichnung
data.Bezeichnung = 'VACODYM 238 TP';

% Dichte des Leiterwerkstoffs rho_Mag [kg/m^3]
data.rho_Mag = 7500;

% Magnetische Remanenzinduktion [T]
data.B_r = 1.37;

% Relative Permeabilitaet [-]
data.mu_r = 1.03;

% Spezifischer elektrischer Widerstand [Ohm*mm^2/m]
% V/A: Widerstand parallel zu magnetischer Vorzugsrichtung
data.rho_el = 1.4; % zwischen 1.4 und 1.6

% Koerzitivfeldstaerke [A/m]
data.H_c = 1058e3;

% Temperaturkoeffizient der magnetischen Remanenzinduktion [%/°C]
% V/A: RT = 100°C
data.TK_B_r = -0.111;

%}

%% Datei erzeugen
path = ['5_Materialien/3_Magnet/', data.Bezeichnung];
save(path,'data');
