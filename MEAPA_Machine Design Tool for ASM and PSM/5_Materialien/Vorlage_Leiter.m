% -------------------------------------------------------------------------
% TU Muenchen - Lehrstuhl fuer Fahrzeugtechnik (FTM)
% -------------------------------------------------------------------------
% Modell fuer den Entwurf und die Analyse einer PMSM oder ASM (MEAPA)
% -------------------------------------------------------------------------
% Autor: Svenja Kalt (kalt@ftm.mw.tum.de)
%        Jonathan Erhard
% -------------------------------------------------------------------------

% Hinweis: Um neue Dateien fuer Leiter zu erzeugen muss ein struct mit den
%          Eintraegen
%           - Bezeichnung
%           - rho_Fe
%           - rho_20
%           - alpha
%          angelegt werden (siehe Kupfer). Anschliessend muss das Skript 
%          nur noch ausgefuehrt werden.

clear data

%% Kupfer
%{%

% Bezeichnung
data.Bezeichnung = 'Kupfer';

% Dichte des Leiterwerkstoffs rho_Le [kg/m^3]
data.rho_Le = 8940;

% Richtwert spezifischer Widerstand rho_20 bei 20°C fuer Kupferdraht [mm^2/S*m]
% Quelle: [Mueller08, S.435 - Tabelle 6.3.1]
data.rho_20 = 1/58;

% Temperaturbeiwert alpha fuer Kupferdraht [1/K]
% Quelle: [Mueller08, S.435]
data.alpha = 0.392 * 10e-3;

%}

%% Aluminiumdraht
%{

% Bezeichnung
data.Bezeichnung = 'Aluminiumdraht';

% Dichte des Leiterwerkstoffs rho_Fe [kg/m^3]
data.rho_Le = 2700;

% Richtwert spezifischer Widerstand rho_20 bei 20°C fuer Aluminiumdraht [mm^2/S*m]
% Quelle: [Mueller08, S.435 - Tabelle 6.3.1]
data.rho_20 = 1/37;

% Temperaturbeiwert alpha fuer Aluminiumdraht [1/K]
% Quelle: [Mueller08, S.435]
data.alpha = 0.4 * 10e-3;

%}

%% Aluminiumguss
%{

% Bezeichnung
data.Bezeichnung = 'Aluminiumguss';

% Dichte des Leiterwerkstoffs rho_Fe [kg/m^3]
data.rho_Le = 2700;

% Richtwert spezifischer Widerstand rho_20 bei 20°C fuer Aluminiumdraht [mm^2/S*m]
% Quelle: [Mueller08, S.435 - Tabelle 6.3.1]
data.rho_20 = 1/30;

% Temperaturbeiwert alpha fuer Aluminiumdraht [1/K]
% Quelle: [Mueller08, S.435]
data.alpha = 0.4 * 10e-3;

%}

%% Datei erzeugen
path = ['5_Materialien/2_Leiter/', data.Bezeichnung];
save(path,'data');
