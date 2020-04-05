% -------------------------------------------------------------------------
% TU Muenchen - Lehrstuhl fuer Fahrzeugtechnik (FTM)
% -------------------------------------------------------------------------
% Modell fuer den Entwurf und die Analyse einer PMSM oder ASM (MEAPA)
% -------------------------------------------------------------------------
% Autor: Svenja Kalt (kalt@ftm.mw.tum.de)
%        Jonathan Erhard
% -------------------------------------------------------------------------

%% Literaturverzeichnis
% [Binder12] - Andraeas Binder, Elektrische Maschinen und Antriebe - Grundlagen, Betriebsverhalten, 978-3-540-71849-9
% [Meyer09] - Wolfgang Meyer, Automatisierter Entwurf elektromechanischer Wandler, 978-3-897-91406-3
% [Meyer18] - Wolfgang Meyer, Entwurf elektrischer Maschinen, Vorlesungsskriptum
% [Mueller08] - Germar Mueller, Berechnung elektrischer Maschinen, 3-527-40525-9
% [Mueller14] - Germar Mueller, Grundlagen elektrischer Maschinen, 3-527-41205-1
% [Pyr14] - Juha Pyrhoenen, Design of rotating electrical machines, 978-1-118-58157-5

%% Abkuerzungsverzeichnis
% V/A = Vereinfachung / Annahme
% Index 1: Stator
% Index 2: Rotor

%% Initialisierung Skript
clear all, close all force, home

% Pfade wiederherstellen
restoredefaultpath;
path(pathdef)
addpath(genpath(pwd));

% Auswahl mit oder ohne GUI
Enable_GUI = 1;

%% Start
if(Enable_GUI)
    GUI_Entwurf;
else
    % Nennleistung P_N [W]
    rated.P_N   = 3000;
    % Nenndrehzahl n_N [U/min]
    rated.n_N   = 1500;
    % Nennspannung U_N [V]
    rated.U_N   = 400;
    % Polpaarzahl p [-]
    rated.p     = 2;
    % Nennfrequenz f_N [-]
    rated.f_N   = (rated.p * rated.n_N) / 60;
    % Strangzahl m [-]
    rated.m     = 3;
    
    [Entwurf, Analyse] = MEAPA_Skript(rated);
end

clear RESTOREDEFAULTPATH_EXECUTED Enable_GUI rated
