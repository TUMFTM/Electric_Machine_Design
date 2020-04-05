% -------------------------------------------------------------------------
% TU Muenchen - Lehrstuhl fuer Fahrzeugtechnik (FTM)
% -------------------------------------------------------------------------
% Modell fuer den Entwurf und die Analyse einer PMSM oder ASM (MEAPA)
% -------------------------------------------------------------------------
% Autor: Svenja Kalt (kalt@ftm.mw.tum.de)
%        Jonathan Erhard
% -------------------------------------------------------------------------

% Hinweis: Die Richtwerte bestehen immer aus einem max-,min- und Standard-
%          Wert je Zeile. Bei mehreren Zeilen sind die unterschiedlichen 
%          Bedingungen anzugeben.
%          Beispiel: |4.0 1.0 0.8| <- für p==1
%                    |2.5 0.5 2.0| <- für p>1
%                      ^   ^   ^
%                      |   |   |
%                     max min Standard

% Richtwert fuer die relative Ankerlaenge lambda [-]
% Quelle: [Mueller08, S.577 - Tabelle 9.1.3], [Meyer18, S.117]
% richt.lambda = [p==1; p>1]
richt.lambda = [4.0 1.0 2.5; 2.5 0.5 2.0];
richt.lambda_opt = [1; Inf];

% Richtwert fuer die Kanalbreite der Ventilationskanaele l_v [m]
% Quelle: [Meyer09, S.41], [Mueller08, S.585]
richt.l_v = [0.01 0.006 0.01];

% Richtwert fuer die Amplitude der Luftspaltinduktion B_p [T]
% Quelle: [Mueller08, S.582 - Tabelle 9.1.5]
richt.B_p = [1.05 0.75 0.85];

% Richtwert fuer den Strombelag A [A/mm]
% Quelle: [Mueller08, S.580 - Tabelle 9.1.4
% richt.A = [Luft (indirekt); Wasser (direkt)]
richt.A = [120.0 30.0 50.0; 300.0 160.0 180.0];
richt.A_opt = [{'Luft (indirekt)'}; {'Wasser (direkt)'}];

% Richtwert fuer die Stromdichte (Stator) S_1 [A/mm^2]
% Quelle: [Mueller08, S.580 - Tabelle 9.1.4]
% richt.S_1 = [Luft (indirekt); Wasser (direkt)]
richt.S_1 = [7.0 3.0 7.0; 18.0 13.0 15.0];
richt.S_1_opt = [{'Luft (indirekt)'}; {'Wasser (direkt)'}];

% Richtwert fuer die max. zulaessige Induktion im Ruecken (Stator) B_1r_max [T]
% Quelle: [Mueller08, S.582 - Tabelle 9.1.5]
richt.B_1r_max = [1.5 1.0 1.4];

% Richtwert fuer die max. zulaessige Induktion in den Zaehnen (Stator) B_1z_max [T]
% Quelle: [Mueller08, S.582 - Tabelle 9.1.5]
richt.B_1z_max = [2.0 1.6 1.8];

% Richtwert fuer den Nutfuellfaktor (Stator) phi_1n [-]
% Quelle: [Mueller08, S.586 - Tabelle 9.1.6]
% V/A: Niederspannung
% richt.phi_1n = [Runddrahtwicklungen; Formspulen- oder Stabwicklungen]
richt.phi_1n = [0.5 0.3 0.5; 0.6 0.35 0.6];
richt.phi_1n_opt = [{'Runddraht'}; {'Formspule- oder Stab'}];

% Richtwert fuer den Wicklungsfaktor (Stator) xi_1p [-]
% Quelle: [Mueller08, S.603]
richt.xi_1p = [0.96 0.96 0.96];

% Richtwert fuer die minimale Nutteilung (Stator) tau_1n_min [m]
% Quelle: [Meyer09, S.46], [Pyr14]
richt.tau_1n_min = [0.07 0.007 0.007];

% Richtwert fuer den Eisenfuellfaktor (Stator) phi_1Fe [-]
% Quelle: [Mueller08, S.599]
richt.phi_1Fe = [1.0 0.9 0.95];

% Richtwert fuer die max. zulaessige Induktion im Ruecken (Rotor) B_2r_max [T]
% Quelle: [Mueller08, S.582 - Tabelle 9.1.5]
richt.B_2r_max = [1.5 1.0 1.4];

% Richtwert fuer den Eisenfuellfaktor (Rotor) phi_2Fe [-]
% Quelle: [Mueller08, S.599]
richt.phi_2Fe = [1.0 0.9 0.95];
