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
%          Beispiel: |1.0 0.6 0.8| <- für p==1
%                    |4.0 1.0 2.0| <- für p>1
%                      ^   ^   ^
%                      |   |   |
%                     max min Standard

% Richtwert fuer die relative Ankerlaenge lambda [-]
% Quelle: [Mueller08, S.577 - Tabelle 9.1.3], [Meyer18, S.117]
% richt.lambda = [p==1; p>1]
richt.lambda = [1.0 0.6 0.8; 4.0 1.0 2.5];
richt.lambda_opt = [1; Inf];

% Richtwert fuer die Kanalbreite der Ventilationskanaele l_v [m]
% Quelle: [Meyer09, S.41], [Mueller08, S.585]
richt.l_v = [0.01 0.006 0.01];

% Richtwert fuer den Mittelwert der Luftspaltinduktion B_m [T]
% Quelle: [Mueller08, S.582 - Tabelle 9.1.5]
richt.B_m = [0.65 0.4 0.58];

% Richtwert fuer den Strombelag A [A/mm]
% Quelle: [Mueller08, S.580 - Tabelle 9.1.4]
richt.A = [120.0 20.0 50.0];

% Richtwert fuer die Stromdichte (Stator) S_1 [A/mm^2]
% Quelle: [Mueller08, S.580 - Tabelle 9.1.4]
richt.S_1 = [8.0 3.0 7.0];

% Richtwert fuer die max. zulaessige Induktion im Ruecken (Stator) B_1r_max [T]
% Quelle: [Mueller08, S.582 - Tabelle 9.1.5]
richt.B_1r_max = [1.65 1.3 1.4];

% Richtwert fuer die max. zulaessige Induktion in den Zaehnen (Stator) B_1z_max [T]
% Quelle: [Mueller08, S.582 - Tabelle 9.1.5]
richt.B_1z_max = [2.1 1.4 1.8];

% Richtwert fuer den Nutfuellfaktor (Stator) phi_1n [-]
% Quelle: [Mueller08, S.586 - Tabelle 9.1.6]
% V/A: Niederspannung
% richt.phi_1n = [Runddrahtwicklungen; Formspulen- oder Stabwicklungen]
richt.phi_1n = [0.5 0.3 0.5; 0.6 0.35 0.6];
richt.phi_1n_opt = [{'Runddraht'}; {'Formspule- oder Stab'}];

% Richtwert fuer den Wicklungsfaktor (Stator) xi_1p [-]
% Quelle: [Mueller08, S.596]
% V/A: keine Unterscheidung zwischen Einschicht- und Zweischichtwicklung,
% da lediglich fuer initiale Abschaetzung benoetigt (wird im weiteren
% Entwurfsverlauf genau berechnet)
richt.xi_1p = [0.96 0.92 0.96];

% Richtwert fuer die minimale Nutteilung (Stator) tau_1n_min [m]
% Quelle: [Meyer09, S.46], [Pyr14]
richt.tau_1n_min = [0.07 0.007 0.007];

% Richtwert fuer den Eisenfuellfaktor (Stator) phi_1Fe [-]
% Quelle: [Mueller08, S.599]
richt.phi_1Fe = [1.0 0.9 0.95];

% Richtwert fuer die Stromdichte (Rotor) S_2 [A/mm^2]
% Quelle: [Mueller08, S.580 - Tabelle 9.1.4]
% richt.S_2 = [Kupfer; Aluminiumguss; Aluminiumdraht]
richt.S_2 = [8.0 3.0 5.0; 6.5 3.0 5.0; 6.5 3.0 5.0];
richt.S_2_opt = [{'Kupfer'}; {'Aluminiumguss'}; {'Aluminiumdraht'}];

% Richtwert fuer die Stabstromdichte (Rotor) S_2s [A/mm^2]
% Quelle: [Mueller08, S.580 - Tabelle 9.1.4]
% richt.S_2s = [Kupfer; Aluminiumguss; Aluminiumdraht]
richt.S_2s = [8.0 3.0 4.0; 6.5 3.0 3.8; 6.5 3.0 3.8;];
richt.S_2s_opt = [{'Kupfer'}; {'Aluminiumguss'}; {'Aluminiumdraht'}];

% Richtwert fuer die Ringstromdichte (Rotor) S_2r [A/mm^2]
% Quelle: [Mueller08, S.580 - Tabelle 9.1.4]
% richt.S_2r = [Kupfer; Aluminiumguss; Aluminiumdraht]
richt.S_2r = [8.0 3.0 5.0; 6.5 3.0 4.0; 6.5 3.0 4.0];
richt.S_2r_opt = [{'Kupfer'}; {'Aluminiumguss'}; {'Aluminiumdraht'}];

% Richtwert fuer die max. zulaessige Induktion im Ruecken (Rotor) B_2r_max [T]
% Quelle: [Mueller08, S.582 - Tabelle 9.1.5]
richt.B_2r_max = [1.6 0.4 1.5];

% Richtwert fuer die max. zulaessige Induktion in den Zaehnen (Rotor) B_2z_max [T]
% Quelle: [Mueller08, S.582 - Tabelle 9.1.5]
richt.B_2z_max = [2.2 1.5 1.9];

% Richtwert fuer den Nutfuellfaktor (Rotor) phi_2n [-]
% Quelle: [Mueller08, S.586 - Tabelle 9.1.6]
% V/A: Niederspannung
% richt.phi_2n = [Runddrahtwicklungen; Formspulen- oder Stabwicklungen]
richt.phi_2n = [0.5 0.3 0.5; 0.6 0.35 0.6];
richt.phi_2n_opt = [{'Runddraht'}; {'Formspule- oder Stab'}];

% Richtwert fuer den Wicklungsfaktor (Rotor) xi_2p [-]
% Quelle: [Mueller08, S.596]
% V/A: keine Unterscheidung zwischen Einschicht- und Zweischichtwicklung,
% da lediglich fuer initiale Abschaetzung benoetigt (wird im weiteren
% Entwurfsverlauf genau berechnet)
richt.xi_2p = [0.96 0.92 0.96];

% Richtwert fuer die minimale Nutteilung (Rotor) tau_2n_min [m]
% Quelle: [Meyer09, S.46], [Pyr14]
richt.tau_2n_min = [0.07 0.007 0.007];

% Richtwert fuer den Eisenfuellfaktor (Rotor) phi_2Fe [-]
% Quelle: [Mueller08, S.599]
richt.phi_2Fe = [1.0 0.9 0.95];
