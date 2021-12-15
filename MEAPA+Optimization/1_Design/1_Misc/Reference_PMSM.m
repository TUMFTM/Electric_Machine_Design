% -------------------------------------------------------------------------
% TU Munich - Institute of Automotive Technology
% -------------------------------------------------------------------------
% Modell for the design and analysis of PMSM or ASM (MEAPA)
% -------------------------------------------------------------------------
% Autor:    Svenja Kalt (svenja.kalt@tum.de), 
%           Jonathan Erhard 
% -------------------------------------------------------------------------


% Note: The reference values always consist of a max, min and standard % value per line.
% value per line. If there are several lines, the different 
% conditions must be specified.
% Example: |4.0 1.0 0.8| <- for p==1
% |2.5 0.5 2.0| <- for p>1
% ^ ^ ^
% | | |
% max min Standard

% Standard value for relative anchor length lambda [-]
% Source: [Mueller08, p.577 - Table 9.1.3], [Meyer18, p.117]
% richt.lambda = [p==1; p>1]
richt.lambda = [4.0 1.0 2.5; 2.5 0.5 2.0];
richt.lambda_opt = [1; Inf];

% Guideline value for the duct width of ventilation ducts l_v [m].
% Source: [Meyer09, p.41], [Mueller08, p.585]
richt.l_v = [0.01 0.006 0.01];

% Guide value for the amplitude of the air gap induction B_p [T].
% Source: [Mueller08, p.582 - Table 9.1.5]
richt.B_p = [1.05 0.75 0.85];

% Reference value for the current coating A [A/mm].
% Source: [Mueller08, p.580 - Table 9.1.4
% richt.A = [air (indirect); water (direct)]
richt.A = [120.0 30.0 50.0; 300.0 160.0 180.0];
richt.A_opt = [{'Wasser (direkt)'}; {'Luft (indirekt)'}];

% Guide value for current density (stator) S_1 [A/mm^2]
% Source: [Mueller08, p.580 - Table 9.1.4]
% richt.S_1 = [air (indirect); water (direct)]
richt.S_1 = [7.0 3.0 7.0; 18.0 13.0 15.0];
richt.S_1_opt = [{'Wasser (direkt)'}; {'Luft (indirekt)'}];

% Reference value for the max. permissible induction in the back (stator) B_1r_max [T].
% Source: [Mueller08, p.582 - table 9.1.5]
richt.B_1r_max = [1.5 1.0 1.4];

% Reference value for the max. permissible induction in the teeth (stator) B_1z_max [T].
% Source: [Mueller08, p.582 - table 9.1.5]
richt.B_1z_max = [2.0 1.6 1.8];

% Guide value for the slot filling factor (stator) phi_1n [-]
% Source: [Mueller08, p.586 - Table 9.1.6]
% V/A: Low voltage
% richt.phi_1n = [round wire windings; shaped coil or bar windings]
richt.phi_1n = [0.5 0.3 0.5; 0.6 0.35 0.6];
richt.phi_1n_opt = [{'Runddraht'}; {'Formspule- oder Stab'}];

% Guide value for the winding factor (stator) xi_1p [-]
% Source: [Mueller08, p.603]
richt.xi_1p = [0.96 0.96 0.96];

% Richtwert fuer die minimale Nutteilung (Stator) tau_1n_min [m]
% Quelle: [Meyer09, S.46], [Pyr14]
richt.tau_1n_min = [0.07 0.007 0.007];

% Guide value for the iron fill factor (stator) phi_1Fe [-]
% Source: [Mueller08, p.599]
richt.phi_1Fe = [1.0 0.9 0.95];

% Guide value for the max. permissible induction in the back (rotor) B_2r_max [T].
% Source: [Mueller08, p.582 - table 9.1.5]
richt.B_2r_max = [1.5 1.0 1.4];

% Guide value for the iron fill factor (rotor) phi_2Fe [-]
% Source: [Mueller08, p.599]
richt.phi_2Fe = [1.0 0.9 0.95];