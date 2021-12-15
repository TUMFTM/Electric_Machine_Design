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
% Example: |1.0 0.6 0.8| <- for p==1
% |4.0 1.0 2.0| <- for p>1
% ^ ^ ^
% | | |
% max min Standard

% Standard value for relative anchor length lambda [-]
% Source: [Mueller08, p.577 - Table 9.1.3], [Meyer18, p.117]
% richt.lambda = [p==1; p>1]
richt.lambda = [1.0 0.6 0.8; 4.0 1.0 2.5];
richt.lambda_opt = [1; Inf];

% Guideline value for the duct width of ventilation ducts l_v [m].
% Source: [Meyer09, p.41], [Mueller08, p.585]
richt.l_v = [0.01 0.006 0.01];

% Guide value for the mean value of the air gap induction B_m [T].
% Source: [Mueller08, p.582 - Table 9.1.5]
richt.B_m = [0.65 0.4 0.58];

% Guide value for the current coating A [A/mm].
% Source: [Mueller08, p.580 - Table 9.1.4]
richt.A = [120.0 20.0 50.0];

% Guide value for current density (stator) S_1 [A/mm^2]
% Source: [Mueller08, p.580 - Table 9.1.4]
richt.S_1 = [8.0 3.0 7.0];

% Reference value for the max. permissible induction in the back (stator) B_1r_max [T].
% Source: [Mueller08, p.582 - table 9.1.5]
richt.B_1r_max = [1.65 1.3 1.4];

% Reference value for the max. permissible induction in the teeth (stator) B_1z_max [T].
% Source: [Mueller08, p.582 - table 9.1.5]
richt.B_1z_max = [2.1 1.4 1.8];

% Guide value for the slot filling factor (stator) phi_1n [-]
% Source: [Mueller08, p.586 - Table 9.1.6]
% V/A: Low voltage
% richt.phi_1n = [round wire windings; shaped coil or bar windings]
richt.phi_1n = [0.5 0.3 0.5; 0.6 0.35 0.6];
richt.phi_1n_opt = [{'Runddraht'}; {'Formspule- oder Stab'}];

% Guide value for the winding factor (stator) xi_1p [-]
% Source: [Mueller08, p.596]
% V/A: no distinction between single and double layer winding,
% only needed for initial estimation (will be calculated in detail in further
% design process)
richt.xi_1p = [0.96 0.92 0.96];

% Guide value for the minimum slot pitch (stator) tau_1n_min [m].
% Source: [Meyer09, p.46], [Pyr14]
richt.tau_1n_min = [0.07 0.007 0.007];

% Guide value for the iron fill factor (stator) phi_1Fe [-]
% Source: [Mueller08, p.599]
richt.phi_1Fe = [1.0 0.9 0.95];

% Guide value for current density (rotor) S_2 [A/mm^2]
% Source: [Mueller08, p.580 - Table 9.1.4]
% richt.S_2 = [copper; cast aluminum; aluminum wire]
richt.S_2 = [8.0 3.0 5.0; 6.5 3.0 5.0; 6.5 3.0 5.0];
richt.S_2_opt = [{'Copper'}; {'Aluminumcasting'}; {'Aluminumwire'}];

% Guide value for bar current density (rotor) S_2s [A/mm^2].
% Source: [Mueller08, p.580 - Table 9.1.4]
% richt.S_2s = [copper; cast aluminum; aluminum wire]
richt.S_2s = [8.0 3.0 4.0; 6.5 3.0 3.8; 6.5 3.0 3.8;];
richt.S_2s_opt = [{'Copper'}; {'Aluminumcasting'}; {'Aluminumwire'}];

% Guide value for the ring current density (rotor) S_2r [A/mm^2]
% Source: [Mueller08, p.580 - Table 9.1.4]
% richt.S_2r = [copper; cast aluminum; aluminum wire]
richt.S_2r = [8.0 3.0 5.0; 6.5 3.0 4.0; 6.5 3.0 4.0];
richt.S_2r_opt = [{'Copper'}; {'Aluminumcasting'}; {'Aluminumwire'}];

% Guide value for the max. permissible induction in the back (rotor) B_2r_max [T].
% Source: [Mueller08, p.582 - table 9.1.5]
richt.B_2r_max = [1.6 0.4 1.5];

% Reference value for the max. permissible induction in the teeth (rotor) B_2z_max [T].
% Source: [Mueller08, p.582 - table 9.1.5]
richt.B_2z_max = [2.2 1.5 1.8];

% Guide value for the slot fill factor (rotor) phi_2n [-]
% Source: [Mueller08, p.586 - Table 9.1.6]
% V/A: Low voltage
% richt.phi_2n = [round wire windings; shaped coil or bar windings]
richt.phi_2n = [0.5 0.3 0.5; 0.6 0.35 0.6];
richt.phi_2n_opt = [{'Runddraht'}; {'Formspule- oder Stab'}];

% Guide value for the winding factor (rotor) xi_2p [-]
% Source: [Mueller08, p.596]
% V/A: no distinction between single and double layer winding,
% only needed for initial estimation (will be calculated in detail in further
% design process)
richt.xi_2p = [0.96 0.92 0.96];

% Guide value for the minimum slot pitch (rotor) tau_2n_min [m].
% Source: [Meyer09, p.46], [Pyr14]
richt.tau_2n_min = [0.07 0.007 0.007];

% Guide value for the iron fill factor (rotor) phi_2Fe [-]
% Source: [Mueller08, p.599]
richt.phi_2Fe = [1.0 0.9 0.95];