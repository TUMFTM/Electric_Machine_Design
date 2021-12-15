% -------------------------------------------------------------------------
% TU Munich - Institute of Automotive Technology
% -------------------------------------------------------------------------
% Modell for the design and analysis of PMSM or ASM (MEAPA)
% -------------------------------------------------------------------------
% Autor:    Svenja Kalt (svenja.kalt@tum.de), 
%           Jonathan Erhard 
% -------------------------------------------------------------------------

% Note: In order to create new files for magnets, a struct with the
% entries
% - name
% - rho_Mag
% - B_r
% - mu_r
% - rho_el
% - H_c
% - TK_B_r
% must be created (see VACOMAX 225 TP). Afterwards the
% script must be executed.

clear data

%% VACOMAX 225 TP (Samarium-cobalt alloy)
%{%
% Source: [https://www.vacuumschmelze.de/de/produkte/dauermagnete-systeme/dauermagnete/sm-co/vacomax/vacomax-225-tp.html]

% Description
data.Bezeichnung = 'VACOMAX 225 TP';

% Density of the conductor material rho_Mag [kg/m^3].
data.rho_Mag = 7500;

% Magnetic remanence induction [T]
data.B_r = 1.07;

% Relative permeability [-]
data.mu_r = 1.03;

% Electrical resistivity [Ohm*mm^2/m]
data.rho_el = 0.65; % between 0.65 and 0.95

% Coercivity [A/m]
data.H_c = 879e3;

% Temperature coefficient of magnetic remanence induction [%/°C]
% V/A: RT = 100°C
data.TK_B_r = -0.030;


%% Data
path = ['5_Material/3_Magnet/', data.Bezeichnung];
save(path,'data');