% -------------------------------------------------------------------------
% TU Munich - Institute of Automotive Technology
% -------------------------------------------------------------------------
% Modell for the design and analysis of PMSM or ASM (MEAPA)
% -------------------------------------------------------------------------
% Autor:    Svenja Kalt (svenja.kalt@tum.de), 
%           Jonathan Erhard 
% -------------------------------------------------------------------------

% Note: To create new files for leaders you must create a struct with the
% entries
% - designation
% - rho_Fe
% - rho_20
% - alpha
% must be created (see Copper). Afterwards the script 
% only has to be executed.

clear data

%% Copper

% Description
data.Bezeichnung = 'Copper';

% Density of the conductor material rho_Le [kg/m^3].
data.rho_Le = 8940;

% Reference value of resistivity rho_20 at 20°C for copper wire [mm^2/S*m]
% Source: [Mueller08, S.435 - Tabelle 6.3.1]
data.rho_20 = 1/58;

% Temperature coefficient alpha for copper wire [1/K]
% Source: [Mueller08, S.435]
data.alpha = 0.392 * 10e-3;


%% Data
path = ['5_Material/2_Conductor/', data.Bezeichnung];
save(path,'data');