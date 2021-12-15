% -------------------------------------------------------------------------
% TU Munich - Institute of Automotive Technology
% -------------------------------------------------------------------------
% Modell for the design and analysis of PMSM or ASM (MEAPA)
% -------------------------------------------------------------------------
% Autor:    Svenja Kalt (svenja.kalt@tum.de), 
%           Jonathan Erhard 
% -------------------------------------------------------------------------

% Note: To create new files for electrical sheets a struct 
% with the entries
% - designation
% - rho_Fe
% - H
% - B
% - mu_r
% - p_vFe_B_vec
% - p_vFe_f_vec
% - p_vFe(columns: p_vFe_f_vec, rows: p_vFe_B_vec)
% must be created (see VACOFLUX 48). Afterwards the
% script must be executed.

clear data

%% Source: [xxx]

% Description
data.Bezeichnung = '<Elektroblech Name>';

% Density of the sheet material rho_Fe [kg/m^3]
data.rho_Fe = 0;

% electric steel sheet
% Field strength H [A/m]
data.H = 0;

% Magnetische Induktion [T]
data.B = 0;

% Relative permeability [-]
data.mu_r = 0;

% Frequency [Hz]
data.p_vFe_f_vec = 0;

% Magnetic induction [T]
data.p_vFe_B_vec = 0;

% Spec. iron losses [W/kg]
data.p_vFe = 0;


%% Data
path = ['5_Material/1_ElectricalSheet/', data.Bezeichnung];
save(path,'data');