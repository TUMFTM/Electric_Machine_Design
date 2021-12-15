% -------------------------------------------------------------------------
% TU Munich - Institute of Automotive Technology
% -------------------------------------------------------------------------
% Modell for the design and analysis of PMSM or ASM (MEAPA)
% -------------------------------------------------------------------------
% Autor:    Svenja Kalt (svenja.kalt@tum.de), 
%           Jonathan Erhard 
% -------------------------------------------------------------------------

%% Bibliography
% [Binder12] - Andreas Binder, Electrical Machines and Drives - Fundamentals, Operational Behavior, 978-3-540-71849-9.
% [Meyer09] - Wolfgang Meyer, Automated design of electromechanical converters, 978-3-897-91406-3
% [Meyer18] - Wolfgang Meyer, Design of electrical machines, Lecture notes
% [Mueller08] - Germar Mueller, Calculation of electrical machines, 3-527-40525-9
% [Mueller14] - Germar Mueller, Fundamentals of Electrical Machines, 3-527-41205-1
% [Pyr14] - Juha Pyrhoenen, Design of rotating electrical machines, 978-1-118-58157-5

%% List of abbreviations
% V/A = simplification / assumption
% Index 1: Stator
% Index 2: Rotor

%% Initialisation Script
clear all, close all force, home

% Restore paths
restoredefaultpath;
path(pathdef)
addpath(genpath(pwd));

% Selection with or without GUI
Enable_GUI = 1;


%% Start
if(Enable_GUI)
   
    GUI_Design;
else
    %Maschine Design
        % Maschine topology - ASM oder PMSM
        rated.type = 'PMSM';
        % Rated power P_N [W]
        rated.P_N   = 75000;
        % Rated rotational speed n_N [U/min]
        rated.n_N   = 4800;
        % Rated voltage U_N [V]
        rated.U_N   = 360;
        % Max. rotational speed n_max [U/min]
        rated.nmax = 12000;
        % Number of pole pairs p [-]
        rated.p     = 6;
        % Rated frequency f_N [-]
        rated.f_N   = (rated.p * rated.n_N) / 60;
        % Number of strands m [-]
        rated.m     = 3;
        % load factor [-]
        rated.cos_phi_N = 0.95;
    
    % Set to default BMW i3 vehicle
        rated.LDS.fz_m.String = num2str(1640); 
        rated.LDS.cW.String = num2str(0.3);
        rated.LDS.A.String = num2str(2.38);
        rated.LDS.tyre_r.String = num2str(0.3498);
        rated.LDS.battcap.String = num2str(22);
        rated.LDS.aux.String = num2str(1500);
        rated.LDS.GearRatio.String = num2str(9.7);
          
        rated.LDS.visual_LDS = 0; %visual of LDS desired: 1, visual not desired: 0
        
    [Entwurf, Analyse] = MEAPA_Script(rated);
end

clear RESTOREDEFAULTPATH_EXECUTED Enable_GUI rated

