% -------------------------------------------------------------------------
% TU Munich - Institute of Automotive Technology
% -------------------------------------------------------------------------
% Model for the design of a permanent magnet excited synchronous machine and
% subsequent efficiency map calculation
% -------------------------------------------------------------------------
% Autor:    Svenja Kalt (kalt@ftm.mw.tum.de)
%           Jonathan Erhard
%           Prof. Markus Lienkamp
% -------------------------------------------------------------------------

%% Literature
% [Binder12]    - Andreas Binder, Elektrische Maschinen und Antriebe - Grundlagen, Betriebsverhalten, 978-3-540-71849-9
% [Ionel98]     - D.M. Ionel, Design considerations for permanent magnet synchronous motors for flux weakening applications, http://dx.doi.org/10.1049/ip-epa:19982171
% [Meyer09]     - Wolfgang Meyer, Automatisierter Entwurf elektromechanischer Wandler, 978-3-897-91406-3
% [Meyer18]     - Wolfgang Meyer, Entwurf elektrischer Maschinen, Vorlesungsskriptum
% [Mueller08]   - Germar Mueller, Berechnung elektrischer Maschinen, 3-527-40525-9
% [Mueller14]   - Germar Mueller, Grundlagen elektrischer Maschinen, 3-527-41205-1
% [Pyr14]       - Juha Pyrhoenen, Design of rotating electrical machines, 978-1-118-58157-5

%% Model Structure
% 1 - Preprocessing
% 2 - Design process PMSM
% 3 - Efficiency map calculation PMSM

%% 1 - Preprocessing 
% The input data of the PMSM is entered via a GUI. Starting from 
% this data set, the entire design of the machine and the efficiency map calculation is carried out.
    
    % Initialization Script
    clear all, close all force, home
    path(pathdef)
    addpath(genpath(pwd));

    % Start GUI
    handles = GUI_Auslegung_PMSM;

    % Termination by user
    if(~isstruct(handles))
        error('Operation aborted by user!');
    else
        % Write struct in machine data
        Maschinendaten.Bemessungsgroessen.Primaerparameter = handles.Primaerparameter;
        Maschinendaten.Bemessungsgroessen.Sekundaerparameter = handles.Sekundaerparameter;
        Maschinendaten.Richtwerte = handles.Richtwerte;
        Maschinendaten.Optionen = handles.Optionen;
        Maschinendaten.Optionen.axes_Animate_Nut = handles.axes_Animate_Nut;
    end

    % Reservation of storage space
    Maschinendaten.Entwurf = struct;
    Maschinendaten.Kennfeld = struct;

%% 2 - Design process PMSM
% This function executes the design of the PMSM. More information about the calculations can be found in the function.

    [Maschinendaten] = Auslegung_PMSM(Maschinendaten, handles);

%% 3 - Efficiency map calculation PMSM
% This function executes the Efficiency map calculation of the PMSM. More information about the calculations can be found in the function.

    [Maschinendaten, handles_Kennfeld] = Kennfeld_PMSM(Maschinendaten);
    