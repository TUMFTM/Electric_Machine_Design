% -------------------------------------------------------------------------
% TU Munich - Institute of Automotive Technology
% -------------------------------------------------------------------------
% Modell for the design and analysis of PMSM or ASM (MEAPA)
% -------------------------------------------------------------------------
% Autor:    Svenja Kalt (svenja.kalt@tum.de), 
%           Jonathan Erhard 
% -------------------------------------------------------------------------

% Notes on the use of this script:
% (1) The input parameters can be set arbitrarily, e.g.
% rated can be passed to the script, accordingly the rated part
% must be commented out below. However, it is not mandatory to pass something
% has to be passed.
% (2) If the script is used without the main file, care must be taken to % ensure that the paths are correct. 
% that the paths are set correctly.
% (3) Since no user interface is called to select the winding, only the % classic design can be used.
% is called, only the classical design can be executed (no user
% user interaction necessary)

function [Entwurf, Analyse, Optimierer] = MEAPA_Script(rated)

%% Select machine type
% opt.machine type = rated.type; % 'ASM', 'PMSM'.
% opt.cooling type = rated.cooling type; % 1: air, 2: liquid
% opt.magnet arrangement = rated.magnet arrangement; % 1: 'SPMSM', 2: 'IPMSM (recessed)', 3: 'IPMSM (tangential)', 4: 'IPMSM (V-shape)';


    
%     if opt.Maschinentyp == 1
        opt.Maschinentyp                    = 'PMSM';
%     elseif opt.Maschinentyp == 2
%         opt.Maschinentyp                    = 'ASM';
%     end

%% Input parameters design
if(strcmp(opt.Maschinentyp,'ASM'))
%% RATED VALUES ASM
% % rated power P_N [W]
% rated.P_N = 3000;
%     
% % rated speed n_N [rpm]
% rated.n_N = 1500;
%     
% % rated voltage U_N [V]
% rated.U_N = 400;
%     
% % number of pole pairs p [-]
% rated.p = 2;
%     
% rated.frequency f_N [-]
% rated.f_N = (rated.p * rated.n_N) / 60;
%     
% % number of strands m [-]
% rated.m = 3;

    %% OPTIONS ASM
    % machine version
    opt.Maschinenausfuehrung            = 'Kaefiglaeufer';                      % 'Kaefiglaeufer'
    
    % Circuit
    opt.Schaltung                       = 'Dreieck';                            % 'Stern', 'Dreieck'
    
    % Coil form Stator
    opt.Spulenform_Stator               = 'Runddraht';                          % 'Runddraht'
    
    % Coil form Rotor
    opt.Spulenform_Rotor                = 'Runddraht';                          % 'Runddraht'
    
    % Slot form Stator
    opt.Nutform_Stator                  = 'Trapezform (eckig)';                 % 'Trapezform (eckig)'
    
    % Slot form Rotor
    opt.Nutform_Rotor                   = 'Trapezform (eckig)';                 % 'Trapezform (eckig)'
    
    % Cooling Stator
    %opt.Kuehlungsart                    = 'Innen- oder Kreislaufkuehlung';      % 'Oberflaechenkuehlung', 'Innen- oder Kreislaufkuehlung' 
    
%     if opt.Kuehlungsart == 1
%         opt.Kuehlungsart                    = 'Oberflaechenkuehlung';
%     else
        opt.Kuehlungsart                    = 'Innen- oder Kreislaufkuehlung';
%     end
    
    % iron material Stator
    opt.Stator_Eisenmaterial.String     = 'VACOFLUX 50';                        % 'M250-35A', 'M800-50A', 'VACOFLUX 48', 'VACOFLUX 50'
    opt.Stator_Eisenmaterial = loadMaterial(opt.Stator_Eisenmaterial,'Elektroblech');
    
    % Conductor material Stator
    opt.Stator_Leitermaterial.String    = 'Copper';                        % 'Aluminumwire', 'Aluminumcasting', 'Copper'
    opt.Stator_Leitermaterial = loadMaterial(opt.Stator_Leitermaterial,'Leiter');
    
    % Temperatur conductor Stator [°C]
    opt.theta_1                         = 90;
    
    % Iron material Rotor
    opt.Rotor_Eisenmaterial.String      = 'VACOFLUX 50';                   % 'M250-35A', 'M800-50A', 'VACOFLUX 48', 'VACOFLUX 50'
    opt.Rotor_Eisenmaterial = loadMaterial(opt.Stator_Eisenmaterial,'Elektroblech');
    
    % Iron material Rotor
    opt.Rotor_Leitermaterial.String     = 'Copper';                 % 'Aluminumwire', 'Aluminumcasting', 'Copper'
    opt.Rotor_Leitermaterial = loadMaterial(opt.Rotor_Leitermaterial,'Leiter');
    
    % Temperatur conductor Rotor [°C]
    opt.theta_2                         = 115;
    
    % Modus Winding design
    opt.Mode_Wicklung                   = 'Klassisch';                     % 'Klassisch'
    
    % Optimization goal winding
    opt.Wicklungstyp                    = 'A';                             % 'A','B','C'

     %% GUIDELINE VALUES ASM
    % Guide value for relative anchor length lambda [-]
    % Source: [Mueller08, p.577 - Table 9.1.3], [Meyer18, p.117]
    richt.lambda = 2.5;                                                   % between 0.6 and 1.0 for p=1, between 1.0 and 4.0 for p>1

    % Guideline value for the duct width of ventilation ducts l_v [m].
    % Source: [Meyer09, p.41], [Mueller08, p.585]
    richt.l_v = 0.01;                                                      % between 0.006 and 0.01

    % Guide value for the mean value of the air gap induction B_m [T].
    % Source: [Mueller08, p.582 - Table 9.1.5]
    richt.B_m = 0.58;                                                      % between 0.4 and 0.65

    % Standard value for the current coating A [A/mm]
    % Source: [Mueller08, p.580 - Table 9.1.4]
    % ATTENTION: specify either B_m or A
    %richt.A = 50.0;                                                       % between 20.0 and 120.0

    % Reference value for current density (stator) S_1 [A/mm^2]
    % Source: [Mueller08, p.580 - Table 9.1.4]
    richt.S_1 = 7.0;                                                       % between 3.0 and 8.0

    % Guide value for the max. permissible induction in the back (stator) B_1r_max [T]
    % Source: [Mueller08, p.582 - Table 9.1.5]
    richt.B_1r_max = 1.4;                                                  % between 1.3 and 1.65

    % Reference value for the max. permissible induction in the teeth (stator) B_1z_max [T].
    % Source: [Mueller08, p.582 - Table 9.1.5]
    richt.B_1z_max = 1.8;                                                 % between 1.4 and 2.1

    % Guide value for the slot filling factor (stator) phi_1n [-]
    % Source: [Mueller08, p.586 - Table 9.1.6]
    % V/A: Low voltage
    richt.phi_1n = 0.5;                                                    % between 0.3 and 0.5 for Round wire, between 0.35 and 0.6 for shaped coil or rod

    % Guide value for the winding factor (stator) xi_1p [-]
    % Source: [Mueller08, p.596]
    % V/A: no distinction between single and double layer winding,
    % only needed for initial estimation (will be calculated in detail in further
    % design process)
    richt.xi_1p = 0.96;                                                    % between 0.92 and 0.96

    % Guide value for minimum slot pitch (stator) tau_1n_min [m]
    % Source: [Meyer09, p.46], [Pyr14]
    richt.tau_1n_min = 0.007;                                              % between 0.007 and 0.07

    % Reference value for the iron fill factor (stator) phi_1Fe [-]
    % Source: [Mueller08, p.599]
    richt.phi_1Fe = 0.95;                                                  % between 0.9 and 1.0

    % Reference value for current density (rotor) S_2 [A/mm^2]
    % Source: [Mueller08, p.580 - Table 9.1.4]
    richt.S_2 = 5.0;                                                       % between 3.0 and 8.0 for Copper, between 3.0 and 6.5 for Aluminumcasting, between 3.0 and 6.5 for Aluminumwire

    % Guide value for bar current density (rotor) S_2s [A/mm^2]
    % Source: [Mueller08, p.580 - Table 9.1.4]
    richt.S_2s = 4.0;                                                      % between 3.0 and 8.0 for Copper, between 3.0 and 6.5 for Aluminumcasting, between 3.0 and 6.5 for Aluminumwire

    % Reference value for the ring current density (rotor) S_2r [A/mm^2]
    % Source: [Mueller08, p.580 - Table 9.1.4]
    richt.S_2r = 5.0;                                                      % between 3.0 and 8.0 for Copper, between 3.0 and 6.5 for Aluminumcasting, between 3.0 and 6.5 for Aluminumwire

    % Guide value for the max. permissible induction in the back (rotor) B_2r_max [T]
    % Source: [Mueller08, p.582 - table 9.1.5]
    richt.B_2r_max = 1.5;                                                  % between 0.4 and 1.6

    % Reference value for the max. permissible induction in the teeth (rotor) B_2z_max [T].
    % Source: [Mueller08, p.582 - table 9.1.5]
    richt.B_2z_max = 1.9;                                                  % between 1.5 and 2.2

    % Guide value for the slot filling factor (rotor) phi_2n [-]
    % Source: [Mueller08, p.586 - Table 9.1.6]
    % V/A: Low voltage
    richt.phi_2n = 0.5;                                                    % between 0.3 and 0.5 for Runddraht, between 0.35 and 0.6 for Shaped coil or rod

    % Guide value for the winding factor (rotor) xi_2p [-]
    % Source: [Mueller08, p.596]
    % V/A: no distinction between single and double layer winding,
    % only needed for initial estimation (will be calculated in detail in further
    % design process)
    richt.xi_2p = 0.96;                                                    % between 0.92 and 0.96

    % Guide value for minimum slot pitch (rotor) tau_2n_min [m]
    % Source: [Meyer09, p.46], [Pyr14]
    richt.tau_2n_min = 0.007;                                              % between 0.007 and 0.07

    % Guide value for the iron fill factor (rotor) phi_2Fe [-]
    % Source: [Mueller08, p.599]
    richt.phi_2Fe = 0.95;                                                  % between 0.9 and 1.0

elseif(strcmp(opt.Maschinentyp,'PMSM'))
%% RATED VALUES PMSM
% % rated power P_N [W]
% rated.P_N = 10000;
%     
% % rated speed n_N [rpm]
% rated.n_N = 3000;
%     
% % rated voltage U_N [V]
% rated.U_N = 400;
%     
% % number of pole pairs p [-]
% rated.p = 2;
%     
% rated.frequency f_N [-]
% rated.f_N = (rated.p * rated.n_N) / 60;
%     
% % power factor cos_phi_N [-]
% rated.cos_phi_N = 0.9;
%     
% % number of strands m [-]
% rated.m = 3;

    %% OPTIONS PMSM
    % Machine topology
    %opt.Maschinenausfuehrung            = 'IPMSM (eingelassen)';          % 'SPMSM', 'IPMSM (eingelassen)', 'IPMSM (tangential)', 'IPMSM (V-Form)'
    
%     if opt.Magnetanordnung == 1
%         opt.Maschinenausfuehrung                    = 'SPMSM';
%     elseif opt.Magnetanordnung == 2
%         opt.Maschinenausfuehrung                    = 'IPMSM (eingelassen)';
%     elseif opt.Magnetanordnung == 3
%         opt.Maschinenausfuehrung                    = 'IPMSM (tangential)';
%     elseif opt.Magnetanordnung == 4
        opt.Maschinenausfuehrung                    = 'IPMSM (V-Form)';
%     end
    
    % Circuit
    opt.Schaltung                       = 'Stern';                         % 'Stern', 'Dreieck'
    
    % Coil form Stator
    opt.Spulenform_Stator               = 'Runddraht';                     % 'Runddraht'
    
    % Slot form Stator
    opt.Nutform_Stator                  = 'Trapezform (eckig)';            % 'Trapezform (eckig)'
    
    % Cooling Stator
    %opt.Kuehlungsart                    = 'Luft (indirekt)';               % 'Luft (indirekt)', 'Wasser (direkt)'
    
%      if opt.Kuehlungsart == 1
%         opt.Kuehlungsart                    = 'Luft (indirekt)';
%     else
        opt.Kuehlungsart                    = 'Wasser (direkt)';
%     end
    
    % Iron material Stator
    opt.Stator_Eisenmaterial.String     = 'VACOFLUX 50';                   % 'M250-35A', 'M800-50A', 'VACOFLUX 48', 'VACOFLUX 50'
    opt.Stator_Eisenmaterial = loadMaterial(opt.Stator_Eisenmaterial,'Elektroblech');
    
    % Conductor material Stator
    opt.Stator_Leitermaterial.String    = 'Copper';                        % 'Aluminumwire', 'Aluminumcasting', 'Copper'
    opt.Stator_Leitermaterial = loadMaterial(opt.Stator_Leitermaterial,'Leiter');
    
    % Temperatur Conductormaterial Stator [°C]
    opt.theta_1                         = 90;
    
    % Iron material Rotor
    opt.Rotor_Eisenmaterial.String      = 'VACOFLUX 50';                   % 'M250-35A', 'M800-50A', 'VACOFLUX 48', 'VACOFLUX 50'
    opt.Rotor_Eisenmaterial = loadMaterial(opt.Stator_Eisenmaterial,'Elektroblech');
    
    % Magnetmaterial Rotor
    opt.Rotor_Magnetmaterial.String     = 'VACODYM 238 TP';                % 'VACODYM 238 TP', 'VACODYM 225 TP'
    opt.Rotor_Magnetmaterial = loadMaterial(opt.Rotor_Magnetmaterial,'Magnet');
    
    % Modus winding design
    opt.Mode_Wicklung                   = 'Klassisch';                     % 'Klassisch'
    
    % Optimization goal winding
    opt.Wicklungstyp                    = 'A';                             % 'A','B','C'

    %% GUIDELINE VALUES PMSM
    % Guideline value for relative anchor length lambda [-]
    % Source: [Mueller08, p.577 - Table 9.1.3], [Meyer18, p.117]
    richt.lambda = 2.5;                                                    % between 0.6 and 1.0 for p=1, between 1.0 and 4.0 for p>1

    % Guideline value for the duct width of ventilation ducts l_v [m].
    % Source: [Meyer09, p.41], [Mueller08, p.585]
    richt.l_v = 0.01;                                                      % between 0.006 and 0.01

    % Reference value for the amplitude of the air gap induction B_p [T].
    % Source: [Mueller08, p.582 - Table 9.1.5]
    richt.B_p = 0.85;                                                      % between 0.75 and 1.05

    % Standard value for the current coating A [A/mm]
    % Source: [Mueller08, p.580 - Table 9.1.4]
    % ATTENTION: specify either B_p or A
    %richt.A = 50.0;                                                       % between 30.0 and 120.0 for Luft (indirekt), between 160.0 and 300.0 for Wasser (direkt)

    % Reference value for current density (stator) S_1 [A/mm^2]
    % Source: [Mueller08, p.580 - Table 9.1.4]
    richt.S_1 = 7.0;                                                       % between 3.0 und 7.0 fuer Luft (indirekt), zwischen 13.0 und 18.0 fuer Wasser (direkt)

    % Guide value for the max. permissible induction in the back (stator) B_1r_max [T]
    % Source: [Mueller08, p.582 - Table 9.1.5]
    richt.B_1r_max = 1.4;                                                  % between 1.0 and 1.5

    % Reference value for the max. permissible induction in the teeth (stator) B_1z_max [T].
    % Source: [Mueller08, p.582 - Table 9.1.5]
    richt.B_1z_max = 1.8;                                                  % between 1.6 and 2.0

    % Guide value for the slot filling factor (stator) phi_1n [-]
    % Source: [Mueller08, p.586 - Table 9.1.6]
    % V/A: Low voltage
    richt.phi_1n = 0.5;                                                    % between 0.3 and 0.5 for conductor

    % Reference value for the winding factor (stator) xi_1p [-]
    % Source: [Mueller08, p.596]
    % V/A: no distinction between single and double layer winding,
    % only needed for initial estimation (will be calculated in detail in further
    % design process)
    richt.xi_1p = 0.96;                                                    % between 0.92 and 0.96

    % Guide value for minimum slot pitch (stator) tau_1n_min [m]
    % Source: [Meyer09, p.46], [Pyr14]
    richt.tau_1n_min = 0.007;                                              % between 0.007 and 0.07

    % Reference value for the iron fill factor (stator) phi_1Fe [-]
    % Source: [Mueller08, p.599]
    richt.phi_1Fe = 0.95;                                                  % between 0.9 and 1.0

    % Guide value for the max. permissible induction in the back (rotor) B_2r_max [T]
    % Source: [Mueller08, p.582 - table 9.1.5]
    richt.B_2r_max = 1.4;                                                  % between 1.0 and 1.5

    % Guide value for the iron fill factor (rotor) phi_2Fe [-]
    % Source: [Mueller08, p.599]
    richt.phi_2Fe = 0.95;                                                  % between 0.9 and 1.0
else
    error('Ungueltige Eingabe bei Variable "opt.Maschinentyp"');
end

%% Restore
handles.rated = rated;
handles.richt = richt;
handles.opt = opt;
clear rated richt opt

%% Start design
if(strcmp(handles.opt.Maschinentyp,'ASM'))
    [handles.Entwurf] = Design_ASM(handles);
elseif(strcmp(handles.opt.Maschinentyp,'PMSM'))
    [handles.Entwurf] = Design_PMSM(handles);
else
    error('Ungueltige Eingabe bei Variable "Entwurf.Optionen.Maschinentyp"');
end
disp('Entwurf abgeschlossen');

%% Inputparameter design
    %% OPTIONS LOSSES
    % winding losses
    handles.opt.P_vw = 1;                                                  % 0, 1
    
    % Remagnetization losses
    handles.opt.P_vu = 1;                                                  % 0, 1
    
    % Iron loss model
    handles.opt.P_vu_Modell = 'Modellansatz Jordan';                       % 'Modellansatz Jordan', 'Abschaetzung ueber verketteten Fluss'

    
    % mech. losses
    handles.opt.P_vme = 1;                                                 % 0, 1
    
    % additional losses
    handles.opt.P_vzus = 1;                                                % 0, 1

    %% OPTIONS CALCULATION
    % Generator
    handles.opt.Generator = 0;                                             % 0, 1
    
    % max. rotational speed [U/min]
    handles.opt.n_max = handles.rated.nmax;
    
    % Resolution rotational speed
    handles.opt.n_tics = 160;
    
    % Resolution torque
    handles.opt.M_tics = 160;
    
    % Control max. voltage [V]
    handles.opt.u_1max = handles.Entwurf.EMAG.U_1Str * sqrt(2);
    
    % Control max. Current [A]
    handles.opt.i_1max = handles.Entwurf.EMAG.I_1Str * sqrt(2);

%% Start analysis
if(strcmp(handles.Entwurf.Optionen.Maschinentyp,'ASM'))
    [handles.Analyse] = Analysis_ASM(handles);
elseif(strcmp(handles.Entwurf.Optionen.Maschinentyp,'PMSM'))
    [handles.Analyse] = Analysis_PMSM(handles);
else
    error('Ungueltige Eingabe bei Variable "handles.Entwurf.Optionen.Maschinentyp"');
end
disp('Analyse abgeschlossen');

Entwurf = handles.Entwurf;
Analyse = handles.Analyse;

%% Save
% 
handles.opt.unique_id = datestr(now,'yyyymmdd_HHMMSS');
handles.opt.folder_id = [handles.opt.Maschinentyp '_data_' handles.opt.unique_id];
handles.opt.file_id = [handles.opt.Maschinentyp,'_',handles.opt.unique_id];
handles.text_ID.String = ['ID: Entwurf_' handles.opt.file_id];
% mkdir('3_Ergebnisse',handles.opt.folder_id);
% mkdir(['3_Ergebnisse/',handles.opt.folder_id],'1_Entwurf');
% mkdir(['3_Ergebnisse/',handles.opt.folder_id],'2_Analyse');
% 
% % Save struct Design
% save(['3_Ergebnisse/',handles.opt.folder_id,'/1_Entwurf/Entwurf_',handles.opt.file_id,'.mat'],'Entwurf');
% % Save DXF file
% 
% saveDXF(handles,handles.Entwurf.Geometrie)

% Save Excel file
%	saveExcel(handles, Entwurf)
    
% Export Excel
% function export_table = export_excel(input)
%     var = fieldnames(input);
%     
%     for i = 1:length(var)
%         export_struct{i,1} = input.(var{i});
%     end
%     
%     export_table = table(export_struct,'VariableNames',{'var'},'RowNames',var);
% end
% Save Excel
% function saveExcel(handles, Entwurf)
%     copyfile('3_Ergebnisse/1_Misc/Vorlage_Ergebnisse.xlsx',['3_Ergebnisse/',handles.opt.folder_id,'/1_Entwurf/Entwurf_',handles.opt.file_id,'.xlsx']);
%     
%     rated_export = export_excel(Entwurf.Bemessungswerte);
%     richt_export = export_excel(Entwurf.Richtwerte);
%     
% % Entwurf.Optionen = rmfield(Entwurf.Optionen,{'Content','Locked','Saved','file_id','folder_id','unique_id'});
%     var1 = fieldnames(Entwurf.Optionen);
%     for i = 1:length(var1)
%         if(isstruct(Entwurf.Optionen.(var1{i})))
%             Entwurf.Optionen.(var1{i}) = Entwurf.Optionen.(var1{i}).String;
%         end
%     end
%     opt_export = export_excel(Entwurf.Optionen);
%     
%     range_oben = 2;
%     range_unten = range_oben + height(rated_export) - 1;
%     range = ['B' num2str(range_oben) ':C' num2str(range_unten)];
%     writetable(rated_export,['3_Ergebnisse/',handles.opt.folder_id,'/1_Entwurf/Entwurf_',handles.opt.file_id,'.xlsx'],'Sheet','Bemessungswerte','Range',range,'WriteVariableNames',0,'WriteRowNames',1)
%     
%     range_oben = 2;
%     range_unten = range_oben + height(richt_export) - 1;
%     range = ['B' num2str(range_oben) ':C' num2str(range_unten)];
%     writetable(richt_export,['3_Ergebnisse/',handles.opt.folder_id,'/1_Entwurf/Entwurf_',handles.opt.file_id,'.xlsx'],'Sheet','Richtwerte','Range',range,'WriteVariableNames',0,'WriteRowNames',1)
%     
%     range_oben = 2;
%     range_unten = range_oben + height(opt_export) - 1;
%     range = ['B' num2str(range_oben) ':C' num2str(range_unten)];
%     writetable(opt_export,['3_Ergebnisse/',handles.opt.folder_id,'/1_Entwurf/Entwurf_',handles.opt.file_id,'.xlsx'],'Sheet','Optionen','Range',range,'WriteVariableNames',0,'WriteRowNames',1)
%     
%         wick_export = export_excel(Entwurf.Wicklung);
%         
%         Entwurf.Geometrie.Nut_1 = rmfield(Entwurf.Geometrie.Nut_1,{'iter','maxIter'});
%         geo_Nut_1_export = export_excel(Entwurf.Geometrie.Nut_1);
%         
%          if(strcmp(handles.Entwurf.Optionen.Maschinentyp,'ASM'))
%              Entwurf.Geometrie.Nut_2 = rmfield(Entwurf.Geometrie.Nut_2,{'iter','maxIter'});
%              geo_Nut_2_export = export_excel(Entwurf.Geometrie.Nut_2);
%              Entwurf.Geometrie = rmfield(Entwurf.Geometrie,{'Nut_1','Nut_2'});
%          elseif(strcmp(handles.Entwurf.Optionen.Maschinentyp,'PMSM'))
%             Entwurf.Geometrie = rmfield(Entwurf.Geometrie,{'Nut_1'});
%         else
%             error('Ungueltige Eingabe bei Variable "Entwurf.Optionen.Maschinentyp"');
%         end
%         
%         var2 = [fieldnames(Entwurf.Geometrie)' fieldnames(Entwurf.Geometrie.misc)'; struct2cell(Entwurf.Geometrie)' struct2cell(Entwurf.Geometrie.misc)'];
%         Entwurf.Geometrie = struct(var2{:});
%         Entwurf.Geometrie = orderfields(Entwurf.Geometrie);
%         Entwurf.Geometrie = rmfield(Entwurf.Geometrie,{'misc'});
%         var1 = fieldnames(Entwurf.Geometrie);
%         for i = 1:length(var1)
%             if(length(Entwurf.Geometrie.(var1{i}))>1)
%                 Entwurf.Geometrie = rmfield(Entwurf.Geometrie,var1(i));
%             end
%         end
%         geo_export = export_excel(Entwurf.Geometrie);
%         
%         var2 = [fieldnames(Entwurf.EMAG)' fieldnames(Entwurf.EMAG.misc)'; struct2cell(Entwurf.EMAG)' struct2cell(Entwurf.EMAG.misc)'];
%         Entwurf.EMAG = struct(var2{:});
%         Entwurf.EMAG = orderfields(Entwurf.EMAG);
%         if(strcmp(handles.Entwurf.Optionen.Maschinentyp,'ASM'))
%         elseif(strcmp(handles.Entwurf.Optionen.Maschinentyp,'PMSM'))
%             Entwurf.EMAG = rmfield(Entwurf.EMAG,{'misc'});
%         else
%             error('Ungueltige Eingabe bei Variable "Entwurf.Optionen.Maschinentyp"');
%         end
%         emag_export = export_excel(Entwurf.EMAG);
%         
%         range_oben = 2;
%         range_unten = range_oben + height(wick_export) - 1;
%         range = ['B' num2str(range_oben) ':C' num2str(range_unten)];
%         writetable(wick_export,['3_Ergebnisse/',handles.opt.folder_id,'/1_Entwurf/Entwurf_',handles.opt.file_id,'.xlsx'],'Sheet','Wicklung','Range',range,'WriteVariableNames',0,'WriteRowNames',1)
% 
%         range_oben = 2;
%         range_unten = range_oben + height(geo_export) - 1;
%         range = ['B' num2str(range_oben) ':C' num2str(range_unten)];
%         writetable(geo_export,['3_Ergebnisse/',handles.opt.folder_id,'/1_Entwurf/Entwurf_',handles.opt.file_id,'.xlsx'],'Sheet','Geometrie','Range',range,'WriteVariableNames',0,'WriteRowNames',1)
% 
%         range_oben = range_unten + 2;
%         range_unten = range_oben + height(geo_Nut_1_export) - 1;
%         range = ['B' num2str(range_oben) ':C' num2str(range_unten)];
%         writetable(geo_Nut_1_export,['3_Ergebnisse/',handles.opt.folder_id,'/1_Entwurf/Entwurf_',handles.opt.file_id,'.xlsx'],'Sheet','Geometrie','Range',range,'WriteVariableNames',0,'WriteRowNames',1)
% 
%         if(strcmp(handles.Entwurf.Optionen.Maschinentyp,'ASM'))
%             range_oben = range_unten + 2;
%             range_unten = range_oben + height(geo_Nut_2_export) - 1;
%             range = ['B' num2str(range_oben) ':C' num2str(range_unten)];
%             writetable(geo_Nut_2_export,['3_Ergebnisse/',handles.opt.folder_id,'/1_Entwurf/Entwurf_',handles.opt.file_id,'.xlsx'],'Sheet','Geometrie','Range',range,'WriteVariableNames',0,'WriteRowNames',1) 
%         elseif(strcmp(handles.Entwurf.Optionen.Maschinentyp,'PMSM'))
%         else
%             error('Ungueltige Eingabe bei Variable "Entwurf.Optionen.Maschinentyp"');
%         end
%         
%         range_oben = 2;
%         range_unten = range_oben + height(emag_export) - 1;
%         range = ['B' num2str(range_oben) ':C' num2str(range_unten)];
%         writetable(emag_export,['3_Ergebnisse/',handles.opt.folder_id,'/1_Entwurf/Entwurf_',handles.opt.file_id,'.xlsx'],'Sheet','EMAG','Range',range,'WriteVariableNames',0,'WriteRowNames',1)
%     end


% % Save Geometrie in DXF file
% function saveDXF(handles,var)
%     % Open DXF file
%     FID = dxf_open(['3_Ergebnisse/',handles.opt.folder_id,'/1_Entwurf/Geometrie_',handles.opt.file_id,'.dxf']);
% 
%     % Produce Stator aussen
%     dxf_polyline(FID,var.Stator_aussen_x,var.Stator_aussen_y,zeros(length(var.Stator_aussen_x),1));
% 
%     % Produce Stator innen
%     dxf_polyline(FID,var.Stator_innen_x,var.Stator_innen_y,zeros(length(var.Stator_innen_x),1));
% 
%     % Produce Rotor aussen
%     dxf_polyline(FID,var.Rotor_aussen_x,var.Rotor_aussen_y,zeros(length(var.Rotor_aussen_x),1));
% 
%     % Produce Rotor innen
%     dxf_polyline(FID,var.Rotor_innen_x,var.Rotor_innen_y,zeros(length(var.Rotor_innen_x),1));
%     
%     % Produce Stator Nutfuellung
%     for i = 1:handles.Entwurf.Wicklung.N_1
%         dxf_polyline(FID,var.Fuellung_Stator_x((length(var.Fuellung_Stator_x)/handles.Entwurf.Wicklung.N_1)*(i-1)+1:(length(var.Fuellung_Stator_x)/handles.Entwurf.Wicklung.N_1)*i,1),var.Fuellung_Stator_y((length(var.Fuellung_Stator_y)/handles.Entwurf.Wicklung.N_1)*(i-1)+1:(length(var.Fuellung_Stator_y)/handles.Entwurf.Wicklung.N_1)*i,1),zeros((length(var.Fuellung_Stator_y)/handles.Entwurf.Wicklung.N_1),1));
%     end
%     
%     if(strcmp(handles.opt.Maschinentyp,'ASM'))
%         % Produce Rotor Nutfuellung
%         for i = 1:handles.Entwurf.Wicklung.N_2
%             dxf_polyline(FID,var.Fuellung_Rotor_x((length(var.Fuellung_Rotor_x)/handles.Entwurf.Wicklung.N_2)*(i-1)+1:(length(var.Fuellung_Rotor_x)/handles.Entwurf.Wicklung.N_2)*i,1),var.Fuellung_Rotor_y((length(var.Fuellung_Rotor_y)/handles.Entwurf.Wicklung.N_2)*(i-1)+1:(length(var.Fuellung_Rotor_y)/handles.Entwurf.Wicklung.N_2)*i,1),zeros((length(var.Fuellung_Rotor_y)/handles.Entwurf.Wicklung.N_2),1));
%         end
%     elseif(strcmp(handles.opt.Maschinentyp,'PMSM'))
%         % Produce Magnets
%         if(strcmp(handles.opt.Maschinenausfuehrung,'SPMSM'))
%             for i = 1:(2*handles.Entwurf.Bemessungswerte.p)
%                 dxf_polyline(FID,var.Magnet_Rotor_x(var.points_Magnet(i):var.points_Magnet(i+1)-1,1),var.Magnet_Rotor_y(var.points_Magnet(i):var.points_Magnet(i+1)-1,1),zeros(length(var.points_Magnet(i):var.points_Magnet(i+1)-1),1));
%             end
%         elseif(strcmp(handles.opt.Maschinenausfuehrung,'IPMSM (eingelassen)'))
%             for i = 1:(2*handles.Entwurf.Bemessungswerte.p)
%                 dxf_polyline(FID,var.Magnet_Rotor_x(var.points_Magnet(i):var.points_Magnet(i+1)-1,1),var.Magnet_Rotor_y(var.points_Magnet(i):var.points_Magnet(i+1)-1,1),zeros(length(var.points_Magnet(i):var.points_Magnet(i+1)-1),1));
%             end
%         elseif(strcmp(handles.opt.Maschinenausfuehrung,'IPMSM (tangential)'))
%             for i = 1:(2*handles.Entwurf.Bemessungswerte.p)
%                 dxf_polyline(FID,var.Magnet_Rotor_x((length(var.Magnet_Rotor_x)/(2*handles.Entwurf.Bemessungswerte.p))*(i-1)+1:(length(var.Magnet_Rotor_x)/(2*handles.Entwurf.Bemessungswerte.p))*i,1),var.Magnet_Rotor_y((length(var.Magnet_Rotor_y)/(2*handles.Entwurf.Bemessungswerte.p))*(i-1)+1:(length(var.Magnet_Rotor_y)/(2*handles.Entwurf.Bemessungswerte.p))*i,1),zeros((length(var.Magnet_Rotor_y)/(2*handles.Entwurf.Bemessungswerte.p)),1));
%             end    
%         elseif(strcmp(handles.opt.Maschinenausfuehrung,'IPMSM (V-Form)'))
%             for i = 1:(4*handles.Entwurf.Bemessungswerte.p)
%                 dxf_polyline(FID,var.Magnet_Rotor_x((length(var.Magnet_Rotor_x)/(4*handles.Entwurf.Bemessungswerte.p))*(i-1)+1:(length(var.Magnet_Rotor_x)/(4*handles.Entwurf.Bemessungswerte.p))*i,1),var.Magnet_Rotor_y((length(var.Magnet_Rotor_y)/(4*handles.Entwurf.Bemessungswerte.p))*(i-1)+1:(length(var.Magnet_Rotor_y)/(4*handles.Entwurf.Bemessungswerte.p))*i,1),zeros((length(var.Magnet_Rotor_y)/(4*handles.Entwurf.Bemessungswerte.p)),1));
%             end
%         else
%             error('Ungueltige Eingabe bei Variable "opt.Maschinenausfuehrung"')
%         end
%     else
%         error('Ungueltige Eingabe bei Variable "handles.opt.Maschinentyp"');
%     end
%     
%     % Close DXF file
%     dxf_close(FID);
% end

%% Longitudinal dynamics simulation Initialisation 

% Handover to LDS
    Maschinendaten = struct;
    Maschinendaten.Entwurf = handles.Entwurf;
    Maschinendaten.Analyse = handles.Analyse;
    Maschinendaten.Entwurf.Optionen.folder_id = handles.opt.folder_id;
    Maschinendaten.Entwurf.Optionen.file_id = handles.opt.file_id;
    
    handles_GUI_LDS = handles.rated.LDS;

    %% Interface LDS
    LDS_values = Interface_LDS(Maschinendaten); 

    %LDS GUI
%     handles_GUI_LDS = GUI_LDS;
% 
%     uiwait(handles_GUI_LDS.output,5)

    FZD_LDS=struct;

    FZD_LDS.fz_m = str2num(handles_GUI_LDS.fz_m.String); 
    FZD_LDS.cW = str2num(handles_GUI_LDS.cW.String); 
    FZD_LDS.A = str2num(handles_GUI_LDS.A.String);
    FZD_LDS.tyre_r = str2num(handles_GUI_LDS.tyre_r.String);
    FZD_LDS.battcap = str2num(handles_GUI_LDS.battcap.String);
    FZD_LDS.aux = str2num(handles_GUI_LDS.aux.String);
    FZD_LDS.GearRatio = str2num(handles_GUI_LDS.GearRatio.String);
    FZD_LDS.AnzMasch = 'GM_X'; %str2num(handles_GUI_LDS.AnzMasch.String);

    %% Save vehicle data for LDS

    folder_id = [Maschinendaten.Entwurf.Optionen.folder_id]; %file_id, '_data'];
    file_id = [Maschinendaten.Entwurf.Optionen.file_id];
    %LD_export = export_excel(struct2table(FZD_LDS));
    %copyfile('Ergebnisse/Vorlage_Ergebnisse.xlsx',['Ergebnisse/',folder_id,'/Ergebnisse_',file_id,'.xlsx']);
    %writetable(LD_export,['Ergebnisse/',folder_id,'/Ergebnisse_',file_id,'.xlsx'],'Sheet','LDS','Range','C3:C10','WriteVariableNames',0,'WriteRowNames',0)
    %xlswrite(['3_Ergebnisse/',folder_id, '/1_Entwurf','/Entwurf_',file_id,'.xlsx'],(struct2array(FZD_LDS))','LDS','C3:C9')

    %% Start LDS
    % Convert to strings
    FZD_LDS.fz_m = handles_GUI_LDS.fz_m.String; 
    FZD_LDS.cW = handles_GUI_LDS.cW.String; 
    FZD_LDS.A = handles_GUI_LDS.A.String;
    FZD_LDS.tyre_r = handles_GUI_LDS.tyre_r.String;
    FZD_LDS.battcap = handles_GUI_LDS.battcap.String;
    FZD_LDS.aux = handles_GUI_LDS.aux.String;
    FZD_LDS.GearRatio = handles_GUI_LDS.GearRatio.String;

    [FZG_LDS1, Optimierer] = RUN_simulation(LDS_values, handles_GUI_LDS, FZD_LDS,Maschinendaten);


% Save map neu
        % Save map
        KennfeldMaschine = struct;
        KennfeldMaschineSpeichern = Maschinendaten;
        handles = Maschinendaten;
        KennfeldMaschine.n = KennfeldMaschineSpeichern.Analyse.Betriebsdaten.n_m_mesh;
        KennfeldMaschine.M = KennfeldMaschineSpeichern.Analyse.Betriebsdaten.M_max_mesh;
        KennfeldMaschine.etages = KennfeldMaschineSpeichern.Analyse.Wirkungsgrad.eta_ges_mesh;
        %save (['3_Ergebnisse/',handles.Entwurf.Optionen.folder_id,'/','2_Analyse/', 'KennfeldMaschine_Gesamtwirkungsgrad','.mat'],'KennfeldMaschine');
        folder_id = [KennfeldMaschineSpeichern.Entwurf.Optionen.folder_id]; %file_id, '_data'];
        file_id = [KennfeldMaschineSpeichern.Entwurf.Optionen.file_id  ];
        
        % Convert data to correct form
        % Design sizes
        Kennfeld_Gesamtwirkungsgrad=struct;
        Kennfeld_Gesamtwirkungsgrad.type = KennfeldMaschineSpeichern.Entwurf.Optionen.Maschinentyp;
        Kennfeld_Gesamtwirkungsgrad.power = KennfeldMaschineSpeichern.Entwurf.Bemessungswerte.P_N/1000;
        Kennfeld_Gesamtwirkungsgrad.n_n = KennfeldMaschineSpeichern.Entwurf.Bemessungswerte.n_N;
        Kennfeld_Gesamtwirkungsgrad.n_max = KennfeldMaschineSpeichern.Analyse.Optionen.n_max;
        Kennfeld_Gesamtwirkungsgrad.U_n = KennfeldMaschineSpeichern.Entwurf.Bemessungswerte.U_N;
        
        % Speed axis - X-axis map
        Kennfeld_Gesamtwirkungsgrad.eff_n_axis = transpose(KennfeldMaschineSpeichern.Analyse.Betriebsdaten.mot.n_m_vec);

        % Torque - Y-Achse
        KennfeldMaschineSpeichern.Analyse.Betriebsdaten.mot.M_vec=transpose(fliplr(KennfeldMaschineSpeichern.Analyse.Betriebsdaten.mot.M_vec));
        Kennfeld_Gesamtwirkungsgrad.eff_T_axis = [(fliplr(transpose(KennfeldMaschineSpeichern.Analyse.Betriebsdaten.mot.M_vec(2:end)))) * (-1), transpose(KennfeldMaschineSpeichern.Analyse.Betriebsdaten.mot.M_vec)];
        %Kennfeld_Gesamtwirkungsgrad.eff_T_axis = [fliplr(transpose(KennfeldMaschineSpeichern.Analyse.Betriebsdaten.mot.M_vec)),(transpose(KennfeldMaschineSpeichern.Analyse.Betriebsdaten.mot.M_vec(2:end))) * (-1)];

       %Diagram                                                             
        % Manipulate efficiency map
        KennfeldMaschineSpeichern.Analyse.Wirkungsgrad.eta_ges_mesh=flipud(KennfeldMaschineSpeichern.Analyse.Wirkungsgrad.eta_ges_mesh); %hier Unterschied Verlustarten
        KennfeldMaschineSpeichern.Analyse.Wirkungsgrad.eta_ges_mesh(KennfeldMaschineSpeichern.Analyse.Wirkungsgrad.eta_ges_mesh<0.1) = 0.1; %Hier Unterschied Verlustarten %set small efficiency values to 0.1
        KennfeldMaschineSpeichern.Analyse.Wirkungsgrad.eta_ges_mesh(1,:) = 0; %first colomn is zero %hier Unterschied Verlustarten
        KennfeldMaschineSpeichern.Analyse.Wirkungsgrad.eta_ges_mesh(:,1) = 0.1; %first row is zero %hier Unterschied Verlustarten
        % Create one "larger" eff. map that includes acceleration and recuperation
        Kennfeld_Gesamtwirkungsgrad.eff_recu = 1./KennfeldMaschineSpeichern.Analyse.Wirkungsgrad.eta_ges_mesh; %recuperation efficiency %hier Unterschied Verlustarten
        Kennfeld_Gesamtwirkungsgrad.eff_recu(Kennfeld_Gesamtwirkungsgrad.eff_recu==inf)=100;
        Kennfeld_Gesamtwirkungsgrad.eff_recu = Kennfeld_Gesamtwirkungsgrad.eff_recu(2:end,:);
        Kennfeld_Gesamtwirkungsgrad.eff = [flipud(Kennfeld_Gesamtwirkungsgrad.eff_recu); KennfeldMaschineSpeichern.Analyse.Wirkungsgrad.eta_ges_mesh]; %hier Unterschied Verlustarten
    
        %Full load characteristic
        %motor.T_max = transpose (Maschinendaten.Analyse.Betriebsdaten.mot.M_max_vec);                                                            %noch auf 2000 interpolieren?
        KennfeldMaschineSpeichern.Analyse.Betriebsdaten.mot.M_max_vec=transpose(KennfeldMaschineSpeichern.Analyse.Betriebsdaten.mot.M_max_vec);
        Kennfeld_Gesamtwirkungsgrad.T_max = transpose(KennfeldMaschineSpeichern.Analyse.Betriebsdaten.mot.M_max_vec);  
        
        % Skale torque
        Kennfeld_Gesamtwirkungsgrad.T_max_n_axis = transpose(KennfeldMaschineSpeichern.Analyse.Betriebsdaten.mot.M_vec);                                                                          %noch auf 2000 interpolieren?

%         save (['3_Ergebnisse/',handles.Entwurf.Optionen.folder_id,'/','2_Analyse/', KennfeldMaschineSpeichern.Entwurf.Optionen.Maschinentyp,'_',num2str(KennfeldMaschineSpeichern.Entwurf.Bemessungswerte.P_N/1000),'_',num2str(KennfeldMaschineSpeichern.Entwurf.Bemessungswerte.n_N),'_',num2str(KennfeldMaschineSpeichern.Analyse.Optionen.n_max),'_',num2str(KennfeldMaschineSpeichern.Entwurf.Bemessungswerte.U_N),'_','Kennfeld_Gesamtwirkungsgrad','.mat'],'Kennfeld_Gesamtwirkungsgrad');
%         save (['LDS/Vehicle/Para_Powertrain/Motor_efficiency/MOT_memory//', KennfeldMaschineSpeichern.Entwurf.Optionen.Maschinentyp,'_',num2str(KennfeldMaschineSpeichern.Entwurf.Bemessungswerte.P_N/1000),'_',num2str(KennfeldMaschineSpeichern.Entwurf.Bemessungswerte.n_N),'_',num2str(KennfeldMaschineSpeichern.Analyse.Optionen.n_max),'_',num2str(KennfeldMaschineSpeichern.Entwurf.Bemessungswerte.U_N),'_', 'Kennfeld_Gesamtwirkungsgrad' ,'.mat'],'Kennfeld_Gesamtwirkungsgrad');
                 
        % Save map
        KennfeldMaschine_Verluste_gesamt = struct;
        KennfeldMaschine_Verluste_gesamt_Speichern = handles;
        KennfeldMaschine_Verluste_gesamt.n = KennfeldMaschine_Verluste_gesamt_Speichern.Analyse.Betriebsdaten.n_m_mesh;
        KennfeldMaschine_Verluste_gesamt.M = KennfeldMaschine_Verluste_gesamt_Speichern.Analyse.Betriebsdaten.M_max_mesh;
        KennfeldMaschine_Verluste_gesamt.verlges = KennfeldMaschine_Verluste_gesamt_Speichern.Analyse.Verluste.P_vges_mesh;
        %save (['3_Ergebnisse/',handles.Entwurf.Optionen.folder_id,'/','2_Analyse/', 'KennfeldMaschine_Verluste_gesamt','.mat'],'KennfeldMaschine_Verluste_gesamt');
        folder_id = [KennfeldMaschine_Verluste_gesamt_Speichern.Entwurf.Optionen.folder_id]; %file_id, '_data'];
        file_id = [KennfeldMaschine_Verluste_gesamt_Speichern.Entwurf.Optionen.file_id  ];
        
        % Convert data to correct form
        % Design sizes
        Kennfeld_Gesamtverluste=struct;
        Kennfeld_Gesamtverluste.type = KennfeldMaschine_Verluste_gesamt_Speichern.Entwurf.Optionen.Maschinentyp;
        Kennfeld_Gesamtverluste.power = KennfeldMaschine_Verluste_gesamt_Speichern.Entwurf.Bemessungswerte.P_N/1000;
        Kennfeld_Gesamtverluste.n_n = KennfeldMaschine_Verluste_gesamt_Speichern.Entwurf.Bemessungswerte.n_N;
        Kennfeld_Gesamtverluste.n_max = KennfeldMaschine_Verluste_gesamt_Speichern.Analyse.Optionen.n_max;
        Kennfeld_Gesamtverluste.U_n = KennfeldMaschine_Verluste_gesamt_Speichern.Entwurf.Bemessungswerte.U_N;
        
        % Speed axis - X-axis map
        Kennfeld_Gesamtverluste.eff_n_axis = transpose(KennfeldMaschine_Verluste_gesamt_Speichern.Analyse.Betriebsdaten.mot.n_m_vec);

        % Torque - Y-Achse
        KennfeldMaschine_Verluste_gesamt_Speichern.Analyse.Betriebsdaten.mot.M_vec=transpose(fliplr(KennfeldMaschine_Verluste_gesamt_Speichern.Analyse.Betriebsdaten.mot.M_vec));
        Kennfeld_Gesamtverluste.eff_T_axis = [(fliplr(transpose(KennfeldMaschine_Verluste_gesamt_Speichern.Analyse.Betriebsdaten.mot.M_vec(2:end)))) * (-1), transpose(KennfeldMaschine_Verluste_gesamt_Speichern.Analyse.Betriebsdaten.mot.M_vec)];
        %Kennfeld_Gesamtwirkungsgrad.eff_T_axis = [fliplr(transpose(KennfeldMaschineSpeichern.Analyse.Betriebsdaten.mot.M_vec)),(transpose(KennfeldMaschineSpeichern.Analyse.Betriebsdaten.mot.M_vec(2:end))) * (-1)];

        % Diagram                                                            
        % Manipulate efficiency map
        KennfeldMaschine_Verluste_gesamt_Speichern.Analyse.Verluste.P_vges_mesh=flipud(KennfeldMaschine_Verluste_gesamt_Speichern.Analyse.Verluste.P_vges_mesh); %hier Unterschied Verlustarten
        KennfeldMaschine_Verluste_gesamt_Speichern.Analyse.Verluste.P_vges_mesh(KennfeldMaschine_Verluste_gesamt_Speichern.Analyse.Verluste.P_vges_mesh<0.1) = 0.1; %Hier Unterschied Verlustarten %set small efficiency values to 0.1
        KennfeldMaschine_Verluste_gesamt_Speichern.Analyse.Verluste.P_vges_mesh(1,:) = 0; %first colomn is zero %hier Unterschied Verlustarten
        KennfeldMaschine_Verluste_gesamt_Speichern.Analyse.Verluste.P_vges_mesh(:,1) = 0.1; %first row is zero %hier Unterschied Verlustarten
        % Create one "larger" eff. map that includes acceleration and recuperation
        Kennfeld_Gesamtverluste.eff_recu = 1./KennfeldMaschine_Verluste_gesamt_Speichern.Analyse.Verluste.P_vges_mesh; %recuperation efficiency %hier Unterschied Verlustarten
        Kennfeld_Gesamtverluste.eff_recu(Kennfeld_Gesamtverluste.eff_recu==inf)=100;
        Kennfeld_Gesamtverluste.eff_recu = Kennfeld_Gesamtverluste.eff_recu(2:end,:);
        Kennfeld_Gesamtverluste.eff = [flipud(Kennfeld_Gesamtverluste.eff_recu); KennfeldMaschine_Verluste_gesamt_Speichern.Analyse.Verluste.P_vges_mesh]; %hier Unterschied Verlustarten
    
        % Full load characteristic
        %motor.T_max = transpose (Maschinendaten.Analyse.Betriebsdaten.mot.M_max_vec);                                                            %noch auf 2000 interpolieren?
        KennfeldMaschine_Verluste_gesamt_Speichern.Analyse.Betriebsdaten.mot.M_max_vec=transpose(KennfeldMaschine_Verluste_gesamt_Speichern.Analyse.Betriebsdaten.mot.M_max_vec);
        Kennfeld_Gesamtverluste.T_max = transpose(KennfeldMaschine_Verluste_gesamt_Speichern.Analyse.Betriebsdaten.mot.M_max_vec);  
        
        % Speed axis scaled
        Kennfeld_Gesamtverluste.T_max_n_axis = transpose(KennfeldMaschine_Verluste_gesamt_Speichern.Analyse.Betriebsdaten.mot.M_vec);                                                                          %noch auf 2000 interpolieren?
        
%         save (['3_Ergebnisse/',handles.Entwurf.Optionen.folder_id,'/','2_Analyse/', KennfeldMaschine_Verluste_gesamt_Speichern.Entwurf.Optionen.Maschinentyp,'_',num2str(KennfeldMaschine_Verluste_gesamt_Speichern.Entwurf.Bemessungswerte.P_N/1000),'_',num2str(KennfeldMaschine_Verluste_gesamt_Speichern.Entwurf.Bemessungswerte.n_N),'_',num2str(KennfeldMaschine_Verluste_gesamt_Speichern.Analyse.Optionen.n_max),'_',num2str(KennfeldMaschine_Verluste_gesamt_Speichern.Entwurf.Bemessungswerte.U_N),'_','Kennfeld_Gesamtverluste','.mat'],'Kennfeld_Gesamtverluste');
%         save (['LDS/Vehicle/Para_Powertrain/Motor_efficiency/MOT_memory//', KennfeldMaschine_Verluste_gesamt_Speichern.Entwurf.Optionen.Maschinentyp,'_',num2str(KennfeldMaschine_Verluste_gesamt_Speichern.Entwurf.Bemessungswerte.P_N/1000),'_',num2str(KennfeldMaschine_Verluste_gesamt_Speichern.Entwurf.Bemessungswerte.n_N),'_',num2str(KennfeldMaschine_Verluste_gesamt_Speichern.Analyse.Optionen.n_max),'_',num2str(KennfeldMaschine_Verluste_gesamt_Speichern.Entwurf.Bemessungswerte.U_N),'_', 'Kennfeld_Gesamtverluste' ,'.mat'],'Kennfeld_Gesamtverluste');
        
        % Save map
        KennfeldMaschine_Wirkungsgrad_Wicklungsverluste = struct;
        KennfeldMaschine_Wirkungsgrad_Wicklungsverluste_Speichern = handles;
        KennfeldMaschine_Wirkungsgrad_Wicklungsverluste.n = KennfeldMaschine_Wirkungsgrad_Wicklungsverluste_Speichern.Analyse.Betriebsdaten.n_m_mesh;
        KennfeldMaschine_Wirkungsgrad_Wicklungsverluste.M = KennfeldMaschine_Wirkungsgrad_Wicklungsverluste_Speichern.Analyse.Betriebsdaten.M_max_mesh;
        KennfeldMaschine_Wirkungsgrad_Wicklungsverluste.wirkwick = KennfeldMaschine_Wirkungsgrad_Wicklungsverluste_Speichern.Analyse.Wirkungsgrad.eta_vw_mesh;
        %save (['3_Ergebnisse/',handles.Entwurf.Optionen.folder_id,'/','2_Analyse/', 'KennfeldMaschine_Wirkungsgrad_Wicklungsverluste','.mat'],'KennfeldMaschine_Wirkungsgrad_Wicklungsverluste');
        folder_id = [KennfeldMaschine_Wirkungsgrad_Wicklungsverluste_Speichern.Entwurf.Optionen.folder_id]; %file_id, '_data'];
        file_id = [KennfeldMaschine_Wirkungsgrad_Wicklungsverluste_Speichern.Entwurf.Optionen.file_id];
        
        % Convert data to correct form
        % Design sizes
        Kennfeld_Wirkungsgrad_Wicklungsverluste = struct;
        Kennfeld_Wirkungsgrad_Wicklungsverluste.type = KennfeldMaschine_Wirkungsgrad_Wicklungsverluste_Speichern.Entwurf.Optionen.Maschinentyp;
        Kennfeld_Wirkungsgrad_Wicklungsverluste.power = KennfeldMaschine_Wirkungsgrad_Wicklungsverluste_Speichern.Entwurf.Bemessungswerte.P_N/1000;
        Kennfeld_Wirkungsgrad_Wicklungsverluste.n_n = KennfeldMaschine_Wirkungsgrad_Wicklungsverluste_Speichern.Entwurf.Bemessungswerte.n_N;
        Kennfeld_Wirkungsgrad_Wicklungsverluste.n_max = KennfeldMaschine_Wirkungsgrad_Wicklungsverluste_Speichern.Analyse.Optionen.n_max;
        Kennfeld_Wirkungsgrad_Wicklungsverluste.U_n = KennfeldMaschine_Wirkungsgrad_Wicklungsverluste_Speichern.Entwurf.Bemessungswerte.U_N;
        
        % Speed axis - X-axis map
        Kennfeld_Wirkungsgrad_Wicklungsverluste.eff_n_axis = transpose(KennfeldMaschine_Wirkungsgrad_Wicklungsverluste_Speichern.Analyse.Betriebsdaten.mot.n_m_vec);

        % Torque - Y-Achse
        KennfeldMaschine_Wirkungsgrad_Wicklungsverluste_Speichern.Analyse.Betriebsdaten.mot.M_vec=transpose(fliplr(KennfeldMaschine_Wirkungsgrad_Wicklungsverluste_Speichern.Analyse.Betriebsdaten.mot.M_vec));
        Kennfeld_Wirkungsgrad_Wicklungsverluste.eff_T_axis = [(fliplr(transpose(KennfeldMaschine_Wirkungsgrad_Wicklungsverluste_Speichern.Analyse.Betriebsdaten.mot.M_vec(2:end)))) * (-1), transpose(KennfeldMaschine_Wirkungsgrad_Wicklungsverluste_Speichern.Analyse.Betriebsdaten.mot.M_vec)];
        
        % Diagram area                                                             
        % Manipulate efficiency map
        KennfeldMaschine_Wirkungsgrad_Wicklungsverluste_Speichern.Analyse.Wirkungsgrad.eta_vw_mesh=flipud(KennfeldMaschine_Wirkungsgrad_Wicklungsverluste_Speichern.Analyse.Wirkungsgrad.eta_vw_mesh); %hier Unterschied Verlustarten
        KennfeldMaschine_Wirkungsgrad_Wicklungsverluste_Speichern.Analyse.Wirkungsgrad.eta_vw_mesh(KennfeldMaschine_Wirkungsgrad_Wicklungsverluste_Speichern.Analyse.Wirkungsgrad.eta_vw_mesh<0.1) = 0.1; %Hier Unterschied Verlustarten %set small efficiency values to 0.1
        KennfeldMaschine_Wirkungsgrad_Wicklungsverluste_Speichern.Analyse.Wirkungsgrad.eta_vw_mesh(1,:) = 0; %first colomn is zero %hier Unterschied Verlustarten
        KennfeldMaschine_Wirkungsgrad_Wicklungsverluste_Speichern.Analyse.Wirkungsgrad.eta_vw_mesh(:,1) = 0.1; %first row is zero %hier Unterschied Verlustarten
        % Create one "larger" eff. map that includes acceleration and recuperation
        Kennfeld_Wirkungsgrad_Wicklungsverluste.eff_recu = 1./KennfeldMaschine_Wirkungsgrad_Wicklungsverluste_Speichern.Analyse.Wirkungsgrad.eta_vw_mesh; %recuperation efficiency %hier Unterschied Verlustarten
        Kennfeld_Wirkungsgrad_Wicklungsverluste.eff_recu(Kennfeld_Wirkungsgrad_Wicklungsverluste.eff_recu==inf)=100;
        Kennfeld_Wirkungsgrad_Wicklungsverluste.eff_recu = Kennfeld_Wirkungsgrad_Wicklungsverluste.eff_recu(2:end,:);
        Kennfeld_Wirkungsgrad_Wicklungsverluste.eff = [flipud(Kennfeld_Wirkungsgrad_Wicklungsverluste.eff_recu); KennfeldMaschine_Wirkungsgrad_Wicklungsverluste_Speichern.Analyse.Wirkungsgrad.eta_vw_mesh]; %hier Unterschied Verlustarten
    
        % Full load characteristic
        KennfeldMaschine_Wirkungsgrad_Wicklungsverluste_Speichern.Analyse.Betriebsdaten.mot.M_max_vec=transpose(KennfeldMaschine_Wirkungsgrad_Wicklungsverluste_Speichern.Analyse.Betriebsdaten.mot.M_max_vec);
        Kennfeld_Wirkungsgrad_Wicklungsverluste.T_max = transpose(KennfeldMaschine_Wirkungsgrad_Wicklungsverluste_Speichern.Analyse.Betriebsdaten.mot.M_max_vec);  
        
        % Speed axis scaled
        Kennfeld_Wirkungsgrad_Wicklungsverluste.T_max_n_axis = transpose(KennfeldMaschine_Wirkungsgrad_Wicklungsverluste_Speichern.Analyse.Betriebsdaten.mot.M_vec);                                                                          %noch auf 2000 interpolieren?
        
%         save (['3_Ergebnisse/',handles.Entwurf.Optionen.folder_id,'/','2_Analyse/', KennfeldMaschine_Wirkungsgrad_Wicklungsverluste_Speichern.Entwurf.Optionen.Maschinentyp,'_',num2str(KennfeldMaschine_Wirkungsgrad_Wicklungsverluste_Speichern.Entwurf.Bemessungswerte.P_N/1000),'_',num2str(KennfeldMaschine_Wirkungsgrad_Wicklungsverluste_Speichern.Entwurf.Bemessungswerte.n_N),'_',num2str(KennfeldMaschine_Wirkungsgrad_Wicklungsverluste_Speichern.Analyse.Optionen.n_max),'_',num2str(KennfeldMaschine_Wirkungsgrad_Wicklungsverluste_Speichern.Entwurf.Bemessungswerte.U_N),'_','Kennfeld_Wirkungsgrad_Wicklungsverluste','.mat'],'Kennfeld_Wirkungsgrad_Wicklungsverluste');
%         save (['LDS/Vehicle/Para_Powertrain/Motor_efficiency/MOT_memory//', KennfeldMaschine_Wirkungsgrad_Wicklungsverluste_Speichern.Entwurf.Optionen.Maschinentyp,'_',num2str(KennfeldMaschine_Wirkungsgrad_Wicklungsverluste_Speichern.Entwurf.Bemessungswerte.P_N/1000),'_',num2str(KennfeldMaschine_Wirkungsgrad_Wicklungsverluste_Speichern.Entwurf.Bemessungswerte.n_N),'_',num2str(KennfeldMaschine_Wirkungsgrad_Wicklungsverluste_Speichern.Analyse.Optionen.n_max),'_',num2str(KennfeldMaschine_Wirkungsgrad_Wicklungsverluste_Speichern.Entwurf.Bemessungswerte.U_N),'_', 'Kennfeld_Wirkungsgrad_Wicklungsverluste' ,'.mat'],'Kennfeld_Wirkungsgrad_Wicklungsverluste');

        % Save map
        KennfeldMaschine_Wicklungsverluste = struct;
        KennfeldMaschine_Wicklungsverluste_Speichern = handles;
        KennfeldMaschine_Wicklungsverluste.n = KennfeldMaschine_Wicklungsverluste_Speichern.Analyse.Betriebsdaten.n_m_mesh;
        KennfeldMaschine_Wicklungsverluste.M = KennfeldMaschine_Wicklungsverluste_Speichern.Analyse.Betriebsdaten.M_max_mesh;
        KennfeldMaschine_Wicklungsverluste.verlwick = KennfeldMaschine_Wicklungsverluste_Speichern.Analyse.Verluste.P_vw_mesh;
        %save (['3_Ergebnisse/',handles.Entwurf.Optionen.folder_id,'/','2_Analyse/', 'KennfeldMaschine_Wicklungsverluste','.mat'],'KennfeldMaschine_Wicklungsverluste');
        folder_id = [KennfeldMaschine_Wicklungsverluste_Speichern.Entwurf.Optionen.folder_id]; %file_id, '_data'];
        file_id = [KennfeldMaschine_Wicklungsverluste_Speichern.Entwurf.Optionen.file_id];
        
        % Convert data to correct form
        % Design sizes
        Kennfeld_Wicklungsverluste = struct;
        Kennfeld_Wicklungsverluste.type = KennfeldMaschine_Wicklungsverluste_Speichern.Entwurf.Optionen.Maschinentyp;
        Kennfeld_Wicklungsverluste.power = KennfeldMaschine_Wicklungsverluste_Speichern.Entwurf.Bemessungswerte.P_N/1000;
        Kennfeld_Wicklungsverluste.n_n = KennfeldMaschine_Wicklungsverluste_Speichern.Entwurf.Bemessungswerte.n_N;
        Kennfeld_Wicklungsverluste.n_max = KennfeldMaschine_Wicklungsverluste_Speichern.Analyse.Optionen.n_max;
        Kennfeld_Wicklungsverluste.U_n = KennfeldMaschine_Wicklungsverluste_Speichern.Entwurf.Bemessungswerte.U_N;
        
        % Speed axis - X-axis map
        Kennfeld_Wicklungsverluste.eff_n_axis = transpose(KennfeldMaschine_Wicklungsverluste_Speichern.Analyse.Betriebsdaten.mot.n_m_vec);

        % Torque - Y-Achse
        KennfeldMaschine_Wicklungsverluste_Speichern.Analyse.Betriebsdaten.mot.M_vec=transpose(fliplr(KennfeldMaschine_Wicklungsverluste_Speichern.Analyse.Betriebsdaten.mot.M_vec));
        Kennfeld_Wicklungsverluste.eff_T_axis = [(fliplr(transpose(KennfeldMaschine_Wicklungsverluste_Speichern.Analyse.Betriebsdaten.mot.M_vec(2:end)))) * (-1), transpose(KennfeldMaschine_Wicklungsverluste_Speichern.Analyse.Betriebsdaten.mot.M_vec)];
        
        % Diagrammbereich                                                             
        % Manipulieren der Effizienzkarte
        KennfeldMaschine_Wicklungsverluste_Speichern.Analyse.Verluste.P_vw_mesh=flipud(KennfeldMaschine_Wicklungsverluste_Speichern.Analyse.Verluste.P_vw_mesh); %hier Unterschied Verlustarten
        KennfeldMaschine_Wicklungsverluste_Speichern.Analyse.Verluste.P_vw_mesh(KennfeldMaschine_Wicklungsverluste_Speichern.Analyse.Verluste.P_vw_mesh<0.1) = 0.1; %Hier Unterschied Verlustarten %set small efficiency values to 0.1
        KennfeldMaschine_Wicklungsverluste_Speichern.Analyse.Verluste.P_vw_mesh(1,:) = 0; %first colomn is zero %hier Unterschied Verlustarten
        KennfeldMaschine_Wicklungsverluste_Speichern.Analyse.Verluste.P_vw_mesh(:,1) = 0.1; %first row is zero %hier Unterschied Verlustarten
        
        % Create one "larger" eff. map that includes acceleration and recuperation
        Kennfeld_Wicklungsverluste.eff_recu = 1./KennfeldMaschine_Wicklungsverluste_Speichern.Analyse.Verluste.P_vw_mesh; %recuperation efficiency %hier Unterschied Verlustarten
        Kennfeld_Wicklungsverluste.eff_recu(Kennfeld_Wicklungsverluste.eff_recu==inf)=100;
        Kennfeld_Wicklungsverluste.eff_recu = Kennfeld_Wicklungsverluste.eff_recu(2:end,:);
        Kennfeld_Wicklungsverluste.eff = [flipud(Kennfeld_Wicklungsverluste.eff_recu); KennfeldMaschine_Wicklungsverluste_Speichern.Analyse.Verluste.P_vw_mesh]; %hier Unterschied Verlustarten
    
        % Full load characteristic
        KennfeldMaschine_Wicklungsverluste_Speichern.Analyse.Betriebsdaten.mot.M_max_vec=transpose(KennfeldMaschine_Wicklungsverluste_Speichern.Analyse.Betriebsdaten.mot.M_max_vec);
        Kennfeld_Wicklungsverluste.T_max = transpose(KennfeldMaschine_Wicklungsverluste_Speichern.Analyse.Betriebsdaten.mot.M_max_vec);  
        
        % Speed axis scaled
        Kennfeld_Wicklungsverluste.T_max_n_axis = transpose(KennfeldMaschine_Wicklungsverluste_Speichern.Analyse.Betriebsdaten.mot.M_vec);                                                          %noch auf 2000 interpolieren?
        
%         save (['3_Ergebnisse/',handles.Entwurf.Optionen.folder_id,'/', '2_Analyse/', KennfeldMaschine_Wicklungsverluste_Speichern.Entwurf.Optionen.Maschinentyp,'_',num2str(KennfeldMaschine_Wicklungsverluste_Speichern.Entwurf.Bemessungswerte.P_N/1000),'_',num2str(KennfeldMaschine_Wicklungsverluste_Speichern.Entwurf.Bemessungswerte.n_N),'_',num2str(KennfeldMaschine_Wicklungsverluste_Speichern.Analyse.Optionen.n_max),'_',num2str(KennfeldMaschine_Wicklungsverluste_Speichern.Entwurf.Bemessungswerte.U_N),'_','Kennfeld_Wicklungsverluste','.mat'],'Kennfeld_Wicklungsverluste');
%         save (['LDS/Vehicle/Para_Powertrain/Motor_efficiency/MOT_memory//', KennfeldMaschine_Wicklungsverluste_Speichern.Entwurf.Optionen.Maschinentyp,'_',num2str(KennfeldMaschine_Wicklungsverluste_Speichern.Entwurf.Bemessungswerte.P_N/1000),'_',num2str(KennfeldMaschine_Wicklungsverluste_Speichern.Entwurf.Bemessungswerte.n_N),'_',num2str(KennfeldMaschine_Wicklungsverluste_Speichern.Analyse.Optionen.n_max),'_',num2str(KennfeldMaschine_Wicklungsverluste_Speichern.Entwurf.Bemessungswerte.U_N),'_', 'Kennfeld_Wicklungsverluste' ,'.mat'],'Kennfeld_Wicklungsverluste');

        % Save map
        KennfeldMaschine_Wirkungsgrad_Ummagnetisierungsverluste = struct;
        KennfeldMaschine_WirkUmmagnetisierungsverluste_Speichern = handles;
        KennfeldMaschine_Wirkungsgrad_Ummagnetisierungsverluste.n = KennfeldMaschine_WirkUmmagnetisierungsverluste_Speichern.Analyse.Betriebsdaten.n_m_mesh;
        KennfeldMaschine_Wirkungsgrad_Ummagnetisierungsverluste.M = KennfeldMaschine_WirkUmmagnetisierungsverluste_Speichern.Analyse.Betriebsdaten.M_max_mesh;
        KennfeldMaschine_Wirkungsgrad_Ummagnetisierungsverluste.wirkummagn = KennfeldMaschine_WirkUmmagnetisierungsverluste_Speichern.Analyse.Wirkungsgrad.eta_fe_mesh;
        %save (['3_Ergebnisse/',handles.Entwurf.Optionen.folder_id,'/','2_Analyse/','KennfeldMaschine_Wirkungsgrad_Ummagnetisierungsverluste','.mat'],'KennfeldMaschine_Wirkungsgrad_Ummagnetisierungsverluste');
        folder_id = [KennfeldMaschine_WirkUmmagnetisierungsverluste_Speichern.Entwurf.Optionen.folder_id]; %file_id, '_data'];
        file_id = [KennfeldMaschine_WirkUmmagnetisierungsverluste_Speichern.Entwurf.Optionen.file_id];
        
        % Convert data into the correct form
        % Design sizes
        Kennfeld_Wirkungsgrad_Ummagnetisierungsverluste = struct;
        Kennfeld_Wirkungsgrad_Ummagnetisierungsverluste.type = KennfeldMaschine_WirkUmmagnetisierungsverluste_Speichern.Entwurf.Optionen.Maschinentyp;
        Kennfeld_Wirkungsgrad_Ummagnetisierungsverluste.power = KennfeldMaschine_WirkUmmagnetisierungsverluste_Speichern.Entwurf.Bemessungswerte.P_N/1000;
        Kennfeld_Wirkungsgrad_Ummagnetisierungsverluste.n_n = KennfeldMaschine_WirkUmmagnetisierungsverluste_Speichern.Entwurf.Bemessungswerte.n_N;
        Kennfeld_Wirkungsgrad_Ummagnetisierungsverluste.n_max = KennfeldMaschine_WirkUmmagnetisierungsverluste_Speichern.Analyse.Optionen.n_max;
        Kennfeld_Wirkungsgrad_Ummagnetisierungsverluste.U_n = KennfeldMaschine_WirkUmmagnetisierungsverluste_Speichern.Entwurf.Bemessungswerte.U_N;
        
        % Speed axis - X-axis map
        Kennfeld_Wirkungsgrad_Ummagnetisierungsverluste.eff_n_axis = transpose(KennfeldMaschine_WirkUmmagnetisierungsverluste_Speichern.Analyse.Betriebsdaten.mot.n_m_vec);

        % Torque - Y-Achse
        KennfeldMaschine_WirkUmmagnetisierungsverluste_Speichern.Analyse.Betriebsdaten.mot.M_vec=transpose(fliplr(KennfeldMaschine_WirkUmmagnetisierungsverluste_Speichern.Analyse.Betriebsdaten.mot.M_vec));
        Kennfeld_Wirkungsgrad_Ummagnetisierungsverluste.eff_T_axis = [(fliplr(transpose(KennfeldMaschine_WirkUmmagnetisierungsverluste_Speichern.Analyse.Betriebsdaten.mot.M_vec(2:end)))) * (-1), transpose(KennfeldMaschine_WirkUmmagnetisierungsverluste_Speichern.Analyse.Betriebsdaten.mot.M_vec)];
        
        % Diagram                                                             
        % Manipulate efficiency map
        KennfeldMaschine_WirkUmmagnetisierungsverluste_Speichern.Analyse.Wirkungsgrad.eta_fe_mesh=flipud(KennfeldMaschine_WirkUmmagnetisierungsverluste_Speichern.Analyse.Wirkungsgrad.eta_fe_mesh); %hier Unterschied Verlustarten
        KennfeldMaschine_WirkUmmagnetisierungsverluste_Speichern.Analyse.Wirkungsgrad.eta_fe_mesh(KennfeldMaschine_WirkUmmagnetisierungsverluste_Speichern.Analyse.Wirkungsgrad.eta_fe_mesh<0.1) = 0.1; %Hier Unterschied Verlustarten %set small efficiency values to 0.1
        KennfeldMaschine_WirkUmmagnetisierungsverluste_Speichern.Analyse.Wirkungsgrad.eta_fe_mesh(1,:) = 0; %first colomn is zero %hier Unterschied Verlustarten
        KennfeldMaschine_WirkUmmagnetisierungsverluste_Speichern.Analyse.Wirkungsgrad.eta_fe_mesh(:,1) = 0.1; %first row is zero %hier Unterschied Verlustarten
        
        % Create one "larger" eff. map that includes acceleration and recuperation
        Kennfeld_Wirkungsgrad_Ummagnetisierungsverluste.eff_recu = 1./KennfeldMaschine_WirkUmmagnetisierungsverluste_Speichern.Analyse.Wirkungsgrad.eta_fe_mesh; %recuperation efficiency %hier Unterschied Verlustarten
        Kennfeld_Wirkungsgrad_Ummagnetisierungsverluste.eff_recu(Kennfeld_Wirkungsgrad_Ummagnetisierungsverluste.eff_recu==inf)=100;
        Kennfeld_Wirkungsgrad_Ummagnetisierungsverluste.eff_recu = Kennfeld_Wirkungsgrad_Ummagnetisierungsverluste.eff_recu(2:end,:);
        Kennfeld_Wirkungsgrad_Ummagnetisierungsverluste.eff = [flipud(Kennfeld_Wirkungsgrad_Ummagnetisierungsverluste.eff_recu); KennfeldMaschine_WirkUmmagnetisierungsverluste_Speichern.Analyse.Wirkungsgrad.eta_fe_mesh]; %hier Unterschied Verlustarten
    
        % Full load characteristic
        KennfeldMaschine_WirkUmmagnetisierungsverluste_Speichern.Analyse.Betriebsdaten.mot.M_max_vec=transpose(KennfeldMaschine_WirkUmmagnetisierungsverluste_Speichern.Analyse.Betriebsdaten.mot.M_max_vec);
        Kennfeld_Wirkungsgrad_Ummagnetisierungsverluste.T_max = transpose(KennfeldMaschine_WirkUmmagnetisierungsverluste_Speichern.Analyse.Betriebsdaten.mot.M_max_vec);  
        
        % Speed axis scaled
        Kennfeld_Wirkungsgrad_Ummagnetisierungsverluste.T_max_n_axis = transpose(KennfeldMaschine_WirkUmmagnetisierungsverluste_Speichern.Analyse.Betriebsdaten.mot.M_vec);                                                          %noch auf 2000 interpolieren?
        
%         save (['3_Ergebnisse/',handles.Entwurf.Optionen.folder_id,'/','2_Analyse/',KennfeldMaschine_WirkUmmagnetisierungsverluste_Speichern.Entwurf.Optionen.Maschinentyp,'_',num2str(KennfeldMaschine_WirkUmmagnetisierungsverluste_Speichern.Entwurf.Bemessungswerte.P_N/1000),'_',num2str(KennfeldMaschine_WirkUmmagnetisierungsverluste_Speichern.Entwurf.Bemessungswerte.n_N),'_',num2str(KennfeldMaschine_WirkUmmagnetisierungsverluste_Speichern.Analyse.Optionen.n_max),'_',num2str(KennfeldMaschine_WirkUmmagnetisierungsverluste_Speichern.Entwurf.Bemessungswerte.U_N),'_','Kennfeld_Wirkungsgrad_Ummagnetisierungsverluste','.mat'],'Kennfeld_Wirkungsgrad_Ummagnetisierungsverluste');
%         save (['LDS/Vehicle/Para_Powertrain/Motor_efficiency/MOT_memory//', KennfeldMaschine_WirkUmmagnetisierungsverluste_Speichern.Entwurf.Optionen.Maschinentyp,'_',num2str(KennfeldMaschine_WirkUmmagnetisierungsverluste_Speichern.Entwurf.Bemessungswerte.P_N/1000),'_',num2str(KennfeldMaschine_WirkUmmagnetisierungsverluste_Speichern.Entwurf.Bemessungswerte.n_N),'_',num2str(KennfeldMaschine_WirkUmmagnetisierungsverluste_Speichern.Analyse.Optionen.n_max),'_',num2str(KennfeldMaschine_WirkUmmagnetisierungsverluste_Speichern.Entwurf.Bemessungswerte.U_N),'_', 'Kennfeld_Wirkungsgrad_Ummagnetisierungsverluste' ,'.mat'],'Kennfeld_Wirkungsgrad_Ummagnetisierungsverluste');
 
        % Save map
        KennfeldMaschine_Ummagnetisierungsverluste = struct;
        KennfeldMaschine_Ummagnetisierungsverluste_Speichern = handles;
        KennfeldMaschine_Ummagnetisierungsverluste.n = KennfeldMaschine_Ummagnetisierungsverluste_Speichern.Analyse.Betriebsdaten.n_m_mesh;
        KennfeldMaschine_Ummagnetisierungsverluste.M = KennfeldMaschine_Ummagnetisierungsverluste_Speichern.Analyse.Betriebsdaten.M_max_mesh;
        KennfeldMaschine_Ummagnetisierungsverluste.verlummagn = KennfeldMaschine_Ummagnetisierungsverluste_Speichern.Analyse.Verluste.P_vu_mesh;
        %save (['3_Ergebnisse/', handles.Entwurf.Optionen.folder_id,'/','2_Analyse/', 'KennfeldMaschine_Ummagnetisierungsverluste','.mat'],'KennfeldMaschine_Ummagnetisierungsverluste');
        folder_id = [KennfeldMaschine_Ummagnetisierungsverluste_Speichern.Entwurf.Optionen.folder_id]; %file_id, '_data'];
        file_id = [KennfeldMaschine_Ummagnetisierungsverluste_Speichern.Entwurf.Optionen.file_id];
        
        % Convert data to correct form
        % Design sizes
        Kennfeld_Ummagnetisierungsverluste = struct;
        Kennfeld_Ummagnetisierungsverluste.type = KennfeldMaschine_Ummagnetisierungsverluste_Speichern.Entwurf.Optionen.Maschinentyp;
        Kennfeld_Ummagnetisierungsverluste.power = KennfeldMaschine_Ummagnetisierungsverluste_Speichern.Entwurf.Bemessungswerte.P_N/1000;
        Kennfeld_Ummagnetisierungsverluste.n_n = KennfeldMaschine_Ummagnetisierungsverluste_Speichern.Entwurf.Bemessungswerte.n_N;
        Kennfeld_Ummagnetisierungsverluste.n_max = KennfeldMaschine_Ummagnetisierungsverluste_Speichern.Analyse.Optionen.n_max;
        Kennfeld_Ummagnetisierungsverluste.U_n = KennfeldMaschine_Ummagnetisierungsverluste_Speichern.Entwurf.Bemessungswerte.U_N;
        
        % Speed axis - X-axis map
        Kennfeld_Ummagnetisierungsverluste.eff_n_axis = transpose(KennfeldMaschine_Ummagnetisierungsverluste_Speichern.Analyse.Betriebsdaten.mot.n_m_vec);

        % Torque - Y-axis
        KennfeldMaschine_Ummagnetisierungsverluste_Speichern.Analyse.Betriebsdaten.mot.M_vec=transpose(fliplr(KennfeldMaschine_Ummagnetisierungsverluste_Speichern.Analyse.Betriebsdaten.mot.M_vec));
        Kennfeld_Ummagnetisierungsverluste.eff_T_axis = [(fliplr(transpose(KennfeldMaschine_Ummagnetisierungsverluste_Speichern.Analyse.Betriebsdaten.mot.M_vec(2:end)))) * (-1), transpose(KennfeldMaschine_Ummagnetisierungsverluste_Speichern.Analyse.Betriebsdaten.mot.M_vec)];
        
        % Diagram                                                             
        % Manipulate efficiency map
        KennfeldMaschine_Ummagnetisierungsverluste_Speichern.Analyse.Verluste.P_vu_mesh=flipud(KennfeldMaschine_Ummagnetisierungsverluste_Speichern.Analyse.Verluste.P_vu_mesh); %hier Unterschied Verlustarten
        KennfeldMaschine_Ummagnetisierungsverluste_Speichern.Analyse.Verluste.P_vu_mesh(KennfeldMaschine_Ummagnetisierungsverluste_Speichern.Analyse.Verluste.P_vu_mesh<0.1) = 0.1; %Hier Unterschied Verlustarten %set small efficiency values to 0.1
        KennfeldMaschine_Ummagnetisierungsverluste_Speichern.Analyse.Verluste.P_vu_mesh(1,:) = 0; %first colomn is zero %hier Unterschied Verlustarten
        KennfeldMaschine_Ummagnetisierungsverluste_Speichern.Analyse.Verluste.P_vu_mesh(:,1) = 0.1; %first row is zero %hier Unterschied Verlustarten
        
        % Create one "larger" eff. map that includes acceleration and recuperation
        Kennfeld_Ummagnetisierungsverluste.eff_recu = 1./KennfeldMaschine_Ummagnetisierungsverluste_Speichern.Analyse.Verluste.P_vu_mesh; %recuperation efficiency %hier Unterschied Verlustarten
        Kennfeld_Ummagnetisierungsverluste.eff_recu(Kennfeld_Ummagnetisierungsverluste.eff_recu==inf)=100;
        Kennfeld_Ummagnetisierungsverluste.eff_recu = Kennfeld_Ummagnetisierungsverluste.eff_recu(2:end,:);
        Kennfeld_Ummagnetisierungsverluste.eff = [flipud(Kennfeld_Ummagnetisierungsverluste.eff_recu); KennfeldMaschine_Ummagnetisierungsverluste_Speichern.Analyse.Verluste.P_vu_mesh]; %hier Unterschied Verlustarten
    
        % Full load characteristic
        KennfeldMaschine_Ummagnetisierungsverluste_Speichern.Analyse.Betriebsdaten.mot.M_max_vec=transpose(KennfeldMaschine_Ummagnetisierungsverluste_Speichern.Analyse.Betriebsdaten.mot.M_max_vec);
        Kennfeld_Ummagnetisierungsverluste.T_max = transpose(KennfeldMaschine_Ummagnetisierungsverluste_Speichern.Analyse.Betriebsdaten.mot.M_max_vec);  
        
        % Speed axis scaled
        Kennfeld_Ummagnetisierungsverluste.T_max_n_axis = transpose(KennfeldMaschine_Ummagnetisierungsverluste_Speichern.Analyse.Betriebsdaten.mot.M_vec);                                                          %noch auf 2000 interpolieren?

%         save (['3_Ergebnisse/',handles.Entwurf.Optionen.folder_id,'/', '2_Analyse/', KennfeldMaschine_Ummagnetisierungsverluste_Speichern.Entwurf.Optionen.Maschinentyp,'_',num2str(KennfeldMaschine_Ummagnetisierungsverluste_Speichern.Entwurf.Bemessungswerte.P_N/1000),'_',num2str(KennfeldMaschine_Ummagnetisierungsverluste_Speichern.Entwurf.Bemessungswerte.n_N),'_',num2str(KennfeldMaschine_Ummagnetisierungsverluste_Speichern.Analyse.Optionen.n_max),'_',num2str(KennfeldMaschine_Ummagnetisierungsverluste_Speichern.Entwurf.Bemessungswerte.U_N),'_','Kennfeld_Ummagnetisierungsverluste','.mat'],'Kennfeld_Ummagnetisierungsverluste');
%         save (['LDS/Vehicle/Para_Powertrain/Motor_efficiency/MOT_memory//', KennfeldMaschine_Ummagnetisierungsverluste_Speichern.Entwurf.Optionen.Maschinentyp,'_',num2str(KennfeldMaschine_Ummagnetisierungsverluste_Speichern.Entwurf.Bemessungswerte.P_N/1000),'_',num2str(KennfeldMaschine_Ummagnetisierungsverluste_Speichern.Entwurf.Bemessungswerte.n_N),'_',num2str(KennfeldMaschine_Ummagnetisierungsverluste_Speichern.Analyse.Optionen.n_max),'_',num2str(KennfeldMaschine_Ummagnetisierungsverluste_Speichern.Entwurf.Bemessungswerte.U_N),'_', 'Kennfeld_Ummagnetisierungsverluste' ,'.mat'],'Kennfeld_Ummagnetisierungsverluste');

        % Save map
        KennfeldMaschine_WirkmechVerluste = struct;
        KennfeldMaschine_WirkmechVerluste_Speichern = handles;
        KennfeldMaschine_WirkmechVerluste.n = KennfeldMaschine_WirkmechVerluste_Speichern.Analyse.Betriebsdaten.n_m_mesh;
        KennfeldMaschine_WirkmechVerluste.M = KennfeldMaschine_WirkmechVerluste_Speichern.Analyse.Betriebsdaten.M_max_mesh;
        KennfeldMaschine_WirkmechVerluste.wirkmechverl = KennfeldMaschine_WirkmechVerluste_Speichern.Analyse.Wirkungsgrad.eta_vme_mesh;
        %save (['3_Ergebnisse/',handles.Entwurf.Optionen.folder_id,'/','2_Analyse/', 'KennfeldMaschine_WirkmechVerluste','.mat'],'KennfeldMaschine_WirkmechVerluste');
        folder_id = [KennfeldMaschine_WirkmechVerluste_Speichern.Entwurf.Optionen.folder_id]; %file_id, '_data'];
        file_id = [KennfeldMaschine_WirkmechVerluste_Speichern.Entwurf.Optionen.file_id];
        
        % Convert data to correct form
        % Design sizes
        Kennfeld_WirkmechVerluste = struct;
        Kennfeld_WirkmechVerluste.type = KennfeldMaschine_WirkmechVerluste_Speichern.Entwurf.Optionen.Maschinentyp;
        Kennfeld_WirkmechVerluste.power = KennfeldMaschine_WirkmechVerluste_Speichern.Entwurf.Bemessungswerte.P_N/1000;
        Kennfeld_WirkmechVerluste.n_n = KennfeldMaschine_WirkmechVerluste_Speichern.Entwurf.Bemessungswerte.n_N;
        Kennfeld_WirkmechVerluste.n_max = KennfeldMaschine_WirkmechVerluste_Speichern.Analyse.Optionen.n_max;
        Kennfeld_WirkmechVerluste.U_n = KennfeldMaschine_WirkmechVerluste_Speichern.Entwurf.Bemessungswerte.U_N;
        
        % Speed axis - X-axis map
        Kennfeld_WirkmechVerluste.eff_n_axis = transpose(KennfeldMaschine_WirkmechVerluste_Speichern.Analyse.Betriebsdaten.mot.n_m_vec);

        % Torque - Y-axis
        KennfeldMaschine_WirkmechVerluste_Speichern.Analyse.Betriebsdaten.mot.M_vec=transpose(fliplr(KennfeldMaschine_WirkmechVerluste_Speichern.Analyse.Betriebsdaten.mot.M_vec));
        Kennfeld_WirkmechVerluste.eff_T_axis = [(fliplr(transpose(KennfeldMaschine_WirkmechVerluste_Speichern.Analyse.Betriebsdaten.mot.M_vec(2:end)))) * (-1), transpose(KennfeldMaschine_WirkmechVerluste_Speichern.Analyse.Betriebsdaten.mot.M_vec)];
        
        % Diagram                                                            
        % Manipulate efficiency map
        KennfeldMaschine_WirkmechVerluste_Speichern.Analyse.Wirkungsgrad.eta_vme_mesh=flipud(KennfeldMaschine_WirkmechVerluste_Speichern.Analyse.Wirkungsgrad.eta_vme_mesh); %hier Unterschied Verlustarten
        KennfeldMaschine_WirkmechVerluste_Speichern.Analyse.Wirkungsgrad.eta_vme_mesh(KennfeldMaschine_WirkmechVerluste_Speichern.Analyse.Wirkungsgrad.eta_vme_mesh<0.1) = 0.1; %Hier Unterschied Verlustarten %set small efficiency values to 0.1
        KennfeldMaschine_WirkmechVerluste_Speichern.Analyse.Wirkungsgrad.eta_vme_mesh(1,:) = 0; %first colomn is zero %hier Unterschied Verlustarten
        KennfeldMaschine_WirkmechVerluste_Speichern.Analyse.Wirkungsgrad.eta_vme_mesh(:,1) = 0.1; %first row is zero %hier Unterschied Verlustarten
        
        % Create one "larger" eff. map that includes acceleration and recuperation
        Kennfeld_WirkmechVerluste.eff_recu = 1./KennfeldMaschine_WirkmechVerluste_Speichern.Analyse.Wirkungsgrad.eta_vme_mesh; %recuperation efficiency %hier Unterschied Verlustarten
        Kennfeld_WirkmechVerluste.eff_recu(Kennfeld_WirkmechVerluste.eff_recu==inf)=100;
        Kennfeld_WirkmechVerluste.eff_recu = Kennfeld_WirkmechVerluste.eff_recu(2:end,:);
        Kennfeld_WirkmechVerluste.eff = [flipud(Kennfeld_WirkmechVerluste.eff_recu); KennfeldMaschine_WirkmechVerluste_Speichern.Analyse.Wirkungsgrad.eta_vme_mesh]; %hier Unterschied Verlustarten
    
        % Full load characteristic

        KennfeldMaschine_WirkmechVerluste_Speichern.Analyse.Betriebsdaten.mot.M_max_vec=transpose(KennfeldMaschine_WirkmechVerluste_Speichern.Analyse.Betriebsdaten.mot.M_max_vec);
        Kennfeld_WirkmechVerluste.T_max = transpose(KennfeldMaschine_WirkmechVerluste_Speichern.Analyse.Betriebsdaten.mot.M_max_vec);  
        
        % Speed axis scaled
        Kennfeld_WirkmechVerluste.T_max_n_axis = transpose(KennfeldMaschine_WirkmechVerluste_Speichern.Analyse.Betriebsdaten.mot.M_vec);                                                          %noch auf 2000 interpolieren?
        
%         save (['3_Ergebnisse/',handles.Entwurf.Optionen.folder_id,'/', '2_Analyse/', KennfeldMaschine_WirkmechVerluste_Speichern.Entwurf.Optionen.Maschinentyp,'_',num2str(KennfeldMaschine_WirkmechVerluste_Speichern.Entwurf.Bemessungswerte.P_N/1000),'_',num2str(KennfeldMaschine_WirkmechVerluste_Speichern.Entwurf.Bemessungswerte.n_N),'_',num2str(KennfeldMaschine_WirkmechVerluste_Speichern.Analyse.Optionen.n_max),'_',num2str(KennfeldMaschine_WirkmechVerluste_Speichern.Entwurf.Bemessungswerte.U_N),'_','Kennfeld_WirkmechVerluste','.mat'],'Kennfeld_WirkmechVerluste');
%         save (['LDS/Vehicle/Para_Powertrain/Motor_efficiency/MOT_memory//', KennfeldMaschine_WirkmechVerluste_Speichern.Entwurf.Optionen.Maschinentyp,'_',num2str(KennfeldMaschine_WirkmechVerluste_Speichern.Entwurf.Bemessungswerte.P_N/1000),'_',num2str(KennfeldMaschine_WirkmechVerluste_Speichern.Entwurf.Bemessungswerte.n_N),'_',num2str(KennfeldMaschine_WirkmechVerluste_Speichern.Analyse.Optionen.n_max),'_',num2str(KennfeldMaschine_WirkmechVerluste_Speichern.Entwurf.Bemessungswerte.U_N),'_', 'Kennfeld_WirkmechVerluste' ,'.mat'],'Kennfeld_WirkmechVerluste');

        % Save map
        KennfeldMaschine_mechVerluste = struct;
        KennfeldMaschine_mechVerluste_Speichern = handles;
        KennfeldMaschine_mechVerluste.n = KennfeldMaschine_mechVerluste_Speichern.Analyse.Betriebsdaten.n_m_mesh;
        KennfeldMaschine_mechVerluste.M = KennfeldMaschine_mechVerluste_Speichern.Analyse.Betriebsdaten.M_max_mesh;
        KennfeldMaschine_mechVerluste.mechverl = KennfeldMaschine_mechVerluste_Speichern.Analyse.Verluste.P_vme_mesh;
        %save (['3_Ergebnisse/', handles.Entwurf.Optionen.folder_id,'/','2_Analyse/','KennfeldMaschine_mechVerluste','.mat'],'KennfeldMaschine_mechVerluste');
        folder_id = [KennfeldMaschine_mechVerluste_Speichern.Entwurf.Optionen.folder_id]; %file_id, '_data'];
        file_id = [KennfeldMaschine_mechVerluste_Speichern.Entwurf.Optionen.file_id];
        
        % Convert data to correct form
        % Design sizes
        Kennfeld_mechVerluste = struct;
        Kennfeld_mechVerluste.type = KennfeldMaschine_mechVerluste_Speichern.Entwurf.Optionen.Maschinentyp;
        Kennfeld_mechVerluste.power = KennfeldMaschine_mechVerluste_Speichern.Entwurf.Bemessungswerte.P_N/1000;
        Kennfeld_mechVerluste.n_n = KennfeldMaschine_mechVerluste_Speichern.Entwurf.Bemessungswerte.n_N;
        Kennfeld_mechVerluste.n_max = KennfeldMaschine_mechVerluste_Speichern.Analyse.Optionen.n_max;
        Kennfeld_mechVerluste.U_n = KennfeldMaschine_mechVerluste_Speichern.Entwurf.Bemessungswerte.U_N;
        
        % Speed axis - X-axis map
        Kennfeld_mechVerluste.eff_n_axis = transpose(KennfeldMaschine_mechVerluste_Speichern.Analyse.Betriebsdaten.mot.n_m_vec);

        % Torque - Y-axis
        KennfeldMaschine_mechVerluste_Speichern.Analyse.Betriebsdaten.mot.M_vec=transpose(fliplr(KennfeldMaschine_mechVerluste_Speichern.Analyse.Betriebsdaten.mot.M_vec));
        Kennfeld_mechVerluste.eff_T_axis = [(fliplr(transpose(KennfeldMaschine_mechVerluste_Speichern.Analyse.Betriebsdaten.mot.M_vec(2:end)))) * (-1), transpose(KennfeldMaschine_mechVerluste_Speichern.Analyse.Betriebsdaten.mot.M_vec)];
        
        % Diagram                                                            
        % Manipulate efficiency map
        KennfeldMaschine_mechVerluste_Speichern.Analyse.Verluste.P_vme_mesh=flipud(KennfeldMaschine_mechVerluste_Speichern.Analyse.Verluste.P_vme_mesh); %hier Unterschied Verlustarten
        KennfeldMaschine_mechVerluste_Speichern.Analyse.Verluste.P_vme_mesh(KennfeldMaschine_mechVerluste_Speichern.Analyse.Verluste.P_vme_mesh<0.1) = 0.1; %Hier Unterschied Verlustarten %set small efficiency values to 0.1
        KennfeldMaschine_mechVerluste_Speichern.Analyse.Verluste.P_vme_mesh(1,:) = 0; %first colomn is zero %hier Unterschied Verlustarten
        KennfeldMaschine_mechVerluste_Speichern.Analyse.Verluste.P_vme_mesh(:,1) = 0.1; %first row is zero %hier Unterschied Verlustarten
        
        % Create one "larger" eff. map that includes acceleration and recuperation
        Kennfeld_mechVerluste.eff_recu = 1./KennfeldMaschine_mechVerluste_Speichern.Analyse.Verluste.P_vme_mesh; %recuperation efficiency %hier Unterschied Verlustarten
        Kennfeld_mechVerluste.eff_recu(Kennfeld_mechVerluste.eff_recu==inf)=100;
        Kennfeld_mechVerluste.eff_recu = Kennfeld_mechVerluste.eff_recu(2:end,:);
        Kennfeld_mechVerluste.eff = [flipud(Kennfeld_mechVerluste.eff_recu); KennfeldMaschine_mechVerluste_Speichern.Analyse.Verluste.P_vme_mesh]; %hier Unterschied Verlustarten
    
        % Full load characteristic
        KennfeldMaschine_mechVerluste_Speichern.Analyse.Betriebsdaten.mot.M_max_vec=transpose(KennfeldMaschine_mechVerluste_Speichern.Analyse.Betriebsdaten.mot.M_max_vec);
        Kennfeld_mechVerluste.T_max = transpose(KennfeldMaschine_mechVerluste_Speichern.Analyse.Betriebsdaten.mot.M_max_vec);  
        
        % Speed axis scaled
        Kennfeld_mechVerluste.T_max_n_axis = transpose(KennfeldMaschine_mechVerluste_Speichern.Analyse.Betriebsdaten.mot.M_vec);                                                          %noch auf 2000 interpolieren?

%         save (['3_Ergebnisse/',handles.Entwurf.Optionen.folder_id,'/', '2_Analyse/', KennfeldMaschine_mechVerluste_Speichern.Entwurf.Optionen.Maschinentyp,'_',num2str(KennfeldMaschine_mechVerluste_Speichern.Entwurf.Bemessungswerte.P_N/1000),'_',num2str(KennfeldMaschine_mechVerluste_Speichern.Entwurf.Bemessungswerte.n_N),'_',num2str(KennfeldMaschine_mechVerluste_Speichern.Analyse.Optionen.n_max),'_',num2str(KennfeldMaschine_mechVerluste_Speichern.Entwurf.Bemessungswerte.U_N),'_','Kennfeld_mechVerluste','.mat'],'Kennfeld_mechVerluste');
%         save (['LDS/Vehicle/Para_Powertrain/Motor_efficiency/MOT_memory//', KennfeldMaschine_mechVerluste_Speichern.Entwurf.Optionen.Maschinentyp,'_',num2str(KennfeldMaschine_mechVerluste_Speichern.Entwurf.Bemessungswerte.P_N/1000),'_',num2str(KennfeldMaschine_mechVerluste_Speichern.Entwurf.Bemessungswerte.n_N),'_',num2str(KennfeldMaschine_mechVerluste_Speichern.Analyse.Optionen.n_max),'_',num2str(KennfeldMaschine_mechVerluste_Speichern.Entwurf.Bemessungswerte.U_N),'_', 'Kennfeld_mechVerluste' ,'.mat'],'Kennfeld_mechVerluste');
        
        % Save map
        KennfeldMaschine_WirkZusVerl = struct;
        KennfeldMaschine_WirkZusVerl_Speichern = handles;
        KennfeldMaschine_WirkZusVerl.n = KennfeldMaschine_WirkZusVerl_Speichern.Analyse.Betriebsdaten.n_m_mesh;
        KennfeldMaschine_WirkZusVerl.M = KennfeldMaschine_WirkZusVerl_Speichern.Analyse.Betriebsdaten.M_max_mesh;
        KennfeldMaschine_WirkZusVerl.wirkzusverl = KennfeldMaschine_WirkZusVerl_Speichern.Analyse.Wirkungsgrad.eta_zus_mesh;
        %save (['3_Ergebnisse/',handles.Entwurf.Optionen.folder_id,'/','2_Analyse/', 'KennfeldMaschine_WirkZusVerl','.mat'],'KennfeldMaschine_WirkZusVerl');
        folder_id = [KennfeldMaschine_WirkZusVerl_Speichern.Entwurf.Optionen.folder_id]; %file_id, '_data'];
        file_id = [KennfeldMaschine_WirkZusVerl_Speichern.Entwurf.Optionen.file_id];
        
        % Convert data to correct form
        % Design sizes
        Kennfeld_WirkZusVerl = struct;
        Kennfeld_WirkZusVerl.type = KennfeldMaschine_WirkZusVerl_Speichern.Entwurf.Optionen.Maschinentyp;
        Kennfeld_WirkZusVerl.power = KennfeldMaschine_WirkZusVerl_Speichern.Entwurf.Bemessungswerte.P_N/1000;
        Kennfeld_WirkZusVerl.n_n = KennfeldMaschine_WirkZusVerl_Speichern.Entwurf.Bemessungswerte.n_N;
        Kennfeld_WirkZusVerl.n_max = KennfeldMaschine_WirkZusVerl_Speichern.Analyse.Optionen.n_max;
        Kennfeld_WirkZusVerl.U_n = KennfeldMaschine_WirkZusVerl_Speichern.Entwurf.Bemessungswerte.U_N;
        
        % Speed axis - X-axis map
        Kennfeld_WirkZusVerl.eff_n_axis = transpose(KennfeldMaschine_WirkZusVerl_Speichern.Analyse.Betriebsdaten.mot.n_m_vec);

        % Torque - Y-axis
        KennfeldMaschine_WirkZusVerl_Speichern.Analyse.Betriebsdaten.mot.M_vec=transpose(fliplr(KennfeldMaschine_WirkZusVerl_Speichern.Analyse.Betriebsdaten.mot.M_vec));
        Kennfeld_WirkZusVerl.eff_T_axis = [(fliplr(transpose(KennfeldMaschine_WirkZusVerl_Speichern.Analyse.Betriebsdaten.mot.M_vec(2:end)))) * (-1), transpose(KennfeldMaschine_WirkZusVerl_Speichern.Analyse.Betriebsdaten.mot.M_vec)];
        
        % Diagram                                                            
        % Manipulate efficiency map
        KennfeldMaschine_WirkZusVerl_Speichern.Analyse.Wirkungsgrad.eta_zus_mesh=flipud(KennfeldMaschine_WirkZusVerl_Speichern.Analyse.Wirkungsgrad.eta_zus_mesh); %hier Unterschied Verlustarten
        KennfeldMaschine_WirkZusVerl_Speichern.Analyse.Wirkungsgrad.eta_zus_mesh(KennfeldMaschine_WirkZusVerl_Speichern.Analyse.Wirkungsgrad.eta_zus_mesh<0.1) = 0.1; %Hier Unterschied Verlustarten %set small efficiency values to 0.1
        KennfeldMaschine_WirkZusVerl_Speichern.Analyse.Wirkungsgrad.eta_zus_mesh(1,:) = 0; %first colomn is zero %hier Unterschied Verlustarten
        KennfeldMaschine_WirkZusVerl_Speichern.Analyse.Wirkungsgrad.eta_zus_mesh(:,1) = 0.1; %first row is zero %hier Unterschied Verlustarten
        
        % Create one "larger" eff. map that includes acceleration and recuperation
        Kennfeld_WirkZusVerl.eff_recu = 1./KennfeldMaschine_WirkZusVerl_Speichern.Analyse.Wirkungsgrad.eta_zus_mesh; %recuperation efficiency %hier Unterschied Verlustarten
        Kennfeld_WirkZusVerl.eff_recu(Kennfeld_WirkZusVerl.eff_recu==inf)=100;
        Kennfeld_WirkZusVerl.eff_recu = Kennfeld_WirkZusVerl.eff_recu(2:end,:);
        Kennfeld_WirkZusVerl.eff = [flipud(Kennfeld_WirkZusVerl.eff_recu); KennfeldMaschine_WirkZusVerl_Speichern.Analyse.Wirkungsgrad.eta_zus_mesh]; %hier Unterschied Verlustarten
    
        % Full load characteristic
        KennfeldMaschine_WirkZusVerl_Speichern.Analyse.Betriebsdaten.mot.M_max_vec=transpose(KennfeldMaschine_WirkZusVerl_Speichern.Analyse.Betriebsdaten.mot.M_max_vec);
        Kennfeld_WirkZusVerl.T_max = transpose(KennfeldMaschine_WirkZusVerl_Speichern.Analyse.Betriebsdaten.mot.M_max_vec);  
        
        % Speed axis scaled
        Kennfeld_WirkZusVerl.T_max_n_axis = transpose(KennfeldMaschine_WirkZusVerl_Speichern.Analyse.Betriebsdaten.mot.M_vec);                                                          %noch auf 2000 interpolieren?
        
%         save (['3_Ergebnisse/',handles.Entwurf.Optionen.folder_id,'/','2_Analyse/', KennfeldMaschine_WirkZusVerl_Speichern.Entwurf.Optionen.Maschinentyp,'_',num2str(KennfeldMaschine_WirkZusVerl_Speichern.Entwurf.Bemessungswerte.P_N/1000),'_',num2str(KennfeldMaschine_WirkZusVerl_Speichern.Entwurf.Bemessungswerte.n_N),'_',num2str(KennfeldMaschine_WirkZusVerl_Speichern.Analyse.Optionen.n_max),'_',num2str(KennfeldMaschine_WirkZusVerl_Speichern.Entwurf.Bemessungswerte.U_N),'_','Kennfeld_WirkZusVerl','.mat'],'Kennfeld_WirkZusVerl');
%         save (['LDS/Vehicle/Para_Powertrain/Motor_efficiency/MOT_memory//', KennfeldMaschine_WirkZusVerl_Speichern.Entwurf.Optionen.Maschinentyp,'_',num2str(KennfeldMaschine_WirkZusVerl_Speichern.Entwurf.Bemessungswerte.P_N/1000),'_',num2str(KennfeldMaschine_WirkZusVerl_Speichern.Entwurf.Bemessungswerte.n_N),'_',num2str(KennfeldMaschine_WirkZusVerl_Speichern.Analyse.Optionen.n_max),'_',num2str(KennfeldMaschine_WirkZusVerl_Speichern.Entwurf.Bemessungswerte.U_N),'_', 'Kennfeld_WirkZusVerl' ,'.mat'],'Kennfeld_WirkZusVerl');
                        
        % Save map
        KennfeldMaschine_ZusVerl = struct;
        KennfeldMaschine_ZusVerl_Speichern = handles;
        KennfeldMaschine_ZusVerl.n = KennfeldMaschine_ZusVerl_Speichern.Analyse.Betriebsdaten.n_m_mesh;
        KennfeldMaschine_ZusVerl.M = KennfeldMaschine_ZusVerl_Speichern.Analyse.Betriebsdaten.M_max_mesh;
        KennfeldMaschine_ZusVerl.zusverl = KennfeldMaschine_ZusVerl_Speichern.Analyse.Verluste.P_vzus_mesh;
        %save (['3_Ergebnisse/', handles.Entwurf.Optionen.folder_id,'/','2_Analyse/', 'KennfeldMaschine_ZusVerl','.mat'],'KennfeldMaschine_ZusVerl');
        folder_id = [KennfeldMaschine_ZusVerl_Speichern.Entwurf.Optionen.folder_id]; %file_id, '_data'];
        file_id = [KennfeldMaschine_ZusVerl_Speichern.Entwurf.Optionen.file_id];
        
        % Convert data to correct form
        % Design sizes
        Kennfeld_ZusVerl = struct;
        Kennfeld_ZusVerl.type = KennfeldMaschine_ZusVerl_Speichern.Entwurf.Optionen.Maschinentyp;
        Kennfeld_ZusVerl.power = KennfeldMaschine_ZusVerl_Speichern.Entwurf.Bemessungswerte.P_N/1000;
        Kennfeld_ZusVerl.n_n = KennfeldMaschine_ZusVerl_Speichern.Entwurf.Bemessungswerte.n_N;
        Kennfeld_ZusVerl.n_max = KennfeldMaschine_ZusVerl_Speichern.Analyse.Optionen.n_max;
        Kennfeld_ZusVerl.U_n = KennfeldMaschine_ZusVerl_Speichern.Entwurf.Bemessungswerte.U_N;
        
        % Speed axis - X-axis map
        Kennfeld_ZusVerl.eff_n_axis = transpose(KennfeldMaschine_ZusVerl_Speichern.Analyse.Betriebsdaten.mot.n_m_vec);

        % Torque - Y-axis
        KennfeldMaschine_ZusVerl_Speichern.Analyse.Betriebsdaten.mot.M_vec=transpose(fliplr(KennfeldMaschine_ZusVerl_Speichern.Analyse.Betriebsdaten.mot.M_vec));
        Kennfeld_ZusVerl.eff_T_axis = [(fliplr(transpose(KennfeldMaschine_ZusVerl_Speichern.Analyse.Betriebsdaten.mot.M_vec(2:end)))) * (-1), transpose(KennfeldMaschine_ZusVerl_Speichern.Analyse.Betriebsdaten.mot.M_vec)];
        
        % Diagram                                                             
        % Manipulate efficiency map
        KennfeldMaschine_ZusVerl_Speichern.Analyse.Verluste.P_vzus_mesh=flipud(KennfeldMaschine_ZusVerl_Speichern.Analyse.Verluste.P_vzus_mesh); %hier Unterschied Verlustarten
        KennfeldMaschine_ZusVerl_Speichern.Analyse.Verluste.P_vzus_mesh(KennfeldMaschine_ZusVerl_Speichern.Analyse.Verluste.P_vzus_mesh<0.1) = 0.1; %Hier Unterschied Verlustarten %set small efficiency values to 0.1
        KennfeldMaschine_ZusVerl_Speichern.Analyse.Verluste.P_vzus_mesh(1,:) = 0; %first colomn is zero %hier Unterschied Verlustarten
        KennfeldMaschine_ZusVerl_Speichern.Analyse.Verluste.P_vzus_mesh(:,1) = 0.1; %first row is zero %hier Unterschied Verlustarten
        
        % Create one "larger" eff. map that includes acceleration and recuperation
        Kennfeld_ZusVerl.eff_recu = 1./KennfeldMaschine_ZusVerl_Speichern.Analyse.Verluste.P_vzus_mesh; %recuperation efficiency %hier Unterschied Verlustarten
        Kennfeld_ZusVerl.eff_recu(Kennfeld_ZusVerl.eff_recu==inf)=100;
        Kennfeld_ZusVerl.eff_recu = Kennfeld_ZusVerl.eff_recu(2:end,:);
        Kennfeld_ZusVerl.eff = [flipud(Kennfeld_ZusVerl.eff_recu); KennfeldMaschine_ZusVerl_Speichern.Analyse.Verluste.P_vzus_mesh]; %hier Unterschied Verlustarten
    
        % Full load characteristic
        KennfeldMaschine_ZusVerl_Speichern.Analyse.Betriebsdaten.mot.M_max_vec=transpose(KennfeldMaschine_ZusVerl_Speichern.Analyse.Betriebsdaten.mot.M_max_vec);
        Kennfeld_ZusVerl.T_max = transpose(KennfeldMaschine_ZusVerl_Speichern.Analyse.Betriebsdaten.mot.M_max_vec);  
        
        % Speed axis scaled
        Kennfeld_ZusVerl.T_max_n_axis = transpose(KennfeldMaschine_ZusVerl_Speichern.Analyse.Betriebsdaten.mot.M_vec);                                                          %noch auf 2000 interpolieren?

%         save (['3_Ergebnisse/',handles.Entwurf.Optionen.folder_id,'/','2_Analyse/', KennfeldMaschine_ZusVerl_Speichern.Entwurf.Optionen.Maschinentyp,'_',num2str(KennfeldMaschine_ZusVerl_Speichern.Entwurf.Bemessungswerte.P_N/1000),'_',num2str(KennfeldMaschine_ZusVerl_Speichern.Entwurf.Bemessungswerte.n_N),'_',num2str(KennfeldMaschine_ZusVerl_Speichern.Analyse.Optionen.n_max),'_',num2str(KennfeldMaschine_ZusVerl_Speichern.Entwurf.Bemessungswerte.U_N),'_','Kennfeld_ZusVerl','.mat'],'Kennfeld_ZusVerl');
%         save (['LDS/Vehicle/Para_Powertrain/Motor_efficiency/MOT_memory//', KennfeldMaschine_ZusVerl_Speichern.Entwurf.Optionen.Maschinentyp,'_',num2str(KennfeldMaschine_ZusVerl_Speichern.Entwurf.Bemessungswerte.P_N/1000),'_',num2str(KennfeldMaschine_ZusVerl_Speichern.Entwurf.Bemessungswerte.n_N),'_',num2str(KennfeldMaschine_ZusVerl_Speichern.Analyse.Optionen.n_max),'_',num2str(KennfeldMaschine_ZusVerl_Speichern.Entwurf.Bemessungswerte.U_N),'_', 'Kennfeld_ZusVerl' ,'.mat'],'Kennfeld_ZusVerl');
        
        % Save map
        KennfeldMaschine_i_ld = handles.Analyse.Momentensteuerung.i_1d_mesh;
%         save (['3_Ergebnisse/',handles.Entwurf.Optionen.folder_id,'/','2_Analyse/', 'KennfeldMaschine_i_ld','.mat'],'KennfeldMaschine_i_ld');

        % Save map
        KennfeldMaschine_i_lq = handles.Analyse.Momentensteuerung.i_1q_mesh;
%         save (['3_Ergebnisse/',handles.Entwurf.Optionen.folder_id,'/','2_Analyse/', 'KennfeldMaschine_i_lq','.mat'],'KennfeldMaschine_i_lq');
    
        % Save map
        KennfeldMaschine_u_ld = handles.Analyse.Momentensteuerung.u_1d_mesh;
%         save (['3_Ergebnisse/',handles.Entwurf.Optionen.folder_id,'/','2_Analyse/', 'KennfeldMaschine_u_ld','.mat'],'KennfeldMaschine_u_ld');
    
        % Save map
        KennfeldMaschine_u_lq = handles.Analyse.Momentensteuerung.u_1q_mesh;
%         save (['3_Ergebnisse/',handles.Entwurf.Optionen.folder_id,'/','2_Analyse/', 'KennfeldMaschine_u_lq','.mat'],'KennfeldMaschine_u_lq');
        


end


%% Additional function
% Load material
function data = loadMaterial(var,typ)

    if(strcmp(var.String,'-'))
        data.String = var.String;
        return
    end
    switch typ
        case 'Elektroblech'
            filepath = '5_Material/1_ElectricalSheet/';
        case 'Leiter'
            filepath = '5_Material/2_Conductor/';
        case 'Magnet'
            filepath = '5_Material/3_Magnet/';
        otherwise
            error('Undefined material type')
    end

    load([filepath var.String '.mat']);
    data.String = var.String;
end
